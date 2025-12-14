#include "simulation.hpp"
#include "simulation_helpers.hpp"
#include "utils.hpp"
#include "logger.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <functional>
#include <cstdint>
#include <algorithm>
#include <random>
#include <set>

// Structure to hold simulation state
struct SimulationState {
    std::map<int, DiskInfo> disks;
    std::map<std::string, bool> failed_nodes_and_enclosures;
    std::map<int, ECGroupFailureInfo> failure_info_per_ec_group;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> failed_events;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> repair_events;
    std::unordered_map<NodeFailureKey, GraphStructure, NodeFailureKeyHash> failed_hardware_graph_table;
    std::unordered_map<NodeFailureKey, DisconnectedStatus, NodeFailureKeyHash> disconnected_table;
    std::unordered_map<FailureStateKey, FlowsAndSpeedEntry, FailureStateKeyHash> flows_and_speed_table;
    std::unordered_map<std::string, int> node_index_map;
    int node_count = 0;
    int total_group_count = 0;

    double up_time = 0.0;
    double perf_up_time_read = 0.0;
    double perf_up_time_write = 0.0;
    double perf_up_time_min = 0.0;
    double timestamp = 0.0;
    double prev_time = 0.0;

    double total_mttdl = 0.0;
    int data_loss_events = 0;
    int disk_failure_events = 0;
    double total_data_loss_ratio = 0.0;
    double prev_up_time = 0.0;
    bool prev_had_data_loss = false;

    // Year-based durability tracking (per EC group-year)
    int current_year = 0;                          // Current year in simulation
    std::set<int> groups_with_loss_this_year;      // EC groups that had data loss this year
    int total_group_years = 0;                     // Total group-years simulated
    int group_years_with_data_loss = 0;            // Group-years that had data loss

    double total_time_for_rebuilding = 0.0;
    int count_for_rebuilding = 0;

    double total_read_bandwidth = 0.0;
    double total_write_bandwidth = 0.0;
    double bandwidth_time = 0.0;

    // Rebuild-specific metrics
    double total_rebuild_speed_time = 0.0;      // Sum of (rebuild_speed * duration)
    double total_time_in_rebuild = 0.0;         // Total time any disk is rebuilding
    double total_host_read_bw_during_rebuild = 0.0;   // Sum of (read_bw * duration) during rebuild
    double total_host_write_bw_during_rebuild = 0.0;  // Sum of (write_bw * duration) during rebuild

    DisconnectedStatus last_disconnected;
};

// Initialize simulation state
SimulationState initialize_simulation_state(
    const SimulationParams& params,
    const GraphStructure& hardware_graph,
    const ErasureCodingScheme& ec_scheme,
    const std::map<std::string, std::string>& node_to_module_map,
    const DiskIOModuleManager& disk_io_manager,
    const nlohmann::json& options,
    const EmpiricalCDF* disk_failure_cdf = nullptr,
    double disk_arr = 0.0) {

    SimulationState state;

    // Initialize disks
    for (int index = 0; index < params.total_disks; ++index) {
        std::string io_module = disk_io_manager.get_primary_io_module(index);
        state.disks[index] = DiskInfo(index, io_module, params.disk_read_bw, params.disk_write_bw);
    }

    // Initialize failed nodes and build node index map
    auto nodes = hardware_graph.get_nodes();
    state.node_count = static_cast<int>(nodes.size());
    for (int index = 0; index < state.node_count; ++index) {
        const std::string& node = nodes[index];
        if (!is_disk_node(node)) {
            state.failed_nodes_and_enclosures[node] = false;
            state.node_index_map[node] = index;
        }
    }

    // Add enclosures to node index map
    for (const auto& [enclosure, _] : hardware_graph.enclosures_) {
        state.failed_nodes_and_enclosures[enclosure] = false;
        state.node_index_map[enclosure] = state.node_count++;
    }

    // Initialize failure info per group
    state.total_group_count = ec_scheme.get_total_groups(params.total_disks);
    for (int group_index = 0; group_index < state.total_group_count; ++group_index) {
        state.failure_info_per_ec_group[group_index] = ECGroupFailureInfo();
    }

    // Generate initial failure events
    state.failed_events = generate_first_failure_events(
        hardware_graph, node_to_module_map, params.total_disks, params.disk_mttf,
        disk_failure_cdf, disk_arr);

    return state;
}

// Process a single simulation event
void process_simulation_event(
    const Event& event,
    SimulationState& state,
    const SimulationParams& params,
    const ErasureCodingScheme& ec_scheme,
    const GraphStructure& hardware_graph,
    const std::map<std::string, std::string>& node_to_module_map,
    const std::map<std::string, std::vector<std::string>>& enclosure_to_node_map,
    double max_read_performance_without_any_failure,
    double max_write_performance_without_any_failure,
    const nlohmann::json& options,
    const DiskIOModuleManager* disk_io_manager,
    const EmpiricalCDF* disk_failure_cdf = nullptr,
    double disk_arr = 0.0) {

    std::string event_type = event.event_type;
    std::string event_node = event.event_node;
    double event_time = event.time;
    double repair_start_time = event.start_time;
    bool software_repair_only = event.software_repair_only;

    // Track rebuild duration per disk (failure or reconnect resync).
    if ((event_type == "repair" || event_type == "rebuild_complete") && is_disk_node(event_node)) {
        int disk_index = get_disk_index(event_node);
        auto it = state.disks.find(disk_index);
        if (it != state.disks.end()) {
            DiskInfo& disk = it->second;
            if (disk.rebuild_start_timestamp > 0 && event_time >= disk.rebuild_start_timestamp) {
                state.total_time_for_rebuilding += (event_time - disk.rebuild_start_timestamp);
                state.count_for_rebuilding++;
            }
            disk.rebuild_start_timestamp = 0.0;
        }
    }

    // Update failure info
    update_failure_info(event_type, event_node, state.failure_info_per_ec_group,
                       state.failed_nodes_and_enclosures, ec_scheme);

    // Recalculate hardware and flow state
    NodeFailureKey node_key = build_node_failure_key(
        state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);

    calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, node_key,
                             enclosure_to_node_map, options, state.failed_hardware_graph_table,
                             state.disconnected_table, disk_io_manager, params.total_disks);
    DisconnectedStatus disconnected = state.disconnected_table[node_key];

    // Update disk state
    if (is_disk_node(event_node)) {
        int disk_index = get_disk_index(event_node);
        update_disk_state(disk_index, state.failure_info_per_ec_group, state.disks,
                         params.disk_capacity, event_type, event_time, params.disk_replace_time,
                         ec_scheme, disconnected);
    }

    // Handle disconnection state changes
    bool disconnection_changed = false;
    for (int i = 0; i < params.total_disks; ++i) {
        bool was_disconnected = state.last_disconnected.is_disk_disconnected(i);
        bool is_disconnected = disconnected.is_disk_disconnected(i);
        if (was_disconnected != is_disconnected) {
            disconnection_changed = true;

            // Handle disconnect event
            if (is_disconnected && !state.disks[i].is_failed) {
                state.disks[i].is_disconnected = true;
                state.disks[i].disconnect_timestamp = event_time;
            }
            // Handle reconnect event
            else if (!is_disconnected && state.disks[i].is_disconnected && !state.disks[i].is_failed) {
                double disconnect_duration = event_time - state.disks[i].disconnect_timestamp;
                state.disks[i].remaining_capacity_to_rebuild = disconnect_duration * params.disk_write_bw * 3600.0;
                state.disks[i].is_disconnected = false;
                state.disks[i].disconnect_timestamp = 0.0;
                state.disks[i].rebuild_start_timestamp = event_time;
                if (state.disks[i].remaining_capacity_to_rebuild > 0) {
                    // Re-sync rebuild after reconnect should NOT create a "repair" event that reschedules failures.
                    // It is a rebuild completion event for remaining_capacity_to_rebuild only.
                    Event ev{event_time + BIG_NUMBER, "rebuild_complete", get_disk_name(i), event_time, event_time, false};
                    state.repair_events.push(ev);
                }
            }
        }
    }

    if (disconnection_changed) {
        update_all_disk_states(state.failure_info_per_ec_group, state.disks, ec_scheme, disconnected, event_time);
        state.last_disconnected = disconnected;
    }

    // Handle failure events
    if (event_type == "fail") {
        push_repair_event(state.repair_events, event_node, event_time, node_to_module_map,
                         hardware_graph, software_repair_only, &options);
        // Count disk failures
        if (is_disk_node(event_node)) {
            state.disk_failure_events++;
        }
    }

    // Handle repair events
    if (event_type == "repair") {
        push_failed_event(state.failed_events, event_node, event_time, node_to_module_map,
                         hardware_graph, params.disk_mttf, disk_failure_cdf, disk_arr);
    }

    // Handle rebuild completion after reconnect (no new failure scheduling)
    if (event_type == "rebuild_complete" && is_disk_node(event_node)) {
        int disk_index = get_disk_index(event_node);
        DiskInfo& disk_info = state.disks[disk_index];
        disk_info.remaining_capacity_to_rebuild = 0;
        disk_info.rebuild_speed = 0;
        disk_info.state = DiskState::NORMAL;
    }

    // Update flow entry after disk/disconnection state changes so rebuild speeds reflect this event.
    FailureStateKey flow_key = build_failure_state_key(
        state.failed_nodes_and_enclosures,
        state.disks,
        state.node_index_map, state.node_count, params.total_disks);

    calculate_flows_and_speed(state.failed_hardware_graph_table[node_key], state.failure_info_per_ec_group,
                              ec_scheme, params, state.flows_and_speed_table,
                              max_read_performance_without_any_failure, max_write_performance_without_any_failure,
                              options.value("start_module", "switch"), options.value("end_module", "io_module"),
                              disconnected, flow_key, disk_io_manager, &state.disks);
    const FlowsAndSpeedEntry& flow_entry = state.flows_and_speed_table[flow_key];

    // Update repair events for disks
    update_repair_event_for_disks(state.repair_events, event_time, state.disks, flow_entry,
                                  state.failure_info_per_ec_group, ec_scheme);
}

// Single-core simulation function
SimulationResult simulation_per_core(
    int simulation_idx,
    const std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int num_iterations,
    const nlohmann::json& options) {

    // Extract parameters
    SimulationParams params;

    // EC configuration
    params.ec_config.type = parse_ec_type(options.value("ec_type", "standard"));
    params.ec_config.m = params_and_results.at("m");
    params.ec_config.k = params_and_results.at("k");
    params.ec_config.n = options.value("n", params.ec_config.m + params.ec_config.k);
    params.ec_config.local_m = options.value("local_m", 0);
    params.ec_config.local_k = options.value("local_k", 0);
    params.ec_config.local_n = options.value("local_n", 0);
    params.ec_config.local_type = parse_local_type(options.value("local_type", "ec"));
    params.ec_config.validate();

    // Disk parameters
    params.total_disks = params_and_results.at("total_ssds");
    params.disk_capacity = params_and_results.at("capacity");
    params.disk_read_bw = params_and_results.at("ssd_read_bw");
    params.disk_write_bw = params_and_results.at("ssd_write_bw");

    // Calculate disk MTTF
    double dwpd = params_and_results.at("dwpd");
    double dwpd_limit = params_and_results.at("dwpd_limit");
    int guaranteed_years = params_and_results.at("guaranteed_years");
    params.disk_mttf = guaranteed_years * 365.0 * 24.0 * dwpd_limit / dwpd;

    params.disk_replace_time = options.value("prep_time_for_rebuilding", 1.0);

    // Rebuild parameters
    params.rebuild_bw_ratio = options.value("rebuild_bw_ratio", 0.2);
    params.degraded_ratio = options.value("degraded_ratio", 0.2);
    params.data_loss_rebuild_bw = options.value("data_loss_rebuild_bw", 0.0);

    // Performance availability threshold
    params.target_performance_ratio = options.value("target_performance_ratio", 0.0);

    // Simulation control
    params.nprocs = params_and_results.at("nprocs");
    params.simulation_years = options.value("simulation_years", 10);

    // EC encoding speed (based on 256KB chunk encoding latency)
    const double EC_CHUNK_SIZE = 256.0 * 1024.0;  // 256KB per chunk
    double ec_latency_sec = Utils::get_encoding_latency_sec(params.ec_config.m, params.ec_config.k);
    if (ec_latency_sec > 0) {
        // Speed = chunk_size / latency_per_chunk (bytes/sec)
        params.ec_encoding_speed = EC_CHUNK_SIZE / ec_latency_sec;
    } else {
        params.ec_encoding_speed = 1e15;  // Very large number (no bottleneck)
    }

    // Create EC scheme
    ErasureCodingScheme ec_scheme;

    // Check if custom ec_groups are defined in the options
    if (options.contains("ec_groups") && options["ec_groups"].is_array()) {
        std::vector<std::vector<int>> ec_groups;
        for (const auto& group : options["ec_groups"]) {
            if (group.is_array()) {
                std::vector<int> disk_indices;
                for (const auto& disk : group) {
                    if (disk.is_number_integer()) {
                        disk_indices.push_back(disk.get<int>());
                    }
                }
                ec_groups.push_back(disk_indices);
            }
        }
        if (!ec_groups.empty()) {
            if (simulation_idx == 0) {
                Logger::getInstance().info("Using custom EC groups from config (" +
                    std::to_string(ec_groups.size()) + " groups)");
            }
            ec_scheme.initialize_with_custom_groups(params.ec_config, ec_groups);
        } else {
            ec_scheme.initialize(params.ec_config);
        }
    } else {
        ec_scheme.initialize(params.ec_config);
    }

    // Set number of EC groups for durability calculation
    params.num_ec_groups = ec_scheme.get_total_groups(params.total_disks);

    // Load empirical CDF for disk failure times if trace file is provided
    EmpiricalCDF disk_failure_cdf;
    const EmpiricalCDF* disk_failure_cdf_ptr = nullptr;
    std::string trace_file = options.value("ssd_failure_trace", "");
    if (!trace_file.empty()) {
        if (disk_failure_cdf.load_from_file(trace_file)) {
            disk_failure_cdf_ptr = &disk_failure_cdf;
            if (simulation_idx == 0) {
                Logger::getInstance().info("Using empirical disk failure distribution from trace");
            }
        } else {
            Logger::getInstance().warning("Failed to load trace file, using exponential distribution");
        }
    }

    // Load disk ARR (Annual Replacement Rate) for spliced distribution
    double disk_arr = options.value("ssd_arr", 0.0);
    if (disk_failure_cdf_ptr != nullptr && disk_arr > 0 && simulation_idx == 0) {
        Logger::getInstance().info("Using spliced distribution with ARR=" + std::to_string(disk_arr));
    }

    // Initialize simulation
    GraphStructure hardware_graph = graph_structure_origin.clone();
    auto [node_to_module_map, enclosure_to_node_map, max_read_performance_without_any_failure, max_write_performance_without_any_failure] =
        initialize_simulation(hardware_graph, params.disk_read_bw, params.total_disks, options);

    // Initialize Disk-IO Module manager
    DiskIOModuleManager disk_io_manager;
    std::vector<DiskIOModuleMapping> disk_io_mappings;

    // Build disk-IO module mapping from enclosures
    // enclosure_to_node_map contains: enclosure_name -> [io_module0, io_module1, ...]
    if (!enclosure_to_node_map.empty()) {
        // Sort enclosures by name to ensure consistent ordering
        std::vector<std::string> enclosure_names;
        for (const auto& [enc_name, io_modules] : enclosure_to_node_map) {
            enclosure_names.push_back(enc_name);
        }
        std::sort(enclosure_names.begin(), enclosure_names.end());

        // Calculate SSDs per enclosure
        int num_enclosures = static_cast<int>(enclosure_names.size());
        int ssds_per_enclosure = (num_enclosures > 0) ? (params.total_disks / num_enclosures) : params.total_disks;

        for (const auto& enc_name : enclosure_names) {
            const auto& io_modules = enclosure_to_node_map.at(enc_name);
            DiskIOModuleMapping mapping;
            mapping.disk_count = ssds_per_enclosure;
            mapping.io_modules = io_modules;
            mapping.start_disk_index = 0;  // Will be computed by initialize()
            disk_io_mappings.push_back(mapping);
        }

        if (simulation_idx == 0 && !disk_io_mappings.empty()) {
            Logger::getInstance().info("Using enclosure-based disk-IO module mapping (" +
                                       std::to_string(disk_io_mappings.size()) + " enclosures, " +
                                       std::to_string(ssds_per_enclosure) + " SSDs each)");
        }
    }

    // Collect all io_modules from hardware graph
    std::set<std::string> all_io_modules;
    for (const auto& node : hardware_graph.get_nodes()) {
        if (node.find("io_module") != std::string::npos) {
            all_io_modules.insert(node);
        }
    }

    disk_io_manager.initialize(disk_io_mappings, params.total_disks, all_io_modules);

    // Simulation parameters (per iteration)
    double simulation_hours = params.simulation_years * 365.0 * 24.0;

    SimulationResult aggregated;
    aggregated.runs = num_iterations;

    int percentage = 1;
    if (simulation_idx == 0) {
        Logger::getInstance().info("Simulation " + std::to_string(simulation_idx) + " started");
        Utils::progress_bar(percentage, 100, 50);
    }

    for (int iter = 0; iter < num_iterations; ++iter) {
        // Initialize state (fresh, independent simulation)
        SimulationState state = initialize_simulation_state(params, hardware_graph, ec_scheme,
                                                            node_to_module_map, disk_io_manager,
                                                            options, disk_failure_cdf_ptr, disk_arr);

        // Initial state calculation
        NodeFailureKey initial_node_key = build_node_failure_key(
            state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);
        calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, initial_node_key,
                                 enclosure_to_node_map, options, state.failed_hardware_graph_table,
                                 state.disconnected_table, &disk_io_manager, params.total_disks);

        FailureStateKey flow_key = build_failure_state_key(
            state.failed_nodes_and_enclosures,
            state.disks,
            state.node_index_map, state.node_count, params.total_disks);

        calculate_flows_and_speed(state.failed_hardware_graph_table[initial_node_key], state.failure_info_per_ec_group,
                                  ec_scheme, params, state.flows_and_speed_table,
                                  max_read_performance_without_any_failure, max_write_performance_without_any_failure,
                                  options.value("start_module", "switch"), options.value("end_module", "io_module"),
                                  state.disconnected_table[initial_node_key], flow_key, &disk_io_manager, &state.disks);

        // Main event loop (run for fixed horizon)
        while (true) {
            Event event = pop_event(state.failed_events, state.repair_events);
            double event_time = event.time;

            bool break_flag = false;
            if (event_time > simulation_hours) {
                event_time = simulation_hours;
                break_flag = true;
            }

            double time_diff = event_time - state.prev_time;

            // Calculate availability for this time period
            const FlowsAndSpeedEntry& flows_and_speed_entry = state.flows_and_speed_table[flow_key];
            double availability_ratio = flows_and_speed_entry.availability_ratio;
            double data_loss_ratio = flows_and_speed_entry.data_loss_ratio;

            // Track MTTDL (time to data loss)
            if (!state.prev_had_data_loss && data_loss_ratio > 0) {
                // Data loss just occurred
                state.total_mttdl += (state.prev_time - state.prev_up_time);
                state.data_loss_events++;
                state.prev_had_data_loss = true;
            } else if (state.prev_had_data_loss && data_loss_ratio == 0) {
                // Recovered from data loss
                state.prev_up_time = state.prev_time;
                state.prev_had_data_loss = false;
            }

            // Update availability and data loss tracking
            state.up_time += time_diff * availability_ratio;
            state.total_data_loss_ratio += time_diff * data_loss_ratio;
            state.timestamp += time_diff;

            // Year-based durability tracking
            constexpr double hours_per_year = 8760.0;
            int new_year = static_cast<int>(event_time / hours_per_year);

            // Check if year changed - finalize previous year(s)
            while (state.current_year < new_year) {
                // Finalize the current year
                state.group_years_with_data_loss += static_cast<int>(state.groups_with_loss_this_year.size());
                state.total_group_years += params.num_ec_groups;
                state.groups_with_loss_this_year.clear();
                state.current_year++;
            }

            // Track which EC groups had data loss this year
            for (int group_idx : flows_and_speed_entry.data_loss_group_indices) {
                state.groups_with_loss_this_year.insert(group_idx);
            }

            // Update performance-based availability
            if (params.target_performance_ratio > 0) {
                if (flows_and_speed_entry.performance_ratio_read >= params.target_performance_ratio) {
                    state.perf_up_time_read += time_diff;
                }
                if (flows_and_speed_entry.performance_ratio_write >= params.target_performance_ratio) {
                    state.perf_up_time_write += time_diff;
                }
                if (flows_and_speed_entry.performance_ratio_min >= params.target_performance_ratio) {
                    state.perf_up_time_min += time_diff;
                }
            }

            // Track bandwidth
            state.total_read_bandwidth += time_diff * flows_and_speed_entry.read_bandwidth;
            state.total_write_bandwidth += time_diff * flows_and_speed_entry.write_bandwidth;
            state.bandwidth_time += time_diff;

            // Track rebuild-specific metrics: check if any disk is rebuilding
            bool any_rebuilding = false;
            double max_rebuild_speed = 0.0;
            for (const auto& [disk_idx, disk_info] : state.disks) {
                if (disk_info.needs_rebuild() && !disk_info.is_disconnected && disk_info.rebuild_speed > 0) {
                    any_rebuilding = true;
                    max_rebuild_speed = std::max(max_rebuild_speed, disk_info.rebuild_speed);
                }
            }
            if (any_rebuilding) {
                state.total_time_in_rebuild += time_diff;
                state.total_rebuild_speed_time += time_diff * max_rebuild_speed;
                // Track both read and write bandwidth during rebuild
                state.total_host_read_bw_during_rebuild += time_diff * flows_and_speed_entry.read_bandwidth;
                state.total_host_write_bw_during_rebuild += time_diff * flows_and_speed_entry.write_bandwidth;
            }

            if (break_flag) break;

            // Process event
            process_simulation_event(event, state, params, ec_scheme, hardware_graph,
                                     node_to_module_map, enclosure_to_node_map,
                                     max_read_performance_without_any_failure, max_write_performance_without_any_failure,
                                     options, &disk_io_manager, disk_failure_cdf_ptr, disk_arr);

            state.prev_time = event_time;

            // Update flow key after event processing
            NodeFailureKey node_key = build_node_failure_key(
                state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);
            calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, node_key,
                                     enclosure_to_node_map, options, state.failed_hardware_graph_table,
                                     state.disconnected_table, &disk_io_manager, params.total_disks);

            flow_key = build_failure_state_key(
                state.failed_nodes_and_enclosures,
                state.disks,
                state.node_index_map, state.node_count, params.total_disks);

            calculate_flows_and_speed(state.failed_hardware_graph_table[node_key], state.failure_info_per_ec_group,
                                      ec_scheme, params, state.flows_and_speed_table,
                                      max_read_performance_without_any_failure, max_write_performance_without_any_failure,
                                      options.value("start_module", "switch"), options.value("end_module", "io_module"),
                                      state.disconnected_table[node_key], flow_key, &disk_io_manager, &state.disks);
        }

        // Finalize the last year for durability tracking
        state.group_years_with_data_loss += static_cast<int>(state.groups_with_loss_this_year.size());
        state.total_group_years += params.num_ec_groups;

        // Aggregate per-iteration outputs.
        aggregated.up_time += state.up_time;
        aggregated.simulation_time += state.timestamp;
        aggregated.perf_up_time_read += state.perf_up_time_read;
        aggregated.perf_up_time_write += state.perf_up_time_write;
        aggregated.perf_up_time_min += state.perf_up_time_min;
        aggregated.mttdl += state.total_mttdl;
        aggregated.data_loss_events += state.data_loss_events;
        aggregated.disk_failure_events += state.disk_failure_events;
        aggregated.total_data_loss_ratio += (state.total_data_loss_ratio / (state.timestamp > 0 ? state.timestamp : 1.0));
        aggregated.time_for_rebuilding += state.total_time_for_rebuilding;
        aggregated.count_for_rebuilding += state.count_for_rebuilding;
        aggregated.total_group_years += state.total_group_years;
        aggregated.group_years_with_data_loss += state.group_years_with_data_loss;

        if (state.bandwidth_time > 0) {
            aggregated.average_read_bandwidth += state.total_read_bandwidth / state.bandwidth_time;
            aggregated.average_write_bandwidth += state.total_write_bandwidth / state.bandwidth_time;
        }

        // Aggregate rebuild-specific metrics
        aggregated.total_time_in_rebuild += state.total_time_in_rebuild;
        aggregated.total_rebuild_speed_time += state.total_rebuild_speed_time;
        aggregated.total_host_read_bw_during_rebuild += state.total_host_read_bw_during_rebuild;
        aggregated.total_host_write_bw_during_rebuild += state.total_host_write_bw_during_rebuild;

        // Progress bar update (best-effort; shows only thread0's portion)
        if (simulation_idx == 0) {
            int new_percentage = static_cast<int>((iter + 1) * 100.0 / num_iterations);
            if (new_percentage > percentage) {
                percentage = new_percentage;
                Utils::progress_bar(percentage, 100, 50);
            }
        }
    }

    return aggregated;
}
