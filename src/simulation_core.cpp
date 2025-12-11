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
    double perf_up_time = 0.0;  // Time with performance >= target
    double timestamp = 0.0;
    double prev_time = 0.0;

    double total_mttdl = 0.0;
    int data_loss_events = 0;
    double total_data_loss_ratio = 0.0;
    double prev_up_time = 0.0;
    bool prev_had_data_loss = false;

    double total_time_for_rebuilding = 0.0;
    int count_for_rebuilding = 0;
    double last_timestamp_for_rebuilding = 0.0;

    double total_read_bandwidth = 0.0;
    double total_write_bandwidth = 0.0;
    double bandwidth_time = 0.0;

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
    const nlohmann::json& options,
    const DiskIOModuleManager* disk_io_manager,
    const EmpiricalCDF* disk_failure_cdf = nullptr,
    double disk_arr = 0.0) {

    std::string event_type = event.event_type;
    std::string event_node = event.event_node;
    double event_time = event.time;
    double repair_start_time = event.start_time;
    bool software_repair_only = event.software_repair_only;

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

    FailureStateKey flow_key = build_failure_state_key(
        state.failed_nodes_and_enclosures, state.failure_info_per_ec_group,
        state.node_index_map, state.node_count, state.total_group_count);

    calculate_flows_and_speed(state.failed_hardware_graph_table[node_key], state.failure_info_per_ec_group,
                              ec_scheme, params, state.flows_and_speed_table,
                              max_read_performance_without_any_failure, disconnected, flow_key, disk_io_manager);
    const FlowsAndSpeedEntry& flow_entry = state.flows_and_speed_table[flow_key];

    // Update disk state
    if (is_disk_node(event_node)) {
        int disk_index = get_disk_index(event_node);
        update_disk_state(disk_index, state.failure_info_per_ec_group, state.disks,
                         params.disk_capacity, event_type, params.disk_replace_time,
                         ec_scheme, disconnected);

        if (event_type == "fail") {
            state.last_timestamp_for_rebuilding = event_time;
        }
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
    }

    // Handle repair events
    if (event_type == "repair") {
        if (is_disk_node(event_node)) {
            state.total_time_for_rebuilding += event_time - state.last_timestamp_for_rebuilding;
            state.count_for_rebuilding++;
        }

        push_failed_event(state.failed_events, event_node, event_time, node_to_module_map,
                         hardware_graph, params.disk_mttf, disk_failure_cdf, disk_arr);
    }

    // Update repair events for disks
    update_repair_event_for_disks(state.repair_events, event_time, state.disks, flow_entry,
                                  state.failure_info_per_ec_group, ec_scheme);
}

// Single-core simulation function
SimulationResult simulation_per_core(
    int simulation_idx,
    const std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int batch_size,
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
    ec_scheme.initialize(params.ec_config);

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
    auto [node_to_module_map, enclosure_to_node_map, max_read_performance_without_any_failure] =
        initialize_simulation(hardware_graph, params.disk_read_bw, params.total_disks, options);

    // Initialize Disk-IO Module manager
    DiskIOModuleManager disk_io_manager;
    std::vector<DiskIOModuleMapping> disk_io_mappings;

    // Parse disk_io_modules from options if present
    if (options.contains("ssd_io_modules") && options["ssd_io_modules"].is_array()) {
        for (const auto& entry : options["ssd_io_modules"]) {
            DiskIOModuleMapping mapping;
            mapping.disk_count = entry.value("ssds", 0);
            mapping.io_modules = entry.value("io_modules", std::vector<std::string>{});
            mapping.start_disk_index = 0;  // Will be computed by initialize()
            disk_io_mappings.push_back(mapping);
        }
        if (simulation_idx == 0 && !disk_io_mappings.empty()) {
            Logger::getInstance().info("Using disk-IO module mapping from config (" +
                                       std::to_string(disk_io_mappings.size()) + " mapping groups)");
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

    // Initialize state
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
        state.failed_nodes_and_enclosures, state.failure_info_per_ec_group,
        state.node_index_map, state.node_count, state.total_group_count);

    calculate_flows_and_speed(state.failed_hardware_graph_table[initial_node_key], state.failure_info_per_ec_group,
                              ec_scheme, params, state.flows_and_speed_table,
                              max_read_performance_without_any_failure, state.disconnected_table[initial_node_key],
                              flow_key, &disk_io_manager);

    // Simulation parameters
    double simulation_hours = params.simulation_years * 365.0 * 24.0;
    int percentage = 1;

    if (simulation_idx == 0) {
        Logger::getInstance().info("Simulation " + std::to_string(simulation_idx) + " started");
        Utils::progress_bar(percentage, 100, 50);
    }

    // Main event loop
    while (true) {
        Event event = pop_event(state.failed_events, state.repair_events);
        double event_time = event.time;

        // Check if simulation should end
        bool break_flag = false;
        if (event_time > simulation_hours * batch_size) {
            event_time = simulation_hours * batch_size;
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

        // Update performance-based availability
        double perf_ratio = flows_and_speed_entry.performance_ratio;
        if (params.target_performance_ratio > 0 && perf_ratio >= params.target_performance_ratio) {
            state.perf_up_time += time_diff;
        }

        // Track bandwidth
        state.total_read_bandwidth += time_diff * flows_and_speed_entry.read_bandwidth;
        state.total_write_bandwidth += time_diff * flows_and_speed_entry.write_bandwidth;
        state.bandwidth_time += time_diff;

        if (break_flag) break;

        // Process event
        process_simulation_event(event, state, params, ec_scheme, hardware_graph,
                                node_to_module_map, enclosure_to_node_map,
                                max_read_performance_without_any_failure, options,
                                &disk_io_manager, disk_failure_cdf_ptr, disk_arr);

        state.prev_time = event_time;

        // Update flow key after event processing
        NodeFailureKey node_key = build_node_failure_key(
            state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);
        calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, node_key,
                                 enclosure_to_node_map, options, state.failed_hardware_graph_table,
                                 state.disconnected_table, &disk_io_manager, params.total_disks);

        flow_key = build_failure_state_key(
            state.failed_nodes_and_enclosures, state.failure_info_per_ec_group,
            state.node_index_map, state.node_count, state.total_group_count);

        calculate_flows_and_speed(state.failed_hardware_graph_table[node_key], state.failure_info_per_ec_group,
                                  ec_scheme, params, state.flows_and_speed_table,
                                  max_read_performance_without_any_failure, state.disconnected_table[node_key],
                                  flow_key, &disk_io_manager);

        // Progress bar update
        if (simulation_idx == 0) {
            int new_percentage = static_cast<int>(event_time * 100.0 / (simulation_hours * batch_size));
            if (new_percentage > percentage) {
                percentage = new_percentage;
                Utils::progress_bar(percentage, 100, 50);
            }
        }
    }

    // Return result
    SimulationResult result;
    result.up_time = state.up_time;
    result.perf_up_time = state.perf_up_time;
    result.simulation_time = state.timestamp;
    result.availability = (state.timestamp > 0) ? (state.up_time / state.timestamp) : 1.0;
    result.perf_availability = (state.timestamp > 0) ? (state.perf_up_time / state.timestamp) : 1.0;
    result.mttdl = state.total_mttdl;
    result.data_loss_events = state.data_loss_events;
    result.total_data_loss_ratio = state.total_data_loss_ratio / (state.timestamp > 0 ? state.timestamp : 1.0);
    result.time_for_rebuilding = state.total_time_for_rebuilding;
    result.count_for_rebuilding = state.count_for_rebuilding;

    // Average bandwidth
    if (state.bandwidth_time > 0) {
        result.average_read_bandwidth = state.total_read_bandwidth / state.bandwidth_time;
        result.average_write_bandwidth = state.total_write_bandwidth / state.bandwidth_time;
    }

    return result;
}
