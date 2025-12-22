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

// Structure to hold simulation parameters
struct SimulationCoreParams {
    int m, k, l;
    int cached_m, cached_k, cached_l;
    int cached_network_m, cached_network_k, cached_network_l;
    int total_ssds, cached_ssds;
    int network_m, network_k, network_l;
    int inter_replicas, intra_replicas;
    double dwpd, dwpd_limit, cached_dwpd_limit;
    double cached_write_ratio;
    int guaranteed_years;
    double read_bw, write_bw;
    double cached_read_bw, cached_write_bw;
    double capacity;
    double box_mttf, io_module_mttr;
    bool qlc, qlc_cache;
    double mttf, cached_mttf;
    double cached_ssd_cost, uncached_ssd_cost;
    double rebuild_bw_ratio;
    double target_perf_ratio;
    bool active_active;
    bool single_port_ssd;
};

// Extract and calculate simulation parameters
SimulationCoreParams extract_simulation_params(
    const std::map<std::string, nlohmann::json>& params_and_results,
    const nlohmann::json& options) {

    SimulationCoreParams params;

    // Extract basic parameters
    params.m = params_and_results.at("m");
    params.k = params_and_results.at("k");
    params.l = params_and_results.at("l");
    params.total_ssds = params_and_results.at("total_ssds");
    params.cached_m = params_and_results.at("cached_m");
    params.cached_k = params_and_results.at("cached_k");
    params.cached_l = params_and_results.at("cached_l");
    params.cached_network_m = params_and_results.at("cached_network_m");
    params.cached_network_k = params_and_results.at("cached_network_k");
    params.cached_network_l = params_and_results.at("cached_network_l");
    params.cached_ssds = params_and_results.at("cached_ssds");
    params.dwpd = params_and_results.at("dwpd");
    params.cached_dwpd_limit = params_and_results.at("cached_dwpd_limit");
    params.dwpd_limit = params_and_results.at("dwpd_limit");
    params.cached_write_ratio = params_and_results.at("cached_write_ratio");
    params.guaranteed_years = params_and_results.at("guaranteed_years");
    params.read_bw = params_and_results.at("ssd_read_bw");
    params.write_bw = params_and_results.at("ssd_write_bw");
    params.cached_read_bw = params_and_results.at("cached_ssd_read_bw");
    params.cached_write_bw = params_and_results.at("cached_ssd_write_bw");
    params.capacity = params_and_results.at("capacity");
    params.network_m = params_and_results.at("network_m");
    params.network_k = params_and_results.at("network_k");
    params.network_l = params_and_results.at("network_l");
    params.inter_replicas = params_and_results.at("inter_replicas");
    params.intra_replicas = params_and_results.at("intra_replicas");
    params.box_mttf = params_and_results.at("box_mttf");
    params.io_module_mttr = params_and_results.at("io_module_mttr");
    params.qlc = params_and_results.at("qlc");
    params.qlc_cache = params_and_results.at("qlc_cache");
    params.rebuild_bw_ratio = params_and_results.at("rebuild_bw_ratio");
    params.target_perf_ratio = params_and_results.at("target_perf_ratio");
    params.active_active = params_and_results.at("active_active");
    params.single_port_ssd = params_and_results.at("single_port_ssd");

    // Adjust cached DWPD limit first (before MTTF calculation)
    params.cached_dwpd_limit = params.cached_dwpd_limit * Utils::get_waf_from_op(0.07);
    params.cached_dwpd_limit /= 2.0;

    // Calculate MTTF (using adjusted cached_dwpd_limit)
    params.mttf = params.guaranteed_years * 365.0 * 24.0 * params.dwpd_limit / params.dwpd;
    params.cached_mttf = params.guaranteed_years * 365.0 * 24.0 * params.cached_dwpd_limit / params.dwpd;

    // Recalculate MTTF for tiered storage
    if (params.cached_m > 0) {
        double total_tbwpd = params.capacity * (params.total_ssds - params.cached_ssds) * params.dwpd / 1e12;
        if (params.cached_write_ratio > 0) {
            double cached_tbwpd = total_tbwpd / params.cached_ssds;
            double uncached_tbwpd = total_tbwpd * (1.0 - params.cached_write_ratio) / (params.total_ssds - params.cached_ssds);
            params.cached_mttf = params.guaranteed_years * 365.0 * 24.0 *
                                (params.cached_dwpd_limit * (params.capacity / 2.0) / 1e12) / cached_tbwpd;
            params.mttf = params.guaranteed_years * 365.0 * 24.0 *
                         (params.dwpd_limit * params.capacity / 1e12) / uncached_tbwpd;
        }
    }

    // Calculate SSD costs
    double tlc_cost_per_gb = options.value("tlc_cost_per_gb", 0.1);
    double qlc_cost_per_gb = options.value("qlc_cost_per_gb", 0.05);
    params.cached_ssd_cost = (params.qlc_cache ? qlc_cost_per_gb : tlc_cost_per_gb) / 1e9 * params.capacity;
    params.uncached_ssd_cost = (params.qlc ? qlc_cost_per_gb : tlc_cost_per_gb) / 1e9 * params.capacity;

    return params;
}

// Structure to hold simulation state
struct SimulationState {
    std::map<int, SSDInfo> SSDs;
    std::map<std::string, bool> failed_nodes_and_enclosures;
    std::map<int, FailureInfo> failure_info_per_ssd_group;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> failed_events;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> repair_events;
    std::map<NodeFailureKey, GraphStructure> failed_hardware_graph_table;
    std::map<NodeFailureKey, DisconnectedStatus> disconnected_table;
    std::map<FailureStateKey, FlowsAndSpeedEntry> flows_and_speed_table;
    std::unordered_map<std::string, int> node_index_map;
    int node_count = 0;
    int total_group_count = 0;

    double up_time = 0.0;
    double cached_up_time = 0.0;
    double effective_up_time = 0.0;
    double credit_up_time = 0.0;
    double credit2_up_time = 0.0;
    double timestamp = 0.0;
    double prev_time = 0.0;

    double total_mttdl = 0.0;
    int total_mttdl_count = 0;
    double prev_up_time = 0.0;
    int prev_up_first_group = 1;

    double total_cached_ssd_repair_cost = 0.0;
    double total_uncached_ssd_repair_cost = 0.0;
    double total_time_for_rebuilding_ssd0 = 0.0;
    int count_for_rebuilding_ssd0 = 0;
    double last_timestamp_for_rebuilding_ssd0 = 0.0;

    double total_cost = 0.0;
    double initial_cost = 0.0;
    double cached_initial_cost = 0.0;
    double uncached_initial_cost = 0.0;

    std::map<std::pair<double, double>, double> effective_availabilities;  // (eff_avail, avail_ratio) -> time
    DisconnectedStatus last_disconnected;

    double prev_event_time = 0.0;
    double prev_credit_uptime = 0.0;
    int simulation_years_idx = 0;
    double penalty_ratio = 0.0;
};

// Initialize simulation state
SimulationState initialize_simulation_state(
    const SimulationCoreParams& params,
    const GraphStructure& hardware_graph,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const std::map<std::string, std::string>& node_to_module_map,
    const std::map<std::string, double>& costs,
    const nlohmann::json& options,
    const EmpiricalCDF* ssd_failure_cdf = nullptr,
    double ssd_arr_qlc = 0.0,
    double ssd_arr_tlc = 0.0) {

    SimulationState state;

    // Initialize SSDs
    for (int index = 0; index < params.total_ssds; ++index) {
        state.SSDs[index] = SSDInfo();
    }

    // Initialize failed nodes and build node index map
    auto nodes = hardware_graph.get_nodes();
    state.node_count = static_cast<int>(nodes.size());
    for (int index = 0; index < state.node_count; ++index) {
        const std::string& node = nodes[index];
        state.failed_nodes_and_enclosures[node] = false;
        state.node_index_map[node] = index;
    }

    // Initialize failure info per group
    state.total_group_count = ssd_redun_scheme.get_total_group_count();
    for (int group_index = 0; group_index < state.total_group_count; ++group_index) {
        state.failure_info_per_ssd_group[group_index] = FailureInfo();
    }

    // Generate initial failure events
    state.failed_events = generate_first_failure_events(
        hardware_graph, node_to_module_map, params.total_ssds, ssd_redun_scheme, ssd_failure_cdf, ssd_arr_qlc, ssd_arr_tlc);

    // Calculate initial cost
    auto [initial_cost, cached_initial_cost, uncached_initial_cost] = get_initial_cost(
        hardware_graph, node_to_module_map, params.total_ssds, ssd_redun_scheme,
        params.cached_ssd_cost, params.uncached_ssd_cost, costs, options);

    state.initial_cost = initial_cost;
    state.cached_initial_cost = cached_initial_cost;
    state.uncached_initial_cost = uncached_initial_cost;
    state.total_cost = initial_cost;

    return state;
}

// Load network availability from file
std::tuple<double, double> load_network_availability(const nlohmann::json& options, int network_k) {
    double availability_without_network_parity = 0.0;
    double availability_without_network_parity_for_cached_ssds = 0.0;

    std::string avail_file = options.value("avail_file", "last_avail.txt");
    std::ifstream avail_in(avail_file);
    if (avail_in.is_open()) {
        if (!(avail_in >> availability_without_network_parity)) {
            availability_without_network_parity = 0.0;
        }

        if (!(avail_in >> availability_without_network_parity_for_cached_ssds)) {
            availability_without_network_parity_for_cached_ssds = 0.0;
        }
    } else if (network_k > 0) {
        Logger::getInstance().error("Failed to open availability file '" + avail_file +
                                    "' required for network parity simulation");
        throw std::runtime_error("Missing availability file for network parity simulation");
    }

    return {availability_without_network_parity, availability_without_network_parity_for_cached_ssds};
}

// Process a single simulation event
void process_simulation_event(
    const Event& event,
    SimulationState& state,
    const SimulationCoreParams& params,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const GraphStructure& hardware_graph,
    const std::map<std::string, std::string>& node_to_module_map,
    const std::map<std::string, std::vector<std::string>>& enclosure_to_node_map,
    const NetworkAvailabilityTable& network_availability_table,
    double max_read_performance_without_any_failure,
    const nlohmann::json& options,
    const std::map<std::string, double>& costs,
    const SSDIOModuleManager* ssd_io_manager,
    const EmpiricalCDF* ssd_failure_cdf = nullptr,
    double ssd_arr_qlc = 0.0,
    double ssd_arr_tlc = 0.0) {

    std::string event_type = event.event_type;
    std::string event_node = event.event_node;
    double event_time = event.time;
    double repair_start_time = event.start_time;
    bool software_repair_only = event.software_repair_only;

    // Update failure info
    update_failure_info(event_type, event_node, state.failure_info_per_ssd_group,
                       state.failed_nodes_and_enclosures, ssd_redun_scheme);

    bool network_state_changed = false;
    if (event_type == "network") {
        network_state_changed = update_network_state(
            state.failure_info_per_ssd_group, ssd_redun_scheme, network_availability_table);
        Event network_event{event_time + 1000.0, "network", "network", event_time, event_time, false};
        state.failed_events.push(network_event);
    }

    // Recalculate hardware and flow state
    NodeFailureKey node_key = build_node_failure_key(
        state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);

    calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, node_key,
                             enclosure_to_node_map, options, state.failed_hardware_graph_table,
                             state.disconnected_table, ssd_io_manager, params.total_ssds);
    DisconnectedStatus disconnected = state.disconnected_table[node_key];

    FailureStateKey flow_key = build_failure_state_key(
        state.failed_nodes_and_enclosures, state.failure_info_per_ssd_group,
        state.node_index_map, state.node_count, state.total_group_count);

    calculate_flows_and_speed(state.failed_hardware_graph_table[node_key], state.failure_info_per_ssd_group,
                              ssd_redun_scheme, options, state.flows_and_speed_table,
                              max_read_performance_without_any_failure, disconnected, flow_key, ssd_io_manager);
    const FlowsAndSpeedEntry& flow_entry = state.flows_and_speed_table[flow_key];

    // Update SSD state
    if (is_event_node_ssd(event_node)) {
        double prep_time = options.value("prep_time_for_rebuilding", 1.0);
        update_ssd_state(event_node, state.failure_info_per_ssd_group, state.SSDs, params.capacity,
                        event_type, prep_time, ssd_redun_scheme, disconnected);
    }

    if ((state.last_disconnected.local_module != disconnected.local_module ||
         state.last_disconnected.common_module != disconnected.common_module) || network_state_changed) {
        update_all_ssd_states(state.failure_info_per_ssd_group, state.SSDs, ssd_redun_scheme, disconnected, event_time);
        state.last_disconnected = disconnected;
    }

    // Handle failure events
    if (event_type == "fail") {
        if (is_event_node_ssd(event_node) && get_ssd_index(event_node) == 0) {
            state.last_timestamp_for_rebuilding_ssd0 = event_time;
        }
        push_repair_event(state.repair_events, event_node, event_time, node_to_module_map,
                         hardware_graph, software_repair_only, &options);
    }

    // Handle repair events
    if (event_type == "repair") {
        if (is_event_node_ssd(event_node) && get_ssd_index(event_node) == 0) {
            state.total_time_for_rebuilding_ssd0 += event_time - state.last_timestamp_for_rebuilding_ssd0;
            state.count_for_rebuilding_ssd0++;
        }

        // Calculate repair cost
        double module_cost = 0.0;
        if (event_time - repair_start_time > SOFT_ERROR_THRESHOLD_HOURS && !software_repair_only) {
            module_cost = calculate_module_cost(event_node, node_to_module_map, costs,
                                               ssd_redun_scheme, params.cached_ssd_cost,
                                               params.uncached_ssd_cost, options);
        }

        if (is_event_node_ssd(event_node)) {
            if (ssd_redun_scheme.is_ssd_index_cached(get_ssd_index(event_node))) {
                state.total_cached_ssd_repair_cost += module_cost;
            } else {
                state.total_uncached_ssd_repair_cost += module_cost;
            }
        }
        state.total_cost += module_cost;

        push_failed_event(state.failed_events, event_node, event_time, node_to_module_map,
                         hardware_graph, ssd_redun_scheme, ssd_failure_cdf, ssd_arr_qlc, ssd_arr_tlc);
    }

    // Update repair events for SSDs
    update_repair_event_for_SSDs(state.repair_events, event_time, state.SSDs, flow_entry,
                                 state.failure_info_per_ssd_group, ssd_redun_scheme);
}

// Single-core simulation function
SimulationResult simulation_per_core(
    int simulation_idx,
    const std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int batch_size,
    const nlohmann::json& base_options,
    const std::map<std::string, double>& costs) {

    // Seed the random number generator for this thread
    uint64_t seed = base_options.contains("seed") ? base_options["seed"].get<uint64_t>() : 0;
    seed_rng(seed, simulation_idx);

    // Extract parameters
    SimulationCoreParams params = extract_simulation_params(params_and_results, base_options);

    // Merge CLI/runtime overrides into options
    nlohmann::json options = base_options;
    if (params.rebuild_bw_ratio > 0.0) {
        options["rebuild_bw_ratio"] = params.rebuild_bw_ratio;
    }
    options["active_active"] = params.active_active;
    options["single_port_ssd"] = params.single_port_ssd;
    options["target_perf_ratio"] = params.target_perf_ratio;

    // Create SSD redundancy scheme
    SSDRedundancyScheme ssd_redun_scheme(
        params.write_bw, params.read_bw, params.mttf,
        params.cached_write_ratio, params.cached_write_bw, params.cached_read_bw, params.cached_mttf,
        params.m, params.k, params.l,
        params.cached_m, params.cached_k, params.cached_l,
        params.network_m, params.network_k, params.network_l,
        params.cached_network_m, params.cached_network_k, params.cached_network_l,
        params.cached_ssds, params.total_ssds,
        params.inter_replicas, params.intra_replicas
    );

    // Load empirical CDF for SSD failure times if trace file is provided
    EmpiricalCDF ssd_failure_cdf;
    const EmpiricalCDF* ssd_failure_cdf_ptr = nullptr;
    std::string trace_file = options.value("ssd_failure_trace", "");
    if (!trace_file.empty()) {
        if (ssd_failure_cdf.load_from_file(trace_file)) {
            ssd_failure_cdf_ptr = &ssd_failure_cdf;
            if (simulation_idx == 0) {
                Logger::getInstance().info("Using empirical SSD failure distribution from trace");
            }
        } else {
            Logger::getInstance().warning("Failed to load trace file, using exponential distribution");
        }
    }

    // Load SSD ARR (Annual Replacement Rate) for spliced distribution
    // Support separate ARR for QLC (uncached) and TLC (cached) SSDs
    double ssd_arr_qlc = options.value("ssd_arr", options.value("ssd_arr_qlc", 0.0));
    double ssd_arr_tlc = options.value("ssd_arr_tlc", ssd_arr_qlc);  // Default to QLC ARR if not specified
    if (ssd_failure_cdf_ptr != nullptr && simulation_idx == 0) {
        if (ssd_arr_qlc > 0 || ssd_arr_tlc > 0) {
            Logger::getInstance().info("Using spliced distribution with ARR: QLC=" +
                std::to_string(ssd_arr_qlc) + ", TLC=" + std::to_string(ssd_arr_tlc));
        }
    }

    // Initialize simulation
    GraphStructure hardware_graph = graph_structure_origin;
    auto [node_to_module_map, enclosure_to_node_map, max_read_performance_without_any_failure] =
        initialize_simulation(hardware_graph, params.read_bw, params.total_ssds, options,
                            params.box_mttf, params.io_module_mttr);

    // Initialize SSD-IO Module manager
    SSDIOModuleManager ssd_io_manager;
    std::vector<SSDIOModuleMapping> ssd_io_mappings;

    // Parse ssd_io_modules from options if present
    if (options.contains("ssd_io_modules") && options["ssd_io_modules"].is_array()) {
        for (const auto& entry : options["ssd_io_modules"]) {
            SSDIOModuleMapping mapping;
            mapping.ssd_count = entry.value("ssds", 0);
            mapping.io_modules = entry.value("io_modules", std::vector<std::string>{});
            mapping.start_ssd_index = 0;  // Will be computed by initialize()
            ssd_io_mappings.push_back(mapping);
        }
        if (simulation_idx == 0 && !ssd_io_mappings.empty()) {
            Logger::getInstance().info("Using SSD-IO module mapping from config (" +
                                       std::to_string(ssd_io_mappings.size()) + " mapping groups)");
        }
    }

    // Collect all io_modules from hardware graph
    std::set<std::string> all_io_modules;
    for (const auto& node : hardware_graph.get_nodes()) {
        if (node.find("io_module") != std::string::npos) {
            all_io_modules.insert(node);
        }
    }

    ssd_io_manager.initialize(ssd_io_mappings, params.total_ssds, all_io_modules);

    // Initialize state
    SimulationState state = initialize_simulation_state(params, hardware_graph, ssd_redun_scheme,
                                                        node_to_module_map, costs, options, ssd_failure_cdf_ptr, ssd_arr_qlc, ssd_arr_tlc);

    // Load network availability
    auto [avail_no_parity, avail_no_parity_cached] = load_network_availability(options, params.network_k);
    NetworkAvailabilityTable network_availability_table;
    generate_network_failure_table(
        params.cached_network_m + params.cached_network_k,
        params.network_m + params.network_k,
        avail_no_parity, avail_no_parity_cached,
        network_availability_table
    );

    // Initial state calculation
    NodeFailureKey initial_node_key = build_node_failure_key(
        state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);
    calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, initial_node_key,
                             enclosure_to_node_map, options, state.failed_hardware_graph_table,
                             state.disconnected_table, &ssd_io_manager, params.total_ssds);

    FailureStateKey flow_key = build_failure_state_key(
        state.failed_nodes_and_enclosures, state.failure_info_per_ssd_group,
        state.node_index_map, state.node_count, state.total_group_count);

    calculate_flows_and_speed(state.failed_hardware_graph_table[initial_node_key], state.failure_info_per_ssd_group,
                              ssd_redun_scheme, options, state.flows_and_speed_table,
                              max_read_performance_without_any_failure, state.disconnected_table[initial_node_key], flow_key, &ssd_io_manager);

    // Simulation parameters
    int simulation_years = options.value("simulation_years", 10);
    double simulation_hours = simulation_years * 365.0 * 24.0;
    int percentage = 1;

    Logger::getInstance().info("Simulation " + std::to_string(simulation_idx) + " started");

    if (simulation_idx == 0) {
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
        double availability_ratio = flows_and_speed_entry.availability_ratio.availability;
        double cached_availability_ratio = flows_and_speed_entry.availability_ratio.cached_availability;
        double effective_availability_ratio = flows_and_speed_entry.eff_availability_ratio;
        double credit_availability_ratio = flows_and_speed_entry.credit_availability_ratio;
        double credit2_availability_ratio = flows_and_speed_entry.credit2_availability_ratio;
        int up_first_group = flows_and_speed_entry.up_first_group;

        // Track MTTDL
        if (state.prev_up_first_group == 1 && up_first_group == 0) {
            state.total_mttdl += (state.prev_time - state.prev_up_time);
            state.total_mttdl_count++;
            state.prev_up_first_group = up_first_group;
        } else if (state.prev_up_first_group == 0 && up_first_group == 1) {
            state.prev_up_time = state.prev_time;
            state.prev_up_first_group = up_first_group;
        }

        // Update availability times
        state.effective_up_time += time_diff * effective_availability_ratio;
        state.up_time += time_diff * availability_ratio;
        state.cached_up_time += time_diff * cached_availability_ratio;
        state.credit_up_time += time_diff * credit_availability_ratio;
        state.credit2_up_time += time_diff * credit2_availability_ratio;
        state.timestamp += time_diff;
        state.effective_availabilities[{effective_availability_ratio, availability_ratio * cached_availability_ratio}] += time_diff;

        if (break_flag) break;

        // Process event
        process_simulation_event(event, state, params, ssd_redun_scheme, hardware_graph,
                                node_to_module_map, enclosure_to_node_map, network_availability_table,
                                max_read_performance_without_any_failure, options, costs,
                                &ssd_io_manager, ssd_failure_cdf_ptr, ssd_arr_qlc, ssd_arr_tlc);

        state.prev_time = event_time;

        // Update flow key after event processing
        NodeFailureKey node_key = build_node_failure_key(
            state.failed_nodes_and_enclosures, state.node_index_map, state.node_count);
        calculate_hardware_graph(hardware_graph, state.failed_nodes_and_enclosures, node_key,
                                 enclosure_to_node_map, options, state.failed_hardware_graph_table,
                                 state.disconnected_table, &ssd_io_manager, params.total_ssds);

        flow_key = build_failure_state_key(
            state.failed_nodes_and_enclosures, state.failure_info_per_ssd_group,
            state.node_index_map, state.node_count, state.total_group_count);

        calculate_flows_and_speed(state.failed_hardware_graph_table[node_key], state.failure_info_per_ssd_group,
                                  ssd_redun_scheme, options, state.flows_and_speed_table,
                                  max_read_performance_without_any_failure, state.disconnected_table[node_key], flow_key);

        // Progress bar update
        if (simulation_idx == 0) {
            int new_percentage = static_cast<int>(event_time * 100.0 / (simulation_hours * batch_size));
            if (new_percentage > percentage) {
                percentage = new_percentage;
                Utils::progress_bar(percentage, 100, 50);
            }
        }

        // Calculate penalty ratio per year
        if (event_time > simulation_hours * (state.simulation_years_idx + 1)) {
            double year_duration = event_time - state.prev_event_time;
            double year_credit_uptime = state.credit_up_time - state.prev_credit_uptime;
            double year_avail = year_credit_uptime / year_duration;
            state.penalty_ratio += (year_duration / simulation_hours) * nine_to_credit(Utils::get_nines(year_avail));
            state.simulation_years_idx++;
            state.prev_credit_uptime = state.credit_up_time;
            state.prev_event_time = event_time;
        }
    }

    // Calculate final metrics
    state.penalty_ratio /= (state.simulation_years_idx > 0 ? state.simulation_years_idx : 1);

    // Return result
    SimulationResult result;
    result.up_time = state.up_time;
    result.cached_up_time = state.cached_up_time;
    result.credit_up_time = state.credit_up_time;
    result.credit2_up_time = state.credit2_up_time;
    result.simulation_time = state.timestamp;
    result.effective_up_time = state.effective_up_time;
    result.effective_availabilities = state.effective_availabilities;
    result.initial_cost = state.initial_cost;
    result.total_cost = state.total_cost;
    result.time_for_rebuilding = state.total_time_for_rebuilding_ssd0;
    result.count_for_rebuilding = state.count_for_rebuilding_ssd0;
    result.cached_mttf = params.cached_mttf;
    result.mttf = params.mttf;
    result.cached_ssd_repair_cost = state.total_cached_ssd_repair_cost;
    result.uncached_ssd_repair_cost = state.total_uncached_ssd_repair_cost;
    result.cached_initial_cost = state.cached_initial_cost;
    result.uncached_initial_cost = state.uncached_initial_cost;
    result.penalty_ratio = state.penalty_ratio;
    result.total_mttdl = state.total_mttdl;
    result.total_mttdl_count = state.total_mttdl_count;

    return result;
}
