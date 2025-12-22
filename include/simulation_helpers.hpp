#pragma once

#include "simulation.hpp"
#include "graph_structure.hpp"
#include <map>
#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

// Helper function declarations from simulation_helpers.cpp

// Build keys from current state
NodeFailureKey build_node_failure_key(
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const std::unordered_map<std::string, int>& node_index_map,
    int node_count);

FailureStateKey build_failure_state_key(
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const std::unordered_map<std::string, int>& node_index_map,
    int node_count,
    int total_group_count);

// Calculate hardware graph with failures
void calculate_hardware_graph(
    const GraphStructure& hardware_graph,
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const NodeFailureKey& key,
    const std::map<std::string, std::vector<std::string>>& enclosure_to_node_map,
    const nlohmann::json& options,
    std::map<NodeFailureKey, GraphStructure>& failed_hardware_graph_table,
    std::map<NodeFailureKey, DisconnectedStatus>& disconnected_table,
    const SSDIOModuleManager* ssd_io_manager = nullptr,
    int total_ssds = 0);

// Initialize simulation
std::tuple<std::map<std::string, std::string>,
           std::map<std::string, std::vector<std::string>>,
           double>
initialize_simulation(
    GraphStructure& hardware_graph,
    double ssd_read_bw,
    int total_ssd_count,
    const nlohmann::json& options,
    double box_mttf,
    double io_module_mttr);

// Simulation result structure
struct SimulationResult {
    double up_time;
    double cached_up_time;
    double credit_up_time;
    double credit2_up_time;
    double simulation_time;
    double effective_up_time;
    std::map<std::pair<double, double>, double> effective_availabilities;  // (eff_avail, avail_ratio) -> time
    double initial_cost;
    double total_cost;
    double time_for_rebuilding;
    int count_for_rebuilding;
    double cached_mttf;
    double mttf;
    double cached_ssd_repair_cost;
    double uncached_ssd_repair_cost;
    double cached_initial_cost;
    double uncached_initial_cost;
    double penalty_ratio;
    double total_mttdl;
    int total_mttdl_count;
};

// Single-core simulation function (from simulation_core.cpp)
SimulationResult simulation_per_core(
    int simulation_idx,
    const std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int batch_size,
    const nlohmann::json& options,
    const std::map<std::string, double>& costs);
