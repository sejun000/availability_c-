#pragma once

#include "simulation.hpp"
#include "graph_structure.hpp"
#include "disk.hpp"
#include "erasure_coding.hpp"
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
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
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
    std::unordered_map<NodeFailureKey, GraphStructure, NodeFailureKeyHash>& failed_hardware_graph_table,
    std::unordered_map<NodeFailureKey, DisconnectedStatus, NodeFailureKeyHash>& disconnected_table,
    const DiskIOModuleManager* disk_io_manager = nullptr,
    int total_disks = 0);

// Initialize simulation
std::tuple<std::map<std::string, std::string>,
           std::map<std::string, std::vector<std::string>>,
           double>
initialize_simulation(
    GraphStructure& hardware_graph,
    double disk_read_bw,
    int total_disk_count,
    const nlohmann::json& options);

// Single-core simulation function (from simulation_core.cpp)
SimulationResult simulation_per_core(
    int simulation_idx,
    const std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int batch_size,
    const nlohmann::json& options);

// Data loss calculation for EC schemes
double calculate_data_loss_ratio(
    const ECGroupFailureInfo& failure_info,
    const ECConfig& ec_config);

// Calculate durability from simulation results
double calculate_durability(
    double total_data_loss_ratio,
    double simulation_time_hours,
    int simulation_years);
