#include "simulation_helpers.hpp"
#include "logger.hpp"
#include <algorithm>
#include <iostream>
#include <set>

// Build node failure key from current state
NodeFailureKey build_node_failure_key(
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const std::unordered_map<std::string, int>& node_index_map,
    int node_count) {

    NodeFailureKey key;
    key.failed_flags.resize(node_count, 0);

    for (const auto& [node, failed] : failed_nodes_and_enclosures) {
        if (failed) {
            auto it = node_index_map.find(node);
            if (it != node_index_map.end()) {
                key.failed_flags[it->second] = 1;
            }
        }
    }

    return key;
}

// Build failure state key from current state
FailureStateKey build_failure_state_key(
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const std::map<int, DiskInfo>& disks,
    const std::unordered_map<std::string, int>& node_index_map,
    int node_count,
    int total_disks) {

    FailureStateKey key;
    key.node_key = build_node_failure_key(failed_nodes_and_enclosures, node_index_map, node_count);
    key.disk_state.resize(std::max(0, total_disks), 0);

    for (int disk_index = 0; disk_index < total_disks; ++disk_index) {
        auto it = disks.find(disk_index);
        if (it == disks.end()) {
            continue;
        }
        const DiskInfo& disk = it->second;

        if (disk.is_failed || disk.is_disconnected) {
            key.disk_state[disk_index] = 1;
        } else if (disk.needs_rebuild()) {
            key.disk_state[disk_index] = 2;
        } else {
            key.disk_state[disk_index] = 0;
        }
    }

    return key;
}

// Helper function to recursively disable nodes in nested enclosures
void disable_enclosure_recursively(
    const std::string& enclosure_name,
    const std::map<std::string, std::vector<std::string>>& enclosure_to_node_map,
    GraphStructure& graph_copy,
    std::set<std::string>& visited) {

    // Prevent infinite recursion
    if (visited.count(enclosure_name)) return;
    visited.insert(enclosure_name);

    auto encl_it = enclosure_to_node_map.find(enclosure_name);
    if (encl_it == enclosure_to_node_map.end()) return;

    for (const auto& enclosed_node : encl_it->second) {
        graph_copy.disable_node(enclosed_node);

        // Check if this node is also an enclosure (nested enclosure)
        if (enclosure_to_node_map.find(enclosed_node) != enclosure_to_node_map.end()) {
            disable_enclosure_recursively(enclosed_node, enclosure_to_node_map, graph_copy, visited);
        }
    }
}

// Calculate hardware graph with failures applied
void calculate_hardware_graph(
    const GraphStructure& hardware_graph,
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const NodeFailureKey& key,
    const std::map<std::string, std::vector<std::string>>& enclosure_to_node_map,
    const nlohmann::json& options,
    std::unordered_map<NodeFailureKey, GraphStructure, NodeFailureKeyHash>& failed_hardware_graph_table,
    std::unordered_map<NodeFailureKey, DisconnectedStatus, NodeFailureKeyHash>& disconnected_table,
    const DiskIOModuleManager* disk_io_manager,
    int total_disks) {

    // Check if already computed
    if (failed_hardware_graph_table.find(key) != failed_hardware_graph_table.end()) {
        return;
    }

    // Clone the hardware graph
    GraphStructure graph_copy = hardware_graph.clone();
    DisconnectedStatus disconnected;
    disconnected.failed_nodes = failed_nodes_and_enclosures;

    // Disable failed nodes
    for (const auto& [node, failed] : failed_nodes_and_enclosures) {
        if (failed) {
            graph_copy.disable_node(node);

            // Also disable nodes in enclosures (recursively for nested enclosures)
            std::set<std::string> visited;
            disable_enclosure_recursively(node, enclosure_to_node_map, graph_copy, visited);
        }
    }

    // Update disk disconnection status
    if (disk_io_manager != nullptr) {
        for (int disk_idx = 0; disk_idx < total_disks; ++disk_idx) {
            bool is_disconnected = disk_io_manager->is_disk_disconnected(disk_idx, failed_nodes_and_enclosures);
            disconnected.set_disk_disconnected(disk_idx, is_disconnected);
        }
    }

    failed_hardware_graph_table[key] = graph_copy;
    disconnected_table[key] = disconnected;
}

// Initialize simulation
std::tuple<std::map<std::string, std::string>,
           std::map<std::string, std::vector<std::string>>,
           double,
           double>
initialize_simulation(
    GraphStructure& hardware_graph,
    double disk_read_bw,
    int total_disk_count,
    const nlohmann::json& options) {

    std::map<std::string, std::string> node_to_module_map;
    std::map<std::string, std::vector<std::string>> enclosure_to_node_map;

    std::string start_module = options["start_module"];
    std::string end_module = options["end_module"];

    // Build node to module map
    for (const auto& node : hardware_graph.nodes_) {
        // Extract module name (prefix before underscore or number)
        size_t pos = node.find_first_of("0123456789");
        if (pos != std::string::npos && pos > 0) {
            // Remove trailing underscore if present
            std::string module = node.substr(0, pos);
            if (!module.empty() && module.back() == '_') {
                module.pop_back();
            }
            node_to_module_map[node] = module;
        } else {
            node_to_module_map[node] = node;
        }
    }

    // Build enclosure to node map
    for (const auto& [enclosure, nodes] : hardware_graph.enclosures_) {
        enclosure_to_node_map[enclosure] = nodes;

        // Add enclosure to node_to_module_map
        size_t pos = enclosure.find_first_of("0123456789");
        if (pos != std::string::npos && pos > 0) {
            std::string module = enclosure.substr(0, pos);
            if (!module.empty() && module.back() == '_') {
                module.pop_back();
            }
            node_to_module_map[enclosure] = module;
        } else {
            node_to_module_map[enclosure] = enclosure;
        }
    }

    // Build parent map for LCA calculation
    hardware_graph.build_parent_map(start_module);

    // Calculate baseline host bandwidths without any failure.
    GraphStructure hardware_graph_copy = hardware_graph.clone();
    double max_read = hardware_graph_copy.calculate_max_flow(end_module, start_module, FlowDirection::UPSTREAM);
    double max_write = hardware_graph_copy.calculate_max_flow(start_module, end_module, FlowDirection::DOWNSTREAM);

    Logger::getInstance().info("Max host READ bw without failure: " + std::to_string(max_read / 1e9) + " GB/s");
    Logger::getInstance().info("Max host WRITE bw without failure: " + std::to_string(max_write / 1e9) + " GB/s");

    (void)disk_read_bw;
    (void)total_disk_count;

    return {node_to_module_map, enclosure_to_node_map, max_read, max_write};
}

// Calculate data loss ratio for EC schemes
double calculate_data_loss_ratio(
    const ECGroupFailureInfo& failure_info,
    const ECConfig& ec_config) {

    int unavailable = failure_info.get_unavailable_count();

    if (unavailable <= ec_config.k) {
        return 0.0;  // Can recover
    }

    // Data loss: fraction of group's data that is lost
    // In declustered parity, each disk holds 1/n of data
    // With k+x failures (x > 0), we lose approximately x/n of data
    int excess_failures = unavailable - ec_config.k;
    int group_size = ec_config.n > 0 ? ec_config.n : (ec_config.m + ec_config.k);

    return static_cast<double>(excess_failures) / group_size;
}

// Calculate durability from simulation results
double calculate_durability(
    double total_data_loss_ratio,
    double simulation_time_hours,
    int simulation_years) {

    if (simulation_time_hours <= 0 || simulation_years <= 0) {
        return 1.0;
    }

    double hours_per_year = 365.0 * 24.0;
    double expected_years = simulation_time_hours / hours_per_year;

    if (expected_years <= 0) {
        return 1.0;
    }

    // Durability = 1 - (data loss ratio per year)
    double data_loss_per_year = total_data_loss_ratio / expected_years;

    // Clamp to valid range
    if (data_loss_per_year >= 1.0) {
        return 0.0;
    }

    return 1.0 - data_loss_per_year;
}
