#include "simulation.hpp"
#include "utils.hpp"
#include <algorithm>
#include <unordered_map>

NodeFailureKey build_node_failure_key(
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const std::unordered_map<std::string, int>& node_index_map,
    int node_count) {

    NodeFailureKey key;
    key.failed_flags.assign(static_cast<size_t>(node_count), 0);

    for (const auto& [node, failed] : failed_nodes_and_enclosures) {
        auto it = node_index_map.find(node);
        if (it == node_index_map.end()) continue;
        key.failed_flags[static_cast<size_t>(it->second)] = failed ? 1 : 0;
    }

    return key;
}

FailureStateKey build_failure_state_key(
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const std::unordered_map<std::string, int>& node_index_map,
    int node_count,
    int total_group_count) {

    FailureStateKey key;
    key.node_key = build_node_failure_key(failed_nodes_and_enclosures, node_index_map, node_count);
    key.failure_counts.assign(static_cast<size_t>(total_group_count), 0);
    key.network_failure_counts.assign(static_cast<size_t>(total_group_count), 0);

    for (const auto& [group_index, failure_info] : failure_info_per_ssd_group) {
        if (group_index < 0 || group_index >= total_group_count) continue;
        key.failure_counts[static_cast<size_t>(group_index)] =
            static_cast<uint8_t>(std::clamp(failure_info.failure_count, 0, 255));
        key.network_failure_counts[static_cast<size_t>(group_index)] =
            static_cast<uint8_t>(std::clamp(failure_info.network_failure_count, 0, 255));
    }

    return key;
}

void calculate_hardware_graph(
    const GraphStructure& hardware_graph,
    const std::map<std::string, bool>& failed_nodes_and_enclosures,
    const NodeFailureKey& key,
    const std::map<std::string, std::vector<std::string>>& enclosure_to_node_map,
    const nlohmann::json& options,
    std::map<NodeFailureKey, GraphStructure>& failed_hardware_graph_table,
    std::map<NodeFailureKey, DisconnectedStatus>& disconnected_table) {

    if (failed_hardware_graph_table.find(key) != failed_hardware_graph_table.end()) {
        return;
    }

    GraphStructure hardware_graph_copy = hardware_graph;

    for (const auto& [failed_node_or_enclosure, is_failed] : failed_nodes_and_enclosures) {
        if (!is_failed) continue;

        if (hardware_graph_copy.has_node(failed_node_or_enclosure)) {
            hardware_graph_copy.remove_node(failed_node_or_enclosure);
        } else {
            auto it = enclosure_to_node_map.find(failed_node_or_enclosure);
            if (it != enclosure_to_node_map.end()) {
                for (const std::string& node : it->second) {
                    if (hardware_graph_copy.has_node(node)) {
                        hardware_graph_copy.remove_node(node);
                    }
                }
            }
        }
    }

    std::string start_module = options["start_module"];
    std::string end_module = options["end_module"];

    hardware_graph_copy.add_virtual_nodes(start_module, end_module);
    bool connected = hardware_graph_copy.has_path("virtual_source", "virtual_sink");
    hardware_graph_copy.remove_virtual_nodes();

    DisconnectedStatus disconnected;
    disconnected.local_module = !connected;
    disconnected.common_module = false;

    disconnected_table[key] = disconnected;
    failed_hardware_graph_table[key] = hardware_graph_copy;
}

std::tuple<std::map<std::string, std::string>,
           std::map<std::string, std::vector<std::string>>,
           double>
initialize_simulation(
    GraphStructure& hardware_graph,
    double ssd_read_bw,
    int total_ssd_count,
    const nlohmann::json& options,
    double box_mttf,
    double io_module_mttr) {

    std::map<std::string, std::string> node_to_module_map;
    std::map<std::string, std::vector<std::string>> enclosure_to_node_map;

    for (const auto& [enclosure, nodes] : hardware_graph.enclosures_) {
        for (const auto& [module, mttf] : hardware_graph.mttfs_) {
            if (enclosure.find(module) != std::string::npos) {
                node_to_module_map[enclosure] = module;
                if (box_mttf > 0) {
                    hardware_graph.mttfs_[module] = box_mttf;
                }
                break;
            }
        }
        enclosure_to_node_map[enclosure] = nodes;
    }

    for (const std::string& node : hardware_graph.get_nodes()) {
        for (const auto& [module, mttf] : hardware_graph.mttfs_) {
            if (node.find(module) != std::string::npos) {
                if (node.find("io_module") != std::string::npos ||
                    node.find("backend_module") != std::string::npos ||
                    node.find("host_module") != std::string::npos) {
                    double origin_mtr = hardware_graph.mtrs_[module];
                    if (io_module_mttr > origin_mtr) {
                        hardware_graph.mtrs_[module] = io_module_mttr;
                    }
                }
                node_to_module_map[node] = module;
                break;
            }
        }
    }

    GraphStructure hardware_graph_copy = hardware_graph;
    std::string start_module = options["start_module"];
    std::string end_module = options["end_module"];
    std::string lowest_common_module = options["lowest_common_module"];

    hardware_graph_copy.add_virtual_nodes(start_module, end_module);
    double max_flow_hardware = hardware_graph_copy.maximum_flow("virtual_source", "virtual_sink");
    hardware_graph_copy.remove_virtual_nodes();

    hardware_graph_copy.add_virtual_nodes(start_module, lowest_common_module);
    double common_module_max_flow = hardware_graph_copy.maximum_flow("virtual_source", "virtual_sink");
    hardware_graph_copy.remove_virtual_nodes();

    double ssd_max_read_performance = std::min(max_flow_hardware, total_ssd_count * ssd_read_bw);
    int network_nodes = options["network_nodes"];
    ssd_max_read_performance = std::min(ssd_max_read_performance, common_module_max_flow / network_nodes);

    double dram_bandwidth = get_dram_bandwidth(options);
    ssd_max_read_performance = std::min(ssd_max_read_performance, dram_bandwidth);

    return {node_to_module_map, enclosure_to_node_map, ssd_max_read_performance};
}
