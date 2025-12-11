#include "graph_structure.hpp"
#include <queue>
#include <algorithm>
#include <limits>
#include <iostream>

GraphStructure::GraphStructure() {}

GraphStructure::GraphStructure(
    const std::vector<std::tuple<std::string, std::string, std::string>>& edges,
    const std::map<std::string, std::vector<std::string>>& enclosures,
    const std::map<std::string, double>& mttfs,
    const std::map<std::string, double>& mtrs
) : enclosures_(enclosures), mttfs_(mttfs), mtrs_(mtrs) {
    add_edges(edges);
}

void GraphStructure::add_edges(const std::vector<std::tuple<std::string, std::string, std::string>>& edges) {
    for (const auto& [start, end, weight] : edges) {
        double capacity = parse_weight(weight);

        Edge edge{start, end, capacity, weight};
        edges_.push_back(edge);

        // Add to adjacency list
        adjacency_list_[start][end] = capacity;
        nodes_.insert(start);
        nodes_.insert(end);
    }
}

double GraphStructure::parse_weight(const std::string& weight) const {
    if (weight.empty()) return 0.0;

    char last_char = weight.back();
    if (last_char == 'G') {
        return std::stod(weight.substr(0, weight.length() - 1)) * 1e9;
    } else if (last_char == 'M') {
        return std::stod(weight.substr(0, weight.length() - 1)) * 1e6;
    } else if (last_char == 'K') {
        return std::stod(weight.substr(0, weight.length() - 1)) * 1e3;
    }

    return std::stod(weight);
}

bool GraphStructure::get_module_degraded(const std::string& end_module, bool active_active) const {
    int cnt = 0;
    for (const auto& node : nodes_) {
        if (node.find(end_module) != std::string::npos) {
            cnt++;
        }
    }

    if (cnt == 1 && active_active) {
        return true;
    }
    return false;
}

double GraphStructure::calculate_max_flow(const std::string& start_node_module,
                                          const std::string& leaf_node_module) {
    add_virtual_nodes(start_node_module, leaf_node_module);
    double flow_value = maximum_flow("virtual_source", "virtual_sink");
    remove_virtual_nodes();
    return flow_value;
}

void GraphStructure::add_virtual_nodes(const std::string& start_node_module,
                                       const std::string& leaf_node_module) {
    nodes_.insert("virtual_source");
    nodes_.insert("virtual_sink");

    for (const auto& node : nodes_) {
        if (node.find(start_node_module) != std::string::npos && node != "virtual_source") {
            adjacency_list_["virtual_source"][node] = MAX_EDGE_VALUE;
        }
        if (node.find(leaf_node_module) != std::string::npos && node != "virtual_sink") {
            adjacency_list_[node]["virtual_sink"] = MAX_EDGE_VALUE;
        }
    }
}

void GraphStructure::add_virtual_source(const std::string& start_node_module) {
    nodes_.insert("virtual_source");

    for (const auto& node : nodes_) {
        if (node.find(start_node_module) != std::string::npos && node != "virtual_source") {
            adjacency_list_["virtual_source"][node] = MAX_EDGE_VALUE;
        }
    }
}

void GraphStructure::add_virtual_sink(const std::string& leaf_node_module,
                                      const std::set<std::string>& exception_nodes) {
    nodes_.insert("virtual_sink");

    for (const auto& node : nodes_) {
        if (node.find(leaf_node_module) != std::string::npos &&
            node != "virtual_sink" &&
            exception_nodes.find(node) == exception_nodes.end()) {
            adjacency_list_[node]["virtual_sink"] = MAX_EDGE_VALUE;
        }
    }
}

void GraphStructure::remove_virtual_source() {
    if (nodes_.find("virtual_source") != nodes_.end()) {
        nodes_.erase("virtual_source");
        adjacency_list_.erase("virtual_source");
        for (auto& [node, neighbors] : adjacency_list_) {
            neighbors.erase("virtual_source");
        }
    }
}

void GraphStructure::remove_virtual_sink() {
    if (nodes_.find("virtual_sink") != nodes_.end()) {
        nodes_.erase("virtual_sink");
        adjacency_list_.erase("virtual_sink");
        for (auto& [node, neighbors] : adjacency_list_) {
            neighbors.erase("virtual_sink");
        }
    }
}

void GraphStructure::remove_virtual_nodes() {
    remove_virtual_source();
    remove_virtual_sink();
}

void GraphStructure::remove_node(const std::string& node) {
    if (nodes_.find(node) == nodes_.end()) return;

    nodes_.erase(node);
    adjacency_list_.erase(node);

    // Remove edges pointing to this node
    for (auto& [n, neighbors] : adjacency_list_) {
        neighbors.erase(node);
    }
}

bool GraphStructure::has_node(const std::string& node) const {
    return nodes_.find(node) != nodes_.end();
}

std::vector<std::string> GraphStructure::get_nodes() const {
    return std::vector<std::string>(nodes_.begin(), nodes_.end());
}

bool GraphStructure::has_path(const std::string& source, const std::string& sink) const {
    if (nodes_.find(source) == nodes_.end() || nodes_.find(sink) == nodes_.end()) {
        return false;
    }

    std::set<std::string> visited;
    std::queue<std::string> q;

    q.push(source);
    visited.insert(source);

    while (!q.empty()) {
        std::string current = q.front();
        q.pop();

        if (current == sink) {
            return true;
        }

        auto it = adjacency_list_.find(current);
        if (it != adjacency_list_.end()) {
            for (const auto& [neighbor, capacity] : it->second) {
                if (visited.find(neighbor) == visited.end() && capacity > 0) {
                    visited.insert(neighbor);
                    q.push(neighbor);
                }
            }
        }
    }

    return false;
}

bool GraphStructure::bfs(const std::string& source, const std::string& sink,
                         std::map<std::string, std::string>& parent,
                         const std::map<std::string, std::map<std::string, double>>& residual_graph) const {
    std::set<std::string> visited;
    std::queue<std::string> q;

    q.push(source);
    visited.insert(source);
    parent[source] = "";

    while (!q.empty()) {
        std::string u = q.front();
        q.pop();

        auto it = residual_graph.find(u);
        if (it == residual_graph.end()) continue;

        for (const auto& [v, capacity] : it->second) {
            if (visited.find(v) == visited.end() && capacity > 1e-9) {
                visited.insert(v);
                parent[v] = u;

                if (v == sink) {
                    return true;
                }

                q.push(v);
            }
        }
    }

    return false;
}

double GraphStructure::maximum_flow(const std::string& source, const std::string& sink) {
    // Create residual graph
    std::map<std::string, std::map<std::string, double>> residual_graph;

    for (const auto& [u, neighbors] : adjacency_list_) {
        for (const auto& [v, capacity] : neighbors) {
            residual_graph[u][v] = capacity;
            if (residual_graph[v].find(u) == residual_graph[v].end()) {
                residual_graph[v][u] = 0.0;
            }
        }
    }

    std::map<std::string, std::string> parent;
    double max_flow = 0.0;

    // Edmonds-Karp algorithm
    while (bfs(source, sink, parent, residual_graph)) {
        // Find minimum capacity along the path
        double path_flow = std::numeric_limits<double>::max();
        std::string v = sink;

        while (v != source) {
            std::string u = parent[v];
            path_flow = std::min(path_flow, residual_graph[u][v]);
            v = u;
        }

        // Update residual capacities
        v = sink;
        while (v != source) {
            std::string u = parent[v];
            residual_graph[u][v] -= path_flow;
            residual_graph[v][u] += path_flow;
            v = u;
        }

        max_flow += path_flow;
        parent.clear();
    }

    return max_flow;
}

// Build parent map for LCA calculation using BFS from root
void GraphStructure::build_parent_map(const std::string& root_module) {
    parent_map_.clear();
    depth_map_.clear();

    // Find the actual root node(s) that match root_module
    std::vector<std::string> root_nodes;
    for (const auto& node : nodes_) {
        if (node.find(root_module) != std::string::npos) {
            root_nodes.push_back(node);
        }
    }

    if (root_nodes.empty()) {
        return;
    }

    // If multiple root nodes, create a virtual root
    if (root_nodes.size() > 1) {
        root_node_ = "virtual_root";
        parent_map_[root_node_] = "";
        depth_map_[root_node_] = 0;

        for (const auto& root : root_nodes) {
            parent_map_[root] = root_node_;
            depth_map_[root] = 1;
        }
    } else {
        root_node_ = root_nodes[0];
        parent_map_[root_node_] = "";
        depth_map_[root_node_] = 0;
    }

    // BFS to build parent relationships (top-down)
    std::queue<std::string> q;
    std::set<std::string> visited;

    for (const auto& root : root_nodes) {
        q.push(root);
        visited.insert(root);
    }

    while (!q.empty()) {
        std::string current = q.front();
        q.pop();

        int current_depth = depth_map_[current];

        // Find children (nodes that current points to)
        auto it = adjacency_list_.find(current);
        if (it != adjacency_list_.end()) {
            for (const auto& [child, capacity] : it->second) {
                if (visited.find(child) == visited.end() && capacity > 0) {
                    visited.insert(child);
                    parent_map_[child] = current;
                    depth_map_[child] = current_depth + 1;
                    q.push(child);
                }
            }
        }
    }
}

// Get path from node to root
std::vector<std::string> GraphStructure::get_path_to_root(const std::string& node) const {
    std::vector<std::string> path;

    auto it = parent_map_.find(node);
    if (it == parent_map_.end()) {
        return path;  // Node not in tree
    }

    std::string current = node;
    while (!current.empty()) {
        path.push_back(current);
        auto parent_it = parent_map_.find(current);
        if (parent_it == parent_map_.end() || parent_it->second.empty()) {
            break;
        }
        current = parent_it->second;
    }

    return path;
}

// Find LCA of two nodes
std::string GraphStructure::find_lca(const std::string& node1, const std::string& node2) const {
    std::vector<std::string> path1 = get_path_to_root(node1);
    std::vector<std::string> path2 = get_path_to_root(node2);

    if (path1.empty() || path2.empty()) {
        return "";
    }

    // Convert path2 to set for O(1) lookup
    std::set<std::string> path2_set(path2.begin(), path2.end());

    // Find first node in path1 that exists in path2
    for (const auto& node : path1) {
        if (path2_set.find(node) != path2_set.end()) {
            return node;
        }
    }

    return "";
}

// Find LCA of multiple nodes
std::string GraphStructure::find_lca_multiple(const std::vector<std::string>& nodes) const {
    if (nodes.empty()) {
        return "";
    }
    if (nodes.size() == 1) {
        return nodes[0];
    }

    std::string lca = nodes[0];
    for (size_t i = 1; i < nodes.size(); ++i) {
        lca = find_lca(lca, nodes[i]);
        if (lca.empty()) {
            return "";
        }
    }

    return lca;
}

// Get all nodes at the same depth as LCA that are ancestors of given leaf nodes
std::set<std::string> GraphStructure::get_lca_level_nodes(const std::vector<std::string>& leaf_nodes) const {
    std::set<std::string> lca_level_nodes;

    if (leaf_nodes.empty()) {
        return lca_level_nodes;
    }

    // Find LCA of all leaf nodes
    std::string lca = find_lca_multiple(leaf_nodes);
    if (lca.empty()) {
        return lca_level_nodes;
    }

    // Get depth of LCA
    auto depth_it = depth_map_.find(lca);
    if (depth_it == depth_map_.end()) {
        return lca_level_nodes;
    }
    int lca_depth = depth_it->second;

    // For each leaf node, traverse up to LCA depth and collect those nodes
    for (const auto& leaf : leaf_nodes) {
        std::string current = leaf;

        while (!current.empty()) {
            auto curr_depth_it = depth_map_.find(current);
            if (curr_depth_it == depth_map_.end()) {
                break;
            }

            if (curr_depth_it->second == lca_depth) {
                lca_level_nodes.insert(current);
                break;
            }

            auto parent_it = parent_map_.find(current);
            if (parent_it == parent_map_.end() || parent_it->second.empty()) {
                break;
            }
            current = parent_it->second;
        }
    }

    return lca_level_nodes;
}

// Calculate max flow from source SSDs to target SSD through LCA level
double GraphStructure::calculate_rebuild_max_flow(const std::vector<std::string>& source_ssds,
                                                   const std::string& target_ssd) {
    if (source_ssds.empty()) {
        return 0.0;
    }

    // Combine source SSDs and target to find their common ancestors
    std::vector<std::string> all_ssds = source_ssds;
    all_ssds.push_back(target_ssd);

    // Get LCA level nodes (these are the nodes at the common ancestor level)
    std::set<std::string> lca_level_nodes = get_lca_level_nodes(all_ssds);

    if (lca_level_nodes.empty()) {
        // Fallback: all SSDs might be under the same parent
        // In this case, use direct connection
        return MAX_EDGE_VALUE;
    }

    // Add virtual source connected to all LCA level nodes
    nodes_.insert("rebuild_virtual_source");
    for (const auto& lca_node : lca_level_nodes) {
        adjacency_list_["rebuild_virtual_source"][lca_node] = MAX_EDGE_VALUE;
    }

    // Add virtual sink connected from target SSD
    nodes_.insert("rebuild_virtual_sink");
    adjacency_list_[target_ssd]["rebuild_virtual_sink"] = MAX_EDGE_VALUE;

    // Calculate max flow
    double flow = maximum_flow("rebuild_virtual_source", "rebuild_virtual_sink");

    // Clean up virtual nodes
    nodes_.erase("rebuild_virtual_source");
    nodes_.erase("rebuild_virtual_sink");
    adjacency_list_.erase("rebuild_virtual_source");
    adjacency_list_.erase("rebuild_virtual_sink");
    for (auto& [node, neighbors] : adjacency_list_) {
        neighbors.erase("rebuild_virtual_source");
        neighbors.erase("rebuild_virtual_sink");
    }

    return flow;
}

// Calculate rebuild max flow using io_module names
// This finds LCA of source io_modules and calculates max flow to target io_module
double GraphStructure::calculate_rebuild_max_flow_via_io_modules(
    const std::vector<std::string>& source_io_modules,
    const std::string& target_io_module) {

    if (source_io_modules.empty()) {
        return 0.0;
    }

    // Combine all io_modules to find their common ancestors
    std::vector<std::string> all_io_modules = source_io_modules;
    // Check if target is already in source (same io_module case)
    bool target_in_source = false;
    for (const auto& src : source_io_modules) {
        if (src == target_io_module) {
            target_in_source = true;
            break;
        }
    }
    if (!target_in_source) {
        all_io_modules.push_back(target_io_module);
    }

    // If all io_modules are the same, no network traversal needed
    std::set<std::string> unique_modules(all_io_modules.begin(), all_io_modules.end());
    if (unique_modules.size() == 1) {
        // All within same io_module, return infinite bandwidth (limited by SSD)
        return MAX_EDGE_VALUE;
    }

    // Get LCA level nodes
    std::set<std::string> lca_level_nodes = get_lca_level_nodes(all_io_modules);

    if (lca_level_nodes.empty()) {
        // Fallback: couldn't find LCA, use direct parent relationship
        return MAX_EDGE_VALUE;
    }

    // If LCA level is the same as io_module level (all under same parent already)
    // check if LCA itself is one of the io_modules
    bool lca_is_io_module = false;
    for (const auto& lca_node : lca_level_nodes) {
        for (const auto& io_mod : all_io_modules) {
            if (lca_node == io_mod) {
                lca_is_io_module = true;
                break;
            }
        }
        if (lca_is_io_module) break;
    }

    if (lca_is_io_module) {
        // LCA is at io_module level, meaning they share the same io_module
        return MAX_EDGE_VALUE;
    }

    // Add virtual source connected to all LCA level nodes
    nodes_.insert("rebuild_virtual_source");
    for (const auto& lca_node : lca_level_nodes) {
        adjacency_list_["rebuild_virtual_source"][lca_node] = MAX_EDGE_VALUE;
    }

    // Add virtual sink connected from target io_module
    nodes_.insert("rebuild_virtual_sink");
    adjacency_list_[target_io_module]["rebuild_virtual_sink"] = MAX_EDGE_VALUE;

    // Calculate max flow
    double flow = maximum_flow("rebuild_virtual_source", "rebuild_virtual_sink");

    // Clean up virtual nodes
    nodes_.erase("rebuild_virtual_source");
    nodes_.erase("rebuild_virtual_sink");
    adjacency_list_.erase("rebuild_virtual_source");
    adjacency_list_.erase("rebuild_virtual_sink");
    for (auto& [node, neighbors] : adjacency_list_) {
        neighbors.erase("rebuild_virtual_source");
        neighbors.erase("rebuild_virtual_sink");
    }

    return flow;
}
