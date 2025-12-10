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
