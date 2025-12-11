#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include "utils.hpp"

constexpr double MAX_EDGE_VALUE = 1.0e15;

// Edge structure
struct Edge {
    std::string start;
    std::string end;
    double capacity;
    std::string label;
};

// Graph structure for hardware topology and max flow calculation
class GraphStructure {
public:
    GraphStructure();
    GraphStructure(
        const std::vector<std::tuple<std::string, std::string, std::string>>& edges,
        const std::map<std::string, std::vector<std::string>>& enclosures,
        const std::map<std::string, double>& mttfs,
        const std::map<std::string, double>& mtrs
    );

    // Add edges to graph
    void add_edges(const std::vector<std::tuple<std::string, std::string, std::string>>& edges);

    // Parse weight string (e.g., "100G", "50M") to double
    double parse_weight(const std::string& weight) const;

    // Check if a module is degraded
    bool get_module_degraded(const std::string& end_module, bool active_active) const;

    // Calculate maximum flow between start and end modules
    double calculate_max_flow(const std::string& start_node_module, const std::string& leaf_node_module);

    // Virtual node management
    void add_virtual_nodes(const std::string& start_node_module, const std::string& leaf_node_module);
    void add_virtual_source(const std::string& start_node_module);
    void add_virtual_sink(const std::string& leaf_node_module, const std::set<std::string>& exception_nodes = {});
    void remove_virtual_source();
    void remove_virtual_sink();
    void remove_virtual_nodes();

    // Node operations
    void remove_node(const std::string& node);
    bool has_node(const std::string& node) const;
    std::vector<std::string> get_nodes() const;

    // Check if path exists between two nodes
    bool has_path(const std::string& source, const std::string& sink) const;

    // Maximum flow calculation using Edmonds-Karp algorithm
    double maximum_flow(const std::string& source, const std::string& sink);

    // LCA (Lowest Common Ancestor) functions
    // Build parent map for LCA calculation (call once after graph construction)
    void build_parent_map(const std::string& root_module);

    // Get path from node to root
    std::vector<std::string> get_path_to_root(const std::string& node) const;

    // Find LCA of two nodes
    std::string find_lca(const std::string& node1, const std::string& node2) const;

    // Find LCA of multiple nodes
    std::string find_lca_multiple(const std::vector<std::string>& nodes) const;

    // Get all nodes at the same level as LCA that are ancestors of the given nodes
    std::set<std::string> get_lca_level_nodes(const std::vector<std::string>& leaf_nodes) const;

    // Calculate max flow from LCA level nodes to a target SSD
    // Uses virtual source connected to all LCA level nodes
    double calculate_rebuild_max_flow(const std::vector<std::string>& source_ssds,
                                       const std::string& target_ssd);

    // Calculate rebuild max flow using io_module names directly
    // source_io_modules: io_modules that source SSDs are connected to
    // target_io_module: io_module that target SSD is connected to
    double calculate_rebuild_max_flow_via_io_modules(
        const std::vector<std::string>& source_io_modules,
        const std::string& target_io_module);

public:
    // Graph data structures
    std::map<std::string, std::map<std::string, double>> adjacency_list_;  // node -> (neighbor -> capacity)
    std::set<std::string> nodes_;
    std::vector<Edge> edges_;
    std::map<std::string, std::vector<std::string>> enclosures_;
    std::map<std::string, double> mttfs_;
    std::map<std::string, double> mtrs_;

    // Parent map for LCA calculation (child -> parent)
    std::map<std::string, std::string> parent_map_;
    // Depth map for LCA calculation (node -> depth from root)
    std::map<std::string, int> depth_map_;
    // Root node for LCA
    std::string root_node_;

private:

    // BFS for finding augmenting path
    bool bfs(const std::string& source, const std::string& sink,
             std::map<std::string, std::string>& parent,
             const std::map<std::string, std::map<std::string, double>>& residual_graph) const;
};
