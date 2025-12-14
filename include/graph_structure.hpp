#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <queue>
#include "utils.hpp"

constexpr double MAX_EDGE_VALUE = 1.0e15;

// Direction for bandwidth in bidirectional graph
enum class FlowDirection {
    UPSTREAM,       // From leaf (disk) towards root (switch) - READ direction
    DOWNSTREAM      // From root (switch) towards leaf (disk) - WRITE direction
};

// Edge structure for bidirectional graph
struct Edge {
    std::string from;
    std::string to;
    double upstream_capacity;    // Capacity from->to (READ direction: disk->switch)
    double downstream_capacity;  // Capacity to->from (WRITE direction: switch->disk)
    std::string label;

    Edge() : upstream_capacity(0), downstream_capacity(0) {}
    Edge(const std::string& f, const std::string& t, double up_cap, double down_cap, const std::string& lbl = "")
        : from(f), to(t), upstream_capacity(up_cap), downstream_capacity(down_cap), label(lbl) {}
};

// Graph structure for hardware topology and max flow calculation
// Supports bidirectional edges with different capacities for read/write
class GraphStructure {
public:
    GraphStructure();

    // Constructor with edges (format: from, to, capacity_string)
    // For bidirectional, use format "100G/50G" for upstream/downstream
    GraphStructure(
        const std::vector<std::tuple<std::string, std::string, std::string>>& edges,
        const std::map<std::string, std::vector<std::string>>& enclosures,
        const std::map<std::string, double>& mttfs,
        const std::map<std::string, double>& mtrs
    );

    // Add edges to graph
    void add_edges(const std::vector<std::tuple<std::string, std::string, std::string>>& edges);

    // Add a single bidirectional edge
    void add_edge(const std::string& from, const std::string& to,
                  double upstream_capacity, double downstream_capacity,
                  const std::string& label = "");

    // Add disk nodes to the graph
    // Connects disks to their io_modules based on mapping
    void add_disk_nodes(int total_disks,
                        const std::map<int, std::string>& disk_to_io_module,
                        double disk_read_bw,
                        double disk_write_bw);

    // Parse weight string (e.g., "100G", "50M", "100G/50G" for up/down)
    std::pair<double, double> parse_weight(const std::string& weight) const;

    // Check if a module is degraded
    bool get_module_degraded(const std::string& end_module, bool active_active) const;

    // Calculate maximum flow between start and end modules
    double calculate_max_flow(const std::string& start_node_module,
                              const std::string& leaf_node_module,
                              FlowDirection direction = FlowDirection::UPSTREAM);

    // Virtual node management
    void add_virtual_nodes(const std::string& start_node_module,
                           const std::string& leaf_node_module,
                           FlowDirection direction = FlowDirection::UPSTREAM);
    void add_virtual_source(const std::string& source_module,
                            FlowDirection direction = FlowDirection::UPSTREAM);
    void add_virtual_sink(const std::string& sink_module,
                          const std::set<std::string>& exception_nodes = {},
                          FlowDirection direction = FlowDirection::UPSTREAM);
    void remove_virtual_source();
    void remove_virtual_sink();
    void remove_virtual_nodes();

    // Node operations
    void remove_node(const std::string& node);
    void disable_node(const std::string& node);  // Temporarily disable without removing
    void enable_node(const std::string& node);   // Re-enable a disabled node
    bool has_node(const std::string& node) const;
    bool is_node_enabled(const std::string& node) const;
    std::vector<std::string> get_nodes() const;
    std::vector<std::string> get_enabled_nodes() const;

    // Check if path exists between two nodes
    bool has_path(const std::string& source, const std::string& sink,
                  FlowDirection direction = FlowDirection::UPSTREAM) const;

    // Maximum flow calculation using Edmonds-Karp algorithm
    double maximum_flow(const std::string& source, const std::string& sink,
                        FlowDirection direction = FlowDirection::UPSTREAM);

    // Maximum flow calculation that also outputs per-edge flow (directed).
    // `flow_out[(u,v)]` is the net flow sent on u->v (only positive entries are meaningful).
    double maximum_flow_with_flows(
        const std::string& source,
        const std::string& sink,
        FlowDirection direction,
        std::map<std::pair<std::string, std::string>, double>& flow_out);

    // LCA (Lowest Common Ancestor) functions
    void build_parent_map(const std::string& root_module);
    std::vector<std::string> get_path_to_root(const std::string& node) const;
    std::string find_lca(const std::string& node1, const std::string& node2) const;
    std::string find_lca_multiple(const std::vector<std::string>& nodes) const;
    std::set<std::string> get_lca_level_nodes(const std::vector<std::string>& leaf_nodes) const;

    // Calculate rebuild max flow from source disks to target disk
    // source_disks: indices of disks that can provide data for rebuild
    // target_disk: index of disk being rebuilt
    // Returns max flow considering network topology
    double calculate_rebuild_max_flow(const std::vector<int>& source_disk_indices,
                                       int target_disk_index);

    // Calculate rebuild max flow using io_module names
    double calculate_rebuild_max_flow_via_io_modules(
        const std::vector<std::string>& source_io_modules,
        const std::string& target_io_module);

    // Calculate max flow from root to all io_modules (for backup store recovery)
    double calculate_max_flow_from_root();

    // Calculate max flow from all io_modules to root (host READ direction).
    double calculate_max_flow_to_root();

    // Best-effort root selection used by calculate_max_flow_* helpers.
    std::string find_root_node() const;

    // Get edge capacity in specified direction
    double get_edge_capacity(const std::string& from, const std::string& to,
                             FlowDirection direction) const;

    // Modify edge capacity (for bandwidth reservation during rebuild)
    void reserve_bandwidth(const std::string& from, const std::string& to,
                           double amount, FlowDirection direction);
    void release_bandwidth(const std::string& from, const std::string& to,
                           double amount, FlowDirection direction);

    // Clone the graph structure
    GraphStructure clone() const;

public:
    // Graph data structures
    // Bidirectional adjacency: node -> (neighbor -> {upstream_cap, downstream_cap})
    std::map<std::string, std::map<std::string, std::pair<double, double>>> adjacency_list_;
    std::set<std::string> nodes_;
    std::set<std::string> disabled_nodes_;  // Nodes that are temporarily disabled
    std::vector<Edge> edges_;
    std::map<std::string, std::vector<std::string>> enclosures_;
    std::map<std::string, double> mttfs_;
    std::map<std::string, double> mtrs_;

    // Disk-related data
    std::map<int, std::string> disk_to_io_module_;  // disk_index -> io_module_name
    int total_disks_ = 0;

    // Parent map for LCA calculation (child -> parent)
    std::map<std::string, std::string> parent_map_;
    // Depth map for LCA calculation (node -> depth from root)
    std::map<std::string, int> depth_map_;
    // Root node for LCA
    std::string root_node_;

    // Bandwidth reservation tracking
    std::map<std::pair<std::string, std::string>, double> reserved_upstream_;
    std::map<std::pair<std::string, std::string>, double> reserved_downstream_;

private:
    // BFS for finding augmenting path
    bool bfs(const std::string& source, const std::string& sink,
             std::map<std::string, std::string>& parent,
             const std::map<std::string, std::map<std::string, double>>& residual_graph,
             FlowDirection direction) const;

    // Build residual graph for max flow in specified direction
    std::map<std::string, std::map<std::string, double>> build_residual_graph(
        FlowDirection direction) const;
};
