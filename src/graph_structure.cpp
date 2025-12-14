#include "graph_structure.hpp"
#include "disk.hpp"
#include <queue>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <limits>

GraphStructure::GraphStructure() = default;

GraphStructure::GraphStructure(
    const std::vector<std::tuple<std::string, std::string, std::string>>& edges,
    const std::map<std::string, std::vector<std::string>>& enclosures,
    const std::map<std::string, double>& mttfs,
    const std::map<std::string, double>& mtrs)
    : enclosures_(enclosures), mttfs_(mttfs), mtrs_(mtrs)
{
    add_edges(edges);
}

void GraphStructure::add_edges(const std::vector<std::tuple<std::string, std::string, std::string>>& edges) {
    for (const auto& [from, to, weight] : edges) {
        auto [upstream, downstream] = parse_weight(weight);
        add_edge(from, to, upstream, downstream);
    }
}

void GraphStructure::add_edge(const std::string& from, const std::string& to,
                               double upstream_capacity, double downstream_capacity,
                               const std::string& label) {
    nodes_.insert(from);
    nodes_.insert(to);

    // Bidirectional adjacency: store both directions
    adjacency_list_[from][to] = {upstream_capacity, downstream_capacity};
    adjacency_list_[to][from] = {downstream_capacity, upstream_capacity};

    edges_.emplace_back(from, to, upstream_capacity, downstream_capacity, label);
}

void GraphStructure::add_disk_nodes(int total_disks,
                                     const std::map<int, std::string>& disk_to_io_module,
                                     double disk_read_bw,
                                     double disk_write_bw) {
    total_disks_ = total_disks;
    disk_to_io_module_ = disk_to_io_module;

    for (int i = 0; i < total_disks; ++i) {
        std::string disk_name = get_disk_name(i);
        auto it = disk_to_io_module.find(i);
        if (it != disk_to_io_module.end()) {
            // Connect disk to its io_module
            // disk -> io_module: READ (upstream)
            // io_module -> disk: WRITE (downstream)
            add_edge(disk_name, it->second, disk_read_bw, disk_write_bw);
        }
    }
}

std::pair<double, double> GraphStructure::parse_weight(const std::string& weight) const {
    // Handle format "100G/50G" for upstream/downstream
    size_t slash_pos = weight.find('/');
    if (slash_pos != std::string::npos) {
        std::string up_str = weight.substr(0, slash_pos);
        std::string down_str = weight.substr(slash_pos + 1);
        double upstream = Utils::KMG_to_bytes(up_str);
        double downstream = Utils::KMG_to_bytes(down_str);
        return {upstream, downstream};
    }

    // Single value: same for both directions
    double capacity = Utils::KMG_to_bytes(weight);
    return {capacity, capacity};
}

bool GraphStructure::get_module_degraded(const std::string& end_module, bool active_active) const {
    if (active_active) {
        return false;
    }

    int total_count = 0;
    int enabled_count = 0;

    for (const auto& node : nodes_) {
        if (node.find(end_module) == 0) {
            total_count++;
            if (disabled_nodes_.find(node) == disabled_nodes_.end()) {
                enabled_count++;
            }
        }
    }

    return enabled_count < total_count;
}

double GraphStructure::calculate_max_flow(const std::string& start_node_module,
                                           const std::string& leaf_node_module,
                                           FlowDirection direction) {
    add_virtual_nodes(start_node_module, leaf_node_module, direction);
    double flow = maximum_flow("virtual_source", "virtual_sink", direction);
    remove_virtual_nodes();
    return flow;
}

void GraphStructure::add_virtual_nodes(const std::string& start_node_module,
                                        const std::string& leaf_node_module,
                                        FlowDirection direction) {
    add_virtual_source(start_node_module, direction);
    add_virtual_sink(leaf_node_module, {}, direction);
}

void GraphStructure::add_virtual_source(const std::string& source_module, FlowDirection direction) {
    nodes_.insert("virtual_source");

    for (const auto& node : nodes_) {
        if (node.find(source_module) == 0 && node != "virtual_source" && node != "virtual_sink") {
            if (disabled_nodes_.find(node) == disabled_nodes_.end()) {
                adjacency_list_["virtual_source"][node] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
                adjacency_list_[node]["virtual_source"] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
            }
        }
    }
}

void GraphStructure::add_virtual_sink(const std::string& sink_module,
                                       const std::set<std::string>& exception_nodes,
                                       FlowDirection direction) {
    nodes_.insert("virtual_sink");

    for (const auto& node : nodes_) {
        if (node.find(sink_module) == 0 && node != "virtual_source" && node != "virtual_sink") {
            if (exception_nodes.find(node) == exception_nodes.end() &&
                disabled_nodes_.find(node) == disabled_nodes_.end()) {
                adjacency_list_[node]["virtual_sink"] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
                adjacency_list_["virtual_sink"][node] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
            }
        }
    }
}

void GraphStructure::remove_virtual_source() {
    nodes_.erase("virtual_source");
    adjacency_list_.erase("virtual_source");
    for (auto& [node, neighbors] : adjacency_list_) {
        neighbors.erase("virtual_source");
    }
}

void GraphStructure::remove_virtual_sink() {
    nodes_.erase("virtual_sink");
    adjacency_list_.erase("virtual_sink");
    for (auto& [node, neighbors] : adjacency_list_) {
        neighbors.erase("virtual_sink");
    }
}

void GraphStructure::remove_virtual_nodes() {
    remove_virtual_source();
    remove_virtual_sink();
}

void GraphStructure::remove_node(const std::string& node) {
    nodes_.erase(node);
    disabled_nodes_.erase(node);
    adjacency_list_.erase(node);

    for (auto& [n, neighbors] : adjacency_list_) {
        neighbors.erase(node);
    }

    // Remove from disk mapping if it's a disk
    if (is_disk_node(node)) {
        int disk_idx = get_disk_index(node);
        disk_to_io_module_.erase(disk_idx);
    }
}

void GraphStructure::disable_node(const std::string& node) {
    if (nodes_.find(node) != nodes_.end()) {
        disabled_nodes_.insert(node);
    }
}

void GraphStructure::enable_node(const std::string& node) {
    disabled_nodes_.erase(node);
}

bool GraphStructure::has_node(const std::string& node) const {
    return nodes_.find(node) != nodes_.end();
}

bool GraphStructure::is_node_enabled(const std::string& node) const {
    return has_node(node) && disabled_nodes_.find(node) == disabled_nodes_.end();
}

std::vector<std::string> GraphStructure::get_nodes() const {
    return std::vector<std::string>(nodes_.begin(), nodes_.end());
}

std::vector<std::string> GraphStructure::get_enabled_nodes() const {
    std::vector<std::string> result;
    for (const auto& node : nodes_) {
        if (disabled_nodes_.find(node) == disabled_nodes_.end()) {
            result.push_back(node);
        }
    }
    return result;
}

bool GraphStructure::has_path(const std::string& source, const std::string& sink,
                               FlowDirection direction) const {
    if (disabled_nodes_.find(source) != disabled_nodes_.end() ||
        disabled_nodes_.find(sink) != disabled_nodes_.end()) {
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
            for (const auto& [neighbor, caps] : it->second) {
                if (visited.find(neighbor) == visited.end() &&
                    disabled_nodes_.find(neighbor) == disabled_nodes_.end()) {
                    double capacity = (direction == FlowDirection::UPSTREAM) ?
                                     caps.first : caps.second;
                    if (capacity > 0) {
                        visited.insert(neighbor);
                        q.push(neighbor);
                    }
                }
            }
        }
    }
    return false;
}

std::map<std::string, std::map<std::string, double>> GraphStructure::build_residual_graph(
    FlowDirection direction) const {

    std::map<std::string, std::map<std::string, double>> residual;

    for (const auto& [node, neighbors] : adjacency_list_) {
        if (disabled_nodes_.find(node) != disabled_nodes_.end()) {
            continue;
        }

        for (const auto& [neighbor, caps] : neighbors) {
            if (disabled_nodes_.find(neighbor) != disabled_nodes_.end()) {
                continue;
            }

            double capacity = (direction == FlowDirection::UPSTREAM) ? caps.first : caps.second;

            // Apply bandwidth reservations
            auto key = std::make_pair(node, neighbor);
            if (direction == FlowDirection::UPSTREAM) {
                auto it = reserved_upstream_.find(key);
                if (it != reserved_upstream_.end()) {
                    capacity -= it->second;
                }
            } else {
                auto it = reserved_downstream_.find(key);
                if (it != reserved_downstream_.end()) {
                    capacity -= it->second;
                }
            }

            if (capacity > 0) {
                residual[node][neighbor] = capacity;
            }
        }
    }

    return residual;
}

bool GraphStructure::bfs(const std::string& source, const std::string& sink,
                          std::map<std::string, std::string>& parent,
                          const std::map<std::string, std::map<std::string, double>>& residual_graph,
                          FlowDirection direction) const {

    std::set<std::string> visited;
    std::queue<std::string> q;

    q.push(source);
    visited.insert(source);

    while (!q.empty()) {
        std::string current = q.front();
        q.pop();

        auto it = residual_graph.find(current);
        if (it == residual_graph.end()) {
            continue;
        }

        for (const auto& [neighbor, capacity] : it->second) {
            if (visited.find(neighbor) == visited.end() && capacity > 0) {
                visited.insert(neighbor);
                parent[neighbor] = current;

                if (neighbor == sink) {
                    return true;
                }
                q.push(neighbor);
            }
        }
    }

    return false;
}

double GraphStructure::maximum_flow(const std::string& source, const std::string& sink,
                                     FlowDirection direction) {
    if (disabled_nodes_.find(source) != disabled_nodes_.end() ||
        disabled_nodes_.find(sink) != disabled_nodes_.end()) {
        return 0.0;
    }

    auto residual = build_residual_graph(direction);
    double max_flow = 0.0;

    std::map<std::string, std::string> parent;

    while (bfs(source, sink, parent, residual, direction)) {
        // Find minimum residual capacity along the path
        double path_flow = std::numeric_limits<double>::max();
        std::string current = sink;

        while (current != source) {
            std::string prev = parent[current];
            path_flow = std::min(path_flow, residual[prev][current]);
            current = prev;
        }

        // Update residual capacities
        current = sink;
        while (current != source) {
            std::string prev = parent[current];
            residual[prev][current] -= path_flow;
            residual[current][prev] += path_flow;
            current = prev;
        }

        max_flow += path_flow;
        parent.clear();
    }

    return max_flow;
}

double GraphStructure::maximum_flow_with_flows(
    const std::string& source,
    const std::string& sink,
    FlowDirection direction,
    std::map<std::pair<std::string, std::string>, double>& flow_out) {

    flow_out.clear();

    if (disabled_nodes_.find(source) != disabled_nodes_.end() ||
        disabled_nodes_.find(sink) != disabled_nodes_.end()) {
        return 0.0;
    }

    auto residual = build_residual_graph(direction);
    double max_flow = 0.0;

    std::map<std::string, std::string> parent;

    while (bfs(source, sink, parent, residual, direction)) {
        double path_flow = std::numeric_limits<double>::max();
        std::string current = sink;

        while (current != source) {
            const std::string& prev = parent[current];
            auto prev_it = residual.find(prev);
            if (prev_it == residual.end()) {
                path_flow = 0.0;
                break;
            }
            auto cap_it = prev_it->second.find(current);
            if (cap_it == prev_it->second.end()) {
                path_flow = 0.0;
                break;
            }
            path_flow = std::min(path_flow, cap_it->second);
            current = prev;
        }

        if (path_flow <= 0.0 || path_flow == std::numeric_limits<double>::max()) {
            break;
        }

        current = sink;
        while (current != source) {
            const std::string prev = parent[current];

            residual[prev][current] -= path_flow;
            residual[current][prev] += path_flow;

            flow_out[{prev, current}] += path_flow;
            flow_out[{current, prev}] -= path_flow;

            current = prev;
        }

        max_flow += path_flow;
        parent.clear();
    }

    for (auto it = flow_out.begin(); it != flow_out.end(); ) {
        if (it->second <= 0.0) {
            it = flow_out.erase(it);
        } else {
            ++it;
        }
    }

    return max_flow;
}

// LCA functions
void GraphStructure::build_parent_map(const std::string& root_module) {
    parent_map_.clear();
    depth_map_.clear();

    // Find root nodes matching the module
    std::vector<std::string> root_nodes;
    for (const auto& node : nodes_) {
        if (node.find(root_module) == 0) {
            root_nodes.push_back(node);
        }
    }

    if (root_nodes.empty()) {
        return;
    }

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

    // BFS to build parent map
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

        auto it = adjacency_list_.find(current);
        if (it != adjacency_list_.end()) {
            for (const auto& [child, caps] : it->second) {
                if (visited.find(child) == visited.end() && caps.first > 0) {
                    visited.insert(child);
                    parent_map_[child] = current;
                    depth_map_[child] = current_depth + 1;
                    q.push(child);
                }
            }
        }
    }
}

std::vector<std::string> GraphStructure::get_path_to_root(const std::string& node) const {
    std::vector<std::string> path;

    auto it = parent_map_.find(node);
    if (it == parent_map_.end()) {
        return path;
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

std::string GraphStructure::find_lca(const std::string& node1, const std::string& node2) const {
    std::vector<std::string> path1 = get_path_to_root(node1);
    std::vector<std::string> path2 = get_path_to_root(node2);

    if (path1.empty() || path2.empty()) {
        return "";
    }

    std::set<std::string> path2_set(path2.begin(), path2.end());

    for (const auto& node : path1) {
        if (path2_set.find(node) != path2_set.end()) {
            return node;
        }
    }

    return "";
}

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

std::set<std::string> GraphStructure::get_lca_level_nodes(const std::vector<std::string>& leaf_nodes) const {
    std::set<std::string> lca_level_nodes;

    if (leaf_nodes.empty()) {
        return lca_level_nodes;
    }

    std::string lca = find_lca_multiple(leaf_nodes);
    if (lca.empty()) {
        return lca_level_nodes;
    }

    auto depth_it = depth_map_.find(lca);
    if (depth_it == depth_map_.end()) {
        return lca_level_nodes;
    }
    int lca_depth = depth_it->second;

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

double GraphStructure::calculate_rebuild_max_flow(const std::vector<int>& source_disk_indices,
                                                   int target_disk_index) {
    if (source_disk_indices.empty()) {
        return 0.0;
    }

    // Get io_modules for source and target disks
    std::vector<std::string> source_io_modules;
    for (int idx : source_disk_indices) {
        auto it = disk_to_io_module_.find(idx);
        if (it != disk_to_io_module_.end()) {
            source_io_modules.push_back(it->second);
        }
    }

    auto target_it = disk_to_io_module_.find(target_disk_index);
    if (target_it == disk_to_io_module_.end()) {
        return MAX_EDGE_VALUE;  // No network constraint
    }

    return calculate_rebuild_max_flow_via_io_modules(source_io_modules, target_it->second);
}

double GraphStructure::calculate_rebuild_max_flow_via_io_modules(
    const std::vector<std::string>& source_io_modules,
    const std::string& target_io_module) {

    if (source_io_modules.empty()) {
        return 0.0;
    }

    // Combine all io_modules
    std::vector<std::string> all_io_modules = source_io_modules;
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

    // If all same io_module, no network constraint
    std::set<std::string> unique_modules(all_io_modules.begin(), all_io_modules.end());
    if (unique_modules.size() == 1) {
        return MAX_EDGE_VALUE;
    }

    // Get LCA level nodes
    std::set<std::string> lca_level_nodes = get_lca_level_nodes(all_io_modules);

    if (lca_level_nodes.empty()) {
        return MAX_EDGE_VALUE;
    }

    // Check if LCA is at io_module level
    for (const auto& lca_node : lca_level_nodes) {
        for (const auto& io_mod : all_io_modules) {
            if (lca_node == io_mod) {
                return MAX_EDGE_VALUE;
            }
        }
    }

    // Add virtual source connected to all LCA level nodes
    nodes_.insert("rebuild_virtual_source");
    for (const auto& lca_node : lca_level_nodes) {
        adjacency_list_["rebuild_virtual_source"][lca_node] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
    }

    // Add virtual sink connected from target io_module
    nodes_.insert("rebuild_virtual_sink");
    adjacency_list_[target_io_module]["rebuild_virtual_sink"] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};

    // Calculate max flow
    double flow = maximum_flow("rebuild_virtual_source", "rebuild_virtual_sink", FlowDirection::DOWNSTREAM);

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

double GraphStructure::calculate_max_flow_from_root() {
    std::string root_node = find_root_node();
    if (root_node.empty()) {
        return MAX_EDGE_VALUE;  // No root found, assume unlimited
    }

    // Collect all io_modules
    std::vector<std::string> io_modules;
    for (const auto& node : nodes_) {
        if (node.find("io_module") != std::string::npos) {
            io_modules.push_back(node);
        }
    }

    if (io_modules.empty()) {
        return MAX_EDGE_VALUE;
    }

    // Add virtual sink connected from all io_modules
    nodes_.insert("backup_virtual_sink");
    for (const auto& io_mod : io_modules) {
        adjacency_list_[io_mod]["backup_virtual_sink"] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
    }

    // Calculate max flow from root to virtual sink
    double flow = maximum_flow(root_node, "backup_virtual_sink", FlowDirection::DOWNSTREAM);

    // Clean up
    nodes_.erase("backup_virtual_sink");
    for (auto& [node, neighbors] : adjacency_list_) {
        neighbors.erase("backup_virtual_sink");
    }

    return flow;
}

double GraphStructure::calculate_max_flow_to_root() {
    std::string root_node = find_root_node();
    if (root_node.empty()) {
        return MAX_EDGE_VALUE;
    }

    std::vector<std::string> io_modules;
    for (const auto& node : nodes_) {
        if (node.find("io_module") != std::string::npos) {
            io_modules.push_back(node);
        }
    }
    if (io_modules.empty()) {
        return MAX_EDGE_VALUE;
    }

    const std::string vsrc = "read_virtual_source";
    nodes_.insert(vsrc);

    for (const auto& io_mod : io_modules) {
        adjacency_list_[vsrc][io_mod] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
        adjacency_list_[io_mod][vsrc] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
    }

    double flow = maximum_flow(vsrc, root_node, FlowDirection::UPSTREAM);

    nodes_.erase(vsrc);
    adjacency_list_.erase(vsrc);
    for (auto& [node, neighbors] : adjacency_list_) {
        neighbors.erase(vsrc);
    }

    return flow;
}

std::string GraphStructure::find_root_node() const {
    if (nodes_.count("root") > 0) {
        return "root";
    }

    for (const auto& node : nodes_) {
        if (node.find("io_module") != std::string::npos ||
            node.find("disk") != std::string::npos ||
            node.find("virtual") != std::string::npos ||
            node.find("__rebuild_") != std::string::npos ||
            node.find("backup_virtual") != std::string::npos ||
            node.find("read_virtual") != std::string::npos) {
            continue;
        }

        bool has_outgoing = adjacency_list_.count(node) > 0 && !adjacency_list_.at(node).empty();
        if (has_outgoing) {
            return node;
        }
    }

    return "";
}

double GraphStructure::get_edge_capacity(const std::string& from, const std::string& to,
                                          FlowDirection direction) const {
    auto it = adjacency_list_.find(from);
    if (it == adjacency_list_.end()) {
        return 0.0;
    }

    auto neighbor_it = it->second.find(to);
    if (neighbor_it == it->second.end()) {
        return 0.0;
    }

    return (direction == FlowDirection::UPSTREAM) ?
           neighbor_it->second.first : neighbor_it->second.second;
}

void GraphStructure::reserve_bandwidth(const std::string& from, const std::string& to,
                                        double amount, FlowDirection direction) {
    auto key = std::make_pair(from, to);
    if (direction == FlowDirection::UPSTREAM) {
        reserved_upstream_[key] += amount;
    } else {
        reserved_downstream_[key] += amount;
    }
}

void GraphStructure::release_bandwidth(const std::string& from, const std::string& to,
                                        double amount, FlowDirection direction) {
    auto key = std::make_pair(from, to);
    if (direction == FlowDirection::UPSTREAM) {
        reserved_upstream_[key] -= amount;
        if (reserved_upstream_[key] <= 0) {
            reserved_upstream_.erase(key);
        }
    } else {
        reserved_downstream_[key] -= amount;
        if (reserved_downstream_[key] <= 0) {
            reserved_downstream_.erase(key);
        }
    }
}

GraphStructure GraphStructure::clone() const {
    GraphStructure copy;
    copy.adjacency_list_ = adjacency_list_;
    copy.nodes_ = nodes_;
    copy.disabled_nodes_ = disabled_nodes_;
    copy.edges_ = edges_;
    copy.enclosures_ = enclosures_;
    copy.mttfs_ = mttfs_;
    copy.mtrs_ = mtrs_;
    copy.disk_to_io_module_ = disk_to_io_module_;
    copy.total_disks_ = total_disks_;
    copy.parent_map_ = parent_map_;
    copy.depth_map_ = depth_map_;
    copy.root_node_ = root_node_;
    copy.reserved_upstream_ = reserved_upstream_;
    copy.reserved_downstream_ = reserved_downstream_;
    return copy;
}
