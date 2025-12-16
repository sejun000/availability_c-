#include "simulation.hpp"
#include "simulation_helpers.hpp"
#include "utils.hpp"
#include "logger.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <thread>
#include <future>
#include <random>
#include <fstream>
#include <limits>

// Calculate probability of stripe failure under x device failures
double pfail(int m, int k, int l, int x) {
    int N = m + k + l;
    int S = m + k;

    if (x > N) return 1.0;

    auto comb = [](int n, int r) -> double {
        if (r > n || r < 0) return 0.0;
        if (r == 0 || r == n) return 1.0;

        double result = 1.0;
        for (int i = 0; i < r; ++i) {
            result *= (n - i);
            result /= (i + 1);
        }
        return result;
    };

    double total = comb(N, S);
    double surviving = 0.0;

    for (int f = 0; f <= std::min(k, x); ++f) {
        surviving += comb(x, f) * comb(N - x, S - f);
    }

    return 1.0 - surviving / total;
}

// Calculate bottleneck speed considering erasure coding and bandwidth
double calculate_bottleneck_speed(int m, int k, const std::vector<double>& other_bws,
                                   double rebuild_bw_ratio, double ec_encoding_speed) {
    double min_speed = ec_encoding_speed;

    for (const auto& other_bw : other_bws) {
        double effective_bw = other_bw * rebuild_bw_ratio;
        if (effective_bw < min_speed) {
            min_speed = effective_bw;
        }
    }

    return min_speed;
}

// Check if data loss occurred
bool is_data_loss(const ECGroupFailureInfo& failure_info, const ErasureCodingScheme& ec_scheme) {
    // Combine failed and disconnected disk indices
    std::set<int> all_unavailable = failure_info.failed_disk_indices;
    all_unavailable.insert(failure_info.disconnected_disk_indices.begin(),
                          failure_info.disconnected_disk_indices.end());
    return ec_scheme.is_data_loss(all_unavailable);
}

// Judge disk state from failure information
std::string judge_state_from_failure_info(const ECGroupFailureInfo& failure_info,
                                          const ErasureCodingScheme& ec_scheme) {
    // Combine failed and disconnected disk indices
    std::set<int> all_unavailable = failure_info.failed_disk_indices;
    all_unavailable.insert(failure_info.disconnected_disk_indices.begin(),
                          failure_info.disconnected_disk_indices.end());

    if (all_unavailable.empty()) {
        return DISK_STATE_NORMAL;
    }

    if (ec_scheme.is_data_loss(all_unavailable)) {
        return DISK_STATE_DATA_LOSS;
    }

    return DISK_STATE_REBUILDING;
}

// Random number generator (thread-local)
thread_local std::mt19937 rng(std::random_device{}());

// Seed the thread-local RNG for reproducibility
void seed_rng(uint64_t seed) {
    rng.seed(static_cast<std::mt19937::result_type>(seed));
}

// Push failure event to priority queue
void push_failed_event(std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
                       const std::string& event_node,
                       double current_time,
                       const std::map<std::string, std::string>& node_to_module_map,
                       const GraphStructure& hardware_graph,
                       double disk_mttf,
                       const EmpiricalCDF* disk_failure_cdf,
                       double disk_arr) {
    double failure_time;

    if (is_disk_node(event_node)) {
        // Use spliced distribution if CDF and ARR are available
        if (disk_failure_cdf != nullptr && disk_failure_cdf->is_initialized() && disk_arr > 0) {
            failure_time = current_time + disk_failure_cdf->sample_spliced(rng, disk_arr);
        } else {
            std::exponential_distribution<double> exp_dist(1.0 / disk_mttf);
            failure_time = current_time + exp_dist(rng);
        }
    } else {
        auto it = node_to_module_map.find(event_node);
        if (it == node_to_module_map.end()) return;
        std::string module = it->second;
        auto mttf_it = hardware_graph.mttfs_.find(module);
        if (mttf_it == hardware_graph.mttfs_.end()) return;
        double mttf = mttf_it->second;
        std::exponential_distribution<double> exp_dist(1.0 / mttf);
        failure_time = current_time + exp_dist(rng);
    }

    Event event{failure_time, "fail", event_node, current_time, current_time, false};
    failed_events.push(event);
}

// Push repair event to priority queue
void push_repair_event(std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events,
                       const std::string& event_node,
                       double current_time,
                       const std::map<std::string, std::string>& node_to_module_map,
                       const GraphStructure& hardware_graph,
                       bool software_repair_only,
                       const nlohmann::json* options) {
    if (is_disk_node(event_node)) {
        // Disk repair: will be updated when rebuild completes
        Event event{current_time + BIG_NUMBER, "repair", event_node, current_time, current_time, software_repair_only};
        repair_events.push(event);
    } else {
        auto it = node_to_module_map.find(event_node);
        if (it == node_to_module_map.end()) return;
        std::string module = it->second;
        auto mtr_it = hardware_graph.mtrs_.find(module);
        if (mtr_it == hardware_graph.mtrs_.end()) return;

        double mtr = mtr_it->second;
        bool is_hardware_fault = false;

        // Check if this is io_module and apply software/hardware fault distinction
        if (module.find("io_module") != std::string::npos && options != nullptr) {
            double software_fault_ratio = options->value("io_module_software_fault_ratio", 1.0);
            double hardware_mtr = options->value("io_module_hardware_mtr", mtr);

            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            double rand_val = uniform(rng);

            if (rand_val >= software_fault_ratio) {
                mtr = hardware_mtr;
                is_hardware_fault = true;
            }
        }

        double repair_time = current_time + mtr;
        Event event{repair_time, "repair", event_node, current_time, current_time, software_repair_only || !is_hardware_fault};
        repair_events.push(event);
    }
}

// Pop event from priority queues
Event pop_event(std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& events,
                std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events) {
    if (!repair_events.empty() && events.empty()) {
        Event e = repair_events.top();
        repair_events.pop();
        return e;
    } else if (repair_events.empty() && !events.empty()) {
        Event e = events.top();
        events.pop();
        return e;
    } else {
        const Event& repair_event = repair_events.top();
        const Event& event = events.top();
        if (repair_event.time < event.time) {
            Event e = repair_events.top();
            repair_events.pop();
            return e;
        } else {
            Event e = events.top();
            events.pop();
            return e;
        }
    }
}

// Generate first failure events for all components
std::priority_queue<Event, std::vector<Event>, std::greater<Event>> generate_first_failure_events(
    const GraphStructure& hardware_graph,
    const std::map<std::string, std::string>& node_to_module_map,
    int total_disks,
    double disk_mttf,
    const EmpiricalCDF* disk_failure_cdf,
    double disk_arr) {

    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> events;

    // Add hardware node events
    for (const auto& node : hardware_graph.nodes_) {
        if (!is_disk_node(node)) {
            push_failed_event(events, node, 0, node_to_module_map, hardware_graph, disk_mttf, disk_failure_cdf, disk_arr);
        }
    }

    // Add enclosure events
    for (const auto& [enclosure, _] : hardware_graph.enclosures_) {
        push_failed_event(events, enclosure, 0, node_to_module_map, hardware_graph, disk_mttf, disk_failure_cdf, disk_arr);
    }

    // Add disk events
    for (int disk_index = 0; disk_index < total_disks; ++disk_index) {
        std::string disk_name = get_disk_name(disk_index);
        push_failed_event(events, disk_name, 0, node_to_module_map, hardware_graph, disk_mttf, disk_failure_cdf, disk_arr);
    }

    return events;
}

// Update failure info when event occurs
void update_failure_info(const std::string& event_type,
                        const std::string& event_node,
                        std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
                        std::map<std::string, bool>& failed_nodes_and_enclosures,
                        const ErasureCodingScheme& ec_scheme) {
    if (event_type != "fail" && event_type != "repair") {
        return;
    }
    if (is_disk_node(event_node)) {
        int disk_index = get_disk_index(event_node);
        int group_index = ec_scheme.get_disk_group(disk_index);

        if (event_type == "fail") {
            failure_info_per_ec_group[group_index].disk_failure_count++;
            failure_info_per_ec_group[group_index].failed_disk_indices.insert(disk_index);
        } else if (event_type == "repair") {
            if (failure_info_per_ec_group[group_index].failed_disk_indices.erase(disk_index) > 0) {
                failure_info_per_ec_group[group_index].disk_failure_count--;
            }
        }
    } else {
        if (event_type == "fail") {
            failed_nodes_and_enclosures[event_node] = true;
        } else if (event_type == "repair") {
            failed_nodes_and_enclosures[event_node] = false;
        }
    }
}

// Calculate flows and speeds for given hardware and failure state
void calculate_flows_and_speed(
    GraphStructure& hardware_graph_copy,
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    const ErasureCodingScheme& ec_scheme,
    const SimulationParams& params,
    std::unordered_map<FailureStateKey, FlowsAndSpeedEntry, FailureStateKeyHash>& flows_and_speed_table,
    double max_read_performance_without_any_failure,
    double max_write_performance_without_any_failure,
    const std::string& start_module,
    const std::string& end_module,
    const DisconnectedStatus& disconnected,
    const FailureStateKey& key,
    const DiskIOModuleManager* disk_io_manager,
    const std::map<int, DiskInfo>* disks) {

    if (flows_and_speed_table.find(key) != flows_and_speed_table.end()) {
        return;
    }

    FlowsAndSpeedEntry entry;
    const ECConfig& ec_config = ec_scheme.get_config();

    int group_size = ec_scheme.get_group_size();

    // IMPORTANT: do not mutate the cached graph in-place.
    // Apply rebuild/degraded reservations on a clone for this flow calculation.
    GraphStructure graph = hardware_graph_copy.clone();

    auto reserve_scaled_flow = [&](const std::map<std::pair<std::string, std::string>, double>& flow_map,
                                   double scale,
                                   FlowDirection direction) {
        if (scale <= 0.0) return;
        for (const auto& [edge, amount] : flow_map) {
            const auto& [from, to] = edge;
            if (from.rfind("__rebuild_src_", 0) == 0 || to.rfind("__rebuild_src_", 0) == 0 ||
                from.rfind("__rebuild_sink_", 0) == 0 || to.rfind("__rebuild_sink_", 0) == 0) {
                continue;
            }
            graph.reserve_bandwidth(from, to, amount * scale, direction);
        }
    };

    uint64_t vcounter = 0;
    auto max_flow_from_sources_to_module = [&](const std::vector<std::string>& sources,
                                               const std::string& sink_module_prefix,
                                               FlowDirection direction,
                                               std::map<std::pair<std::string, std::string>, double>& flow_map) -> double {
        flow_map.clear();
        if (sources.empty() || sink_module_prefix.empty()) {
            return 0.0;
        }

        // Create unique virtual nodes to avoid collisions.
        const std::string vsrc = "__rebuild_src_" + std::to_string(++vcounter);
        const std::string vsink = "__rebuild_sink_" + std::to_string(vcounter);

        graph.nodes_.insert(vsrc);
        graph.nodes_.insert(vsink);

        // Super source -> all sources
        for (const auto& src : sources) {
            if (src.empty()) continue;
            graph.nodes_.insert(src);
            graph.adjacency_list_[vsrc][src] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
            graph.adjacency_list_[src][vsrc] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
        }

        // Connect all nodes in module prefix -> vsink (virtual sink)
        for (const auto& node : graph.nodes_) {
            if (node == vsrc || node == vsink) continue;
            if (node.find(sink_module_prefix) == 0) {
                graph.adjacency_list_[node][vsink] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
                graph.adjacency_list_[vsink][node] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
            }
        }

        double flow = graph.maximum_flow_with_flows(vsrc, vsink, direction, flow_map);

        // Cleanup virtual nodes/edges
        graph.nodes_.erase(vsrc);
        graph.nodes_.erase(vsink);
        graph.adjacency_list_.erase(vsrc);
        graph.adjacency_list_.erase(vsink);
        for (auto& [node, neighbors] : graph.adjacency_list_) {
            neighbors.erase(vsrc);
            neighbors.erase(vsink);
        }

        return flow;
    };

    auto max_flow_from_module_to_target = [&](const std::string& source_module_prefix,
                                              const std::string& target_node,
                                              FlowDirection direction,
                                              std::map<std::pair<std::string, std::string>, double>& flow_map) -> double {
        flow_map.clear();
        if (source_module_prefix.empty() || target_node.empty()) return 0.0;

        const std::string vsrc = "__rebuild_src_" + std::to_string(++vcounter);
        const std::string vsink = "__rebuild_sink_" + std::to_string(vcounter);

        graph.nodes_.insert(vsrc);
        graph.nodes_.insert(vsink);

        for (const auto& node : graph.nodes_) {
            if (node == vsrc || node == vsink) continue;
            if (node.find(source_module_prefix) == 0) {
                graph.adjacency_list_[vsrc][node] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
                graph.adjacency_list_[node][vsrc] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
            }
        }

        graph.adjacency_list_[target_node][vsink] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};
        graph.adjacency_list_[vsink][target_node] = {MAX_EDGE_VALUE, MAX_EDGE_VALUE};

        double flow = graph.maximum_flow_with_flows(vsrc, vsink, direction, flow_map);

        graph.nodes_.erase(vsrc);
        graph.nodes_.erase(vsink);
        graph.adjacency_list_.erase(vsrc);
        graph.adjacency_list_.erase(vsink);
        for (auto& [node, neighbors] : graph.adjacency_list_) {
            neighbors.erase(vsrc);
            neighbors.erase(vsink);
        }

        return flow;
    };

    auto comb = [](int n, int k) -> double {
        if (k < 0 || n < 0 || k > n) return 0.0;
        if (k == 0 || k == n) return 1.0;
        k = std::min(k, n - k);
        double result = 1.0;
        for (int i = 1; i <= k; ++i) {
            result *= (n - (k - i));
            result /= i;
        }
        return result;
    };

    auto declustered_no_degraded_ratio = [&](int n, int stripe_width, int unavailable) -> double {
        if (unavailable <= 0) return 1.0;
        if (stripe_width <= 0) return 1.0;
        if (n <= 0) return 1.0;
        if (stripe_width >= n) return 0.0; // if any unavailable exists, stripe must intersect it
        if (n - unavailable < stripe_width) return 0.0;
        double denom = comb(n, stripe_width);
        if (denom <= 0) return 0.0;
        double numer = comb(n - unavailable, stripe_width);
        double r = numer / denom;
        if (r < 0) r = 0;
        if (r > 1) r = 1;
        return r;
    };

    auto stripe_width_for_ec = [&](const ECConfig& cfg) -> int {
        if (cfg.type == ECType::STANDARD) {
            return cfg.m + cfg.k;
        }
        if (cfg.type == ECType::LRC) {
            if (cfg.local_m <= 0) return cfg.m + cfg.k;
            int num_local_groups = cfg.m / cfg.local_m;
            return cfg.m + cfg.k + cfg.local_k * num_local_groups;
        }
        // MULTI_EC: handled separately (outer + local).
        return cfg.m + cfg.k;
    };

    // Process each EC group
    double total_data_loss = 0.0;
    int total_groups = 0;
    int groups_with_data_loss = 0;  // Count of EC groups with data loss

    // Get total number of EC groups
    int num_ec_groups = ec_scheme.get_total_groups(params.total_disks);

    for (int group_index = 0; group_index < num_ec_groups; ++group_index) {
        total_groups++;
        std::vector<int> group_disks = ec_scheme.get_group_disks(group_index);
        if (group_disks.empty()) {
            // Fallback for sequential mapping
            int start_disk_idx = ec_scheme.get_group_start_disk(group_index);
            for (int i = 0; i < group_size; ++i) {
                group_disks.push_back(start_disk_idx + i);
            }
        }

        // Get failure info for this group (may not exist if no failures)
        ECGroupFailureInfo failure_info;
        auto it = failure_info_per_ec_group.find(group_index);
        if (it != failure_info_per_ec_group.end()) {
            failure_info = it->second;
        }

        // Calculate effective failure count (failed + disconnected)
        std::set<int> all_unavailable = failure_info.failed_disk_indices;
        for (int disk_idx : group_disks) {
            if (disconnected.is_disk_disconnected(disk_idx)) {
                all_unavailable.insert(disk_idx);
            }
        }
        int effective_failure_count = static_cast<int>(all_unavailable.size());

        // Determine state based on EC scheme (properly handles LRC/Multi-EC)
        std::string state;
        if (effective_failure_count == 0) {
            state = DISK_STATE_NORMAL;
        } else if (ec_scheme.is_data_loss(all_unavailable)) {
            state = DISK_STATE_DATA_LOSS;
        } else {
            state = DISK_STATE_REBUILDING;
        }

        // Compute declustered parity relief ratio r:
        // r = P(random stripe uses only healthy disks) = C(n-x, S) / C(n, S)
        // degraded/rebuild impact happens with (1-r).
        double impact_ratio = 1.0;
        if (ec_config.type != ECType::MULTI_EC) {
            const int stripe_width = stripe_width_for_ec(ec_config);
            const double r_no_degraded = declustered_no_degraded_ratio(group_size, stripe_width, effective_failure_count);
            impact_ratio = 1.0 - r_no_degraded;
        } else {
            // MULTI_EC: model declustered relief at both levels.
            // - Local level: local stripe uses (local_m+local_k) out of local_n disks.
            // - Outer level: outer stripe uses (m+k) chunks out of n chunks.
            const int local_n = ec_config.local_n;
            const int local_width = ec_config.local_m + ec_config.local_k;
            const int outer_n = ec_config.n;
            const int outer_width = ec_config.m + ec_config.k;

            double r_local_avg = 1.0;
            int stripe_count = std::max(1, ec_config.m + ec_config.k);
            if (local_n > 0 && local_width > 0) {
                double sum = 0.0;
                for (int stripe_idx = 0; stripe_idx < stripe_count; ++stripe_idx) {
                    int unavailable_in_stripe = 0;
                    int start_pos = stripe_idx * local_n;
                    int end_pos = std::min(static_cast<int>(group_disks.size()), start_pos + local_n);
                    for (int pos = start_pos; pos < end_pos; ++pos) {
                        int d = group_disks[pos];
                        if (all_unavailable.count(d)) {
                            unavailable_in_stripe++;
                        }
                    }
                    sum += declustered_no_degraded_ratio(local_n, local_width, unavailable_in_stripe);
                }
                r_local_avg = sum / stripe_count;
            }

            int failed_chunks = 0;
            if (local_n > 0 && ec_config.local_k >= 0) {
                for (int stripe_idx = 0; stripe_idx < stripe_count; ++stripe_idx) {
                    int unavailable_in_stripe = 0;
                    int start_pos = stripe_idx * local_n;
                    int end_pos = std::min(static_cast<int>(group_disks.size()), start_pos + local_n);
                    for (int pos = start_pos; pos < end_pos; ++pos) {
                        int d = group_disks[pos];
                        if (all_unavailable.count(d)) {
                            unavailable_in_stripe++;
                        }
                    }
                    if (unavailable_in_stripe > ec_config.local_k) {
                        failed_chunks++;
                    }
                }
            }

            double r_outer = 1.0;
            if (outer_n > 0 && outer_width > 0) {
                r_outer = declustered_no_degraded_ratio(outer_n, outer_width, failed_chunks);
            }

            double r_total = std::clamp(r_local_avg * r_outer, 0.0, 1.0);
            impact_ratio = 1.0 - r_total;
        }

        // Reserve degraded read bandwidth (workload assumed 100%, degraded_ratio share).
        // This consumes extra READ bandwidth on the network from (n-k) sources.
        if (impact_ratio > 0 && params.degraded_ratio > 0 && disk_io_manager != nullptr && effective_failure_count > 0) {
            // Choose sources from available disks.
            std::vector<int> available_disks;
            available_disks.reserve(group_disks.size());
            for (int d : group_disks) {
                if (!all_unavailable.count(d)) {
                    available_disks.push_back(d);
                }
            }
            std::sort(available_disks.begin(), available_disks.end());

            // Global degraded read uses (n-k) reads; local degraded read not modeled separately here.
            int read_sources = std::max(0, group_size - ec_config.k);
            read_sources = std::min(read_sources, static_cast<int>(available_disks.size()));

            if (read_sources > 0) {
                // First, find the most common io_module among sources to prefer it
                std::map<std::string, int> io_module_count;
                for (int i = 0; i < read_sources; ++i) {
                    auto src_io_modules = disk_io_manager->get_io_modules_for_disk(available_disks[i]);
                    for (const auto& io_mod : src_io_modules) {
                        auto it = disconnected.failed_nodes.find(io_mod);
                        if (it == disconnected.failed_nodes.end() || !it->second) {
                            io_module_count[io_mod]++;
                        }
                    }
                }
                std::string preferred_io_module;
                int max_count = 0;
                for (const auto& [io_mod, count] : io_module_count) {
                    if (count > max_count) {
                        max_count = count;
                        preferred_io_module = io_mod;
                    }
                }

                std::vector<std::string> source_modules;
                source_modules.reserve(static_cast<size_t>(read_sources));
                for (int i = 0; i < read_sources; ++i) {
                    // Get an active (non-failed) io_module for each source disk
                    // Prefer using the most common io_module to minimize cross-switch traffic
                    std::string src_module;
                    std::string fallback_module;
                    auto src_io_modules = disk_io_manager->get_io_modules_for_disk(available_disks[i]);
                    for (const auto& io_mod : src_io_modules) {
                        auto it = disconnected.failed_nodes.find(io_mod);
                        if (it == disconnected.failed_nodes.end() || !it->second) {
                            if (fallback_module.empty()) {
                                fallback_module = io_mod;
                            }
                            if (io_mod == preferred_io_module) {
                                src_module = io_mod;
                                break;
                            }
                        }
                    }
                    if (src_module.empty()) {
                        src_module = fallback_module.empty() ?
                            disk_io_manager->get_primary_io_module(available_disks[i]) : fallback_module;
                    }
                    source_modules.push_back(src_module);
                }

                // Find LCA to determine if degraded read traffic crosses switch
                // Build parent map if not already built
                if (graph.parent_map_.empty()) {
                    graph.build_parent_map(start_module);
                }

                std::string lca = graph.find_lca_multiple(source_modules);

                // Check if all sources are on the same io_module
                std::set<std::string> unique_modules(source_modules.begin(), source_modules.end());
                bool is_intra_io_module = (unique_modules.size() == 1);

                // Also check if LCA is at io_module level
                if (!is_intra_io_module && !lca.empty()) {
                    auto lca_depth_it = graph.depth_map_.find(lca);
                    if (lca_depth_it != graph.depth_map_.end() && !source_modules.empty()) {
                        auto src_depth_it = graph.depth_map_.find(source_modules[0]);
                        if (src_depth_it != graph.depth_map_.end()) {
                            if (lca_depth_it->second == src_depth_it->second) {
                                is_intra_io_module = true;
                            }
                        }
                    }
                }

                // degraded read speed per spec-ish: min((n-k)*disk_read_bw*degraded_ratio, ec_speed) * impact_ratio
                // For replication (or MULTI_EC with local replication), no EC encoding needed - skip ec_encoding_speed limit
                double degraded_speed = static_cast<double>(read_sources) * params.disk_read_bw * params.degraded_ratio;
                bool skip_ec_encoding = (ec_config.type == ECType::REPLICATION) ||
                    (ec_config.type == ECType::MULTI_EC && ec_config.local_type == LocalType::REPLICATION);
                if (!skip_ec_encoding) {
                    degraded_speed = std::min(degraded_speed, params.ec_encoding_speed);
                }
                degraded_speed *= impact_ratio;

                // Only reserve cross-switch bandwidth if traffic actually crosses io_modules
                if (!is_intra_io_module) {
                    std::string flow_meeting_point = lca.empty() ? start_module : lca;
                    std::map<std::pair<std::string, std::string>, double> flow_map;
                    double max_up = max_flow_from_sources_to_module(source_modules, flow_meeting_point, FlowDirection::UPSTREAM, flow_map);

                    if (max_up > 0 && degraded_speed > 0) {
                        reserve_scaled_flow(flow_map, degraded_speed / max_up, FlowDirection::UPSTREAM);
                    }
                }
            }
        }

        // Reserve rebuild bandwidth by actually consuming edge capacities along rebuild flows.
        // NOTE: This implements the spec's "subtract flow along rebuild path" behavior at module-graph level.
        if (state == DISK_STATE_REBUILDING && disk_io_manager != nullptr && effective_failure_count > 0 && disks != nullptr) {
            std::vector<int> available_disks;
            available_disks.reserve(group_disks.size());
            for (int d : group_disks) {
                if (!all_unavailable.count(d)) {
                    available_disks.push_back(d);
                }
            }
            std::sort(available_disks.begin(), available_disks.end());

            // Rebuild targets are disks that need rebuild and are currently connected (not disconnected).
            std::vector<int> targets;
            for (int d : group_disks) {
                auto dit = disks->find(d);
                if (dit == disks->end()) continue;
                const DiskInfo& info = dit->second;
                if (!info.needs_rebuild()) continue;
                if (info.is_disconnected) continue; // cannot rebuild while disconnected
                targets.push_back(d);
            }

            for (int target_disk : targets) {
                if (disk_io_manager == nullptr) break;

                // Get an active (non-failed) io_module for the target disk
                std::string target_module;
                auto target_io_modules = disk_io_manager->get_io_modules_for_disk(target_disk);
                for (const auto& io_mod : target_io_modules) {
                    auto it = disconnected.failed_nodes.find(io_mod);
                    if (it == disconnected.failed_nodes.end() || !it->second) {
                        target_module = io_mod;
                        break;
                    }
                }
                if (target_module.empty()) continue;

                int required_reads = ec_scheme.get_rebuild_read_count(all_unavailable, target_disk);
                if (required_reads <= 0) continue;
                if (static_cast<int>(available_disks.size()) < required_reads) {
                    continue;
                }

                // For LRC / Multi-EC: if local rebuild is possible, use local sources within local group first.
                std::vector<int> chosen_sources;
                chosen_sources.reserve(static_cast<size_t>(required_reads));

                if (ec_config.type == ECType::LRC || ec_config.type == ECType::MULTI_EC) {
                    if (ec_scheme.can_local_rebuild(all_unavailable, target_disk) && ec_config.local_m > 0) {
                        int relative_target = -1;
                        for (size_t i = 0; i < group_disks.size(); ++i) {
                            if (group_disks[i] == target_disk) {
                                relative_target = static_cast<int>(i);
                                break;
                            }
                        }
                        int local_group_size = (ec_config.type == ECType::MULTI_EC) ? ec_config.local_n : (ec_config.local_m + ec_config.local_k);
                        if (relative_target >= 0 && local_group_size > 0) {
                            int local_group_id = relative_target / local_group_size;
                            int local_start = local_group_id * local_group_size;
                            int local_end = std::min(static_cast<int>(group_disks.size()), local_start + local_group_size);

                            for (int pos = local_start; pos < local_end && static_cast<int>(chosen_sources.size()) < ec_config.local_m; ++pos) {
                                int d = group_disks[pos];
                                if (d == target_disk) continue;
                                if (all_unavailable.count(d)) continue;
                                chosen_sources.push_back(d);
                            }
                        }
                    }
                }

                // MULTI_EC global rebuild approximation:
                // When local rebuild is not possible, treat each other local stripe as a "chunk"
                // and pick one representative disk from each stripe.
                if (ec_config.type == ECType::MULTI_EC && !ec_scheme.can_local_rebuild(all_unavailable, target_disk)) {
                    chosen_sources.clear();
                    int relative_target = -1;
                    for (size_t i = 0; i < group_disks.size(); ++i) {
                        if (group_disks[i] == target_disk) {
                            relative_target = static_cast<int>(i);
                            break;
                        }
                    }
                    int local_n = ec_config.local_n;
                    int stripe_count = std::max(1, ec_config.m + ec_config.k);
                    if (relative_target >= 0 && local_n > 0) {
                        int target_stripe = relative_target / local_n;
                        for (int stripe_idx = 0; stripe_idx < stripe_count && static_cast<int>(chosen_sources.size()) < required_reads; ++stripe_idx) {
                            if (stripe_idx == target_stripe) continue;
                            int start_pos = stripe_idx * local_n;
                            int end_pos = std::min(static_cast<int>(group_disks.size()), start_pos + local_n);
                            for (int pos = start_pos; pos < end_pos; ++pos) {
                                int d = group_disks[pos];
                                if (all_unavailable.count(d)) continue;
                                chosen_sources.push_back(d);
                                break;
                            }
                        }
                    }
                }

                // Fill remaining sources from global pool.
                for (int d : available_disks) {
                    if (static_cast<int>(chosen_sources.size()) >= required_reads) break;
                    if (d == target_disk) continue;
                    if (std::find(chosen_sources.begin(), chosen_sources.end(), d) != chosen_sources.end()) continue;
                    chosen_sources.push_back(d);
                }

                if (static_cast<int>(chosen_sources.size()) < required_reads) {
                    continue;
                }

                std::vector<std::string> source_modules;
                source_modules.reserve(static_cast<size_t>(required_reads));
                for (int i = 0; i < required_reads; ++i) {
                    // Get an active (non-failed) io_module for each source disk
                    // Prefer using the same io_module as target to minimize network traffic
                    std::string src_module;
                    std::string fallback_module;
                    auto src_io_modules = disk_io_manager->get_io_modules_for_disk(chosen_sources[i]);
                    for (const auto& io_mod : src_io_modules) {
                        auto it = disconnected.failed_nodes.find(io_mod);
                        if (it == disconnected.failed_nodes.end() || !it->second) {
                            if (fallback_module.empty()) {
                                fallback_module = io_mod;
                            }
                            // Prefer using target_module if available
                            if (io_mod == target_module) {
                                src_module = io_mod;
                                break;
                            }
                        }
                    }
                    if (src_module.empty()) {
                        src_module = fallback_module.empty() ?
                            disk_io_manager->get_primary_io_module(chosen_sources[i]) : fallback_module;
                    }
                    source_modules.push_back(src_module);
                }

                // Find LCA of all source modules and target module to determine actual network path
                std::vector<std::string> all_modules = source_modules;
                all_modules.push_back(target_module);

                // Build parent map if not already built
                if (graph.parent_map_.empty()) {
                    graph.build_parent_map(start_module);
                }

                std::string lca = graph.find_lca_multiple(all_modules);

                // Check if all modules are at the same io_module (LCA == io_module level)
                // This means traffic is intra-io_module and doesn't need cross-switch bandwidth
                std::set<std::string> unique_modules(all_modules.begin(), all_modules.end());
                bool is_intra_io_module = (unique_modules.size() == 1);

                // Also check if LCA is at io_module level (all nodes share common io_module-level ancestor)
                if (!is_intra_io_module && !lca.empty()) {
                    auto lca_depth_it = graph.depth_map_.find(lca);
                    auto target_depth_it = graph.depth_map_.find(target_module);
                    if (lca_depth_it != graph.depth_map_.end() && target_depth_it != graph.depth_map_.end()) {
                        // If LCA is at io_module level (same depth as io_modules), traffic stays within that level
                        if (lca_depth_it->second == target_depth_it->second) {
                            is_intra_io_module = true;
                        }
                    }
                }

                // For replication (or MULTI_EC with local replication), no EC encoding needed - skip ec_encoding_speed limit
                bool skip_rebuild_encoding = (ec_config.type == ECType::REPLICATION) ||
                    (ec_config.type == ECType::MULTI_EC && ec_config.local_type == LocalType::REPLICATION);
                double rebuild_speed = skip_rebuild_encoding ?
                    std::numeric_limits<double>::max() : params.ec_encoding_speed;
                std::map<std::pair<std::string, std::string>, double> flow_map;
                std::map<std::pair<std::string, std::string>, double> down_flow_map;
                double max_up = 0.0;
                double max_down = 0.0;

                if (is_intra_io_module) {
                    // Intra-io_module: no cross-switch traffic needed
                    // Speed limited only by disk I/O (and encoding for EC)
                    rebuild_speed = std::min(rebuild_speed, static_cast<double>(required_reads) * params.disk_read_bw);
                    rebuild_speed = std::min(rebuild_speed, params.disk_write_bw);
                } else {
                    // Cross-io_module traffic: compute max flow up to LCA level
                    // Use LCA as the meeting point instead of start_module (switch)
                    std::string flow_meeting_point = lca.empty() ? start_module : lca;

                    // Upstream (read) path: sources -> LCA
                    max_up = max_flow_from_sources_to_module(source_modules, flow_meeting_point, FlowDirection::UPSTREAM, flow_map);
                    if (max_up <= 0.0) {
                        continue;
                    }

                    // Downstream (write) path: LCA -> target_io_module
                    max_down = max_flow_from_module_to_target(flow_meeting_point, target_module, FlowDirection::DOWNSTREAM, down_flow_map);
                    if (max_down <= 0.0) {
                        continue;
                    }

                    // Rebuild speed limited by encoding (for EC) + upstream read path + downstream write path
                    rebuild_speed = std::min(rebuild_speed, max_up * params.rebuild_bw_ratio);
                    rebuild_speed = std::min(rebuild_speed, max_down * params.rebuild_bw_ratio);
                    // Also limited by per-disk IO (approx)
                    rebuild_speed = std::min(rebuild_speed, static_cast<double>(required_reads) * params.disk_read_bw);
                    rebuild_speed = std::min(rebuild_speed, params.disk_write_bw);
                }

                entry.rebuild_speed_by_disk[target_disk] = rebuild_speed;

                // Record a representative rebuild bandwidth for compatibility with existing outputs.
                auto it_bw = entry.rebuild_bandwidth.find(effective_failure_count);
                if (it_bw == entry.rebuild_bandwidth.end()) {
                    entry.rebuild_bandwidth[effective_failure_count] = rebuild_speed;
                } else {
                    it_bw->second = std::min(it_bw->second, rebuild_speed);
                }

                // Reserve impact-scaled rebuild traffic (declustered parity relief).
                // Only reserve if there's cross-io_module traffic
                if (!is_intra_io_module) {
                    const double reserved_write_speed = rebuild_speed * impact_ratio;
                    // Upstream needs required_reads * rebuild_speed (reading from m sources)
                    const double reserved_read_speed = reserved_write_speed * required_reads;
                    if (reserved_write_speed > 0 && max_up > 0 && max_down > 0) {
                        reserve_scaled_flow(flow_map, reserved_read_speed / max_up, FlowDirection::UPSTREAM);
                        reserve_scaled_flow(down_flow_map, reserved_write_speed / max_down, FlowDirection::DOWNSTREAM);
                    }
                }
            }
        }

        if (state == DISK_STATE_DATA_LOSS) {
            // Calculate data loss ratio for this group
            double data_loss = ec_scheme.get_data_loss_ratio(all_unavailable);
            total_data_loss += data_loss / total_groups;
            groups_with_data_loss++;  // Count this group as having data loss
            entry.data_loss_group_indices.push_back(group_index);  // Record which group has data loss

            // Even during data loss, disks need rebuild (slower, fixed speed).
            // This allows recovery from data loss state when disks are replaced.
            if (params.data_loss_rebuild_bw > 0 && disk_io_manager != nullptr && disks != nullptr) {
                // Find targets: disks that need rebuild and are not disconnected
                std::vector<int> targets;
                for (int d : group_disks) {
                    auto dit = disks->find(d);
                    if (dit == disks->end()) continue;
                    const DiskInfo& info = dit->second;
                    if (!info.needs_rebuild()) continue;
                    if (info.is_disconnected) continue;
                    targets.push_back(d);
                }

                // Use fixed data_loss_rebuild_bw for all targets
                const double fixed_rebuild_speed = params.data_loss_rebuild_bw;

                for (int target_disk : targets) {
                    entry.rebuild_speed_by_disk[target_disk] = fixed_rebuild_speed;

                    // Reserve write bandwidth (root -> io_module -> disk direction)
                    // Get an active io_module for this disk
                    std::string target_module;
                    auto target_io_modules = disk_io_manager->get_io_modules_for_disk(target_disk);
                    for (const auto& io_mod : target_io_modules) {
                        auto it = disconnected.failed_nodes.find(io_mod);
                        if (it == disconnected.failed_nodes.end() || !it->second) {
                            target_module = io_mod;
                            break;
                        }
                    }
                    if (target_module.empty()) continue;

                    // Reserve downstream bandwidth for data loss rebuild
                    std::map<std::pair<std::string, std::string>, double> down_flow_map;
                    double max_down = max_flow_from_module_to_target(start_module, target_module,
                                                                      FlowDirection::DOWNSTREAM, down_flow_map);
                    if (max_down > 0) {
                        reserve_scaled_flow(down_flow_map, fixed_rebuild_speed / max_down, FlowDirection::DOWNSTREAM);
                    }
                }
            }
        }
    }

    // Host IO on the graph with reservations applied.
    entry.read_bandwidth = std::max(0.0, graph.calculate_max_flow(end_module, start_module, FlowDirection::UPSTREAM));

    // Apply encoding overhead for write bandwidth using LCA-based cross-traffic calculation
    // When IO modules perform encoding, cross-traffic flows through LCA paths between IO module pairs
    if (params.encoding_config.entity == EncodingEntity::IO_MODULE &&
        disk_io_manager != nullptr && ec_config.k > 0) {

        auto io_module_set = disk_io_manager->get_all_io_modules();
        std::vector<std::string> all_io_modules(io_module_set.begin(), io_module_set.end());
        int N = static_cast<int>(all_io_modules.size());

        if (N > 1) {
            // Calculate cross-traffic factor per edge
            // Each IO module pair exchanges data through their LCA path
            // Cross-traffic per pair = 2 / N^2 (bidirectional, per unit host write)
            std::map<std::pair<std::string, std::string>, double> edge_cross_factor;

            for (size_t i = 0; i < all_io_modules.size(); ++i) {
                for (size_t j = i + 1; j < all_io_modules.size(); ++j) {
                    const std::string& io_i = all_io_modules[i];
                    const std::string& io_j = all_io_modules[j];

                    // Find LCA of this pair
                    std::string lca = graph.find_lca(io_i, io_j);
                    if (lca.empty()) continue;

                    // Get paths from each IO module to root
                    auto path_i = graph.get_path_to_root(io_i);
                    auto path_j = graph.get_path_to_root(io_j);

                    // Cross-traffic amount per pair per unit host write
                    // Each module holds 1/N of data, sends (1/N) to each of (N-1) others
                    // Per pair per direction = 1/N^2, bidirectional = 2/N^2
                    double cross_per_pair = 2.0 / (N * N);

                    // Add factor to edges on path io_i -> lca
                    for (size_t k = 0; k + 1 < path_i.size(); ++k) {
                        if (path_i[k] == lca) break;
                        // Edge: path_i[k] <-> path_i[k+1] (both directions for cross-traffic)
                        edge_cross_factor[{path_i[k+1], path_i[k]}] += cross_per_pair;
                        edge_cross_factor[{path_i[k], path_i[k+1]}] += cross_per_pair;
                    }
                    // Add factor to edges on path io_j -> lca
                    for (size_t k = 0; k + 1 < path_j.size(); ++k) {
                        if (path_j[k] == lca) break;
                        edge_cross_factor[{path_j[k+1], path_j[k]}] += cross_per_pair;
                        edge_cross_factor[{path_j[k], path_j[k+1]}] += cross_per_pair;
                    }
                }
            }

            // Debug: log cross-traffic factors
            static int cross_log_count = 0;
            if (cross_log_count < 1) {
                std::cerr << "[CROSS-TRAFFIC] N=" << N << " pairs=" << edge_cross_factor.size() << std::endl;
                for (const auto& [edge, factor] : edge_cross_factor) {
                    std::cerr << "  " << edge.first << " -> " << edge.second << " factor=" << factor << std::endl;
                }
                cross_log_count++;
            }

            // Calculate write bandwidth with cross-traffic overhead
            // For edges with cross-traffic, effective capacity = capacity / (1 + factor)
            // We clone the graph and adjust capacities
            GraphStructure write_graph = graph.clone();
            for (const auto& [edge, factor] : edge_cross_factor) {
                if (factor > 0) {
                    // Reduce both upstream and downstream capacity
                    double current_up = write_graph.get_edge_capacity(edge.first, edge.second, FlowDirection::UPSTREAM);
                    double current_down = write_graph.get_edge_capacity(edge.first, edge.second, FlowDirection::DOWNSTREAM);
                    if (current_up > 0) {
                        write_graph.reserve_bandwidth(edge.first, edge.second,
                            current_up * factor / (1.0 + factor), FlowDirection::UPSTREAM);
                    }
                    if (current_down > 0) {
                        write_graph.reserve_bandwidth(edge.first, edge.second,
                            current_down * factor / (1.0 + factor), FlowDirection::DOWNSTREAM);
                    }
                }
            }
            entry.write_bandwidth = std::max(0.0, write_graph.calculate_max_flow(start_module, end_module, FlowDirection::DOWNSTREAM));
        } else {
            // Single IO module - no cross-traffic
            entry.write_bandwidth = std::max(0.0, graph.calculate_max_flow(start_module, end_module, FlowDirection::DOWNSTREAM));
        }
    } else {
        entry.write_bandwidth = std::max(0.0, graph.calculate_max_flow(start_module, end_module, FlowDirection::DOWNSTREAM));
    }

    entry.data_loss_ratio = total_data_loss;
    entry.groups_with_data_loss = groups_with_data_loss;

    // Debug: log bandwidth reservation impact per edge
    static thread_local int bw_log_count = 0;
    if (bw_log_count < 2 && !entry.rebuild_bandwidth.empty()) {
        std::string edges_str;
        for (const auto& [key, val] : graph.reserved_upstream_) {
            edges_str += key.first + "->" + key.second + ":" + std::to_string(val / 1e9) + "G, ";
        }
        Logger::getInstance().info(
            "BW_DEBUG: read_bw=" + std::to_string(entry.read_bandwidth / 1e9) + "GB/s" +
            " edges=[" + edges_str + "]");
        bw_log_count++;
    }

    // Availability: 1.0 if no EC group has more than k failures (can still read)
    // Data loss doesn't mean unavailable - it means some data is permanently lost
    // Availability drops to 0 only when we can't serve any reads
    // For now, availability = 1 - (fraction of groups with data loss)
    int groups_with_loss = 0;
    int total_failed_disks = 0;
    int total_disconnected_disks = 0;

    for (int group_index = 0; group_index < num_ec_groups; ++group_index) {
        std::vector<int> group_disks = ec_scheme.get_group_disks(group_index);
        if (group_disks.empty()) {
            int start_disk_idx = ec_scheme.get_group_start_disk(group_index);
            for (int i = 0; i < group_size; ++i) {
                group_disks.push_back(start_disk_idx + i);
            }
        }
        std::set<int> unavailable_in_group;
        int failed_in_group = 0;
        int disconnected_in_group = 0;

        auto fi = failure_info_per_ec_group.find(group_index);
        if (fi != failure_info_per_ec_group.end()) {
            unavailable_in_group = fi->second.failed_disk_indices;
            failed_in_group = static_cast<int>(fi->second.failed_disk_indices.size());
        }
        for (int d : group_disks) {
            if (disconnected.is_disk_disconnected(d)) {
                unavailable_in_group.insert(d);
                disconnected_in_group++;
            }
        }

        total_failed_disks += failed_in_group;
        total_disconnected_disks += disconnected_in_group;

        // Use ec_scheme.is_data_loss() to properly handle LRC/Multi-EC
        if (ec_scheme.is_data_loss(unavailable_in_group)) {
            groups_with_loss++;
        }
    }

    // Debug logging for unavailability
    if (groups_with_loss > 0) {
        // Count failed nodes
        int failed_io_modules = 0;
        int failed_switches = 0;
        int failed_enclosures = 0;
        for (const auto& [node, failed] : disconnected.failed_nodes) {
            if (failed) {
                if (node.find("io_module") != std::string::npos) {
                    failed_io_modules++;
                } else if (node.find("switch") != std::string::npos) {
                    failed_switches++;
                } else if (node.find("Enclosure") != std::string::npos) {
                    failed_enclosures++;
                }
            }
        }

        Logger::getInstance().info(
            "UNAVAIL: groups_with_loss=" + std::to_string(groups_with_loss) +
            "/" + std::to_string(num_ec_groups) +
            ", failed_disks=" + std::to_string(total_failed_disks) +
            ", disconnected_disks=" + std::to_string(total_disconnected_disks) +
            ", failed_io_modules=" + std::to_string(failed_io_modules) +
            ", failed_switches=" + std::to_string(failed_switches) +
            ", failed_enclosures=" + std::to_string(failed_enclosures));
    }

    // Availability = fraction of groups that are still readable (not data_loss)
    int available_groups = num_ec_groups - groups_with_loss;
    entry.availability_ratio = static_cast<double>(available_groups) / num_ec_groups;

    if (max_read_performance_without_any_failure > 0) {
        entry.performance_ratio_read = std::clamp(entry.read_bandwidth / max_read_performance_without_any_failure, 0.0, 1.0);
    } else {
        entry.performance_ratio_read = 1.0;
    }

    if (max_write_performance_without_any_failure > 0) {
        entry.performance_ratio_write = std::clamp(entry.write_bandwidth / max_write_performance_without_any_failure, 0.0, 1.0);
    } else {
        entry.performance_ratio_write = 1.0;
    }

    // Data loss인 그룹은 해당 비율만큼 performance도 감소 (perf_availability <= availability 보장)
    entry.performance_ratio_read = std::min(entry.performance_ratio_read, entry.availability_ratio);
    entry.performance_ratio_write = std::min(entry.performance_ratio_write, entry.availability_ratio);

    entry.performance_ratio_min = std::min(entry.performance_ratio_read, entry.performance_ratio_write);

    flows_and_speed_table[key] = entry;
}

// Update all disk states
void update_all_disk_states(
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    std::map<int, DiskInfo>& disks,
    const ErasureCodingScheme& ec_scheme,
    const DisconnectedStatus& disconnected,
    double current_time) {

    for (auto& [disk_index, disk_info] : disks) {
        if (!disk_info.is_failed && !disk_info.is_disconnected) {
            continue;
        }

        int group_index = ec_scheme.get_disk_group(disk_index);
        auto it = failure_info_per_ec_group.find(group_index);
        if (it == failure_info_per_ec_group.end()) {
            continue;
        }

        const ECGroupFailureInfo& failure_info = it->second;
        std::string state = judge_state_from_failure_info(failure_info, ec_scheme);

        if (state == DISK_STATE_REBUILDING) {
            disk_info.state = DiskState::REBUILDING;
        } else if (state == DISK_STATE_DATA_LOSS) {
            disk_info.state = DiskState::FAILED;
        }

        // Handle reconnection after disconnect
        if (!disconnected.is_disk_disconnected(disk_index) &&
            disk_info.disconnect_timestamp != 0 && !disk_info.is_failed) {
            // Calculate how much data needs to be rebuilt based on disconnect time
            double disconnect_duration = current_time - disk_info.disconnect_timestamp;
            // Proportional rebuild: disconnect time * write rate approximation
            disk_info.remaining_capacity_to_rebuild = disconnect_duration * disk_info.write_bandwidth * 3600.0;
            disk_info.disconnect_timestamp = 0;
            disk_info.is_disconnected = false;
        }
    }
}

// Update single disk state
void update_disk_state(
    int disk_index,
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    std::map<int, DiskInfo>& disks,
    double capacity,
    const std::string& event_type,
    double event_time,
    double replace_time,
    const ErasureCodingScheme& ec_scheme,
    const DisconnectedStatus& disconnected) {

    DiskInfo& disk_info = disks[disk_index];

    if (event_type == "fail") {
        disk_info.is_failed = true;
        disk_info.fail_timestamp = event_time;
        disk_info.rebuild_start_timestamp = event_time;

        if (disk_info.disconnect_timestamp != 0) {
            // Already disconnected, minimal rebuild needed
            disk_info.remaining_capacity_to_rebuild = 0;
            disk_info.remaining_replace_time = replace_time;
        } else {
            // Full disk failure: need to replace and rebuild entire disk
            disk_info.remaining_capacity_to_rebuild = capacity;
            std::exponential_distribution<double> exp_dist(1.0 / replace_time);
            disk_info.remaining_replace_time = exp_dist(rng);
        }
        disk_info.rebuild_speed = 0;
    } else if (event_type == "repair") {
        disk_info.is_failed = false;
        disk_info.is_disconnected = false;
        disk_info.state = DiskState::NORMAL;
        disk_info.remaining_capacity_to_rebuild = 0;
        disk_info.remaining_replace_time = 0;
        disk_info.disconnect_timestamp = 0;
        disk_info.rebuild_start_timestamp = 0;
    }

    // Update state based on EC group
    int group_index = ec_scheme.get_disk_group(disk_index);
    auto it = failure_info_per_ec_group.find(group_index);
    if (it != failure_info_per_ec_group.end()) {
        std::string state_str = judge_state_from_failure_info(it->second, ec_scheme);
        if (state_str == DISK_STATE_REBUILDING) {
            disk_info.state = DiskState::REBUILDING;
        } else if (state_str == DISK_STATE_DATA_LOSS) {
            disk_info.state = DiskState::FAILED;
        }
    }
}

// Update failure events for disks (when io_module fails, disks become disconnected)
void update_failure_event_for_disks(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
    double current_time,
    std::map<int, DiskInfo>& disks) {

    std::vector<Event> updated_failed_events;

    while (!failed_events.empty()) {
        Event popped_event = failed_events.top();
        failed_events.pop();

        if (!is_disk_node(popped_event.event_node)) {
            updated_failed_events.push_back(popped_event);
            continue;
        }

        if (popped_event.event_type != "fail") {
            updated_failed_events.push_back(popped_event);
            continue;
        }

        int disk_index = get_disk_index(popped_event.event_node);
        DiskInfo& disk_info = disks[disk_index];

        if (disk_info.is_failed) {
            updated_failed_events.push_back(popped_event);
            continue;
        }

        // Mark as disconnected and schedule immediate failure
        double failure_time = current_time + 1.0 / BIG_NUMBER;
        disk_info.disconnect_timestamp = current_time;
        disk_info.is_disconnected = true;

        Event new_event{failure_time, "fail", popped_event.event_node,
                       popped_event.prev_time, popped_event.start_time, true};
        updated_failed_events.push_back(new_event);
    }

    for (const auto& event : updated_failed_events) {
        failed_events.push(event);
    }
}

// Update repair events for disks
void update_repair_event_for_disks(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events,
    double current_time,
    std::map<int, DiskInfo>& disks,
    const FlowsAndSpeedEntry& flows_and_speed_entry,
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    const ErasureCodingScheme& ec_scheme) {

    std::vector<Event> updated_repair_events;

    while (!repair_events.empty()) {
        Event popped_event = repair_events.top();
        repair_events.pop();

        if (!is_disk_node(popped_event.event_node)) {
            updated_repair_events.push_back(popped_event);
            continue;
        }

        int disk_index = get_disk_index(popped_event.event_node);
        DiskInfo& disk_info = disks[disk_index];
        double remaining_consumed_time = current_time - popped_event.prev_time;

        // Process replace time first
        if (disk_info.remaining_replace_time > 0) {
            double remaining_replace = disk_info.remaining_replace_time - remaining_consumed_time;
            if (remaining_replace <= 0) {
                disk_info.remaining_replace_time = 0;
                remaining_consumed_time = -remaining_replace;
            } else {
                disk_info.remaining_replace_time = remaining_replace;
                remaining_consumed_time = 0;
            }
        }

        // Update remaining capacity with previous rebuild speed
        if (remaining_consumed_time > 0 && disk_info.rebuild_speed > 0) {
            double rebuild_speed_per_hour = disk_info.rebuild_speed * 3600.0;
            double rebuilt = remaining_consumed_time * rebuild_speed_per_hour;
            disk_info.remaining_capacity_to_rebuild -= rebuilt;
            if (disk_info.remaining_capacity_to_rebuild < 0) {
                disk_info.remaining_capacity_to_rebuild = 0;
            }
        }

        // Get new rebuild speed
        double rebuild_speed = 0;
        int group_index = ec_scheme.get_disk_group(disk_index);
        auto per_disk_it = flows_and_speed_entry.rebuild_speed_by_disk.find(disk_index);
        if (per_disk_it != flows_and_speed_entry.rebuild_speed_by_disk.end()) {
            rebuild_speed = per_disk_it->second;
        } else {
            auto group_it = failure_info_per_ec_group.find(group_index);
            if (group_it != failure_info_per_ec_group.end()) {
                int failure_count = static_cast<int>(group_it->second.failed_disk_indices.size());
                auto bw_it = flows_and_speed_entry.rebuild_bandwidth.find(failure_count);
                if (bw_it != flows_and_speed_entry.rebuild_bandwidth.end()) {
                    rebuild_speed = bw_it->second;
                }
            }
        }

        disk_info.rebuild_speed = rebuild_speed;

        // Calculate new repair time
        if (rebuild_speed <= 0) {
            rebuild_speed = 1.0 / BIG_NUMBER;
        }

        double rebuild_speed_per_hour = rebuild_speed * 3600.0;
        double remaining_time = disk_info.remaining_replace_time +
                               disk_info.remaining_capacity_to_rebuild / rebuild_speed_per_hour;

        // Preserve the event type:
        // - "repair" means physical disk replacement + full rebuild completed
        // - "rebuild_complete" means resync after reconnect completed (should NOT reschedule failures)
        Event new_event{current_time + remaining_time, popped_event.event_type, popped_event.event_node,
                       current_time, popped_event.start_time, popped_event.software_repair_only};
        updated_repair_events.push_back(new_event);
    }

    for (const auto& event : updated_repair_events) {
        repair_events.push(event);
    }
}

// Calculate combinations (n choose k)
static double combinations_count(int n, int k) {
    if (n < k) return 0.0;
    if (k == 0 || k == n) return 1.0;

    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}

// Monte Carlo simulation main function
void monte_carlo_simulation(
    std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int num_simulations,
    const nlohmann::json& options
) {
    Logger::getInstance().info("Monte Carlo simulation starting with " + std::to_string(num_simulations) + " iterations");

    int nprocs = params_and_results.at("nprocs");
    int base = (nprocs > 0) ? (num_simulations / nprocs) : 0;
    int rem = (nprocs > 0) ? (num_simulations % nprocs) : 0;

    Logger::getInstance().info("Using " + std::to_string(nprocs) + " parallel processes");

    std::vector<std::future<SimulationResult>> futures;

    // Launch parallel simulations
    for (int i = 0; i < nprocs; ++i) {
        int num_iterations = base + ((i < rem) ? 1 : 0);
        if (num_iterations <= 0) continue;
        futures.push_back(std::async(std::launch::async, simulation_per_core,
                                     i, std::ref(params_and_results),
                                     std::ref(graph_structure_origin),
                                     num_iterations, std::ref(options)));
    }

    Logger::getInstance().info("All simulation threads launched");

    // Collect results
    int total_runs = 0;
    double total_up_time = 0.0;
    double total_simulation_time = 0.0;
    double total_time_for_rebuilding = 0.0;
    int total_count_for_rebuilding = 0;
    double total_mttdl = 0.0;
    int total_data_loss_events = 0;
    int total_disk_failure_events = 0;
    double total_data_loss_ratio = 0.0;
    double total_read_bw = 0.0;
    double total_write_bw = 0.0;
    double total_perf_up_time_read = 0.0;
    double total_perf_up_time_write = 0.0;
    double total_perf_up_time_min = 0.0;
    // Rebuild-specific metrics
    double total_rebuild_speed_time = 0.0;
    double total_time_in_rebuild = 0.0;
    double total_host_read_bw_during_rebuild = 0.0;
    double total_host_write_bw_during_rebuild = 0.0;
    // Year-based durability metrics
    int total_group_years = 0;
    int total_group_years_with_data_loss = 0;
    // Baseline performance (same for all threads, capture from first result)
    double max_read_performance = 0.0;
    double max_write_performance = 0.0;
    bool baseline_captured = false;

    for (auto& future : futures) {
        SimulationResult result = future.get();
        // Capture baseline from first result
        if (!baseline_captured) {
            max_read_performance = result.max_read_performance;
            max_write_performance = result.max_write_performance;
            baseline_captured = true;
        }
        total_runs += result.runs;
        total_up_time += result.up_time;
        total_simulation_time += result.simulation_time;
        total_time_for_rebuilding += result.time_for_rebuilding;
        total_count_for_rebuilding += result.count_for_rebuilding;
        total_mttdl += result.mttdl;
        total_data_loss_events += result.data_loss_events;
        total_disk_failure_events += result.disk_failure_events;
        total_data_loss_ratio += result.total_data_loss_ratio;
        total_read_bw += result.average_read_bandwidth;
        total_write_bw += result.average_write_bandwidth;
        total_perf_up_time_read += result.perf_up_time_read;
        total_perf_up_time_write += result.perf_up_time_write;
        total_perf_up_time_min += result.perf_up_time_min;
        // Rebuild-specific metrics
        total_rebuild_speed_time += result.total_rebuild_speed_time;
        total_time_in_rebuild += result.total_time_in_rebuild;
        total_host_read_bw_during_rebuild += result.total_host_read_bw_during_rebuild;
        total_host_write_bw_during_rebuild += result.total_host_write_bw_during_rebuild;
        // Year-based durability metrics
        total_group_years += result.total_group_years;
        total_group_years_with_data_loss += result.group_years_with_data_loss;
    }

    Logger::getInstance().info("All simulation threads completed, aggregating results");

    // Calculate averages
    if (total_runs <= 0) {
        total_runs = num_simulations;
    }
    double avg_up_time = total_up_time / total_runs;
    double avg_simulation_time = total_simulation_time / total_runs;
    double availability = avg_up_time / avg_simulation_time;

    // Get EC group count for per-group metrics
    int total_disks = options.value("total_disks", 48);
    int ec_n = options.value("n", 16);
    int num_ec_groups = (ec_n > 0) ? (total_disks / ec_n) : 1;

    // Calculate system-wide MTTDL (total simulation time / data loss events)
    // This gives the correct MTTDL = 1 / data_loss_events_per_hour
    double system_mttdl = (total_data_loss_events > 0) ? (total_simulation_time / total_data_loss_events) : 0.0;

    // Per-EC-group MTTDL: If groups fail independently, per-group MTTDL = system_MTTDL * num_groups
    double mttdl_per_ec_group = system_mttdl * num_ec_groups;

    // Calculate durability and PDL
    int simulation_years = options.value("simulation_years", 10);
    double hours_per_year = 365.0 * 24.0;
    double total_simulation_years = (avg_simulation_time / hours_per_year) * total_runs;

    // Data loss events per year (for informational purposes)
    double data_loss_events_per_year = (total_simulation_years > 0) ?
        (static_cast<double>(total_data_loss_events) / total_simulation_years) : 0.0;

    // PDL (Probability of Data Loss ratio) - time-weighted average fraction of data lost
    double avg_data_loss_ratio = total_data_loss_ratio / total_runs;

    // Durability = fraction of group-years without data loss
    // Each EC group × year is an independent sample
    // Durability = (group-years without data loss) / (total group-years)
    double durability = (total_group_years > 0) ?
        (static_cast<double>(total_group_years - total_group_years_with_data_loss) / total_group_years) : 1.0;

    // Store results
    params_and_results["availability"] = availability;
    params_and_results["avail_nines"] = Utils::get_nines(availability);
    params_and_results["simulation_time"] = avg_simulation_time;
    params_and_results["up_time"] = avg_up_time;
    params_and_results["mttdl_per_ec_group"] = mttdl_per_ec_group;
    params_and_results["mttdl_system"] = system_mttdl;
    params_and_results["num_ec_groups"] = num_ec_groups;
    params_and_results["data_loss_events"] = total_data_loss_events;
    params_and_results["data_loss_events_per_year"] = data_loss_events_per_year;
    params_and_results["disk_failure_events"] = total_disk_failure_events;
    params_and_results["pdl_ratio"] = avg_data_loss_ratio;
    params_and_results["durability"] = durability;
    params_and_results["durability_nines"] = Utils::get_nines(durability);
    params_and_results["total_group_years"] = total_group_years;
    params_and_results["group_years_with_data_loss"] = total_group_years_with_data_loss;

    if (total_count_for_rebuilding > 0) {
        params_and_results["avg_rebuilding_time"] = total_time_for_rebuilding / total_count_for_rebuilding;
    } else {
        params_and_results["avg_rebuilding_time"] = 0.0;
    }

    params_and_results["avg_read_bandwidth"] = total_read_bw / total_runs;
    params_and_results["avg_write_bandwidth"] = total_write_bw / total_runs;
    params_and_results["max_read_performance"] = max_read_performance;
    params_and_results["max_write_performance"] = max_write_performance;

    // Rebuild-specific metrics
    if (total_time_in_rebuild > 0) {
        double avg_rebuild_speed = total_rebuild_speed_time / total_time_in_rebuild;
        double avg_host_read_bw_during_rebuild = total_host_read_bw_during_rebuild / total_time_in_rebuild;
        double avg_host_write_bw_during_rebuild = total_host_write_bw_during_rebuild / total_time_in_rebuild;
        params_and_results["avg_rebuild_speed"] = avg_rebuild_speed;
        params_and_results["avg_host_read_bw_during_rebuild"] = avg_host_read_bw_during_rebuild;
        params_and_results["avg_host_write_bw_during_rebuild"] = avg_host_write_bw_during_rebuild;
        params_and_results["total_time_in_rebuild_hours"] = total_time_in_rebuild;
    }

    // Performance availability (if target_performance_ratio specified)
    double target_perf_ratio = options.value("target_performance_ratio", 0.0);
    if (target_perf_ratio > 0) {
        double avg_perf_up_time_read = total_perf_up_time_read / total_runs;
        double avg_perf_up_time_write = total_perf_up_time_write / total_runs;
        double avg_perf_up_time_min = total_perf_up_time_min / total_runs;

        double perf_availability_read = avg_perf_up_time_read / avg_simulation_time;
        double perf_availability_write = avg_perf_up_time_write / avg_simulation_time;
        double perf_availability_min = avg_perf_up_time_min / avg_simulation_time;

        params_and_results["perf_availability_read"] = perf_availability_read;
        params_and_results["perf_avail_nines_read"] = Utils::get_nines(perf_availability_read);
        params_and_results["perf_availability_write"] = perf_availability_write;
        params_and_results["perf_avail_nines_write"] = Utils::get_nines(perf_availability_write);
        params_and_results["perf_availability_min"] = perf_availability_min;
        params_and_results["perf_avail_nines_min"] = Utils::get_nines(perf_availability_min);
    }

    Logger::getInstance().info("Simulation completed");
    Logger::getInstance().info("Availability: " + std::to_string(availability) +
                              " (" + std::to_string(Utils::get_nines(availability)) + " nines)");
    Logger::getInstance().info("Durability: " + std::to_string(durability) +
                              " (" + std::to_string(Utils::get_nines(durability)) + " nines)");
}
