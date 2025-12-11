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
bool is_data_loss(const ECGroupFailureInfo& failure_info, const ECConfig& ec_config) {
    int unavailable = failure_info.get_unavailable_count();
    return unavailable > ec_config.k;
}

// Judge disk state from failure information
std::string judge_state_from_failure_info(const ECGroupFailureInfo& failure_info,
                                          const ECConfig& ec_config) {
    int unavailable = failure_info.get_unavailable_count();

    if (unavailable == 0) {
        return DISK_STATE_NORMAL;
    }

    if (unavailable > ec_config.k) {
        return DISK_STATE_DATA_LOSS;
    }

    return DISK_STATE_REBUILDING;
}

// Random number generator (thread-local)
thread_local std::mt19937 rng(std::random_device{}());

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
    if (is_disk_node(event_node)) {
        int disk_index = get_disk_index(event_node);
        int group_index = ec_scheme.get_disk_group(disk_index);

        if (event_type == "fail") {
            failure_info_per_ec_group[group_index].disk_failure_count++;
            failure_info_per_ec_group[group_index].failed_disk_indices.insert(disk_index);
        } else if (event_type == "repair") {
            failure_info_per_ec_group[group_index].disk_failure_count--;
            failure_info_per_ec_group[group_index].failed_disk_indices.erase(disk_index);
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
    std::map<FailureStateKey, FlowsAndSpeedEntry>& flows_and_speed_table,
    double max_read_performance_without_any_failure,
    const DisconnectedStatus& disconnected,
    const FailureStateKey& key,
    const DiskIOModuleManager* disk_io_manager) {

    if (flows_and_speed_table.find(key) != flows_and_speed_table.end()) {
        return;
    }

    FlowsAndSpeedEntry entry;
    const ECConfig& ec_config = ec_scheme.get_config();

    int group_size = ec_scheme.get_group_size();

    // Process each EC group
    double total_data_loss = 0.0;
    int total_groups = 0;

    // Get total number of EC groups
    int num_ec_groups = ec_scheme.get_total_groups(params.total_disks);

    for (int group_index = 0; group_index < num_ec_groups; ++group_index) {
        total_groups++;
        int start_disk_idx = ec_scheme.get_group_start_disk(group_index);

        // Get failure info for this group (may not exist if no failures)
        ECGroupFailureInfo failure_info;
        auto it = failure_info_per_ec_group.find(group_index);
        if (it != failure_info_per_ec_group.end()) {
            failure_info = it->second;
        }

        // Calculate effective failure count (failed + disconnected)
        std::set<int> all_unavailable = failure_info.failed_disk_indices;
        for (int i = start_disk_idx; i < start_disk_idx + group_size; ++i) {
            if (disconnected.is_disk_disconnected(i)) {
                all_unavailable.insert(i);
            }
        }
        int effective_failure_count = static_cast<int>(all_unavailable.size());

        // Determine state based on effective unavailable count
        std::string state;
        if (effective_failure_count == 0) {
            state = DISK_STATE_NORMAL;
        } else if (effective_failure_count > ec_config.k) {
            state = DISK_STATE_DATA_LOSS;
        } else {
            state = DISK_STATE_REBUILDING;
        }

        if (state == DISK_STATE_REBUILDING) {
            int degraded_disks = group_size - effective_failure_count;
            if (degraded_disks < 0) degraded_disks = 0;

            // Calculate rebuild bandwidth
            std::vector<double> bws = {params.disk_read_bw, params.disk_write_bw};

            // Check if EC group crosses io_module boundary
            if (disk_io_manager != nullptr && !disk_io_manager->is_legacy_mode()) {
                bool crosses_boundary = disk_io_manager->does_ec_group_cross_io_module(start_disk_idx, group_size);

                if (crosses_boundary) {
                    std::set<std::string> ec_group_io_modules =
                        disk_io_manager->get_io_modules_for_ec_group(start_disk_idx, group_size);
                    std::vector<std::string> source_modules(ec_group_io_modules.begin(), ec_group_io_modules.end());
                    std::string target_module = *ec_group_io_modules.begin();

                    double lca_max_flow = hardware_graph_copy.calculate_rebuild_max_flow_via_io_modules(
                        source_modules, target_module);

                    double lca_bw_per_group = lca_max_flow / std::max(1, static_cast<int>(failure_info_per_ec_group.size()));
                    bws.push_back(lca_bw_per_group);
                }
            }

            double rebuild_bw = calculate_bottleneck_speed(
                ec_config.m, effective_failure_count, bws,
                params.rebuild_bw_ratio, params.ec_encoding_speed);

            entry.rebuild_bandwidth[effective_failure_count] = rebuild_bw;
        }
        else if (state == DISK_STATE_DATA_LOSS) {
            // Calculate data loss ratio for this group
            double data_loss = ec_scheme.get_data_loss_ratio(all_unavailable);
            total_data_loss += data_loss / total_groups;
        }
    }

    // Calculate available read/write bandwidth
    // Count available disks
    int available_disks = 0;
    for (int i = 0; i < params.total_disks; ++i) {
        if (!disconnected.is_disk_disconnected(i)) {
            int group_idx = ec_scheme.get_disk_group(i);
            auto it = failure_info_per_ec_group.find(group_idx);
            if (it == failure_info_per_ec_group.end() ||
                it->second.failed_disk_indices.find(i) == it->second.failed_disk_indices.end()) {
                available_disks++;
            }
        }
    }

    entry.read_bandwidth = params.disk_read_bw * available_disks;
    entry.write_bandwidth = params.disk_write_bw * available_disks;
    entry.data_loss_ratio = total_data_loss;

    // Availability: 1.0 if no EC group has more than k failures (can still read)
    // Data loss doesn't mean unavailable - it means some data is permanently lost
    // Availability drops to 0 only when we can't serve any reads
    // For now, availability = 1 - (fraction of groups with data loss)
    int groups_with_loss = 0;
    for (int group_index = 0; group_index < num_ec_groups; ++group_index) {
        int start_disk_idx = ec_scheme.get_group_start_disk(group_index);
        std::set<int> unavailable_in_group;
        auto fi = failure_info_per_ec_group.find(group_index);
        if (fi != failure_info_per_ec_group.end()) {
            unavailable_in_group = fi->second.failed_disk_indices;
        }
        for (int i = start_disk_idx; i < start_disk_idx + group_size; ++i) {
            if (disconnected.is_disk_disconnected(i)) {
                unavailable_in_group.insert(i);
            }
        }
        if (static_cast<int>(unavailable_in_group.size()) > ec_config.k) {
            groups_with_loss++;
        }
    }
    // Availability = fraction of groups that are still readable
    entry.availability_ratio = 1.0 - static_cast<double>(groups_with_loss) / num_ec_groups;

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
        std::string state = judge_state_from_failure_info(failure_info, ec_scheme.get_config());

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
    double replace_time,
    const ErasureCodingScheme& ec_scheme,
    const DisconnectedStatus& disconnected) {

    DiskInfo& disk_info = disks[disk_index];

    if (event_type == "fail") {
        disk_info.is_failed = true;
        disk_info.fail_timestamp = 0;  // Will be set by caller

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
    }

    // Update state based on EC group
    int group_index = ec_scheme.get_disk_group(disk_index);
    auto it = failure_info_per_ec_group.find(group_index);
    if (it != failure_info_per_ec_group.end()) {
        std::string state_str = judge_state_from_failure_info(it->second, ec_scheme.get_config());
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
        auto group_it = failure_info_per_ec_group.find(group_index);

        if (group_it != failure_info_per_ec_group.end()) {
            int failure_count = static_cast<int>(group_it->second.failed_disk_indices.size());
            auto bw_it = flows_and_speed_entry.rebuild_bandwidth.find(failure_count);
            if (bw_it != flows_and_speed_entry.rebuild_bandwidth.end()) {
                rebuild_speed = bw_it->second;
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

        Event new_event{current_time + remaining_time, "repair", popped_event.event_node,
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
    int batch_size = (num_simulations + nprocs - 1) / nprocs;

    Logger::getInstance().info("Using " + std::to_string(nprocs) + " parallel processes");

    std::vector<std::future<SimulationResult>> futures;

    // Launch parallel simulations
    for (int i = 0; i < nprocs; ++i) {
        futures.push_back(std::async(std::launch::async, simulation_per_core,
                                     i, std::ref(params_and_results),
                                     std::ref(graph_structure_origin),
                                     batch_size, std::ref(options)));
    }

    Logger::getInstance().info("All simulation threads launched");

    // Collect results
    double total_up_time = 0.0;
    double total_simulation_time = 0.0;
    double total_time_for_rebuilding = 0.0;
    int total_count_for_rebuilding = 0;
    double total_mttdl = 0.0;
    int total_data_loss_events = 0;
    double total_data_loss_ratio = 0.0;
    double total_read_bw = 0.0;
    double total_write_bw = 0.0;

    for (auto& future : futures) {
        SimulationResult result = future.get();

        total_up_time += result.up_time;
        total_simulation_time += result.simulation_time;
        total_time_for_rebuilding += result.time_for_rebuilding;
        total_count_for_rebuilding += result.count_for_rebuilding;
        total_mttdl += result.mttdl;
        total_data_loss_events += result.data_loss_events;
        total_data_loss_ratio += result.total_data_loss_ratio;
        total_read_bw += result.average_read_bandwidth;
        total_write_bw += result.average_write_bandwidth;
    }

    Logger::getInstance().info("All simulation threads completed, aggregating results");

    // Calculate averages
    double avg_up_time = total_up_time / nprocs;
    double avg_simulation_time = total_simulation_time / nprocs;
    double availability = avg_up_time / avg_simulation_time;

    // Calculate MTTDL
    double avg_mttdl = (total_data_loss_events > 0) ? (total_mttdl / total_data_loss_events) : 0.0;

    // Calculate durability
    int simulation_years = options.value("simulation_years", 10);
    double avg_data_loss_ratio = total_data_loss_ratio / nprocs;
    double durability = calculate_durability(avg_data_loss_ratio, avg_simulation_time, simulation_years);

    // Store results
    params_and_results["availability"] = availability;
    params_and_results["avail_nines"] = Utils::get_nines(availability);
    params_and_results["simulation_time"] = avg_simulation_time;
    params_and_results["up_time"] = avg_up_time;
    params_and_results["mttdl"] = avg_mttdl;
    params_and_results["data_loss_events"] = total_data_loss_events;
    params_and_results["data_loss_ratio"] = avg_data_loss_ratio;
    params_and_results["durability"] = durability;
    params_and_results["durability_nines"] = Utils::get_nines(durability);

    if (total_count_for_rebuilding > 0) {
        params_and_results["avg_rebuilding_time"] = total_time_for_rebuilding / total_count_for_rebuilding;
    } else {
        params_and_results["avg_rebuilding_time"] = 0.0;
    }

    params_and_results["avg_read_bandwidth"] = total_read_bw / nprocs;
    params_and_results["avg_write_bandwidth"] = total_write_bw / nprocs;

    Logger::getInstance().info("Simulation completed");
    Logger::getInstance().info("Availability: " + std::to_string(availability) +
                              " (" + std::to_string(Utils::get_nines(availability)) + " nines)");
    Logger::getInstance().info("Durability: " + std::to_string(durability) +
                              " (" + std::to_string(Utils::get_nines(durability)) + " nines)");
}
