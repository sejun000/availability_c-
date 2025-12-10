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

    // Calculate total ways to choose stripe members
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

    // Cumulative probability for f = 0 ... k failures inside the stripe
    for (int f = 0; f <= std::min(k, x); ++f) {
        surviving += comb(x, f) * comb(N - x, S - f);
    }

    return 1.0 - surviving / total;
}

// Calculate bottleneck speed considering erasure coding and bandwidth
double calculate_bottleneck_speed(int m, int k, const std::vector<double>& other_bws,
                                   const nlohmann::json& options, bool replication) {
    double erasure_coding_latency = 0.00000001;

    if (k > 0 && !replication) {
        erasure_coding_latency = Utils::get_encoding_latency_sec(m, k);
    }

    double erasure_coding_speed = 256000.0 / erasure_coding_latency;
    double erasure_coding_cores = options.value("erasure_coding_cores", 1);
    double min_speed = erasure_coding_speed * erasure_coding_cores;

    double rebuild_bw_ratio = options.value("rebuild_bw_ratio", 0.2);

    for (const auto& other_bw : other_bws) {
        if (other_bw * rebuild_bw_ratio < min_speed) {
            min_speed = other_bw * rebuild_bw_ratio;
        }
    }

    return min_speed;
}

// Get DRAM bandwidth
double get_dram_bandwidth(const nlohmann::json& options) {
    std::string dram_bw_str = options.value("dram_bandwidth", "100G");
    double dram_bandwidth = Utils::KMG_to_bytes(dram_bw_str);

    bool nic_to_ssd_direct = options.value("nic_to_ssd_direct", true);
    if (!nic_to_ssd_direct) {
        dram_bandwidth /= 2.0;
    }

    return dram_bandwidth;
}

// Convert availability nines to credit ratio
double nine_to_credit(double nine) {
    if (nine >= 6.0) return 0.0;
    else if (nine >= 5.0 && nine < 6.0) return 0.03;
    else if (nine >= 4.0 && nine < 5.0) return 0.05;
    else if (nine >= 3.301 && nine < 4.0) return 0.1;
    else if (nine >= 3.0 && nine < 3.301) return 0.25;
    return 1.0;
}

// Check if failure is catastrophic
bool is_catastrophic_failure(const FailureInfo& failure_info, int ssd_k, bool local_module_disconnected) {
    return failure_info.failure_count > ssd_k || local_module_disconnected;
}

// Check if other nodes have catastrophic failure and are recoverable
bool is_other_nodes_catastrophic_failure_and_recoverable(
    const FailureInfo& failure_info,
    int ssd_k,
    int network_k,
    const DisconnectedStatus& disconnected) {

    if (disconnected.common_module == false &&
        failure_info.network_failure_count > 0 &&
        !is_catastrophic_failure(failure_info, ssd_k, disconnected.local_module) &&
        failure_info.network_failure_count <= network_k &&
        network_k > 0) {
        return true;
    }
    return false;
}

// Judge SSD state from failure information
std::string judge_state_from_failure_info(const FailureInfo& failure_info,
                                          const SSDRedundancyScheme& ssd_redun_scheme,
                                          const DisconnectedStatus& disconnected,
                                          bool cached,
                                          bool ignore_disconnected) {
    bool common_disconnected = disconnected.common_module;
    bool local_disconnected = disconnected.local_module;

    if (ignore_disconnected) {
        local_disconnected = false;
        common_disconnected = false;
    }

    int ssd_k = ssd_redun_scheme.get_k(cached);
    int network_k = ssd_redun_scheme.get_network_k(cached);

    if (common_disconnected) {
        return SSD_STATE_DATA_LOSS;
    }

    bool local_module_disconnected = local_disconnected;

    if (failure_info.failure_count > 0 &&
        !is_catastrophic_failure(failure_info, ssd_k, local_module_disconnected)) {
        return SSD_STATE_INTRA_REBUILDING;
    } else if (is_catastrophic_failure(failure_info, ssd_k, local_module_disconnected) &&
               failure_info.network_failure_count <= network_k - 1 && network_k > 0) {
        return SSD_STATE_INTER_REBUILDING;
    } else if ((network_k > 0 && failure_info.network_failure_count > network_k) ||
               (is_catastrophic_failure(failure_info, ssd_k, local_module_disconnected) &&
                failure_info.network_failure_count > network_k - 1)) {
        return SSD_STATE_DATA_LOSS;
    } else {
        return SSD_STATE_NORMAL;
    }
}

// Random number generator (thread-local)
thread_local std::mt19937 rng(std::random_device{}());

// Push failure event to priority queue
void push_failed_event(std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
                       const std::string& event_node,
                       double current_time,
                       const std::map<std::string, std::string>& node_to_module_map,
                       const GraphStructure& hardware_graph,
                       const SSDRedundancyScheme& ssd_redun_scheme,
                       const EmpiricalCDF* ssd_failure_cdf,
                       double ssd_arr) {
    double failure_time;

    if (is_event_node_ssd(event_node)) {
        // Use spliced distribution if CDF and ARR are available
        if (ssd_failure_cdf != nullptr && ssd_failure_cdf->is_initialized() && ssd_arr > 0) {
            failure_time = current_time + ssd_failure_cdf->sample_spliced(rng, ssd_arr);
        } else {
            int ssd_index = get_ssd_index(event_node);
            bool cached = ssd_redun_scheme.is_ssd_index_cached(ssd_index);
            double mttf = ssd_redun_scheme.get_mttf(cached);
            std::exponential_distribution<double> exp_dist(1.0 / mttf);
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
                       bool software_repair_only) {
    if (is_event_node_ssd(event_node)) {
        Event event{current_time + BIG_NUMBER, "repair", event_node, current_time, current_time, software_repair_only};
        repair_events.push(event);
    } else {
        auto it = node_to_module_map.find(event_node);
        if (it == node_to_module_map.end()) return;
        std::string module = it->second;
        auto mtr_it = hardware_graph.mtrs_.find(module);
        if (mtr_it == hardware_graph.mtrs_.end()) return;
        double repair_time = current_time + mtr_it->second;
        Event event{repair_time, "repair", event_node, current_time, current_time, software_repair_only};
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
    int ssd_total_count,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const EmpiricalCDF* ssd_failure_cdf,
    double ssd_arr) {

    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> events;

    // Add network event if needed
    if (ssd_redun_scheme.get_network_k(ssd_redun_scheme.is_ssd_group_index_cached(0)) > 0) {
        Event network_event{0, "network", "network", 0, 0, false};
        events.push(network_event);
    }

    // Add hardware node events
    for (const auto& node : hardware_graph.nodes_) {
        push_failed_event(events, node, 0, node_to_module_map, hardware_graph, ssd_redun_scheme, ssd_failure_cdf, ssd_arr);
    }

    // Add enclosure events
    for (const auto& [enclosure, _] : hardware_graph.enclosures_) {
        push_failed_event(events, enclosure, 0, node_to_module_map, hardware_graph, ssd_redun_scheme, ssd_failure_cdf, ssd_arr);
    }

    // Add SSD events
    for (int ssd_index = 0; ssd_index < ssd_total_count; ++ssd_index) {
        std::string ssd_name = get_ssd_name(ssd_index);
        push_failed_event(events, ssd_name, 0, node_to_module_map, hardware_graph, ssd_redun_scheme, ssd_failure_cdf, ssd_arr);
    }

    return events;
}

// Update failure info when event occurs
void update_failure_info(const std::string& event_type,
                        const std::string& event_node,
                        std::map<int, FailureInfo>& failure_info_per_ssd_group,
                        std::map<std::string, bool>& failed_nodes_and_enclosures,
                        const SSDRedundancyScheme& ssd_redun_scheme) {
    if (event_node == "network") {
        return;
    }

    if (is_event_node_ssd(event_node)) {
        int ssd_index = get_ssd_index(event_node);
        int ssd_group_index = ssd_redun_scheme.get_ssd_group_index(ssd_index);

        if (event_type == "fail") {
            failure_info_per_ssd_group[ssd_group_index].failure_count++;
        } else if (event_type == "repair") {
            failure_info_per_ssd_group[ssd_group_index].failure_count--;
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
    const GraphStructure& hardware_graph_copy,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const nlohmann::json& options,
    std::map<FailureStateKey, FlowsAndSpeedEntry>& flows_and_speed_table,
    double max_read_performance_without_any_failure,
    const DisconnectedStatus& disconnected,
    const FailureStateKey& key) {

    if (flows_and_speed_table.find(key) != flows_and_speed_table.end()) {
        return;
    }

    FlowsAndSpeedEntry entry;

    double dram_bandwidth = get_dram_bandwidth(options);

    // Calculate hardware max flow
    GraphStructure hw_copy = hardware_graph_copy;
    hw_copy.add_virtual_nodes(options["lowest_common_module"], options["end_module"]);
    double hardware_max_flow = hw_copy.maximum_flow("virtual_source", "virtual_sink");
    hw_copy.remove_virtual_nodes();

    hw_copy.add_virtual_nodes(options["start_module"], options["lowest_common_module"]);
    double common_module_max_flow = hw_copy.maximum_flow("virtual_source", "virtual_sink");
    hw_copy.remove_virtual_nodes();

    // Calculate total read BW for SSDs
    int normal_ssds_count = 0;
    double total_read_bw_for_ssds = 0.0;
    bool end_module_degraded = hw_copy.get_module_degraded(options["end_module"], options["active_active"]);

    for (const auto& [group_index, failure_info] : failure_info_per_ssd_group) {
        bool cached = ssd_redun_scheme.is_ssd_group_index_cached(group_index);
        int ssd_m = ssd_redun_scheme.get_m(cached);
        int ssd_k = ssd_redun_scheme.get_k(cached);
        int ssd_l = ssd_redun_scheme.get_l(cached);

        double ssd_read_bw = ssd_redun_scheme.get_read_bw(cached);
        if (end_module_degraded || (!options["active_active"] && !options["single_port_ssd"])) {
            ssd_read_bw /= 2.0;
        }

        int degraded_ssd_count = ssd_m + ssd_k + ssd_l - failure_info.failure_count;
        total_read_bw_for_ssds += ssd_read_bw * degraded_ssd_count;
        normal_ssds_count += degraded_ssd_count;
    }

    // Calculate local read degradation ratio (DRAM bandwidth bottleneck)
    if (total_read_bw_for_ssds == 0) {
        total_read_bw_for_ssds = 1.0 / BIG_NUMBER;
    }
    double local_read_degradation_ratio = dram_bandwidth / total_read_bw_for_ssds;
    if (local_read_degradation_ratio > 1.0) {
        local_read_degradation_ratio = 1.0;
    }

    total_read_bw_for_ssds *= local_read_degradation_ratio;
    double bottleneck_read_bw = std::min(total_read_bw_for_ssds, hardware_max_flow);

    if (normal_ssds_count == 0) {
        normal_ssds_count = 1;
    }
    double bottleneck_read_bw_per_ssd = bottleneck_read_bw / normal_ssds_count;

    // Process each SSD group
    int index = 0;
    for (const auto& [group_index, failure_info] : failure_info_per_ssd_group) {
        bool cached = ssd_redun_scheme.is_ssd_group_index_cached(group_index);
        if (cached) {
            index = 0;
        }

        double ssd_read_bw = ssd_redun_scheme.get_read_bw(cached);
        double ssd_write_bw = ssd_redun_scheme.get_write_bw(cached);
        if (end_module_degraded || (!options["active_active"] && !options["single_port_ssd"])) {
            ssd_read_bw /= 2.0;
            ssd_write_bw /= 2.0;
        }

        int network_m = ssd_redun_scheme.get_network_m(cached);
        int network_k = ssd_redun_scheme.get_network_k(cached);
        double local_ssd_read_bw = ssd_read_bw * local_read_degradation_ratio;
        std::string prefix = ssd_redun_scheme.get_cached_prefix(cached);
        int tiered_ssds = ssd_redun_scheme.get_tiered_ssds(cached);
        int ssd_m = ssd_redun_scheme.get_m(cached);
        int ssd_k = ssd_redun_scheme.get_k(cached);
        int ssd_l = ssd_redun_scheme.get_l(cached);
        int inter_replicas = ssd_redun_scheme.get_inter_replicas(cached);
        int intra_replicas = ssd_redun_scheme.get_intra_replicas(cached);
        bool inter_replication = (inter_replicas > 0);
        bool intra_replication = (intra_replicas > 0);

        std::string state = judge_state_from_failure_info(failure_info, ssd_redun_scheme, disconnected, cached);

        // Intra rebuilding
        if (state == SSD_STATE_INTRA_REBUILDING) {
            int local_failure_count = failure_info.failure_count;
            int degraded_ssds = ssd_m + ssd_k + ssd_l - local_failure_count;
            int rebuild_speed_up = 1;
            if (ssd_l > 0 && failure_info.failure_count <= ssd_l) {
                rebuild_speed_up = degraded_ssds;
            }

            double degraded_bw = calculate_bottleneck_speed(ssd_m, local_failure_count,
                                                           {local_ssd_read_bw}, options, intra_replication);
            double rebuilding_bw = calculate_bottleneck_speed(ssd_m, local_failure_count,
                                                             {local_ssd_read_bw, ssd_write_bw * rebuild_speed_up},
                                                             options, intra_replication);

            if (prefix == "cached_") {
                entry.cached_intra_rebuilding_bw[local_failure_count] = rebuilding_bw;
            } else {
                entry.intra_rebuilding_bw[local_failure_count] = rebuilding_bw;
            }

            if (intra_replicas > 0) {
                degraded_bw /= intra_replicas;
            }

            total_read_bw_for_ssds -= degraded_bw * degraded_ssds;
            total_read_bw_for_ssds += degraded_bw * local_failure_count;
        }
        // Inter rebuilding
        else if (state == SSD_STATE_INTER_REBUILDING) {
            int network_failure_count = failure_info.network_failure_count;
            double rebuilding_bw = calculate_bottleneck_speed(network_m, network_failure_count + 1,
                                                             {bottleneck_read_bw_per_ssd, ssd_write_bw},
                                                             options, inter_replication);
            if (prefix == "cached_") {
                entry.cached_inter_rebuilding_bw[network_failure_count] = rebuilding_bw;
            } else {
                entry.inter_rebuilding_bw[network_failure_count] = rebuilding_bw;
            }
        }
        // Data loss
        else if (state == SSD_STATE_DATA_LOSS) {
            double rebuilding_bw = calculate_bottleneck_speed(network_m, 0,
                                                             {bottleneck_read_bw_per_ssd, ssd_write_bw},
                                                             options, false);
            rebuilding_bw /= 8.0;
            if (prefix == "cached_") {
                entry.cached_backup_rebuild_speed = rebuilding_bw;
            } else {
                entry.backup_rebuild_speed = rebuilding_bw;
            }

            double data_loss_portion = pfail(ssd_m, ssd_k, ssd_l, failure_info.failure_count);
            int degraded_ssd_count = ssd_m + ssd_k + ssd_l - failure_info.failure_count;
            total_read_bw_for_ssds -= local_ssd_read_bw * degraded_ssd_count;

            double avail_reduction = (ssd_m + ssd_k + ssd_l) / static_cast<double>(tiered_ssds) * data_loss_portion;
            if (prefix == "cached_") {
                entry.availability_ratio.cached_availability -= avail_reduction;
            } else {
                entry.availability_ratio.availability -= avail_reduction;
            }
        }

        // Check data loss ignoring disconnection
        if (judge_state_from_failure_info(failure_info, ssd_redun_scheme, disconnected, cached, true) == SSD_STATE_DATA_LOSS) {
            if (index == 0) {
                entry.up_first_group = 0;
            }
        }

        // Check other nodes catastrophic failure
        if (is_other_nodes_catastrophic_failure_and_recoverable(failure_info, ssd_k, network_k, disconnected)) {
            double degraded_bw = calculate_bottleneck_speed(network_m, failure_info.network_failure_count,
                                                           {bottleneck_read_bw_per_ssd, ssd_write_bw},
                                                           options, inter_replication);
            int degraded_ssd_count = ssd_m + ssd_k + ssd_l - failure_info.failure_count;
            if (inter_replicas > 0) {
                degraded_bw /= inter_replicas;
            }
            bottleneck_read_bw -= degraded_bw * degraded_ssd_count;
        }

        if (!cached) {
            index++;
        }
    }

    // Calculate final performance metrics
    double max_read_performance = std::min(bottleneck_read_bw, total_read_bw_for_ssds);
    max_read_performance = std::min(max_read_performance, common_module_max_flow / options["network_nodes"].get<int>());

    if (!options["active_active"]) {
        max_read_performance_without_any_failure /= 2.0;
    }

    if (disconnected.local_module && ssd_redun_scheme.get_network_k(false) > 0) {
        max_read_performance = (1.0 - options["rebuild_bw_ratio"].get<double>()) * max_read_performance_without_any_failure;
    }

    entry.eff_availability_ratio = max_read_performance / max_read_performance_without_any_failure;

    if (max_read_performance < 1e-3) {
        entry.availability_ratio.availability = 0.0;
        if (ssd_redun_scheme.get_tiered_ssds(true) > 0) {
            entry.availability_ratio.cached_availability = 0.0;
        }
    }

    if (entry.eff_availability_ratio > 1.0) {
        entry.eff_availability_ratio = 1.0;
    }

    // Calculate credit availability
    double credit_avail_ratio = 1.0;
    double operational_utilization = options["target_perf_ratio"];
    if (entry.eff_availability_ratio < operational_utilization) {
        credit_avail_ratio = std::min(entry.eff_availability_ratio / operational_utilization,
                                     credit_avail_ratio);
    }

    double avail_ratio = entry.availability_ratio.availability * entry.availability_ratio.cached_availability;
    entry.credit_availability_ratio = std::min(credit_avail_ratio, avail_ratio);

    flows_and_speed_table[key] = entry;
}

// Update all SSD states
void update_all_ssd_states(
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    std::map<int, SSDInfo>& SSDs,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const DisconnectedStatus& disconnected,
    double current_time) {

    for (int ssd_index = 0; ssd_index < ssd_redun_scheme.get_total_ssds(); ++ssd_index) {
        if (!SSDs[ssd_index].failed) {
            continue;
        }

        int ssd_group_index = ssd_redun_scheme.get_ssd_group_index(ssd_index);
        bool cached = ssd_redun_scheme.is_ssd_index_cached(ssd_index);
        const FailureInfo& failure_info = failure_info_per_ssd_group.at(ssd_group_index);

        std::string changed_state = judge_state_from_failure_info(failure_info, ssd_redun_scheme, disconnected, cached);
        SSDs[ssd_index].state = changed_state;

        if (!disconnected.local_module && !disconnected.common_module &&
            SSDs[ssd_index].disconnected_timestamp != 0 && SSDs[ssd_index].failed) {
            // If local module is not disconnected, we can assume that the SSDs are not disconnected
            SSDs[ssd_index].remaining_prep_time_for_rebuilding = current_time - SSDs[ssd_index].disconnected_timestamp;
            SSDs[ssd_index].disconnected_timestamp = 0;
        }
    }
}

// Update single SSD state
void update_ssd_state(
    const std::string& ssd_name,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    std::map<int, SSDInfo>& SSDs,
    double capacity,
    const std::string& event_type,
    double prep_time_for_rebuilding,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const DisconnectedStatus& disconnected) {

    int ssd_index = get_ssd_index(ssd_name);
    bool cached = ssd_redun_scheme.is_ssd_index_cached(ssd_index);

    if (event_type == "fail") {
        SSDs[ssd_index].failed = true;

        if (SSDs[ssd_index].disconnected_timestamp != 0) {
            // We just set prep time to show rebuild time after disconnected
            SSDs[ssd_index].remaining_capacity_to_rebuild = 0;
            // It can be changed when rebuild speed is positive in update_repair_event_for_SSDs
            SSDs[ssd_index].remaining_prep_time_for_rebuilding = BIG_NUMBER + BIG_NUMBER;
        } else {
            SSDs[ssd_index].remaining_capacity_to_rebuild = capacity;
            if (cached) {
                SSDs[ssd_index].remaining_capacity_to_rebuild = capacity / 2.0;
            }
            SSDs[ssd_index].rebuild_speed = 0;
            std::exponential_distribution<double> exp_dist(1.0 / prep_time_for_rebuilding);
            SSDs[ssd_index].remaining_prep_time_for_rebuilding = exp_dist(rng);
        }
    } else if (event_type == "repair") {
        SSDs[ssd_index].failed = false;
    }

    int m = ssd_redun_scheme.get_m(cached);
    int k = ssd_redun_scheme.get_k(cached);
    int l = ssd_redun_scheme.get_l(cached);
    int n = m + k + l;
    int group_index = ssd_redun_scheme.get_ssd_group_index(ssd_index);

    const FailureInfo& failure_info = failure_info_per_ssd_group.at(group_index);
    std::string changed_state = judge_state_from_failure_info(failure_info, ssd_redun_scheme, disconnected, cached);

    int start_ssd_index = ssd_redun_scheme.get_start_ssd_index(group_index);
    for (int i = start_ssd_index; i < start_ssd_index + n; ++i) {
        if (!SSDs[i].failed) {
            continue;
        }
        SSDs[i].state = changed_state;
    }
}

// Update failure events for SSDs
void update_failure_event_for_SSDs(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
    double current_time,
    std::map<int, SSDInfo>& SSDs) {

    std::vector<Event> updated_failed_events;

    while (!failed_events.empty()) {
        Event popped_event = failed_events.top();
        failed_events.pop();

        if (!is_event_node_ssd(popped_event.event_node)) {
            updated_failed_events.push_back(popped_event);
            continue;
        }

        if (popped_event.event_type != "fail") {
            std::cerr << "Error: failed_event is not fail: " << popped_event.event_type << " " << popped_event.event_node << std::endl;
            continue;
        }

        int ssd_index = get_ssd_index(popped_event.event_node);
        if (SSDs[ssd_index].failed) {
            updated_failed_events.push_back(popped_event);
            continue;
        }

        double failure_time = current_time + 1.0 / BIG_NUMBER;
        SSDs[ssd_index].disconnected_timestamp = current_time;

        Event new_event{failure_time, "fail", popped_event.event_node,
                       popped_event.prev_time, popped_event.start_time, true};
        updated_failed_events.push_back(new_event);
    }

    for (const auto& event : updated_failed_events) {
        failed_events.push(event);
    }
}

// Update repair events for SSDs
void update_repair_event_for_SSDs(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events,
    double current_time,
    std::map<int, SSDInfo>& SSDs,
    const FlowsAndSpeedEntry& flows_and_speed_entry,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const SSDRedundancyScheme& ssd_redun_scheme) {

    std::vector<Event> updated_repair_events;

    while (!repair_events.empty()) {
        Event popped_event = repair_events.top();
        repair_events.pop();

        if (!is_event_node_ssd(popped_event.event_node)) {
            updated_repair_events.push_back(popped_event);
            continue;
        }

        int repaired_ssd_index = get_ssd_index(popped_event.event_node);
        SSDInfo& current_ssd_info = SSDs[repaired_ssd_index];
        double remaining_consumed_time = current_time - popped_event.prev_time;

        // Before the repair, SSD shall be switched to standby
        if (current_ssd_info.remaining_prep_time_for_rebuilding > 0) {
            double remaining_prep_time = current_ssd_info.remaining_prep_time_for_rebuilding - remaining_consumed_time;
            if (remaining_prep_time <= 0) {
                current_ssd_info.remaining_prep_time_for_rebuilding = 0;
                remaining_consumed_time = -remaining_prep_time;
            } else {
                current_ssd_info.remaining_prep_time_for_rebuilding = remaining_prep_time;
                remaining_consumed_time = 0;
            }
        }

        // Calculate the remaining capacity with previous rebuild speed
        if (remaining_consumed_time > 0) {
            double rebuild_speed = current_ssd_info.rebuild_speed;
            double rebuild_speed_per_hour = rebuild_speed * 3600.0;
            double remaining_capacity = current_ssd_info.remaining_capacity_to_rebuild -
                                       remaining_consumed_time * rebuild_speed_per_hour;
            current_ssd_info.remaining_capacity_to_rebuild = remaining_capacity;
        }

        // Change rebuild speed to the current state
        double rebuild_speed = 0;
        int ssd_group_index = ssd_redun_scheme.get_ssd_group_index(repaired_ssd_index);
        bool cached = ssd_redun_scheme.is_ssd_group_index_cached(ssd_group_index);
        std::string prefix = ssd_redun_scheme.get_cached_prefix(cached);
        int ssd_failure_count = failure_info_per_ssd_group.at(ssd_group_index).failure_count;
        int network_failure_count = failure_info_per_ssd_group.at(ssd_group_index).network_failure_count;

        if (current_ssd_info.state == SSD_STATE_INTRA_REBUILDING) {
            const auto& bw_map = (prefix == "cached_") ?
                flows_and_speed_entry.cached_intra_rebuilding_bw :
                flows_and_speed_entry.intra_rebuilding_bw;
            auto it = bw_map.find(ssd_failure_count);
            if (it != bw_map.end()) {
                rebuild_speed = it->second;
            }
        } else if (current_ssd_info.state == SSD_STATE_INTER_REBUILDING) {
            const auto& bw_map = (prefix == "cached_") ?
                flows_and_speed_entry.cached_inter_rebuilding_bw :
                flows_and_speed_entry.inter_rebuilding_bw;
            auto it = bw_map.find(network_failure_count);
            if (it != bw_map.end()) {
                rebuild_speed = it->second;
            }
        } else if (current_ssd_info.state == SSD_STATE_DATA_LOSS) {
            rebuild_speed = (prefix == "cached_") ?
                flows_and_speed_entry.cached_backup_rebuild_speed :
                flows_and_speed_entry.backup_rebuild_speed;
        }

        current_ssd_info.rebuild_speed = rebuild_speed;

        if (rebuild_speed == 0) {
            rebuild_speed = 1.0 / BIG_NUMBER;
        }

        double rebuild_speed_per_hour = rebuild_speed * 3600.0;
        if (current_ssd_info.remaining_capacity_to_rebuild < 0) {
            current_ssd_info.remaining_capacity_to_rebuild = 0;
        }

        double updated_repair_event_time = current_ssd_info.remaining_prep_time_for_rebuilding +
                                           current_ssd_info.remaining_capacity_to_rebuild / rebuild_speed_per_hour;

        Event new_event{current_time + updated_repair_event_time, "repair", popped_event.event_node,
                       current_time, popped_event.start_time, popped_event.software_repair_only};
        updated_repair_events.push_back(new_event);
    }

    for (const auto& event : updated_repair_events) {
        repair_events.push(event);
    }
}

// Calculate combinations (n choose k)
double combinations_count(int n, int k) {
    if (n < k) return 0.0;
    if (k == 0 || k == n) return 1.0;

    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}

// Generate network failure probability table
void generate_network_failure_table(
    int cached_network_n,
    int network_n,
    double availability_without_network_parity,
    double availability_without_network_parity_for_cached_ssds,
    NetworkAvailabilityTable& network_availability_table) {

    double probability = 0.0;
    double single_availability = availability_without_network_parity;
    double single_availability_for_cached_ssds = availability_without_network_parity_for_cached_ssds;

    // Reset tables before computing new probabilities
    network_availability_table.availability.clear();
    network_availability_table.cached_availability.clear();

    // Calculate for regular (uncached) SSDs
    for (int failed_count = 0; failed_count < network_n; ++failed_count) {
        probability += combinations_count(network_n - 1, failed_count) *
                      std::pow(1.0 - single_availability, failed_count) *
                      std::pow(single_availability, network_n - 1 - failed_count);
        network_availability_table.availability[failed_count] = probability;
    }
    if (network_n > 0) {
        network_availability_table.availability[network_n - 1] = 1.0;
    }

    // Reset probability for cached SSDs
    probability = 0.0;

    // Calculate for cached SSDs
    for (int failed_count = 0; failed_count < cached_network_n; ++failed_count) {
        probability += combinations_count(cached_network_n - 1, failed_count) *
                      std::pow(1.0 - single_availability_for_cached_ssds, failed_count) *
                      std::pow(single_availability_for_cached_ssds, cached_network_n - 1 - failed_count);
        network_availability_table.cached_availability[failed_count] = probability;
    }
    if (cached_network_n > 0) {
        network_availability_table.cached_availability[cached_network_n - 1] = 1.0;
    }

}

// Update network failure state based on random sampling
bool update_network_state(
    std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const NetworkAvailabilityTable& network_availability_table) {

    bool changed = false;
    int total_group_count = ssd_redun_scheme.get_total_group_count();

    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    for (int group_index = 0; group_index < total_group_count; ++group_index) {
        double random_value = uniform_dist(rng);
        FailureInfo& failure_info = failure_info_per_ssd_group[group_index];
        bool cached = ssd_redun_scheme.is_ssd_group_index_cached(group_index);
        std::string prefix = ssd_redun_scheme.get_cached_prefix(cached);
        int network_n = ssd_redun_scheme.get_network_m(cached) + ssd_redun_scheme.get_network_k(cached);

        // Select appropriate availability table
        const auto& avail_table = (prefix == "cached_") ?
            network_availability_table.cached_availability :
            network_availability_table.availability;

        // Find the failure count based on random value
        for (int failed_count = 0; failed_count < network_n; ++failed_count) {
            auto it = avail_table.find(failed_count);
            if (it != avail_table.end() && random_value <= it->second) {
                if (failed_count != failure_info.network_failure_count) {
                    failure_info.network_failure_count = failed_count;
                    changed = true;
                }
                break;
            }
        }
    }

    return changed;
}

// Get cost coefficient for a module
double get_coefficient_for_cost(const std::string& module, const nlohmann::json& options) {
    int network_nodes = options.value("network_nodes", 1);

    if (module == SSD_MODULE_NAME ||
        module == "io_module" ||
        module == "host_module" ||
        module == "backend_module" ||
        module == "NVMeEnclosure") {
        return static_cast<double>(network_nodes);
    }

    return 1.0;
}

// Calculate cost of a module
double calculate_module_cost(
    const std::string& node,
    const std::map<std::string, std::string>& node_to_module_map,
    const std::map<std::string, double>& costs,
    const SSDRedundancyScheme& ssd_redun_scheme,
    double cached_ssd_cost,
    double uncached_ssd_cost,
    const nlohmann::json& options) {

    // Check if it's an SSD
    if (node.find(SSD_MODULE_NAME) == std::string::npos) {
        // Regular hardware module
        auto it = node_to_module_map.find(node);
        if (it == node_to_module_map.end()) {
            return 0.0;
        }
        std::string module = it->second;
        auto cost_it = costs.find(module);
        if (cost_it == costs.end()) {
            return 0.0;
        }
        return cost_it->second * get_coefficient_for_cost(module, options);
    } else {
        // SSD module
        int ssd_index = get_ssd_index(node);
        bool cached = ssd_redun_scheme.is_ssd_index_cached(ssd_index);
        if (cached) {
            return cached_ssd_cost;
        }
        return uncached_ssd_cost * get_coefficient_for_cost(SSD_MODULE_NAME, options);
    }
}

// Calculate initial hardware cost
std::tuple<double, double, double> get_initial_cost(
    const GraphStructure& hardware_graph,
    const std::map<std::string, std::string>& node_to_module_map,
    int ssd_total_count,
    const SSDRedundancyScheme& ssd_redun_scheme,
    double cached_ssd_cost,
    double uncached_ssd_cost,
    const std::map<std::string, double>& costs,
    const nlohmann::json& options) {

    double initial_cost = 0.0;
    double cached_initial_cost = 0.0;
    double uncached_initial_cost = 0.0;

    // Cost of hardware nodes
    for (const auto& node : hardware_graph.get_nodes()) {
        auto it = node_to_module_map.find(node);
        if (it == node_to_module_map.end()) continue;

        std::string module = it->second;
        auto cost_it = costs.find(module);
        if (cost_it == costs.end()) continue;

        initial_cost += cost_it->second * get_coefficient_for_cost(module, options);
    }

    // Cost of enclosures
    for (const auto& [enclosure, _] : hardware_graph.enclosures_) {
        auto it = node_to_module_map.find(enclosure);
        if (it == node_to_module_map.end()) continue;

        std::string module = it->second;
        auto cost_it = costs.find(module);
        if (cost_it == costs.end()) continue;

        initial_cost += cost_it->second * get_coefficient_for_cost(module, options);
    }

    // Cost of SSDs
    double ssd_coef = get_coefficient_for_cost(SSD_MODULE_NAME, options);
    for (int ssd_index = 0; ssd_index < ssd_total_count; ++ssd_index) {
        bool cached = ssd_redun_scheme.is_ssd_index_cached(ssd_index);
        if (cached) {
            double cost = cached_ssd_cost * ssd_coef;
            initial_cost += cost;
            cached_initial_cost += cost;
        } else {
            double cost = uncached_ssd_cost * ssd_coef;
            initial_cost += cost;
            uncached_initial_cost += cost;
        }
    }

    return {initial_cost, cached_initial_cost, uncached_initial_cost};
}

// Monte Carlo simulation main function
void monte_carlo_simulation(
    std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int num_simulations,
    const nlohmann::json& options,
    const std::map<std::string, double>& costs
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
                                     batch_size, std::ref(options), std::ref(costs)));
    }

    Logger::getInstance().info("All simulation threads launched");

    // Collect results
    double total_up_time = 0.0;
    double total_cached_up_time = 0.0;
    double total_credit_up_time = 0.0;
    double total_simulation_time = 0.0;
    double total_effective_up_time = 0.0;
    std::map<double, double> total_effective_availabilities;
    double total_initial_cost = 0.0;
    double total_cost = 0.0;
    double total_time_for_rebuilding = 0.0;
    int total_count_for_rebuilding = 0;
    double total_cached_ssd_repair_cost = 0.0;
    double total_uncached_ssd_repair_cost = 0.0;
    double total_cached_initial_cost = 0.0;
    double total_uncached_initial_cost = 0.0;
    double total_penalty_ratio = 0.0;
    double total_mttdl = 0.0;
    int total_mttdl_count = 0;

    double cached_mttf = 0.0;
    double mttf = 0.0;

    for (auto& future : futures) {
        SimulationResult result = future.get();

        total_up_time += result.up_time;
        total_cached_up_time += result.cached_up_time;
        total_credit_up_time += result.credit_up_time;
        total_simulation_time += result.simulation_time;
        total_effective_up_time += result.effective_up_time;

        for (const auto& [avail, time] : result.effective_availabilities) {
            total_effective_availabilities[avail] += time;
        }

        total_initial_cost += result.initial_cost;
        total_cost += result.total_cost;
        total_time_for_rebuilding += result.time_for_rebuilding;
        total_count_for_rebuilding += result.count_for_rebuilding;
        total_cached_ssd_repair_cost += result.cached_ssd_repair_cost;
        total_uncached_ssd_repair_cost += result.uncached_ssd_repair_cost;
        total_cached_initial_cost += result.cached_initial_cost;
        total_uncached_initial_cost += result.uncached_initial_cost;
        total_penalty_ratio += result.penalty_ratio;
        total_mttdl += result.total_mttdl;
        total_mttdl_count += result.total_mttdl_count;

        cached_mttf = result.cached_mttf;
        mttf = result.mttf;
    }

    Logger::getInstance().info("All simulation threads completed, aggregating results");

    // Calculate averages
    double avg_up_time = total_up_time / nprocs;
    double avg_cached_up_time = total_cached_up_time / nprocs;
    double avg_credit_up_time = total_credit_up_time / nprocs;
    double avg_simulation_time = total_simulation_time / nprocs;
    double avg_effective_up_time = total_effective_up_time / nprocs;
    double avg_initial_cost = total_initial_cost / nprocs;
    double avg_total_cost = total_cost / nprocs;
    double avg_penalty_ratio = total_penalty_ratio / nprocs;

    // Calculate availability
    double uncached_availability = avg_up_time / avg_simulation_time;
    double cached_availability = avg_cached_up_time / avg_simulation_time;
    double availability = uncached_availability * cached_availability;
    double effective_availability = avg_effective_up_time / avg_simulation_time;
    double credit_availability = avg_credit_up_time / avg_simulation_time;

    // Calculate MTTDL
    double avg_mttdl = (total_mttdl_count > 0) ? (total_mttdl / total_mttdl_count) : 0.0;

    // Calculate percentiles for effective availability
    auto [avg_eff_avail, median_eff_avail, p99, p999, p9999] =
        Utils::get_percentile_value(total_effective_availabilities, false);

    // Store results
    params_and_results["uncached_availability"] = uncached_availability;
    params_and_results["availability"] = availability;
    params_and_results["cached_availability"] = cached_availability;
    params_and_results["effective_availability"] = effective_availability;
    params_and_results["credit_availability"] = credit_availability;
    params_and_results["avail_nines"] = Utils::get_nines(availability);
    params_and_results["cached_avail_nines"] = Utils::get_nines(cached_availability);
    params_and_results["effective_avail_nines"] = Utils::get_nines(effective_availability);
    params_and_results["credit_avail_nines"] = Utils::get_nines(credit_availability);
    params_and_results["simulation_time"] = avg_simulation_time;
    params_and_results["up_time"] = avg_up_time;
    params_and_results["cached_up_time"] = avg_cached_up_time;
    params_and_results["effective_up_time"] = avg_effective_up_time;
    params_and_results["credit_up_time"] = avg_credit_up_time;

    params_and_results["initial_cost"] = avg_initial_cost;
    params_and_results["total_cost"] = avg_total_cost;
    params_and_results["cached_initial_cost"] = total_cached_initial_cost / nprocs;
    params_and_results["uncached_initial_cost"] = total_uncached_initial_cost / nprocs;
    double avg_cached_ssd_repair_cost = total_cached_ssd_repair_cost / nprocs;
    double avg_uncached_ssd_repair_cost = total_uncached_ssd_repair_cost / nprocs;
    params_and_results["cached_ssd_repair_cost"] = avg_cached_ssd_repair_cost;
    params_and_results["uncached_ssd_repair_cost"] = avg_uncached_ssd_repair_cost;
    params_and_results["penalty_ratio"] = avg_penalty_ratio;

    params_and_results["mttf"] = mttf;
    params_and_results["cached_mttf"] = cached_mttf;
    params_and_results["mttdl"] = avg_mttdl;

    if (total_count_for_rebuilding > 0) {
        params_and_results["avg_rebuilding_time"] = total_time_for_rebuilding / total_count_for_rebuilding;
    } else {
        params_and_results["avg_rebuilding_time"] = 0.0;
    }

    params_and_results["avg_effective_availability"] = avg_eff_avail;
    params_and_results["median_effective_availability"] = median_eff_avail;
    params_and_results["p99_effective_availability"] = p99;
    params_and_results["p999_effective_availability"] = p999;
    params_and_results["p9999_effective_availability"] = p9999;

    // Derived cost and efficiency metrics (parity with original Python tool)
    double simulation_years = options.value("simulation_years", 10);
    double simulation_hours = simulation_years * 365.0 * 24.0;
    double batch_multiplier = (simulation_hours > 0.0) ? (avg_simulation_time / simulation_hours) : 0.0;
    double denominator = (simulation_years > 0.0) ? simulation_years * batch_multiplier : 0.0;

    double repair_cost_per_year = 0.0;
    double cached_ssd_repair_cost_per_year = 0.0;
    double uncached_ssd_repair_cost_per_year = 0.0;
    double down_cost_per_year = 0.0;
    double operation_cost_per_year = 0.0;
    double repair_cost_for_10_years = 0.0;
    double down_cost_for_10_years = 0.0;
    double cached_ssd_repair_cost_for_10_years = 0.0;
    double uncached_ssd_repair_cost_for_10_years = 0.0;
    double total_cost_for_10_years = avg_initial_cost;

    double credit_ratio = avg_penalty_ratio;

    if (denominator > 0.0) {
        repair_cost_per_year = (avg_total_cost - avg_initial_cost) / denominator;
        cached_ssd_repair_cost_per_year = avg_cached_ssd_repair_cost / denominator;
        uncached_ssd_repair_cost_per_year = avg_uncached_ssd_repair_cost / denominator;

        down_cost_per_year = credit_ratio * (avg_initial_cost / 10.0 + repair_cost_per_year);
        operation_cost_per_year = repair_cost_per_year + down_cost_per_year;

        repair_cost_for_10_years = repair_cost_per_year * 10.0;
        down_cost_for_10_years = down_cost_per_year * 10.0;
        cached_ssd_repair_cost_for_10_years = cached_ssd_repair_cost_per_year * 10.0;
        uncached_ssd_repair_cost_for_10_years = uncached_ssd_repair_cost_per_year * 10.0;
        total_cost_for_10_years = avg_initial_cost + operation_cost_per_year * 10.0;
    }

    double capacity = params_and_results.at("capacity").get<double>();
    int total_ssds = params_and_results.at("total_ssds").get<int>();
    int cached_ssds = params_and_results.at("cached_ssds").get<int>();
    int m_data = params_and_results.at("m").get<int>();
    int k_parity = params_and_results.at("k").get<int>();
    int l_spares = params_and_results.at("l").get<int>();
    int network_m = params_and_results.at("network_m").get<int>();
    int network_k = params_and_results.at("network_k").get<int>();

    double effective_capacity = 0.0;
    int uncached_ssds = total_ssds - cached_ssds;
    if (uncached_ssds > 0 && (m_data + k_parity) > 0 && (network_m + network_k) > 0) {
        effective_capacity = capacity * static_cast<double>(uncached_ssds) *
                             static_cast<double>(m_data) / static_cast<double>(m_data + k_parity) *
                             static_cast<double>(network_m) / static_cast<double>(network_m + network_k);

        if (l_spares > 0) {
            int min_val = std::min(k_parity, l_spares);
            effective_capacity *= static_cast<double>(m_data) / static_cast<double>(m_data + min_val);
        }
    }

    int network_nodes = options.value("network_nodes", 1);
    double initial_cost_per_gb = 0.0;
    double total_cost_per_gb = 0.0;
    if (effective_capacity > 0.0 && network_nodes > 0) {
        initial_cost_per_gb = avg_initial_cost / effective_capacity * 1.0e9 / static_cast<double>(network_nodes);
        total_cost_per_gb = total_cost_for_10_years / effective_capacity * 1.0e9 / static_cast<double>(network_nodes);
    }

    params_and_results["simulation_year"] = simulation_years;
    params_and_results["total_credit_ratio"] = credit_ratio;
    params_and_results["repair_cost_per_year"] = repair_cost_per_year;
    params_and_results["cached_ssd_repair_cost_per_year"] = cached_ssd_repair_cost_per_year;
    params_and_results["uncached_ssd_repair_cost_per_year"] = uncached_ssd_repair_cost_per_year;
    params_and_results["down_cost_per_year"] = down_cost_per_year;
    params_and_results["operation_cost_per_year"] = operation_cost_per_year;
    params_and_results["repair_cost_for_10_years"] = repair_cost_for_10_years;
    params_and_results["down_cost_for_10_years"] = down_cost_for_10_years;
    params_and_results["cached_ssd_repair_cost_for_10_years"] = cached_ssd_repair_cost_for_10_years;
    params_and_results["uncached_ssd_repair_cost_for_10_years"] = uncached_ssd_repair_cost_for_10_years;
    params_and_results["total_cost_for_10_years"] = total_cost_for_10_years;
    params_and_results["effective_capacity"] = effective_capacity;
    params_and_results["initial_cost_per_gb"] = initial_cost_per_gb;
    params_and_results["cost_per_gb"] = total_cost_per_gb;

    // Write availability to file for next run only when no network parity is configured
    if (network_k == 0) {
        std::string avail_file = options.value("avail_file", "last_avail.txt");
        std::ofstream avail_out(avail_file);
        if (avail_out.is_open()) {
            avail_out << uncached_availability << "\n";
            avail_out << cached_availability << "\n";
            avail_out.close();
        }
    }

    Logger::getInstance().info("Simulation completed");
    Logger::getInstance().info("Availability: " + std::to_string(availability) +
                              " (" + std::to_string(Utils::get_nines(availability)) + " nines)");
    Logger::getInstance().info("Effective Availability: " + std::to_string(effective_availability) +
                              " (" + std::to_string(Utils::get_nines(effective_availability)) + " nines)");
}
