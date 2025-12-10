#pragma once

#include <map>
#include <vector>
#include <string>
#include <queue>
#include <random>
#include <cstdint>
#include <nlohmann/json.hpp>
#include "graph_structure.hpp"
#include "ssd.hpp"
#include "utils.hpp"

constexpr double BIG_NUMBER = 1e15;
constexpr double SOFT_ERROR_THRESHOLD_HOURS = 1.0;

// SSD States
const std::string SSD_STATE_NORMAL = "normal";
const std::string SSD_STATE_INTRA_REBUILDING = "intra_rebuilding";
const std::string SSD_STATE_INTER_REBUILDING = "inter_rebuilding";
const std::string SSD_STATE_INTER_DEGRADED = "inter_degraded";
const std::string SSD_STATE_DATA_LOSS = "data_loss";

// Event structure for priority queue
struct Event {
    double time;
    std::string event_type;  // "fail", "repair", "network"
    std::string event_node;
    double prev_time;
    double start_time;
    bool software_repair_only;

    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// SSD information structure
struct SSDInfo {
    bool failed;
    double remaining_prep_time_for_rebuilding;
    double remaining_capacity_to_rebuild;
    double rebuild_speed;
    double disconnected_timestamp;
    std::string state;

    SSDInfo() : failed(false), remaining_prep_time_for_rebuilding(0),
                remaining_capacity_to_rebuild(0), rebuild_speed(0),
                disconnected_timestamp(0), state(SSD_STATE_NORMAL) {}
};

// Failure information per SSD group
struct FailureInfo {
    int failure_count;
    int network_failure_count;

    FailureInfo() : failure_count(0), network_failure_count(0) {}
};

// Disconnection status
struct DisconnectedStatus {
    bool local_module;
    bool common_module;
    std::map<int, bool> ssd_disconnected;  // SSD index -> disconnected due to io_module failure

    DisconnectedStatus() : local_module(false), common_module(false) {}

    // Get count of disconnected SSDs in a group
    int get_disconnected_count_in_group(int start_ssd_index, int group_size) const {
        int count = 0;
        for (int i = start_ssd_index; i < start_ssd_index + group_size; ++i) {
            auto it = ssd_disconnected.find(i);
            if (it != ssd_disconnected.end() && it->second) {
                count++;
            }
        }
        return count;
    }

    // Check if specific SSD is disconnected
    bool is_ssd_disconnected(int ssd_index) const {
        auto it = ssd_disconnected.find(ssd_index);
        return it != ssd_disconnected.end() && it->second;
    }
};

// Availability ratio
struct AvailabilityRatio {
    double availability;
    double cached_availability;

    AvailabilityRatio() : availability(1.0), cached_availability(1.0) {}
};

// Flow and speed tables
struct FlowsAndSpeedEntry {
    std::map<int, double> intra_rebuilding_bw;
    std::map<int, double> inter_rebuilding_bw;
    double backup_rebuild_speed;
    std::map<int, double> cached_intra_rebuilding_bw;
    std::map<int, double> cached_inter_rebuilding_bw;
    double cached_backup_rebuild_speed;
    int up_first_group;

    AvailabilityRatio availability_ratio;
    double eff_availability_ratio;
    double credit_availability_ratio;

    FlowsAndSpeedEntry() : backup_rebuild_speed(0), cached_backup_rebuild_speed(0),
                           up_first_group(1), eff_availability_ratio(1.0),
                           credit_availability_ratio(1.0) {}
};

// Key structures for caching
struct NodeFailureKey {
    std::vector<uint8_t> failed_flags;

    bool operator<(const NodeFailureKey& other) const {
        return failed_flags < other.failed_flags;
    }
};

struct FailureStateKey {
    NodeFailureKey node_key;
    std::vector<uint8_t> failure_counts;
    std::vector<uint8_t> network_failure_counts;

    bool operator<(const FailureStateKey& other) const {
        if (node_key < other.node_key) return true;
        if (other.node_key < node_key) return false;
        if (failure_counts < other.failure_counts) return true;
        if (other.failure_counts < failure_counts) return false;
        return network_failure_counts < other.network_failure_counts;
    }
};

// Simulation parameters
struct SimulationParams {
    int m, k, l;
    int cached_m, cached_k, cached_l;
    int network_m, network_k, network_l;
    int cached_network_m, cached_network_k, cached_network_l;
    int total_ssds;
    int cached_ssds;
    int inter_replicas, intra_replicas;
    double cached_write_ratio;
    double capacity;
    double dwpd, dwpd_limit;
    double cached_dwpd_limit;
    double guaranteed_years;
    double read_bw, write_bw;
    double cached_read_bw, cached_write_bw;
    bool qlc, qlc_cache;
    double box_mttf, io_module_mttr;
    double rebuild_bw_ratio;
    double target_perf_ratio;
    bool active_active;
    bool single_port_ssd;
    int nprocs;
    int simulation_years;
};

// Network availability table structure
struct NetworkAvailabilityTable {
    std::map<int, double> availability;
    std::map<int, double> cached_availability;
};

// Simulation helper functions
double pfail(int m, int k, int l, int x);
double calculate_bottleneck_speed(int m, int k, const std::vector<double>& other_bws,
                                   const nlohmann::json& options, bool replication = false);
double get_dram_bandwidth(const nlohmann::json& options);
double nine_to_credit(double nine);
double combinations_count(int n, int k);

// Network failure functions
void generate_network_failure_table(
    int cached_network_n,
    int network_n,
    double availability_without_network_parity,
    double availability_without_network_parity_for_cached_ssds,
    NetworkAvailabilityTable& network_availability_table);

bool update_network_state(
    std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const NetworkAvailabilityTable& network_availability_table);

// Cost calculation functions
double get_coefficient_for_cost(const std::string& module, const nlohmann::json& options);

double calculate_module_cost(
    const std::string& node,
    const std::map<std::string, std::string>& node_to_module_map,
    const std::map<std::string, double>& costs,
    const SSDRedundancyScheme& ssd_redun_scheme,
    double cached_ssd_cost,
    double uncached_ssd_cost,
    const nlohmann::json& options);

std::tuple<double, double, double> get_initial_cost(
    const GraphStructure& hardware_graph,
    const std::map<std::string, std::string>& node_to_module_map,
    int ssd_total_count,
    const SSDRedundancyScheme& ssd_redun_scheme,
    double cached_ssd_cost,
    double uncached_ssd_cost,
    const std::map<std::string, double>& costs,
    const nlohmann::json& options);

// State management
std::string judge_state_from_failure_info(const FailureInfo& failure_info,
                                          const SSDRedundancyScheme& ssd_redun_scheme,
                                          const DisconnectedStatus& disconnected,
                                          bool cached,
                                          bool ignore_disconnected = false);

bool is_catastrophic_failure(const FailureInfo& failure_info, int ssd_k, bool local_module_disconnected);

bool is_other_nodes_catastrophic_failure_and_recoverable(
    const FailureInfo& failure_info,
    int ssd_k,
    int network_k,
    const DisconnectedStatus& disconnected);

// Event management
void push_failed_event(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
    const std::string& event_node,
    double current_time,
    const std::map<std::string, std::string>& node_to_module_map,
    const GraphStructure& hardware_graph,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const EmpiricalCDF* ssd_failure_cdf = nullptr,
    double ssd_arr = 0.0);

void push_repair_event(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events,
    const std::string& event_node,
    double current_time,
    const std::map<std::string, std::string>& node_to_module_map,
    const GraphStructure& hardware_graph,
    bool software_repair_only = false,
    const nlohmann::json* options = nullptr);

Event pop_event(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& events,
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events);

std::priority_queue<Event, std::vector<Event>, std::greater<Event>> generate_first_failure_events(
    const GraphStructure& hardware_graph,
    const std::map<std::string, std::string>& node_to_module_map,
    int ssd_total_count,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const EmpiricalCDF* ssd_failure_cdf = nullptr,
    double ssd_arr = 0.0);

void update_failure_info(
    const std::string& event_type,
    const std::string& event_node,
    std::map<int, FailureInfo>& failure_info_per_ssd_group,
    std::map<std::string, bool>& failed_nodes_and_enclosures,
    const SSDRedundancyScheme& ssd_redun_scheme);

// Flow and speed calculation
void calculate_flows_and_speed(
    const GraphStructure& hardware_graph_copy,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const nlohmann::json& options,
    std::map<FailureStateKey, FlowsAndSpeedEntry>& flows_and_speed_table,
    double max_read_performance_without_any_failure,
    const DisconnectedStatus& disconnected,
    const FailureStateKey& key);

// SSD state update functions
void update_all_ssd_states(
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    std::map<int, SSDInfo>& SSDs,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const DisconnectedStatus& disconnected,
    double current_time);

void update_ssd_state(
    const std::string& ssd_name,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    std::map<int, SSDInfo>& SSDs,
    double capacity,
    const std::string& event_type,
    double prep_time_for_rebuilding,
    const SSDRedundancyScheme& ssd_redun_scheme,
    const DisconnectedStatus& disconnected);

void update_failure_event_for_SSDs(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
    double current_time,
    std::map<int, SSDInfo>& SSDs);

void update_repair_event_for_SSDs(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events,
    double current_time,
    std::map<int, SSDInfo>& SSDs,
    const FlowsAndSpeedEntry& flows_and_speed_entry,
    const std::map<int, FailureInfo>& failure_info_per_ssd_group,
    const SSDRedundancyScheme& ssd_redun_scheme);

// Main simulation function
void monte_carlo_simulation(
    std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int num_simulations,
    const nlohmann::json& options,
    const std::map<std::string, double>& costs
);
