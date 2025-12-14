#pragma once

#include <map>
#include <vector>
#include <string>
#include <queue>
#include <random>
#include <cstdint>
#include <nlohmann/json.hpp>
#include "graph_structure.hpp"
#include "disk.hpp"
#include "erasure_coding.hpp"
#include "utils.hpp"

constexpr double BIG_NUMBER = 1e15;
constexpr double SOFT_ERROR_THRESHOLD_HOURS = 1.0;

// Disk States (for state machine)
const std::string DISK_STATE_NORMAL = "normal";
const std::string DISK_STATE_DISCONNECTED = "disconnected";
const std::string DISK_STATE_FAILED = "failed";
const std::string DISK_STATE_REBUILDING = "rebuilding";
const std::string DISK_STATE_DATA_LOSS = "data_loss";

// Event structure for priority queue
struct Event {
    double time;
    std::string event_type;  // "fail", "repair", "rebuild_complete"
    std::string event_node;
    double prev_time;
    double start_time;
    bool software_repair_only;  // For io_module: software fault (quick repair) vs hardware fault

    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// Failure information per EC group
struct ECGroupFailureInfo {
    int disk_failure_count = 0;          // Number of physically failed disks
    int disconnected_count = 0;          // Number of disconnected disks
    std::set<int> failed_disk_indices;   // Indices of failed disks
    std::set<int> disconnected_disk_indices;  // Indices of disconnected disks

    // Get total unavailable disks
    int get_unavailable_count() const {
        std::set<int> all_unavailable = failed_disk_indices;
        all_unavailable.insert(disconnected_disk_indices.begin(), disconnected_disk_indices.end());
        return static_cast<int>(all_unavailable.size());
    }

    // Check if specific disk is unavailable
    bool is_disk_unavailable(int disk_index) const {
        return failed_disk_indices.count(disk_index) > 0 ||
               disconnected_disk_indices.count(disk_index) > 0;
    }
};

// Disconnection status for tracking which components are down
struct DisconnectedStatus {
    std::map<std::string, bool> failed_nodes;  // node_name -> is_failed
    std::map<int, bool> disk_disconnected;     // disk_index -> is_disconnected (due to upstream failure)

    DisconnectedStatus() = default;

    // Get count of disconnected disks in a group
    int get_disconnected_count_in_group(int start_disk_index, int group_size) const {
        int count = 0;
        for (int i = start_disk_index; i < start_disk_index + group_size; ++i) {
            auto it = disk_disconnected.find(i);
            if (it != disk_disconnected.end() && it->second) {
                count++;
            }
        }
        return count;
    }

    // Check if specific disk is disconnected
    bool is_disk_disconnected(int disk_index) const {
        auto it = disk_disconnected.find(disk_index);
        return it != disk_disconnected.end() && it->second;
    }

    // Mark disk as disconnected
    void set_disk_disconnected(int disk_index, bool disconnected) {
        disk_disconnected[disk_index] = disconnected;
    }
};

// Flow and speed calculation result
struct FlowsAndSpeedEntry {
    std::map<int, double> rebuild_bandwidth;  // failure_count -> rebuild_bandwidth
    std::map<int, double> rebuild_speed_by_disk;  // disk_index -> rebuild_speed (bytes/sec)
    double read_bandwidth = 0.0;              // Available read bandwidth
    double write_bandwidth = 0.0;             // Available write bandwidth
    double availability_ratio = 1.0;          // Current availability (0 if data loss)
    double data_loss_ratio = 0.0;             // Fraction of data lost
    double performance_ratio_read = 1.0;      // read_bw / max_read_bw_no_failure
    double performance_ratio_write = 1.0;     // write_bw / max_write_bw_no_failure
    double performance_ratio_min = 1.0;       // min(read,write) ratio
    int groups_with_data_loss = 0;            // Number of EC groups with data loss
    std::vector<int> data_loss_group_indices; // Indices of EC groups with data loss

    FlowsAndSpeedEntry() = default;
};

// Key structures for caching
struct NodeFailureKey {
    std::vector<uint8_t> failed_flags;

    bool operator<(const NodeFailureKey& other) const {
        return failed_flags < other.failed_flags;
    }

    bool operator==(const NodeFailureKey& other) const {
        return failed_flags == other.failed_flags;
    }
};

struct FailureStateKey {
    NodeFailureKey node_key;
    std::vector<uint8_t> disk_state;  // Per disk: 0=normal, 1=unavailable, 2=rebuilding

    bool operator<(const FailureStateKey& other) const {
        if (node_key < other.node_key) return true;
        if (other.node_key < node_key) return false;
        return disk_state < other.disk_state;
    }

    bool operator==(const FailureStateKey& other) const {
        return node_key == other.node_key && disk_state == other.disk_state;
    }
};

// Hash functions for unordered_map
struct NodeFailureKeyHash {
    size_t operator()(const NodeFailureKey& key) const {
        size_t hash = 0;
        for (uint8_t flag : key.failed_flags) {
            hash ^= std::hash<uint8_t>{}(flag) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

struct FailureStateKeyHash {
    size_t operator()(const FailureStateKey& key) const {
        size_t hash = NodeFailureKeyHash{}(key.node_key);
        for (uint8_t state : key.disk_state) {
            hash ^= std::hash<uint8_t>{}(state) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// Simulation parameters
struct SimulationParams {
    // EC configuration
    ECConfig ec_config;

    // Disk parameters
    int total_disks = 0;
    double disk_capacity = 0;        // Capacity per disk in bytes
    double disk_read_bw = 0;         // Read bandwidth per disk
    double disk_write_bw = 0;        // Write bandwidth per disk
    double disk_mttf = 0;            // Mean time to failure for disks
    double disk_replace_time = 0;    // Time to replace a failed disk (hours)

    // Rebuild parameters
    double rebuild_bw_ratio = 0.2;   // Fraction of bandwidth used for rebuild
    double degraded_ratio = 0.2;     // Fraction of bandwidth for degraded read
    double data_loss_rebuild_bw = 0; // Fixed rebuild bandwidth during data loss (bytes/sec), 0 = disabled

    // Performance availability threshold
    double target_performance_ratio = 0.0;  // If > 0, calculate perf_availability

    // Simulation control
    int nprocs = 40;
    int simulation_years = 10;

    // EC encoding speed (bytes/sec for m+k encoding)
    double ec_encoding_speed = 0;

    // Number of EC groups (for durability calculation)
    int num_ec_groups = 1;
};

// Simulation result structure
struct SimulationResult {
    int runs = 0;                           // Number of independent simulations aggregated
    double up_time = 0.0;
    double simulation_time = 0.0;
    double availability = 0.0;
    double perf_up_time_read = 0.0;         // Time with read_ratio >= target
    double perf_up_time_write = 0.0;        // Time with write_ratio >= target
    double perf_up_time_min = 0.0;          // Time with min(read,write) >= target
    double perf_availability_read = 0.0;
    double perf_availability_write = 0.0;
    double perf_availability_min = 0.0;
    double mttdl = 0.0;                     // Mean time to data loss
    int data_loss_events = 0;               // Number of data loss events
    int disk_failure_events = 0;            // Total number of disk failures
    double total_data_loss_ratio = 0.0;     // Total fraction of data lost
    double durability = 0.0;                // 1 - (data_loss_ratio per year)
    double time_for_rebuilding = 0.0;
    int count_for_rebuilding = 0;
    double average_read_bandwidth = 0.0;
    double average_write_bandwidth = 0.0;
    // Rebuild-specific metrics
    double total_rebuild_speed_time = 0.0;      // Sum of (rebuild_speed * duration)
    double total_time_in_rebuild = 0.0;         // Total time any disk is rebuilding
    double total_host_read_bw_during_rebuild = 0.0;   // Sum of (read_bw * duration) during rebuild
    double total_host_write_bw_during_rebuild = 0.0;  // Sum of (write_bw * duration) during rebuild
    // Year-based durability metrics
    int total_group_years = 0;                  // Total EC group-years simulated
    int group_years_with_data_loss = 0;         // EC group-years that had data loss
};

// Simulation helper functions
double pfail(int m, int k, int l, int x);
double calculate_bottleneck_speed(int m, int k, const std::vector<double>& other_bws,
                                   double rebuild_bw_ratio, double ec_encoding_speed);

// State management
std::string judge_state_from_failure_info(const ECGroupFailureInfo& failure_info,
                                          const ErasureCodingScheme& ec_scheme);

bool is_data_loss(const ECGroupFailureInfo& failure_info, const ErasureCodingScheme& ec_scheme);

// Event management
void push_failed_event(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
    const std::string& event_node,
    double current_time,
    const std::map<std::string, std::string>& node_to_module_map,
    const GraphStructure& hardware_graph,
    double disk_mttf,
    const EmpiricalCDF* disk_failure_cdf = nullptr,
    double disk_arr = 0.0);

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
    int total_disks,
    double disk_mttf,
    const EmpiricalCDF* disk_failure_cdf = nullptr,
    double disk_arr = 0.0);

void update_failure_info(
    const std::string& event_type,
    const std::string& event_node,
    std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    std::map<std::string, bool>& failed_nodes_and_enclosures,
    const ErasureCodingScheme& ec_scheme);

// Flow and speed calculation
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
    const DiskIOModuleManager* disk_io_manager = nullptr,
    const std::map<int, DiskInfo>* disks = nullptr);

// Disk state update functions
void update_all_disk_states(
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    std::map<int, DiskInfo>& disks,
    const ErasureCodingScheme& ec_scheme,
    const DisconnectedStatus& disconnected,
    double current_time);

void update_disk_state(
    int disk_index,
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    std::map<int, DiskInfo>& disks,
    double capacity,
    const std::string& event_type,
    double event_time,
    double replace_time,
    const ErasureCodingScheme& ec_scheme,
    const DisconnectedStatus& disconnected);

void update_failure_event_for_disks(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& failed_events,
    double current_time,
    std::map<int, DiskInfo>& disks);

void update_repair_event_for_disks(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& repair_events,
    double current_time,
    std::map<int, DiskInfo>& disks,
    const FlowsAndSpeedEntry& flows_and_speed_entry,
    const std::map<int, ECGroupFailureInfo>& failure_info_per_ec_group,
    const ErasureCodingScheme& ec_scheme);

// Main simulation function
void monte_carlo_simulation(
    std::map<std::string, nlohmann::json>& params_and_results,
    const GraphStructure& graph_structure_origin,
    int num_simulations,
    const nlohmann::json& options
);
