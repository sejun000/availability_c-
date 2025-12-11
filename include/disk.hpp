#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>

const std::string DISK_MODULE_NAME = "disk";

// Disk state enumeration
enum class DiskState {
    NORMAL,             // Disk is healthy and operational
    DISCONNECTED,       // Disk is inaccessible due to upstream failure (io_module, switch, etc.)
    FAILED,             // Disk itself has failed, needs replacement
    REBUILDING,         // Disk is being rebuilt (after failure or reconnection)
    DEGRADED            // Disk is part of a degraded EC group (other disk failed)
};

// Convert DiskState to string
std::string disk_state_to_string(DiskState state);

// Disk information structure
struct DiskInfo {
    int index;                              // Disk index (0-based)
    std::string name;                       // Disk name (e.g., "disk_0")
    std::string io_module;                  // Connected io_module name

    DiskState state = DiskState::NORMAL;

    // Failure/disconnect tracking
    bool is_failed = false;                 // True if disk itself failed
    bool is_disconnected = false;           // True if disconnected due to upstream failure
    double disconnect_timestamp = 0.0;      // When disconnection started
    double fail_timestamp = 0.0;            // When failure occurred

    // Rebuild tracking
    double remaining_capacity_to_rebuild = 0.0;  // Bytes remaining to rebuild
    double rebuild_speed = 0.0;                  // Current rebuild speed (bytes/sec)
    double remaining_replace_time = 0.0;         // Time until replacement disk arrives (for failures)

    // Performance
    double read_bandwidth = 0.0;            // Read bandwidth (bytes/sec)
    double write_bandwidth = 0.0;           // Write bandwidth (bytes/sec)

    DiskInfo() = default;
    DiskInfo(int idx, const std::string& io_mod, double read_bw, double write_bw);

    // Get effective state considering both failure and disconnection
    DiskState get_effective_state() const;

    // Check if disk needs rebuild
    bool needs_rebuild() const;

    // Get remaining rebuild time at current speed
    double get_remaining_rebuild_time() const;
};

// Disk naming utilities
std::string get_disk_name(int disk_index);
int get_disk_index(const std::string& disk_name);
bool is_disk_node(const std::string& node_name);

// Disk to IO Module mapping entry
struct DiskIOModuleMapping {
    int disk_count;                         // Number of disks in this group
    std::vector<std::string> io_modules;    // IO modules these disks are connected to
    int start_disk_index = 0;               // Starting disk index (computed)
};

// Manages Disk to IO Module connectivity
class DiskIOModuleManager {
public:
    DiskIOModuleManager() = default;

    // Initialize from JSON config
    void initialize(const std::vector<DiskIOModuleMapping>& mappings,
                    int total_disks,
                    const std::set<std::string>& all_io_modules);

    // Check if disk is disconnected given failed nodes
    // Returns true if ALL io_modules for this disk have failed
    bool is_disk_disconnected(int disk_index,
                              const std::map<std::string, bool>& failed_nodes) const;

    // Get list of io_modules that this disk is connected to
    std::vector<std::string> get_io_modules_for_disk(int disk_index) const;

    // Get primary io_module for a disk (first in the list)
    std::string get_primary_io_module(int disk_index) const;

    // Check if EC group crosses io_module boundary
    bool does_ec_group_cross_io_module(int start_disk_index, int group_size) const;

    // Get set of io_modules used by an EC group
    std::set<std::string> get_io_modules_for_ec_group(int start_disk_index, int group_size) const;

    // Check if legacy mode (all disks connect to all io_modules)
    bool is_legacy_mode() const { return legacy_mode_; }

    // Get all io_modules in the system
    const std::set<std::string>& get_all_io_modules() const { return all_io_modules_; }

    // Get total disk count
    int get_total_disks() const { return total_disks_; }

    // Get active port ratio for a disk (considering failed io_modules)
    // Returns (active_ports / total_ports), e.g., 0.5 if 1 of 2 ports failed
    double get_active_port_ratio(int disk_index,
                                  const std::map<std::string, bool>& failed_nodes) const;

    // Get effective bandwidth for a disk considering port failures
    double get_effective_bandwidth(int disk_index, double full_bandwidth,
                                    const std::map<std::string, bool>& failed_nodes) const;

private:
    std::vector<DiskIOModuleMapping> mappings_;
    std::set<std::string> all_io_modules_;
    int total_disks_ = 0;
    bool legacy_mode_ = true;
};
