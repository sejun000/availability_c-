#include "disk.hpp"
#include <sstream>
#include <regex>

// Convert DiskState to string
std::string disk_state_to_string(DiskState state) {
    switch (state) {
        case DiskState::NORMAL: return "normal";
        case DiskState::DISCONNECTED: return "disconnected";
        case DiskState::FAILED: return "failed";
        case DiskState::REBUILDING: return "rebuilding";
        case DiskState::DEGRADED: return "degraded";
    }
    return "unknown";
}

// DiskInfo implementation
DiskInfo::DiskInfo(int idx, const std::string& io_mod, double read_bw, double write_bw)
    : index(idx)
    , name(get_disk_name(idx))
    , io_module(io_mod)
    , state(DiskState::NORMAL)
    , is_failed(false)
    , is_disconnected(false)
    , disconnect_timestamp(0.0)
    , fail_timestamp(0.0)
    , remaining_capacity_to_rebuild(0.0)
    , rebuild_speed(0.0)
    , remaining_replace_time(0.0)
    , read_bandwidth(read_bw)
    , write_bandwidth(write_bw)
{}

DiskState DiskInfo::get_effective_state() const {
    if (is_failed) {
        if (remaining_capacity_to_rebuild > 0 || remaining_replace_time > 0) {
            return DiskState::REBUILDING;
        }
        return DiskState::FAILED;
    }
    if (is_disconnected) {
        return DiskState::DISCONNECTED;
    }
    return state;
}

bool DiskInfo::needs_rebuild() const {
    return (is_failed || is_disconnected) && remaining_capacity_to_rebuild > 0;
}

double DiskInfo::get_remaining_rebuild_time() const {
    if (rebuild_speed <= 0) {
        return 1e15;  // Infinite time
    }
    return remaining_replace_time + (remaining_capacity_to_rebuild / rebuild_speed / 3600.0);
}

// Disk naming utilities
std::string get_disk_name(int disk_index) {
    return DISK_MODULE_NAME + "_" + std::to_string(disk_index);
}

int get_disk_index(const std::string& disk_name) {
    // Extract number from "disk_N" format
    size_t underscore_pos = disk_name.rfind('_');
    if (underscore_pos != std::string::npos) {
        try {
            return std::stoi(disk_name.substr(underscore_pos + 1));
        } catch (...) {
            return -1;
        }
    }
    return -1;
}

bool is_disk_node(const std::string& node_name) {
    return node_name.find(DISK_MODULE_NAME + "_") == 0;
}

// DiskIOModuleManager implementation
void DiskIOModuleManager::initialize(
    const std::vector<DiskIOModuleMapping>& mappings,
    int total_disks,
    const std::set<std::string>& all_io_modules) {

    mappings_ = mappings;
    total_disks_ = total_disks;
    all_io_modules_ = all_io_modules;

    if (mappings_.empty()) {
        // Legacy mode: all disks connect to all io_modules
        legacy_mode_ = true;
        return;
    }

    legacy_mode_ = false;

    // Compute start indices for each mapping group
    int current_index = 0;
    for (auto& mapping : mappings_) {
        mapping.start_disk_index = current_index;
        current_index += mapping.disk_count;
    }
}

bool DiskIOModuleManager::is_disk_disconnected(
    int disk_index,
    const std::map<std::string, bool>& failed_nodes) const {

    if (legacy_mode_) {
        // In legacy mode, disk is disconnected only if ALL io_modules have failed
        for (const auto& io_mod : all_io_modules_) {
            auto it = failed_nodes.find(io_mod);
            if (it == failed_nodes.end() || !it->second) {
                return false;  // At least one io_module is up
            }
        }
        return true;  // All io_modules are down
    }

    // Find the mapping for this disk
    std::vector<std::string> io_modules = get_io_modules_for_disk(disk_index);

    if (io_modules.empty()) {
        return true;  // No io_module = disconnected
    }

    // Check if ALL io_modules for this disk have failed
    for (const auto& io_mod : io_modules) {
        auto it = failed_nodes.find(io_mod);
        if (it == failed_nodes.end() || !it->second) {
            return false;  // At least one io_module is up
        }
    }
    return true;  // All io_modules are down
}

std::vector<std::string> DiskIOModuleManager::get_io_modules_for_disk(int disk_index) const {
    if (legacy_mode_) {
        return std::vector<std::string>(all_io_modules_.begin(), all_io_modules_.end());
    }

    for (const auto& mapping : mappings_) {
        if (disk_index >= mapping.start_disk_index &&
            disk_index < mapping.start_disk_index + mapping.disk_count) {
            return mapping.io_modules;
        }
    }
    return {};
}

std::string DiskIOModuleManager::get_primary_io_module(int disk_index) const {
    auto io_modules = get_io_modules_for_disk(disk_index);
    if (io_modules.empty()) {
        return "";
    }
    return io_modules[0];
}

bool DiskIOModuleManager::does_ec_group_cross_io_module(int start_disk_index, int group_size) const {
    if (legacy_mode_) {
        return false;  // All disks connect to all io_modules
    }

    std::set<std::string> io_module_set;

    for (int i = start_disk_index; i < start_disk_index + group_size; ++i) {
        auto io_modules = get_io_modules_for_disk(i);
        std::set<std::string> disk_io_set(io_modules.begin(), io_modules.end());

        if (io_module_set.empty()) {
            io_module_set = disk_io_set;
        } else if (io_module_set != disk_io_set) {
            return true;  // Different io_module sets
        }
    }
    return false;
}

std::set<std::string> DiskIOModuleManager::get_io_modules_for_ec_group(
    int start_disk_index, int group_size) const {

    std::set<std::string> result;

    for (int i = start_disk_index; i < start_disk_index + group_size; ++i) {
        auto io_modules = get_io_modules_for_disk(i);
        result.insert(io_modules.begin(), io_modules.end());
    }
    return result;
}

double DiskIOModuleManager::get_active_port_ratio(
    int disk_index,
    const std::map<std::string, bool>& failed_nodes) const {

    std::vector<std::string> io_modules = get_io_modules_for_disk(disk_index);

    if (io_modules.empty()) {
        return 0.0;  // No ports = no bandwidth
    }

    int total_ports = static_cast<int>(io_modules.size());
    int active_ports = 0;

    for (const auto& io_mod : io_modules) {
        auto it = failed_nodes.find(io_mod);
        if (it == failed_nodes.end() || !it->second) {
            active_ports++;  // This io_module is up
        }
    }

    if (active_ports == 0) {
        return 0.0;  // All ports down = disconnected
    }

    return static_cast<double>(active_ports) / static_cast<double>(total_ports);
}

double DiskIOModuleManager::get_effective_bandwidth(
    int disk_index, double full_bandwidth,
    const std::map<std::string, bool>& failed_nodes) const {

    double port_ratio = get_active_port_ratio(disk_index, failed_nodes);
    return full_bandwidth * port_ratio;
}
