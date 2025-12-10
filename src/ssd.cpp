#include "ssd.hpp"

// SSDIOModuleManager implementation
void SSDIOModuleManager::initialize(const std::vector<SSDIOModuleMapping>& mappings,
                                     int total_ssds,
                                     const std::set<std::string>& all_io_modules) {
    total_ssds_ = total_ssds;
    all_io_modules_ = all_io_modules;

    if (mappings.empty()) {
        // Legacy mode: all SSDs connect to all io_modules
        legacy_mode_ = true;
        mappings_.clear();
        return;
    }

    legacy_mode_ = false;
    mappings_.clear();

    // Compute start_ssd_index for each mapping
    int current_index = 0;
    for (const auto& m : mappings) {
        SSDIOModuleMapping entry;
        entry.ssd_count = m.ssd_count;
        entry.io_modules = m.io_modules;
        entry.start_ssd_index = current_index;
        mappings_.push_back(entry);
        current_index += m.ssd_count;
    }
}

bool SSDIOModuleManager::is_ssd_disconnected(int ssd_index,
                                              const std::map<std::string, bool>& failed_nodes) const {
    if (legacy_mode_) {
        // Legacy: all io_modules must be failed for SSD to be disconnected
        for (const auto& io_module : all_io_modules_) {
            auto it = failed_nodes.find(io_module);
            if (it == failed_nodes.end() || !it->second) {
                // At least one io_module is up
                return false;
            }
        }
        return true;  // All io_modules failed
    }

    // Find which mapping this SSD belongs to
    for (const auto& mapping : mappings_) {
        int end_index = mapping.start_ssd_index + mapping.ssd_count;
        if (ssd_index >= mapping.start_ssd_index && ssd_index < end_index) {
            // Check if ALL io_modules for this SSD are failed
            for (const auto& io_module : mapping.io_modules) {
                auto it = failed_nodes.find(io_module);
                if (it == failed_nodes.end() || !it->second) {
                    // At least one io_module is up
                    return false;
                }
            }
            return true;  // All connected io_modules failed
        }
    }

    // SSD not found in any mapping - assume connected to all (fallback)
    for (const auto& io_module : all_io_modules_) {
        auto it = failed_nodes.find(io_module);
        if (it == failed_nodes.end() || !it->second) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> SSDIOModuleManager::get_io_modules_for_ssd(int ssd_index) const {
    if (legacy_mode_) {
        return std::vector<std::string>(all_io_modules_.begin(), all_io_modules_.end());
    }

    for (const auto& mapping : mappings_) {
        int end_index = mapping.start_ssd_index + mapping.ssd_count;
        if (ssd_index >= mapping.start_ssd_index && ssd_index < end_index) {
            return mapping.io_modules;
        }
    }

    // Fallback: return all io_modules
    return std::vector<std::string>(all_io_modules_.begin(), all_io_modules_.end());
}

std::string get_ssd_name(int ssd_index) {
    return SSD_MODULE_NAME + std::to_string(ssd_index);
}

int get_ssd_index(const std::string& ssd_name) {
    return std::stoi(ssd_name.substr(SSD_MODULE_NAME.length()));
}

std::string get_ssd_group_name(int ssd_group_index) {
    return SSD_MODULE_NAME + "_group_" + std::to_string(ssd_group_index);
}

int get_ssd_group_index_from_group_name(const std::string& ssd_group_name) {
    std::string prefix = SSD_MODULE_NAME + "_group_";
    return std::stoi(ssd_group_name.substr(prefix.length()));
}

bool is_event_node_ssd(const std::string& event_node) {
    return event_node.find(SSD_MODULE_NAME) == 0;
}

SSDRedundancyScheme::SSDRedundancyScheme(
    double write_bw, double read_bw, double mttf,
    double cached_write_ratio, double cached_write_bw, double cached_read_bw, double cached_mttf,
    int m, int k, int l,
    int cached_m, int cached_k, int cached_l,
    int network_m, int network_k, int network_l,
    int cached_network_m, int cached_network_k, int cached_network_l,
    int cached_ssds, int total_ssds,
    int inter_replicas, int intra_replicas
) : write_bw_(write_bw), read_bw_(read_bw), mttf_(mttf),
    cached_write_ratio_(cached_write_ratio),
    cached_write_bw_(cached_write_bw), cached_read_bw_(cached_read_bw), cached_mttf_(cached_mttf),
    m_(m), k_(k), l_(l),
    cached_m_(cached_m), cached_k_(cached_k), cached_l_(cached_l),
    network_m_(network_m), network_k_(network_k), network_l_(network_l),
    cached_network_m_(cached_network_m), cached_network_k_(cached_network_k), cached_network_l_(cached_network_l),
    cached_ssds_(cached_ssds), total_ssds_(total_ssds),
    inter_replicas_(inter_replicas), intra_replicas_(intra_replicas) {

    ssd_group_size_ = m + k + l;
    cached_ssd_group_size_ = cached_m + cached_k + cached_l;
}

bool SSDRedundancyScheme::is_ssd_index_cached(int ssd_index) const {
    return ssd_index >= total_ssds_ - cached_ssds_;
}

bool SSDRedundancyScheme::is_ssd_group_index_cached(int ssd_group_index) const {
    return ssd_group_index >= (total_ssds_ - cached_ssds_) / ssd_group_size_;
}

int SSDRedundancyScheme::get_start_ssd_index(int ssd_group_index) const {
    if (!is_ssd_group_index_cached(ssd_group_index)) {
        return ssd_group_index * ssd_group_size_;
    } else {
        int cached_group_start_index = (total_ssds_ - cached_ssds_) / ssd_group_size_;
        return (total_ssds_ - cached_ssds_) +
               (ssd_group_index - cached_group_start_index) * cached_ssd_group_size_;
    }
}

int SSDRedundancyScheme::get_total_group_count() const {
    if (cached_ssds_ == 0) {
        return total_ssds_ / ssd_group_size_;
    }
    return (total_ssds_ - cached_ssds_) / ssd_group_size_ +
           cached_ssds_ / cached_ssd_group_size_;
}

int SSDRedundancyScheme::get_ssd_group_index(int ssd_index) const {
    if (!is_ssd_index_cached(ssd_index)) {
        return ssd_index / ssd_group_size_;
    } else {
        int cached_ssd_group_start_index = (total_ssds_ - cached_ssds_) / ssd_group_size_;
        return (ssd_index - (total_ssds_ - cached_ssds_)) / cached_ssd_group_size_ +
               cached_ssd_group_start_index;
    }
}

double SSDRedundancyScheme::get_read_bw(bool cached) const {
    return cached ? cached_read_bw_ : read_bw_;
}

double SSDRedundancyScheme::get_write_bw(bool cached) const {
    return cached ? cached_write_bw_ : write_bw_;
}

double SSDRedundancyScheme::get_mttf(bool cached) const {
    return cached ? cached_mttf_ : mttf_;
}

int SSDRedundancyScheme::get_m(bool cached) const {
    return cached ? cached_m_ : m_;
}

int SSDRedundancyScheme::get_k(bool cached) const {
    return cached ? cached_k_ : k_;
}

int SSDRedundancyScheme::get_l(bool cached) const {
    return cached ? cached_l_ : l_;
}

int SSDRedundancyScheme::get_network_m(bool cached) const {
    return cached ? cached_network_m_ : network_m_;
}

int SSDRedundancyScheme::get_network_k(bool cached) const {
    return cached ? cached_network_k_ : network_k_;
}

int SSDRedundancyScheme::get_network_l(bool cached) const {
    return cached ? cached_network_l_ : network_l_;
}

int SSDRedundancyScheme::get_total_ssds() const {
    return total_ssds_;
}

int SSDRedundancyScheme::get_tiered_ssds(bool cached) const {
    return cached ? cached_ssds_ : total_ssds_ - cached_ssds_;
}

std::string SSDRedundancyScheme::get_cached_prefix(bool cached) const {
    return cached ? "cached_" : "";
}

int SSDRedundancyScheme::get_inter_replicas(bool cached) const {
    return cached ? 0 : inter_replicas_;
}

int SSDRedundancyScheme::get_intra_replicas(bool cached) const {
    return cached ? intra_replicas_ : 0;
}
