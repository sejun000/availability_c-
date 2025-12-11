#include "erasure_coding.hpp"
#include <algorithm>
#include <cmath>

// Parse EC type from string
ECType parse_ec_type(const std::string& type_str) {
    if (type_str == "standard" || type_str == "ec" || type_str == "EC") {
        return ECType::STANDARD;
    } else if (type_str == "lrc" || type_str == "LRC") {
        return ECType::LRC;
    } else if (type_str == "multi" || type_str == "multi_ec" || type_str == "MULTI_EC") {
        return ECType::MULTI_EC;
    }
    return ECType::STANDARD;  // Default
}

// Convert EC type to string
std::string ec_type_to_string(ECType type) {
    switch (type) {
        case ECType::STANDARD: return "standard";
        case ECType::LRC: return "lrc";
        case ECType::MULTI_EC: return "multi_ec";
    }
    return "unknown";
}

// Calculate data loss for given failed disk indices
double ECConfig::calculate_data_loss(const std::set<int>& failed_disk_indices) const {
    int failure_count = static_cast<int>(failed_disk_indices.size());

    if (failure_count == 0) {
        return 0.0;
    }

    switch (type) {
        case ECType::STANDARD: {
            // Standard EC: data loss if more than k failures
            if (failure_count <= k) {
                return 0.0;  // Can recover
            }
            // Data loss: fraction of data affected
            // In declustered parity, each disk has 1/n of the data
            // With more than k failures, we lose data proportional to failures
            return static_cast<double>(failure_count) / n;
        }

        case ECType::LRC: {
            // LRC: can tolerate k global + local_k per group
            // First check if we can recover locally
            int num_local_groups = m / local_m;

            // Count failures per local group
            std::vector<int> local_failures(num_local_groups, 0);
            for (int idx : failed_disk_indices) {
                int local_group = idx / (local_m + local_k);
                if (local_group < num_local_groups) {
                    local_failures[local_group]++;
                }
            }

            // Check each local group
            int groups_needing_global = 0;
            for (int i = 0; i < num_local_groups; ++i) {
                if (local_failures[i] > local_k) {
                    groups_needing_global++;
                }
            }

            // If more groups need global recovery than k allows, data loss
            if (groups_needing_global > k) {
                return static_cast<double>(failure_count) / n;
            }
            return 0.0;
        }

        case ECType::MULTI_EC: {
            // Multi-EC: local stripe failures count as one "chunk" failure
            int local_stripe_size = local_m + local_k;
            int num_local_stripes = m + k;

            // Count failed local stripes (stripes with more than local_k failures)
            std::vector<int> stripe_failures(num_local_stripes, 0);
            for (int idx : failed_disk_indices) {
                int stripe_idx = idx / local_stripe_size;
                if (stripe_idx < num_local_stripes) {
                    stripe_failures[stripe_idx]++;
                }
            }

            int failed_stripes = 0;
            for (int i = 0; i < num_local_stripes; ++i) {
                if (stripe_failures[i] > local_k) {
                    failed_stripes++;
                }
            }

            // If more than k stripes fail, data loss
            if (failed_stripes > k) {
                return static_cast<double>(failure_count) / total_disks;
            }
            return 0.0;
        }
    }
    return 0.0;
}

// Get number of disks needed for rebuild read
int ECConfig::get_rebuild_read_disk_count(const std::set<int>& failed_disk_indices, int target_disk) const {
    switch (type) {
        case ECType::STANDARD:
            // Need m disks to rebuild (data chunks)
            return m;

        case ECType::LRC: {
            // Check if local rebuild is possible
            int local_group = target_disk / (local_m + local_k);
            int local_start = local_group * (local_m + local_k);
            int local_end = local_start + local_m + local_k;

            int local_failures = 0;
            for (int idx : failed_disk_indices) {
                if (idx >= local_start && idx < local_end) {
                    local_failures++;
                }
            }

            if (local_failures <= local_k) {
                // Can rebuild locally, need local_m disks
                return local_m;
            }
            // Need global rebuild, need m disks
            return m;
        }

        case ECType::MULTI_EC: {
            // Check if local rebuild is possible within the stripe
            int local_stripe_size = local_m + local_k;
            int stripe_idx = target_disk / local_stripe_size;
            int stripe_start = stripe_idx * local_stripe_size;
            int stripe_end = stripe_start + local_stripe_size;

            int stripe_failures = 0;
            for (int idx : failed_disk_indices) {
                if (idx >= stripe_start && idx < stripe_end) {
                    stripe_failures++;
                }
            }

            if (stripe_failures <= local_k) {
                // Can rebuild locally within stripe
                return local_m;
            }
            // Need to read from other stripes (m stripes worth of data)
            return m * local_m;
        }
    }
    return m;  // Default
}

// ErasureCodingScheme implementation

void ErasureCodingScheme::initialize(const ECConfig& config) {
    config_ = config;
    if (!config_.validate()) {
        throw std::runtime_error("Invalid EC configuration");
    }
}

int ErasureCodingScheme::get_disk_group(int disk_index) const {
    int group_size = get_group_size();
    return disk_index / group_size;
}

std::vector<int> ErasureCodingScheme::get_stripe_disks(int disk_index) const {
    std::vector<int> disks;
    int group_size = get_group_size();
    int group_start = (disk_index / group_size) * group_size;

    for (int i = 0; i < group_size; ++i) {
        disks.push_back(group_start + i);
    }
    return disks;
}

bool ErasureCodingScheme::can_rebuild(const std::set<int>& failed_disks) const {
    // Can rebuild if data loss is 0
    return config_.calculate_data_loss(failed_disks) == 0.0;
}

double ErasureCodingScheme::get_data_loss_ratio(const std::set<int>& failed_disks) const {
    return config_.calculate_data_loss(failed_disks);
}

int ErasureCodingScheme::get_rebuild_read_count(const std::set<int>& failed_disks, int target_disk) const {
    return config_.get_rebuild_read_disk_count(failed_disks, target_disk);
}

bool ErasureCodingScheme::can_local_rebuild(const std::set<int>& failed_disks, int target_disk) const {
    if (config_.type == ECType::LRC) {
        int local_group = target_disk / (config_.local_m + config_.local_k);
        int local_start = local_group * (config_.local_m + config_.local_k);
        int local_end = local_start + config_.local_m + config_.local_k;

        int local_failures = 0;
        for (int idx : failed_disks) {
            if (idx >= local_start && idx < local_end) {
                local_failures++;
            }
        }
        return local_failures <= config_.local_k;
    } else if (config_.type == ECType::MULTI_EC) {
        int local_stripe_size = config_.local_m + config_.local_k;
        int stripe_idx = target_disk / local_stripe_size;
        int stripe_start = stripe_idx * local_stripe_size;
        int stripe_end = stripe_start + local_stripe_size;

        int stripe_failures = 0;
        for (int idx : failed_disks) {
            if (idx >= stripe_start && idx < stripe_end) {
                stripe_failures++;
            }
        }
        return stripe_failures <= config_.local_k;
    }
    return false;  // Standard EC has no local parity
}

int ErasureCodingScheme::get_total_groups(int total_disks) const {
    int group_size = get_group_size();
    return (total_disks + group_size - 1) / group_size;
}

int ErasureCodingScheme::get_group_start_disk(int group_index) const {
    return group_index * get_group_size();
}
