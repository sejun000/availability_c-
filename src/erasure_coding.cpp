#include "erasure_coding.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

// Calculate combinations C(n, k)
static double combinations(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    if (k > n - k) k = n - k;  // Optimization: C(n,k) = C(n,n-k)

    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result *= static_cast<double>(n - i) / static_cast<double>(i + 1);
    }
    return result;
}

// Parse EC type from string
ECType parse_ec_type(const std::string& type_str) {
    if (type_str == "standard" || type_str == "ec" || type_str == "EC") {
        return ECType::STANDARD;
    } else if (type_str == "lrc" || type_str == "LRC") {
        return ECType::LRC;
    } else if (type_str == "multi" || type_str == "multi_ec" || type_str == "MULTI_EC") {
        return ECType::MULTI_EC;
    } else if (type_str == "replication" || type_str == "REPLICATION" || type_str == "rep") {
        return ECType::REPLICATION;
    }
    return ECType::STANDARD;  // Default
}

// Convert EC type to string
std::string ec_type_to_string(ECType type) {
    switch (type) {
        case ECType::STANDARD: return "standard";
        case ECType::LRC: return "lrc";
        case ECType::MULTI_EC: return "multi_ec";
        case ECType::REPLICATION: return "replication";
    }
    return "unknown";
}

// Parse local type from string
LocalType parse_local_type(const std::string& type_str) {
    if (type_str == "ec" || type_str == "EC" || type_str == "erasure_coding") {
        return LocalType::EC;
    } else if (type_str == "replication" || type_str == "REPLICATION" || type_str == "rep") {
        return LocalType::REPLICATION;
    }
    return LocalType::EC;  // Default
}

// Convert local type to string
std::string local_type_to_string(LocalType type) {
    switch (type) {
        case LocalType::EC: return "ec";
        case LocalType::REPLICATION: return "replication";
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

            // Declustered parity: each stripe uses (m+k) disks out of n
            int stripe_size = m + k;
            if (stripe_size >= n) {
                // Clustered parity (no declustering): single stripe per group
                // If f > k, the entire stripe is unrecoverable = 100% data loss
                return 1.0;
            }

            // Calculate fraction of stripes that have >= k+1 failures
            // P(stripe unsafe) = sum_{i=k+1}^{min(f, stripe_size)} C(f,i) * C(n-f, stripe_size-i) / C(n, stripe_size)
            double total_combinations = combinations(n, stripe_size);
            if (total_combinations <= 0) {
                // Fallback: if can't calculate, assume worst case
                return 1.0;
            }

            double unsafe_combinations = 0.0;
            int max_failures_in_stripe = std::min(failure_count, stripe_size);
            for (int i = k + 1; i <= max_failures_in_stripe; ++i) {
                unsafe_combinations += combinations(failure_count, i) *
                                       combinations(n - failure_count, stripe_size - i);
            }

            double unsafe_ratio = unsafe_combinations / total_combinations;
            return unsafe_ratio;  // Fraction of data lost
        }

        case ECType::LRC: {
            // LRC structure:
            // - num_local_groups local groups, each with (local_m + local_k) disks
            // - k global parity disks at the end
            // Total: num_local_groups * (local_m + local_k) + k = n
            int num_local_groups = m / local_m;
            int local_group_size = local_m + local_k;
            int global_parity_start = num_local_groups * local_group_size;

            // Count failures per local group and global parity failures
            std::vector<int> local_failures(num_local_groups, 0);
            int global_parity_failures = 0;

            for (int idx : failed_disk_indices) {
                if (idx >= global_parity_start && idx < global_parity_start + k) {
                    // This is a global parity disk (only count k global parity disks)
                    global_parity_failures++;
                } else if (idx < global_parity_start) {
                    int local_group = idx / local_group_size;
                    if (local_group < num_local_groups) {
                        local_failures[local_group]++;
                    }
                }
                // idx >= global_parity_start + k: spare/unused disk, ignore
            }

            // Available global parity = k - failed global parity
            int available_global_parity = k - global_parity_failures;
            if (available_global_parity < 0) available_global_parity = 0;

            // Check each local group
            int total_global_parity_needed = 0;
            for (int i = 0; i < num_local_groups; ++i) {
                // If local group needs more recovery than available, data loss
                int recovery_needed = local_failures[i] - local_k;
                if (recovery_needed > 0) {
                    total_global_parity_needed += recovery_needed;
                }
            }

            // If total global parity needed exceeds available, data loss
            if (total_global_parity_needed > available_global_parity) {
                return static_cast<double>(failure_count) / n;
            }
            return 0.0;
        }

        case ECType::MULTI_EC: {
            // Multi-EC: local stripe/group failures count as one "chunk" failure
            int num_local_groups = m + k;

            if (local_type == LocalType::REPLICATION) {
                // Local layer uses replication: local_n copies per chunk
                // A local group fails only if ALL copies fail
                std::vector<int> group_failures(num_local_groups, 0);
                for (int idx : failed_disk_indices) {
                    int group_idx = idx / local_n;
                    if (group_idx < num_local_groups) {
                        group_failures[group_idx]++;
                    }
                }

                int failed_groups = 0;
                for (int i = 0; i < num_local_groups; ++i) {
                    // Replication fails only if all local_n copies fail
                    if (group_failures[i] >= local_n) {
                        failed_groups++;
                    }
                }

                // If more than k groups fail, data loss
                if (failed_groups > k) {
                    return static_cast<double>(failure_count) / total_disks;
                }
            } else {
                // Local layer uses EC (original behavior)
                int local_stripe_size = local_m + local_k;

                // Count failed local stripes (stripes with more than local_k failures)
                std::vector<int> stripe_failures(num_local_groups, 0);
                for (int idx : failed_disk_indices) {
                    int stripe_idx = idx / local_stripe_size;
                    if (stripe_idx < num_local_groups) {
                        stripe_failures[stripe_idx]++;
                    }
                }

                int failed_stripes = 0;
                for (int i = 0; i < num_local_groups; ++i) {
                    if (stripe_failures[i] > local_k) {
                        failed_stripes++;
                    }
                }

                // If more than k stripes fail, data loss
                if (failed_stripes > k) {
                    return static_cast<double>(failure_count) / total_disks;
                }
            }
            return 0.0;
        }

        case ECType::REPLICATION: {
            // n-way replication: data loss only if all n copies fail
            // k = n - 1, so tolerates k failures
            if (failure_count > k) {
                return 1.0;  // All copies lost = 100% data loss
            }
            return 0.0;  // At least one copy survives
        }
    }
    return 0.0;
}

// Get number of disks needed for rebuild read
int ECConfig::get_rebuild_read_disk_count(const std::set<int>& failed_disk_indices, int target_disk) const {
    switch (type) {
        case ECType::STANDARD:
            // Rebuild reads from (n - k) disks in the declustered group.
            // If n == m+k, this equals m.
            return std::max(0, n - k);

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
            // Need global rebuild, use (n - k) disks (global parity count = k)
            return std::max(0, n - k);
        }

        case ECType::MULTI_EC: {
            if (local_type == LocalType::REPLICATION) {
                // Local layer uses replication
                int group_idx = target_disk / local_n;
                int group_start = group_idx * local_n;
                int group_end = group_start + local_n;

                int group_failures = 0;
                for (int idx : failed_disk_indices) {
                    if (idx >= group_start && idx < group_end) {
                        group_failures++;
                    }
                }

                if (group_failures < local_n) {
                    // Can rebuild locally: read from any surviving copy
                    return 1;
                }
                // All copies failed - this is a chunk failure, needs global rebuild
                // Read from m surviving chunks (outer EC)
                return std::max(0, n - k);
            } else {
                // Local layer uses EC (original behavior)
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
                // Need global rebuild, treat the local stripe as a chunk and read (n - k) chunks worth.
                return std::max(0, n - k);
            }
        }

        case ECType::REPLICATION:
            // Replication rebuild: read from any single surviving copy
            return 1;
    }
    return m;  // Default
}

// ErasureCodingScheme implementation

void ErasureCodingScheme::initialize(const ECConfig& config) {
    config_ = config;
    if (!config_.validate()) {
        throw std::runtime_error("Invalid EC configuration");
    }
    use_custom_mapping_ = false;
    disk_to_group_.clear();
    group_to_disks_.clear();
}

void ErasureCodingScheme::initialize_with_custom_groups(const ECConfig& config,
                                                         const std::vector<std::vector<int>>& ec_groups) {
    config_ = config;

    // Build custom mappings
    use_custom_mapping_ = true;
    disk_to_group_.clear();
    group_to_disks_ = ec_groups;

    for (size_t group_idx = 0; group_idx < ec_groups.size(); ++group_idx) {
        for (int disk_idx : ec_groups[group_idx]) {
            disk_to_group_[disk_idx] = static_cast<int>(group_idx);
        }
    }

    // Validate that each group has the right size
    int expected_size = get_group_size();
    for (size_t i = 0; i < ec_groups.size(); ++i) {
        if (static_cast<int>(ec_groups[i].size()) != expected_size) {
            throw std::runtime_error("EC group " + std::to_string(i) +
                " has " + std::to_string(ec_groups[i].size()) +
                " disks, expected " + std::to_string(expected_size));
        }
    }
}

int ErasureCodingScheme::get_disk_group(int disk_index) const {
    if (use_custom_mapping_) {
        auto it = disk_to_group_.find(disk_index);
        if (it != disk_to_group_.end()) {
            return it->second;
        }
        return -1;  // Disk not in any group
    }
    int group_size = get_group_size();
    return disk_index / group_size;
}

std::vector<int> ErasureCodingScheme::get_stripe_disks(int disk_index) const {
    if (use_custom_mapping_) {
        int group_idx = get_disk_group(disk_index);
        if (group_idx >= 0 && group_idx < static_cast<int>(group_to_disks_.size())) {
            return group_to_disks_[group_idx];
        }
        return {};  // Empty if not found
    }

    std::vector<int> disks;
    int group_size = get_group_size();
    int group_start = (disk_index / group_size) * group_size;

    for (int i = 0; i < group_size; ++i) {
        disks.push_back(group_start + i);
    }
    return disks;
}

std::vector<int> ErasureCodingScheme::get_group_disks(int group_index) const {
    if (use_custom_mapping_) {
        if (group_index >= 0 && group_index < static_cast<int>(group_to_disks_.size())) {
            return group_to_disks_[group_index];
        }
        return {};
    }

    // Sequential mapping
    std::vector<int> disks;
    int group_size = get_group_size();
    int group_start = group_index * group_size;

    for (int i = 0; i < group_size; ++i) {
        disks.push_back(group_start + i);
    }
    return disks;
}

bool ErasureCodingScheme::can_rebuild(const std::set<int>& failed_disks) const {
    // Can rebuild if data loss is 0
    return get_data_loss_ratio(failed_disks) == 0.0;
}

double ErasureCodingScheme::get_data_loss_ratio(const std::set<int>& failed_disks) const {
    if (use_custom_mapping_ && !failed_disks.empty()) {
        // Convert absolute indices to relative indices within the group
        // All failed disks should be in the same group (this is called per-group)
        int group_index = get_disk_group(*failed_disks.begin());
        if (group_index >= 0) {
            std::set<int> relative = to_relative_indices(group_index, failed_disks);
            return config_.calculate_data_loss(relative);
        }
    }
    return config_.calculate_data_loss(failed_disks);
}

bool ErasureCodingScheme::is_data_loss(const std::set<int>& failed_disks) const {
    return get_data_loss_ratio(failed_disks) > 0.0;
}

int ErasureCodingScheme::get_rebuild_read_count(const std::set<int>& failed_disks, int target_disk) const {
    if (use_custom_mapping_ && !failed_disks.empty()) {
        // Convert to relative indices
        int group_index = get_disk_group(target_disk);
        if (group_index >= 0) {
            std::set<int> relative_failed = to_relative_indices(group_index, failed_disks);
            int relative_target = -1;

            // Find relative index for target_disk
            const auto& group_disks = group_to_disks_[group_index];
            for (size_t i = 0; i < group_disks.size(); ++i) {
                if (group_disks[i] == target_disk) {
                    relative_target = static_cast<int>(i);
                    break;
                }
            }

            if (relative_target >= 0) {
                return config_.get_rebuild_read_disk_count(relative_failed, relative_target);
            }
        }
    }
    return config_.get_rebuild_read_disk_count(failed_disks, target_disk);
}

bool ErasureCodingScheme::can_local_rebuild(const std::set<int>& failed_disks, int target_disk) const {
    // Convert to relative indices for custom mapping
    std::set<int> relative_failed;
    int relative_target;

    if (use_custom_mapping_) {
        int group_index = get_disk_group(target_disk);
        if (group_index < 0) return false;

        relative_failed = to_relative_indices(group_index, failed_disks);

        // Find relative index for target_disk
        const auto& group_disks = group_to_disks_[group_index];
        relative_target = -1;
        for (size_t i = 0; i < group_disks.size(); ++i) {
            if (group_disks[i] == target_disk) {
                relative_target = static_cast<int>(i);
                break;
            }
        }
        if (relative_target < 0) return false;
    } else {
        relative_failed = failed_disks;
        relative_target = target_disk;
    }

    if (config_.type == ECType::LRC) {
        int local_group_size = config_.local_m + config_.local_k;
        int local_group = relative_target / local_group_size;
        int local_start = local_group * local_group_size;
        int local_end = local_start + local_group_size;

        int local_failures = 0;
        for (int idx : relative_failed) {
            if (idx >= local_start && idx < local_end) {
                local_failures++;
            }
        }
        return local_failures <= config_.local_k;
    } else if (config_.type == ECType::MULTI_EC) {
        if (config_.local_type == LocalType::REPLICATION) {
            // Local layer uses replication: local rebuild if at least one copy survives
            int group_idx = relative_target / config_.local_n;
            int group_start = group_idx * config_.local_n;
            int group_end = group_start + config_.local_n;

            int group_failures = 0;
            for (int idx : relative_failed) {
                if (idx >= group_start && idx < group_end) {
                    group_failures++;
                }
            }
            return group_failures < config_.local_n;  // At least one copy survives
        } else {
            // Local layer uses EC (original behavior)
            int local_stripe_size = config_.local_m + config_.local_k;
            int stripe_idx = relative_target / local_stripe_size;
            int stripe_start = stripe_idx * local_stripe_size;
            int stripe_end = stripe_start + local_stripe_size;

            int stripe_failures = 0;
            for (int idx : relative_failed) {
                if (idx >= stripe_start && idx < stripe_end) {
                    stripe_failures++;
                }
            }
            return stripe_failures <= config_.local_k;
        }
    } else if (config_.type == ECType::REPLICATION) {
        // Replication: always "local" rebuild (read from any single surviving copy)
        return true;
    }
    return false;  // Standard EC has no local parity
}

std::set<int> ErasureCodingScheme::to_relative_indices(int group_index, const std::set<int>& absolute_indices) const {
    std::set<int> relative;

    if (use_custom_mapping_) {
        // Build mapping from absolute to relative index
        if (group_index < static_cast<int>(group_to_disks_.size())) {
            const auto& group_disks = group_to_disks_[group_index];
            for (int abs_idx : absolute_indices) {
                for (size_t i = 0; i < group_disks.size(); ++i) {
                    if (group_disks[i] == abs_idx) {
                        relative.insert(static_cast<int>(i));
                        break;
                    }
                }
            }
        }
    } else {
        // Sequential mapping: relative = absolute - start
        int group_size = get_group_size();
        int start = group_index * group_size;
        for (int abs_idx : absolute_indices) {
            if (abs_idx >= start && abs_idx < start + group_size) {
                relative.insert(abs_idx - start);
            }
        }
    }
    return relative;
}

int ErasureCodingScheme::get_total_groups(int total_disks) const {
    if (use_custom_mapping_) {
        return static_cast<int>(group_to_disks_.size());
    }
    int group_size = get_group_size();
    return (total_disks + group_size - 1) / group_size;
}

int ErasureCodingScheme::get_group_start_disk(int group_index) const {
    if (use_custom_mapping_) {
        // For custom mapping, return the first disk in the group
        if (group_index >= 0 && group_index < static_cast<int>(group_to_disks_.size()) &&
            !group_to_disks_[group_index].empty()) {
            return group_to_disks_[group_index][0];
        }
        return -1;
    }
    return group_index * get_group_size();
}
