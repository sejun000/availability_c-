#pragma once

#include <string>
#include <vector>
#include <set>
#include <map>
#include <stdexcept>

// Erasure Coding Types
enum class ECType {
    STANDARD,   // Standard EC: m+k, tolerates k failures
    LRC,        // Local Reconstruction Code: local_m+local_k groups within m+k
    MULTI_EC    // Multi-level EC: local stripe (local_m+local_k) treated as chunk for outer m+k
};

// Erasure Coding Configuration
struct ECConfig {
    ECType type = ECType::STANDARD;

    // Standard EC parameters
    int m = 0;      // Data chunks (or outer data chunks for Multi-EC)
    int k = 0;      // Parity chunks (or outer parity chunks for Multi-EC)

    // LRC / Multi-EC local parameters
    int local_m = 0;    // Local data chunks
    int local_k = 0;    // Local parity chunks

    // Disk group size (for declustered parity placement)
    int n = 0;          // Disk group size for standard EC / LRC
    int local_n = 0;    // Local disk group size for Multi-EC

    // Computed values
    int total_disks = 0;        // Total disks in the EC group
    double effective_ratio = 0; // Effective capacity ratio

    // Validate and compute derived values
    bool validate() {
        switch (type) {
            case ECType::STANDARD:
                if (m <= 0 || k < 0) return false;
                if (n < m + k) return false;
                total_disks = n;
                effective_ratio = static_cast<double>(m) / (m + k);
                break;

            case ECType::LRC:
                if (m <= 0 || k < 0 || local_m <= 0 || local_k < 0) return false;
                if (m % local_m != 0) return false;  // local_m must divide m
                {
                    int num_local_groups = m / local_m;
                    int total_parities = k + local_k * num_local_groups;
                    if (n < m + total_parities) return false;
                    total_disks = n;
                    effective_ratio = static_cast<double>(m) / (m + k + local_k * num_local_groups);
                }
                break;

            case ECType::MULTI_EC:
                if (m <= 0 || k < 0 || local_m <= 0 || local_k < 0) return false;
                if (local_n < local_m + local_k) return false;
                if (n < m + k) return false;  // n is number of local groups for outer EC
                {
                    int num_local_groups = m + k;  // Each local stripe is a chunk
                    total_disks = local_n * num_local_groups;
                    double local_ratio = static_cast<double>(local_m) / (local_m + local_k);
                    double outer_ratio = static_cast<double>(m) / (m + k);
                    effective_ratio = local_ratio * outer_ratio;
                }
                break;
        }
        return true;
    }

    // Get maximum tolerable failures
    int max_tolerable_failures() const {
        switch (type) {
            case ECType::STANDARD:
                return k;
            case ECType::LRC:
                // Can tolerate k + local_k failures (in specific patterns)
                return k + local_k;
            case ECType::MULTI_EC:
                // Complex: depends on failure distribution
                // Worst case: k outer failures, each with local_k local failures
                return k;  // Conservative estimate
        }
        return 0;
    }

    // Check if data loss occurred given failed disk indices within the EC group
    // Returns fraction of data lost (0.0 = no loss, 1.0 = complete loss)
    double calculate_data_loss(const std::set<int>& failed_disk_indices) const;

    // Get number of disks needed for rebuild read (excluding the failed disk)
    int get_rebuild_read_disk_count(const std::set<int>& failed_disk_indices, int target_disk) const;
};

// Erasure Coding Scheme Manager
class ErasureCodingScheme {
public:
    ErasureCodingScheme() = default;

    // Initialize from JSON config
    void initialize(const ECConfig& config);

    // Get EC configuration
    const ECConfig& get_config() const { return config_; }

    // Get disk group for a given disk index (0-based within total disks)
    int get_disk_group(int disk_index) const;

    // Get all disks in the same EC stripe as the given disk
    std::vector<int> get_stripe_disks(int disk_index) const;

    // Check if rebuild is possible given current failures
    bool can_rebuild(const std::set<int>& failed_disks) const;

    // Calculate data loss ratio
    double get_data_loss_ratio(const std::set<int>& failed_disks) const;

    // Get number of read disks needed for rebuild
    int get_rebuild_read_count(const std::set<int>& failed_disks, int target_disk) const;

    // For LRC: check if local rebuild is possible
    bool can_local_rebuild(const std::set<int>& failed_disks, int target_disk) const;

    // Get the disk group size
    int get_group_size() const { return config_.n > 0 ? config_.n : config_.m + config_.k; }

    // Get total number of EC groups
    int get_total_groups(int total_disks) const;

    // Get start disk index for a group
    int get_group_start_disk(int group_index) const;

private:
    ECConfig config_;
};

// Parse EC type from string
ECType parse_ec_type(const std::string& type_str);

// Convert EC type to string
std::string ec_type_to_string(ECType type);
