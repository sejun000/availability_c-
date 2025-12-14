#pragma once

#include <string>
#include <vector>
#include <set>
#include <map>
#include <stdexcept>

// Encoding entity - who performs the EC encoding
enum class EncodingEntity {
    CONTROLLER,  // Controller performs encoding (no cross-IO-module traffic)
    IO_MODULE    // IO modules perform encoding (requires cross-IO-module traffic for parity)
};

// Parse encoding entity from string
EncodingEntity parse_encoding_entity(const std::string& entity_str);

// Convert encoding entity to string
std::string encoding_entity_to_string(EncodingEntity entity);

// Encoding configuration
struct EncodingConfig {
    EncodingEntity entity = EncodingEntity::CONTROLLER;  // Default: controller does encoding
    double chunk_size = 1024 * 1024;  // Default: 1MB chunks for bandwidth calculation

    // Calculate encoding cross-traffic ratio for write operations
    // Returns the ratio of cross-IO-module traffic to host write traffic
    // io_module_count: number of IO modules in the EC stripe
    // data_chunks: m (number of data chunks)
    // parity_chunks: k (number of parity chunks)
    // io_module_data_distribution: map of io_module_name -> number of data chunks on that module
    double calculate_cross_traffic_ratio(
        int io_module_count,
        int data_chunks,
        int parity_chunks,
        const std::map<std::string, int>& io_module_data_distribution,
        const std::map<std::string, int>& io_module_parity_distribution) const;
};

// Erasure Coding Types
enum class ECType {
    STANDARD,    // Standard EC: m+k, tolerates k failures
    LRC,         // Local Reconstruction Code: local_m+local_k groups within m+k
    MULTI_EC,    // Multi-level EC: local stripe (local_m+local_k) treated as chunk for outer m+k
    REPLICATION  // N-way replication: 1 original + (n-1) replicas, rebuild reads 1 disk
};

// Local layer type for MULTI_EC
enum class LocalType {
    EC,          // Local layer uses erasure coding (default, original MULTI_EC behavior)
    REPLICATION  // Local layer uses replication (each local group is n-way replicated)
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

    // Local layer type for MULTI_EC (default: EC for backward compatibility)
    LocalType local_type = LocalType::EC;

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
                if (m <= 0 || k < 0) return false;
                if (n < m + k) return false;  // n is number of local groups for outer EC
                {
                    int num_local_groups = m + k;  // Each local stripe is a chunk
                    double local_ratio;

                    if (local_type == LocalType::REPLICATION) {
                        // Local layer uses replication: local_n copies per chunk
                        if (local_n <= 0) return false;
                        total_disks = local_n * num_local_groups;
                        local_ratio = 1.0 / local_n;  // Replication efficiency
                    } else {
                        // Local layer uses EC (original behavior)
                        if (local_m <= 0 || local_k < 0) return false;
                        if (local_n < local_m + local_k) return false;
                        total_disks = local_n * num_local_groups;
                        local_ratio = static_cast<double>(local_m) / (local_m + local_k);
                    }

                    double outer_ratio = static_cast<double>(m) / (m + k);
                    effective_ratio = local_ratio * outer_ratio;
                }
                break;

            case ECType::REPLICATION:
                // n-way replication: n copies of data
                // k = n - 1 (tolerates n-1 failures)
                // m = 1 (single data unit)
                if (n <= 0) return false;
                m = 1;
                k = n - 1;
                total_disks = n;
                effective_ratio = 1.0 / n;  // 1/n storage efficiency
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
            case ECType::REPLICATION:
                return k;  // n-1 failures tolerable
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

    // Initialize with custom group mapping
    // ec_groups: vector of vectors, each inner vector contains disk indices for that group
    void initialize_with_custom_groups(const ECConfig& config,
                                        const std::vector<std::vector<int>>& ec_groups);

    // Check if using custom mapping
    bool has_custom_mapping() const { return use_custom_mapping_; }

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

    // Check if data loss occurred (considers LRC local parity)
    bool is_data_loss(const std::set<int>& failed_disks) const;

    // Get number of read disks needed for rebuild
    int get_rebuild_read_count(const std::set<int>& failed_disks, int target_disk) const;

    // For LRC: check if local rebuild is possible
    bool can_local_rebuild(const std::set<int>& failed_disks, int target_disk) const;

    // Get the disk group size
    int get_group_size() const {
        if (config_.type == ECType::MULTI_EC) {
            return config_.total_disks;
        }
        return config_.n > 0 ? config_.n : config_.m + config_.k;
    }

    // Get total number of EC groups
    int get_total_groups(int total_disks) const;

    // Get start disk index for a group (only meaningful for sequential mapping)
    int get_group_start_disk(int group_index) const;

    // Get all disks in a specific group (works with both custom and sequential mapping)
    std::vector<int> get_group_disks(int group_index) const;

    // Convert absolute disk indices to relative indices (0-based within EC group)
    // This is needed for LRC local group calculation with custom mapping
    std::set<int> to_relative_indices(int group_index, const std::set<int>& absolute_indices) const;

private:
    ECConfig config_;

    // Custom mapping support
    bool use_custom_mapping_ = false;
    std::map<int, int> disk_to_group_;           // disk_index -> group_index
    std::vector<std::vector<int>> group_to_disks_;  // group_index -> [disk_indices]
};

// Parse EC type from string
ECType parse_ec_type(const std::string& type_str);

// Convert EC type to string
std::string ec_type_to_string(ECType type);

// Parse local type from string
LocalType parse_local_type(const std::string& type_str);

// Convert local type to string
std::string local_type_to_string(LocalType type);
