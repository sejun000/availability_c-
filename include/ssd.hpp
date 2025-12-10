#pragma once

#include <string>

const std::string SSD_MODULE_NAME = "SSD";

// SSD naming and indexing utilities
std::string get_ssd_name(int ssd_index);
int get_ssd_index(const std::string& ssd_name);
std::string get_ssd_group_name(int ssd_group_index);
int get_ssd_group_index_from_group_name(const std::string& ssd_group_name);
bool is_event_node_ssd(const std::string& event_node);

// SSD Redundancy Scheme
class SSDRedundancyScheme {
public:
    SSDRedundancyScheme(
        double write_bw, double read_bw, double mttf,
        double cached_write_ratio, double cached_write_bw, double cached_read_bw, double cached_mttf,
        int m, int k, int l,
        int cached_m, int cached_k, int cached_l,
        int network_m, int network_k, int network_l,
        int cached_network_m, int cached_network_k, int cached_network_l,
        int cached_ssds, int total_ssds,
        int inter_replicas, int intra_replicas
    );

    // Check if SSD is in cached tier
    bool is_ssd_index_cached(int ssd_index) const;
    bool is_ssd_group_index_cached(int ssd_group_index) const;

    // Get SSD group information
    int get_start_ssd_index(int ssd_group_index) const;
    int get_total_group_count() const;
    int get_ssd_group_index(int ssd_index) const;

    // Get parameters based on cached tier
    double get_read_bw(bool cached) const;
    double get_write_bw(bool cached) const;
    double get_mttf(bool cached) const;
    int get_m(bool cached) const;
    int get_k(bool cached) const;
    int get_l(bool cached) const;
    int get_network_m(bool cached) const;
    int get_network_k(bool cached) const;
    int get_network_l(bool cached) const;

    int get_total_ssds() const;
    int get_tiered_ssds(bool cached) const;
    std::string get_cached_prefix(bool cached) const;
    int get_inter_replicas(bool cached) const;
    int get_intra_replicas(bool cached) const;

private:
    double write_bw_;
    double read_bw_;
    double mttf_;
    double cached_write_ratio_;
    double cached_write_bw_;
    double cached_read_bw_;
    double cached_mttf_;
    int m_;
    int k_;
    int l_;
    int cached_m_;
    int cached_k_;
    int cached_l_;
    int network_m_;
    int network_k_;
    int network_l_;
    int cached_network_m_;
    int cached_network_k_;
    int cached_network_l_;
    int cached_ssds_;
    int total_ssds_;
    int ssd_group_size_;
    int cached_ssd_group_size_;
    int inter_replicas_;
    int intra_replicas_;
};
