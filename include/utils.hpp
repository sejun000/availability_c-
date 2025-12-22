#pragma once

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <mutex>
#include <nlohmann/json.hpp>

// Encoding time data structure
struct EncodingTimeEntry {
    int n;
    int k;
    double encoding_time;  // in microseconds
};

// Empirical CDF for trace-driven failure time sampling
class EmpiricalCDF {
public:
    EmpiricalCDF() : initialized_(false), t_max_(0.0) {}

    // Load failure times from simple text file (one value per line, in hours)
    bool load_from_file(const std::string& file_path);

    // Sample failure time using inverse transform sampling (raw empirical)
    // Returns failure time in hours
    double sample(std::mt19937& rng) const;

    // Sample using spliced distribution (empirical + exponential tail)
    // ARR: Annual Replacement Rate (e.g., 0.0074 for 0.74%)
    // Returns failure time in hours
    double sample_spliced(std::mt19937& rng, double arr) const;

    // Check if CDF is loaded and ready
    bool is_initialized() const { return initialized_; }

    // Get statistics
    size_t get_sample_count() const { return failure_times_.size(); }
    double get_mean() const;
    double get_median() const;
    double get_max() const { return t_max_; }

private:
    std::vector<double> failure_times_;  // sorted failure times in hours
    std::vector<double> cdf_;            // cumulative probabilities
    bool initialized_;
    double t_max_;                       // maximum failure time in trace
};

class Utils {
public:
    // Initialize encoding time data
    static void initialize_encoding_time_data();

    // Get encoding latency in microseconds
    static double get_encoding_latency_usec(int m, int k, bool replication = false);

    // Get encoding latency in seconds
    static double get_encoding_latency_sec(int m, int k, bool replication = false);

    // Convert KMG notation to bytes (e.g., "100G" -> 100000000000)
    static double KMG_to_bytes(const std::string& s);

    // Convert time string to seconds (e.g., "135us", "5s", "2m", "1h")
    static double convert_to_seconds(const std::string& time_str);

    // Convert time string to microseconds
    static double convert_to_microseconds(const std::string& time_str);

    // Calculate availability nines (-log10(1 - availability))
    static double get_nines(double availability);

    // Calculate WAF from overprovisioning ratio
    static double get_waf_from_op(double op);

    // Parse input from JSON file
    static bool parse_input_from_json(
        const std::string& file_path,
        std::vector<std::tuple<std::string, std::string, std::string>>& edges,
        std::map<std::string, std::vector<std::string>>& enclosures,
        std::map<std::string, double>& mttfs,
        std::map<std::string, double>& mtrs,
        std::map<std::string, double>& costs,
        nlohmann::json& options
    );

    // Progress bar display
    static void progress_bar(int i, int total, int bar_len = 50);

    // Calculate percentile values from distribution
    static std::tuple<double, double, double, double, double> get_percentile_value(
        const std::map<double, double>& raw_datas,
        bool ascending = true);

    // Write simulation results to CSV file
    static void write_results_to_csv(
        const std::string& output_file,
        const std::map<std::string, nlohmann::json>& params_and_results);

private:
    static std::map<std::pair<int, int>, double> encoding_time_map_;
    static std::once_flag init_flag_;
    static void do_initialize_encoding_time_data();
};
