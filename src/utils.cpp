#include "utils.hpp"
#include "logger.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sys/stat.h>

std::map<std::pair<int, int>, double> Utils::encoding_time_map_;
std::once_flag Utils::init_flag_;

// EmpiricalCDF implementation
bool EmpiricalCDF::load_from_file(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        Logger::getInstance().error("Failed to open failure trace file: " + file_path);
        return false;
    }

    failure_times_.clear();
    std::string line;

    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        try {
            double value = std::stod(line);
            if (value >= 0) {
                failure_times_.push_back(value);
            }
        } catch (...) {
            // Skip invalid lines
            continue;
        }
    }

    file.close();

    if (failure_times_.empty()) {
        Logger::getInstance().error("No valid failure times found in file");
        return false;
    }

    // Sort failure times
    std::sort(failure_times_.begin(), failure_times_.end());

    // Set T_max
    t_max_ = failure_times_.back();

    // Build CDF
    cdf_.resize(failure_times_.size());
    for (size_t i = 0; i < failure_times_.size(); ++i) {
        cdf_[i] = static_cast<double>(i + 1) / static_cast<double>(failure_times_.size());
    }

    initialized_ = true;
    Logger::getInstance().info("Loaded " + std::to_string(failure_times_.size()) +
                               " failure times from trace (mean: " +
                               std::to_string(get_mean()) + " hours, max: " +
                               std::to_string(t_max_) + " hours)");
    return true;
}

double EmpiricalCDF::sample(std::mt19937& rng) const {
    if (!initialized_ || failure_times_.empty()) {
        // Fallback: return a default value (shouldn't happen if properly initialized)
        return 8760.0;  // 1 year in hours
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(rng);

    // Binary search for inverse CDF
    auto it = std::lower_bound(cdf_.begin(), cdf_.end(), u);
    size_t idx = std::distance(cdf_.begin(), it);

    if (idx >= failure_times_.size()) {
        idx = failure_times_.size() - 1;
    }

    // Linear interpolation for smoother sampling
    if (idx > 0 && idx < failure_times_.size()) {
        double cdf_low = cdf_[idx - 1];
        double cdf_high = cdf_[idx];
        double t_low = failure_times_[idx - 1];
        double t_high = failure_times_[idx];

        if (cdf_high > cdf_low) {
            double alpha = (u - cdf_low) / (cdf_high - cdf_low);
            return t_low + alpha * (t_high - t_low);
        }
    }

    return failure_times_[idx];
}

double EmpiricalCDF::sample_spliced(std::mt19937& rng, double arr) const {
    if (!initialized_ || failure_times_.empty()) {
        // Fallback to exponential with MTTF from ARR
        double mttf_hours = (1.0 / arr) * 365.0 * 24.0;
        std::exponential_distribution<double> exp_dist(1.0 / mttf_hours);
        return exp_dist(rng);
    }

    // Calculate MTTF from ARR (in hours)
    // MTTF = 1/ARR years = (1/ARR) * 365 * 24 hours
    double mttf_hours = (1.0 / arr) * 365.0 * 24.0;

    // Calculate probability p = P(Exp(MTTF) <= T_max) = 1 - exp(-T_max/MTTF)
    double p = 1.0 - std::exp(-t_max_ / mttf_hours);

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(rng);

    if (u < p) {
        // Sample from empirical distribution (within [0, T_max])
        return sample(rng);
    } else {
        // Sample from exponential tail: T_max + Exp(MTTF)
        std::exponential_distribution<double> exp_dist(1.0 / mttf_hours);
        return t_max_ + exp_dist(rng);
    }
}

double EmpiricalCDF::get_mean() const {
    if (failure_times_.empty()) return 0.0;
    double sum = 0.0;
    for (double t : failure_times_) {
        sum += t;
    }
    return sum / failure_times_.size();
}

double EmpiricalCDF::get_median() const {
    if (failure_times_.empty()) return 0.0;
    return failure_times_[failure_times_.size() / 2];
}

void Utils::initialize_encoding_time_data() {
    std::call_once(init_flag_, do_initialize_encoding_time_data);
}

void Utils::do_initialize_encoding_time_data() {
    // k = 1 data
    std::vector<EncodingTimeEntry> data = {
        {2, 1, 65}, {3, 1, 91}, {4, 1, 98}, {5, 1, 115}, {6, 1, 138},
        {7, 1, 150}, {8, 1, 127}, {9, 1, 141}, {10, 1, 157}, {11, 1, 173},
        {12, 1, 199}, {13, 1, 206}, {14, 1, 224}, {15, 1, 240}, {16, 1, 256},
        {17, 1, 273}, {18, 1, 290}, {19, 1, 306}, {20, 1, 321}, {21, 1, 363},
        {22, 1, 353}, {23, 1, 369}, {24, 1, 387}, {25, 1, 401}, {26, 1, 420},
        {27, 1, 449}, {28, 1, 449}, {29, 1, 481}, {30, 1, 497}, {31, 1, 551},
        {32, 1, 759}, {33, 1, 808}, {34, 1, 814}, {35, 1, 830}, {36, 1, 914},
        {37, 1, 921}, {38, 1, 963}, {39, 1, 914}, {40, 1, 960}, {41, 1, 1004},
        {42, 1, 1576}, {43, 1, 1054}, {44, 1, 1032}, {45, 1, 1054}, {46, 1, 1154},
        {47, 1, 1103}, {48, 1, 1173},

        // k = 2 data
        {2, 2, 133}, {3, 2, 157}, {4, 2, 177}, {5, 2, 195}, {6, 2, 273},
        {7, 2, 241}, {8, 2, 191}, {9, 2, 219}, {10, 2, 245}, {11, 2, 265},
        {12, 2, 278}, {13, 2, 300}, {14, 2, 332}, {15, 2, 361}, {16, 2, 366},
        {17, 2, 382}, {18, 2, 404}, {19, 2, 434}, {20, 2, 445}, {21, 2, 474},
        {22, 2, 552}, {23, 2, 520}, {24, 2, 576}, {25, 2, 592}, {26, 2, 604},
        {27, 2, 640}, {28, 2, 660}, {29, 2, 770}, {30, 2, 991}, {31, 2, 1071},
        {32, 2, 1151}, {33, 2, 1128}, {34, 2, 1184}, {35, 2, 1259}, {36, 2, 1299},
        {37, 2, 1250}, {38, 2, 1288}, {39, 2, 1316}, {40, 2, 1364}, {41, 2, 1474},
        {42, 2, 1692}, {43, 2, 1485}, {44, 2, 1525}, {45, 2, 1521}, {46, 2, 1589},
        {47, 2, 1627}, {48, 2, 1679},

        // k = 3 data
        {2, 3, 152}, {3, 3, 177}, {4, 3, 189}, {5, 3, 203}, {6, 3, 215},
        {7, 3, 248}, {8, 3, 279}, {9, 3, 307}, {10, 3, 331}, {11, 3, 367},
        {12, 3, 392}, {13, 3, 438}, {14, 3, 459}, {15, 3, 508}, {16, 3, 527},
        {17, 3, 606}, {18, 3, 593}, {19, 3, 632}, {20, 3, 705}, {21, 3, 692},
        {22, 3, 883}, {23, 3, 774}, {24, 3, 856}, {25, 3, 878}, {26, 3, 889},
        {27, 3, 938}, {28, 3, 1242}, {29, 3, 1273}, {30, 3, 1318}, {31, 3, 1384},
        {32, 3, 1488}, {33, 3, 1464}, {34, 3, 1545}, {35, 3, 1646}, {36, 3, 1677},
        {37, 3, 1671}, {38, 3, 1720}, {39, 3, 1770}, {40, 3, 1809}, {41, 3, 1906},
        {42, 3, 2799}, {43, 3, 1995}, {44, 3, 2025}, {45, 3, 2131}, {46, 3, 2109},
        {47, 3, 2224}, {48, 3, 2203},

        // k = 4 data
        {2, 4, 273}, {3, 4, 273}, {4, 4, 359}, {5, 4, 362}, {6, 4, 377},
        {7, 4, 437}, {8, 4, 358}, {9, 4, 417}, {10, 4, 437}, {11, 4, 476},
        {12, 4, 513}, {13, 4, 566}, {14, 4, 609}, {15, 4, 662}, {16, 4, 698},
        {17, 4, 731}, {18, 4, 783}, {19, 4, 857}, {20, 4, 916}, {21, 4, 939},
        {22, 4, 968}, {23, 4, 1023}, {24, 4, 1084}, {25, 4, 1098}, {26, 4, 1228},
        {27, 4, 1246}, {28, 4, 1610}, {29, 4, 1706}, {30, 4, 1782}, {31, 4, 1853},
        {32, 4, 1888}, {33, 4, 1965}, {34, 4, 2056}, {35, 4, 2069}, {36, 4, 2121},
        {37, 4, 2203}, {38, 4, 2261}, {39, 4, 2321}, {40, 4, 2393}, {41, 4, 4645},
        {42, 4, 2614}, {43, 4, 2657}, {44, 4, 2617}, {45, 4, 2775}, {46, 4, 2752},
        {47, 4, 2901}, {48, 4, 2876},

        // k = 5 data
        {2, 5, 341}, {3, 5, 362}, {4, 5, 386}, {5, 5, 453}, {6, 5, 528},
        {7, 5, 589}, {8, 5, 474}, {9, 5, 608}, {10, 5, 589}, {11, 5, 680},
        {12, 5, 716}, {13, 5, 807}, {14, 5, 831}, {15, 5, 908}, {16, 5, 939},
        {17, 5, 1021}, {18, 5, 1115}, {19, 5, 1146}, {20, 5, 1179}, {21, 5, 1243},
        {22, 5, 1362}, {23, 5, 1351}, {24, 5, 1447}, {25, 5, 1475}, {26, 5, 1582},
        {27, 5, 1657}, {28, 5, 1743}, {29, 5, 2201}, {30, 5, 2262}, {31, 5, 2298},
        {32, 5, 2596}, {33, 5, 2619}, {34, 5, 2771}, {35, 5, 2895}, {36, 5, 2983},
        {37, 5, 3074}, {38, 5, 3194}, {39, 5, 3289}, {40, 5, 3446}, {41, 5, 5398},
        {42, 5, 3498}, {43, 5, 3683}, {44, 5, 3681}, {45, 5, 3872}, {46, 5, 3855},
        {47, 5, 4135}, {48, 5, 4156},

        // k = 6 data
        {2, 6, 401}, {3, 6, 392}, {4, 6, 520}, {5, 6, 503}, {6, 6, 563},
        {7, 6, 639}, {8, 6, 525}, {9, 6, 584}, {10, 6, 713}, {11, 6, 808},
        {12, 6, 783}, {13, 6, 950}, {14, 6, 967}, {15, 6, 987}, {16, 6, 1053},
        {17, 6, 1119}, {18, 6, 1182}, {19, 6, 1245}, {20, 6, 1303}, {21, 6, 1408},
        {22, 6, 1506}, {23, 6, 1570}, {24, 6, 1582}, {25, 6, 1721}, {26, 6, 1805},
        {27, 6, 2027}, {28, 6, 2261}, {29, 6, 2770}, {30, 6, 2702}, {31, 6, 2824},
        {32, 6, 3092}, {33, 6, 3030}, {34, 6, 3101}, {35, 6, 3238}, {36, 6, 3452},
        {37, 6, 3539}, {38, 6, 3667}, {39, 6, 3616}, {40, 6, 3741}, {41, 6, 5460},
        {42, 6, 4075}, {43, 6, 4031}, {44, 6, 4336}, {45, 6, 4337}, {46, 6, 4381},
        {47, 6, 4498}, {48, 6, 4788}
    };

    for (const auto& entry : data) {
        encoding_time_map_[{entry.n, entry.k}] = entry.encoding_time;
    }
}

double Utils::get_encoding_latency_usec(int m, int k, bool replication) {
    if (k == 0 || replication) {
        return 0.0;
    }

    initialize_encoding_time_data();

    auto it = encoding_time_map_.find({m, k});
    if (it != encoding_time_map_.end()) {
        return it->second;
    }

    return 0.0;
}

double Utils::get_encoding_latency_sec(int m, int k, bool replication) {
    return get_encoding_latency_usec(m, k, replication) / 1e6;
}

double Utils::KMG_to_bytes(const std::string& s) {
    if (s.empty()) return 0.0;

    char last_char = s.back();
    std::string number_part = s.substr(0, s.length() - 1);
    double value = std::stod(number_part);

    if (last_char == 'K') {
        return value * 1000;
    } else if (last_char == 'M') {
        return value * 1000 * 1000;
    } else if (last_char == 'G') {
        return value * 1000 * 1000 * 1000;
    }

    // If no suffix, parse the whole string as a number
    return std::stod(s);
}

double Utils::convert_to_seconds(const std::string& time_str) {
    static const std::map<std::string, double> time_units = {
        {"us", 1e-6},
        {"ms", 1e-3},
        {"s", 1.0},
        {"m", 60.0},
        {"h", 3600.0},
        {"d", 86400.0}
    };

    std::regex pattern(R"((\d+\.?\d*)([a-zA-Z]+))");
    std::smatch match;

    if (!std::regex_match(time_str, match, pattern)) {
        throw std::invalid_argument("Invalid time string format");
    }

    double value = std::stod(match[1].str());
    std::string unit = match[2].str();

    auto it = time_units.find(unit);
    if (it == time_units.end()) {
        throw std::invalid_argument("Unsupported time unit: " + unit);
    }

    return value * it->second;
}

double Utils::convert_to_microseconds(const std::string& time_str) {
    return convert_to_seconds(time_str) * 1e6;
}

double Utils::get_nines(double availability) {
    if (availability >= 1.0) {
        return 13.0;  // Max value
    }
    return -std::log10(1.0 - availability);
}

double Utils::get_waf_from_op(double op) {
    return 0.5 * (1.0 + op) / op;
}

bool Utils::parse_input_from_json(
    const std::string& file_path,
    std::vector<std::tuple<std::string, std::string, std::string>>& edges,
    std::map<std::string, std::vector<std::string>>& enclosures,
    std::map<std::string, double>& mttfs,
    std::map<std::string, double>& mtrs,
    std::map<std::string, double>& costs,
    nlohmann::json& options
) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << file_path << std::endl;
        return false;
    }

    nlohmann::json data;
    try {
        file >> data;
    } catch (const std::exception& e) {
        std::cerr << "Failed to parse JSON: " << e.what() << std::endl;
        return false;
    }

    // Parse edges
    if (data.contains("edges")) {
        for (const auto& edge : data["edges"]) {
            std::string start = edge["start"];
            std::string end = edge["end"];
            std::string bandwidth = edge["bandwidth"];
            edges.push_back({start, end, bandwidth});
        }
    }

    // Parse enclosures
    if (data.contains("enclosures")) {
        enclosures = data["enclosures"].get<std::map<std::string, std::vector<std::string>>>();
    }

    // Parse mttf
    if (data.contains("mttf")) {
        mttfs = data["mttf"].get<std::map<std::string, double>>();
    }

    // Parse mtr
    if (data.contains("mtr")) {
        mtrs = data["mtr"].get<std::map<std::string, double>>();
    }

    // Parse cost
    if (data.contains("cost")) {
        costs = data["cost"].get<std::map<std::string, double>>();
    }

    // Parse options
    if (data.contains("options")) {
        options = data["options"];
    }

    return true;
}

void Utils::progress_bar(int i, int total, int bar_len) {
    double frac = static_cast<double>(i) / total;
    int filled = static_cast<int>(bar_len * frac);
    int pct = static_cast<int>(frac * 100);

    std::string bar(filled, '=');
    bar += std::string(bar_len - filled, '-');

    char endchar = (i == total) ? '\n' : '\r';
    std::cout << "[" << bar << "] " << std::setw(3) << pct << "% completed" << endchar << std::flush;
}

std::tuple<double, double, double, double, double> Utils::get_percentile_value(
    const std::map<double, double>& raw_datas,
    bool ascending) {

    if (raw_datas.empty()) {
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }

    // Create sorted vector
    std::vector<std::pair<double, double>> sorted_data(raw_datas.begin(), raw_datas.end());

    if (ascending) {
        std::sort(sorted_data.begin(), sorted_data.end());
    } else {
        std::sort(sorted_data.begin(), sorted_data.end(), std::greater<std::pair<double, double>>());
    }

    // Calculate cumulative time
    std::vector<double> cumulative_time;
    double sum = 0.0;
    for (const auto& [value, interval] : sorted_data) {
        sum += interval;
        cumulative_time.push_back(sum);
    }
    double total_time = sum;

    // Calculate percentiles
    auto find_percentile = [&](double percentile) -> double {
        double target_time = total_time * percentile;
        for (size_t i = 0; i < cumulative_time.size(); ++i) {
            if (cumulative_time[i] >= target_time) {
                return sorted_data[i].first;
            }
        }
        return sorted_data.back().first;
    };

    double p99 = find_percentile(0.99);
    double p99_9 = find_percentile(0.999);
    double p99_99 = find_percentile(0.9999);
    double median = find_percentile(0.5);

    // Calculate average
    double weighted_sum = 0.0;
    double interval_sum = 0.0;
    for (const auto& [value, interval] : raw_datas) {
        weighted_sum += value * interval;
        interval_sum += interval;
    }
    double average = (interval_sum > 0) ? (weighted_sum / interval_sum) : 0.0;

    return {average, median, p99, p99_9, p99_99};
}

void Utils::write_results_to_csv(
    const std::string& output_file,
    const std::map<std::string, nlohmann::json>& params_and_results) {

    // Filter out fields that should not be persisted
    std::map<std::string, nlohmann::json> filtered_results;
    for (const auto& [key, value] : params_and_results) {
        if (key != "df") {
            filtered_results[key] = value;
        }
    }

    // Check if file exists and is empty
    struct stat stat_buf;
    bool file_is_empty = (stat(output_file.c_str(), &stat_buf) != 0) || (stat_buf.st_size == 0);

    std::vector<std::string> header;

    if (file_is_empty) {
        // Create header from keys
        for (const auto& [key, value] : filtered_results) {
            header.push_back(key);
        }

        Logger::getInstance().info("Creating new CSV with " + std::to_string(header.size()) + " columns");

        // Write header
        std::ofstream file(output_file);
        if (!file.is_open()) {
            Logger::getInstance().error("Failed to create output file: " + output_file);
            return;
        }

        for (size_t i = 0; i < header.size(); ++i) {
            file << header[i];
            if (i < header.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
        file.close();
    } else {
        // Read existing header
        std::ifstream file(output_file);
        if (!file.is_open()) {
            Logger::getInstance().error("Failed to open output file for reading: " + output_file);
            return;
        }

        std::string header_line;
        if (std::getline(file, header_line)) {
            std::istringstream ss(header_line);
            std::string column;
            while (std::getline(ss, column, ',')) {
                header.push_back(column);
            }
        }
        file.close();

        // Validate that all keys exist in header
        for (const auto& [key, value] : filtered_results) {
            if (std::find(header.begin(), header.end(), key) == header.end()) {
                Logger::getInstance().error("CSV header missing expected field: " + key);
                throw std::runtime_error("CSV header missing expected field: " + key);
            }
        }
    }

    // Append row matching header order
    std::ofstream file(output_file, std::ios::app);
    if (!file.is_open()) {
        Logger::getInstance().error("Failed to open output file for appending: " + output_file);
        return;
    }

    for (size_t i = 0; i < header.size(); ++i) {
        const std::string& key = header[i];
        auto it = filtered_results.find(key);
        if (it != filtered_results.end()) {
            const nlohmann::json& value = it->second;

            // Handle different JSON types
            if (value.is_string()) {
                file << value.get<std::string>();
            } else if (value.is_number_float()) {
                file << std::fixed << std::setprecision(10) << value.get<double>();
            } else if (value.is_number_integer()) {
                file << value.get<int64_t>();
            } else if (value.is_boolean()) {
                file << (value.get<bool>() ? "true" : "false");
            } else {
                // For complex types, dump as JSON string
                file << value.dump();
            }
        }

        if (i < header.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    file.close();

    Logger::getInstance().info("Appended new row to " + output_file);
}
