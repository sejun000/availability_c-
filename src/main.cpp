#include <iostream>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include "utils.hpp"
#include "graph_structure.hpp"
#include "disk.hpp"
#include "erasure_coding.hpp"
#include "simulation.hpp"
#include "logger.hpp"

// Command line argument parser
struct Arguments {
    int total_disks = 48;
    int m = 15;
    int k = 1;
    int n = 0;  // Disk group size for declustered parity (0 = m+k)
    int local_m = 0;  // LRC / Multi-EC local data chunks
    int local_k = 0;  // LRC / Multi-EC local parity chunks
    int local_n = 0;  // Multi-EC local disk group size
    std::string ec_type = "standard";  // EC type: standard, lrc, multi_ec
    uint64_t capacity = 64000000000000ULL;
    bool qlc = false;
    double dwpd = 1.0;
    int guaranteed_years = 5;
    std::string config_file = "config.json";
    std::string output_file = "results.txt";
    std::string disk_failure_trace = "";  // Path to disk failure trace file
    int nprocs = 40;
    double rebuild_bw_ratio = 0.2;
    double degraded_ratio = 0.2;
    double target_performance_ratio = 0.0;  // 0 = disabled
    std::string data_loss_rebuild_bw_str = "100M";  // Fixed rebuild BW during data loss (default: 100M)
    int num_simulations = 400000;  // Number of Monte Carlo iterations
    int simulation_years = 10;     // Virtual years per simulation
    bool no_result = false;
    bool verbose = false;
};

void print_usage() {
    std::cerr << "Usage: availability_sim [options] <config_file>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  --total_disks <int>       Total number of disks (default: 48)" << std::endl;
    std::cerr << "  --m <int>                 Data chunks in erasure coding (default: 15)" << std::endl;
    std::cerr << "  --k <int>                 Parity chunks in erasure coding (default: 1)" << std::endl;
    std::cerr << "  --n <int>                 Disk group size for declustered parity (default: m+k)" << std::endl;
    std::cerr << "  --ec_type <string>        EC type: standard, lrc, multi_ec (default: standard)" << std::endl;
    std::cerr << "  --local_m <int>           Local data chunks for LRC/Multi-EC (default: 0)" << std::endl;
    std::cerr << "  --local_k <int>           Local parity chunks for LRC/Multi-EC (default: 0)" << std::endl;
    std::cerr << "  --local_n <int>           Local disk group size for Multi-EC (default: 0)" << std::endl;
    std::cerr << "  --capacity <uint64>       Disk capacity in bytes (default: 64TB)" << std::endl;
    std::cerr << "  --qlc                     Use QLC disks (default: TLC)" << std::endl;
    std::cerr << "  --dwpd <double>           Drive writes per day (default: 1.0)" << std::endl;
    std::cerr << "  --guaranteed_years <int>  Guaranteed years (default: 5)" << std::endl;
    std::cerr << "  --output_file <path>      Output file path (default: results.txt)" << std::endl;
    std::cerr << "  --disk_failure_trace <path> Disk failure trace file for empirical distribution" << std::endl;
    std::cerr << "  --nprocs <int>            Number of parallel processes (default: 40)" << std::endl;
    std::cerr << "  --rebuild_bw_ratio <double> Rebuild bandwidth ratio (default: 0.2)" << std::endl;
    std::cerr << "  --degraded_ratio <double>  Degraded read ratio (default: 0.2)" << std::endl;
    std::cerr << "  --target_perf <double>    Target performance ratio for perf_availability (default: 0)" << std::endl;
    std::cerr << "  --data_loss_rebuild_bw <string>  Fixed rebuild BW during data loss (default: 100M)" << std::endl;
    std::cerr << "  --num_sim <int>           Number of Monte Carlo iterations (default: 200000)" << std::endl;
    std::cerr << "  --sim_years <int>         Virtual years per simulation (default: 10)" << std::endl;
    std::cerr << "  --no_result               Don't write results to file" << std::endl;
    std::cerr << "  --verbose                 Enable verbose logging" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Example:" << std::endl;
    std::cerr << "  availability_sim --total_disks 96 --m 22 --k 2 config.json" << std::endl;
}

Arguments parse_arguments(int argc, char* argv[]) {
    Arguments args;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            print_usage();
            exit(0);
        }
        else if (arg == "--total_disks" && i + 1 < argc) args.total_disks = std::stoi(argv[++i]);
        else if (arg == "--m" && i + 1 < argc) args.m = std::stoi(argv[++i]);
        else if (arg == "--k" && i + 1 < argc) args.k = std::stoi(argv[++i]);
        else if (arg == "--n" && i + 1 < argc) args.n = std::stoi(argv[++i]);
        else if (arg == "--ec_type" && i + 1 < argc) args.ec_type = argv[++i];
        else if (arg == "--local_m" && i + 1 < argc) args.local_m = std::stoi(argv[++i]);
        else if (arg == "--local_k" && i + 1 < argc) args.local_k = std::stoi(argv[++i]);
        else if (arg == "--local_n" && i + 1 < argc) args.local_n = std::stoi(argv[++i]);
        else if (arg == "--capacity" && i + 1 < argc) args.capacity = std::stoull(argv[++i]);
        else if (arg == "--qlc") args.qlc = true;
        else if (arg == "--dwpd" && i + 1 < argc) args.dwpd = std::stod(argv[++i]);
        else if (arg == "--guaranteed_years" && i + 1 < argc) args.guaranteed_years = std::stoi(argv[++i]);
        else if (arg == "--output_file" && i + 1 < argc) args.output_file = argv[++i];
        else if (arg == "--disk_failure_trace" && i + 1 < argc) args.disk_failure_trace = argv[++i];
        else if (arg == "--nprocs" && i + 1 < argc) args.nprocs = std::stoi(argv[++i]);
        else if (arg == "--rebuild_bw_ratio" && i + 1 < argc) args.rebuild_bw_ratio = std::stod(argv[++i]);
        else if (arg == "--degraded_ratio" && i + 1 < argc) args.degraded_ratio = std::stod(argv[++i]);
        else if (arg == "--target_perf" && i + 1 < argc) args.target_performance_ratio = std::stod(argv[++i]);
        else if (arg == "--data_loss_rebuild_bw" && i + 1 < argc) args.data_loss_rebuild_bw_str = argv[++i];
        else if (arg == "--num_sim" && i + 1 < argc) args.num_simulations = std::stoi(argv[++i]);
        else if (arg == "--sim_years" && i + 1 < argc) args.simulation_years = std::stoi(argv[++i]);
        else if (arg == "--no_result") args.no_result = true;
        else if (arg == "--verbose") args.verbose = true;
        else if (arg[0] != '-') {
            // Positional argument: config file
            args.config_file = arg;
        }
        else {
            std::cerr << "Error: Unknown argument: " << arg << std::endl;
            print_usage();
            exit(1);
        }
    }

    return args;
}

int main(int argc, char* argv[]) {
    Arguments args = parse_arguments(argc, argv);

    // Configure logger: quiet by default, verbose with --verbose flag
    if (!args.verbose) {
        Logger::getInstance().setLevel(Logger::Level::ERROR);
    }

    // Parse configuration file
    std::vector<std::tuple<std::string, std::string, std::string>> edges;
    std::map<std::string, std::vector<std::string>> enclosures;
    std::map<std::string, double> mttfs;
    std::map<std::string, double> mtrs;
    std::map<std::string, double> costs;
    nlohmann::json options;

    if (!Utils::parse_input_from_json(args.config_file, edges, enclosures, mttfs, mtrs, costs, options)) {
        std::cerr << "Failed to parse configuration file: " << args.config_file << std::endl;
        return 1;
    }

    // Create hardware graph
    GraphStructure hardware_graph(edges, enclosures, mttfs, mtrs);

    // Get disk parameters from options
    double qlc_write_bw = Utils::KMG_to_bytes(options.value("qlc_write_bw", "500M"));
    double qlc_read_bw = Utils::KMG_to_bytes(options.value("qlc_read_bw", "2G"));
    double qlc_dwpd = options.value("qlc_dwpd_limit", 0.3);
    double tlc_write_bw = Utils::KMG_to_bytes(options.value("tlc_write_bw", "2G"));
    double tlc_read_bw = Utils::KMG_to_bytes(options.value("tlc_read_bw", "5G"));
    double tlc_dwpd = options.value("tlc_dwpd_limit", 1.0);

    double write_bw = args.qlc ? qlc_write_bw : tlc_write_bw;
    double read_bw = args.qlc ? qlc_read_bw : tlc_read_bw;
    double dwpd_limit = args.qlc ? qlc_dwpd : tlc_dwpd;

    // Set default n if not specified
    if (args.n == 0) {
        args.n = args.m + args.k;
    }

    // Validation
    if (args.n < args.m + args.k) {
        std::cerr << "Error: n must be >= m + k" << std::endl;
        return 1;
    }

    if (args.total_disks % args.n != 0) {
        std::cerr << "Warning: total_disks (" << args.total_disks << ") is not divisible by n (" << args.n << ")" << std::endl;
    }

    // Print configuration
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Total disks: " << args.total_disks << std::endl;
    std::cout << "  EC type: " << args.ec_type << std::endl;
    std::cout << "  m=" << args.m << ", k=" << args.k << ", n=" << args.n << std::endl;
    if (args.local_m > 0 || args.local_k > 0) {
        std::cout << "  local_m=" << args.local_m << ", local_k=" << args.local_k;
        if (args.local_n > 0) {
            std::cout << ", local_n=" << args.local_n;
        }
        std::cout << std::endl;
    }
    std::cout << "  Disk capacity: " << args.capacity / 1e12 << " TB" << std::endl;
    std::cout << "  Disk read BW: " << read_bw / 1e9 << " GB/s" << std::endl;
    std::cout << "  Disk write BW: " << write_bw / 1e9 << " GB/s" << std::endl;
    double disk_mttf_hours = args.guaranteed_years * 365.0 * 24.0 * dwpd_limit / args.dwpd;
    std::cout << "  Disk MTTF: " << disk_mttf_hours / (365.0 * 24.0) << " years (" << disk_mttf_hours << " hours)" << std::endl;

    // Calculate and show EC encoding speed (256KB chunk)
    const double EC_CHUNK_SIZE = 256.0 * 1024.0;
    double ec_latency_sec = Utils::get_encoding_latency_sec(args.m, args.k);
    double ec_encoding_speed = 0;
    if (ec_latency_sec > 0) {
        ec_encoding_speed = EC_CHUNK_SIZE / ec_latency_sec;
        std::cout << "  EC encoding latency: " << ec_latency_sec * 1e6 << " usec/256KB chunk" << std::endl;
        std::cout << "  EC encoding speed: " << ec_encoding_speed / 1e6 << " MB/s" << std::endl;
    }
    std::cout << "  Rebuild BW ratio: " << args.rebuild_bw_ratio << std::endl;
    std::cout << "  Degraded ratio: " << args.degraded_ratio << std::endl;

    // Prepare simulation parameters
    std::map<std::string, nlohmann::json> params_and_results;
    params_and_results["total_ssds"] = args.total_disks;  // Keep 'ssds' name for compatibility
    params_and_results["m"] = args.m;
    params_and_results["k"] = args.k;
    params_and_results["capacity"] = args.capacity;
    params_and_results["qlc"] = args.qlc;
    params_and_results["dwpd"] = args.dwpd;
    params_and_results["guaranteed_years"] = args.guaranteed_years;
    params_and_results["dwpd_limit"] = dwpd_limit;
    params_and_results["ssd_read_bw"] = read_bw;
    params_and_results["ssd_write_bw"] = write_bw;
    params_and_results["config_file"] = args.config_file;
    params_and_results["nprocs"] = args.nprocs;

    // Set EC-related options
    options["ec_type"] = args.ec_type;
    options["n"] = args.n;
    options["total_disks"] = args.total_disks;
    options["local_m"] = args.local_m;
    options["local_k"] = args.local_k;
    options["local_n"] = args.local_n;
    options["rebuild_bw_ratio"] = args.rebuild_bw_ratio;
    options["degraded_ratio"] = args.degraded_ratio;
    options["target_performance_ratio"] = args.target_performance_ratio;
    options["data_loss_rebuild_bw"] = Utils::KMG_to_bytes(args.data_loss_rebuild_bw_str);
    options["simulation_years"] = args.simulation_years;

    // Set disk failure trace file (CLI overrides JSON config)
    if (!args.disk_failure_trace.empty()) {
        options["ssd_failure_trace"] = args.disk_failure_trace;
    }

    // Run simulation
    monte_carlo_simulation(params_and_results, hardware_graph, args.num_simulations, options);

    // Output results
    std::cout << "\n========================================" << std::endl;
    std::cout << "Simulation Results:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Availability: " << params_and_results["availability"] << std::endl;
    std::cout << "  Availability (nines): " << params_and_results["avail_nines"] << std::endl;
    std::cout << "  Durability: " << params_and_results["durability"] << std::endl;
    std::cout << "  Durability (nines): " << params_and_results["durability_nines"] << std::endl;

    if (params_and_results.find("pdl_ratio") != params_and_results.end()) {
        double pdl = params_and_results["pdl_ratio"].get<double>();
        std::cout << "  PDL (data loss ratio): " << pdl << std::endl;
    }

    if (params_and_results.find("mttdl_per_ec_group") != params_and_results.end() &&
        params_and_results["mttdl_per_ec_group"].get<double>() > 0) {
        double mttdl_per_group = params_and_results["mttdl_per_ec_group"].get<double>();
        int num_groups = 1;
        if (params_and_results.find("num_ec_groups") != params_and_results.end()) {
            num_groups = params_and_results["num_ec_groups"].get<int>();
        }
        std::cout << "  MTTDL (per EC group): " << mttdl_per_group << " hours ("
                  << (mttdl_per_group / 8760.0) << " years)" << std::endl;
        std::cout << "  MTTDL (system): " << params_and_results["mttdl_system"].get<double>()
                  << " hours (" << num_groups << " EC groups)" << std::endl;
    }

    if (params_and_results.find("data_loss_events") != params_and_results.end()) {
        std::cout << "  Data loss events: " << params_and_results["data_loss_events"] << std::endl;
        if (params_and_results.find("data_loss_events_per_year") != params_and_results.end()) {
            std::cout << "  Data loss events/year: "
                      << params_and_results["data_loss_events_per_year"].get<double>() << std::endl;
        }
    }

    if (params_and_results.find("disk_failure_events") != params_and_results.end()) {
        std::cout << "  Disk failure events: " << params_and_results["disk_failure_events"] << std::endl;
    }

    if (params_and_results.find("avg_rebuilding_time") != params_and_results.end()) {
        std::cout << "  Avg rebuild time: " << params_and_results["avg_rebuilding_time"].get<double>() << " hours" << std::endl;
    }

    // Rebuild-specific metrics
    if (params_and_results.find("avg_rebuild_speed") != params_and_results.end()) {
        double avg_rebuild_speed = params_and_results["avg_rebuild_speed"].get<double>();
        std::cout << "  Avg rebuild speed: " << (avg_rebuild_speed / 1e6) << " MB/s" << std::endl;
    }
    if (params_and_results.find("avg_host_read_bw_during_rebuild") != params_and_results.end()) {
        double avg_read_bw = params_and_results["avg_host_read_bw_during_rebuild"].get<double>();
        std::cout << "  Avg host READ BW during rebuild: " << (avg_read_bw / 1e9) << " GB/s" << std::endl;
    }
    if (params_and_results.find("avg_host_write_bw_during_rebuild") != params_and_results.end()) {
        double avg_write_bw = params_and_results["avg_host_write_bw_during_rebuild"].get<double>();
        std::cout << "  Avg host WRITE BW during rebuild: " << (avg_write_bw / 1e9) << " GB/s" << std::endl;
    }

    // Performance availability (if target_perf specified)
    if (args.target_performance_ratio > 0) {
        if (params_and_results.find("perf_availability_read") != params_and_results.end()) {
            std::cout << "  Perf Availability READ (>=" << args.target_performance_ratio * 100 << "%): "
                      << params_and_results["perf_availability_read"] << std::endl;
            std::cout << "  Perf Availability READ (nines): " << params_and_results["perf_avail_nines_read"] << std::endl;
        }
        if (params_and_results.find("perf_availability_write") != params_and_results.end()) {
            std::cout << "  Perf Availability WRITE (>=" << args.target_performance_ratio * 100 << "%): "
                      << params_and_results["perf_availability_write"] << std::endl;
            std::cout << "  Perf Availability WRITE (nines): " << params_and_results["perf_avail_nines_write"] << std::endl;
        }
        if (params_and_results.find("perf_availability_min") != params_and_results.end()) {
            std::cout << "  Perf Availability MIN (>=" << args.target_performance_ratio * 100 << "%): "
                      << params_and_results["perf_availability_min"] << std::endl;
            std::cout << "  Perf Availability MIN (nines): " << params_and_results["perf_avail_nines_min"] << std::endl;
        }
    }

    // Write to output file if needed
    if (!args.no_result) {
        Utils::write_results_to_csv(args.output_file, params_and_results);
        std::cout << "  Results written to: " << args.output_file << std::endl;
    }

    return 0;
}
