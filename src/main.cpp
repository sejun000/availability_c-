#include <iostream>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include "utils.hpp"
#include "graph_structure.hpp"
#include "ssd.hpp"
#include "simulation.hpp"

// Command line argument parser
struct Arguments {
    int total_ssds = 48;
    int m = 15;
    int k = 1;
    int l = 0;
    int cached_ssds = 0;
    int cached_m = 0;
    int cached_k = 0;
    int cached_l = 0;
    int cached_network_m = 0;
    int cached_network_k = 0;
    int cached_network_l = 0;
    int inter_replicas = 0;
    int intra_replicas = 0;
    double cached_write_ratio = 0.0;
    int network_m = 8;
    int network_k = 0;
    int network_l = 0;
    uint64_t capacity = 64000000000000ULL;
    bool qlc = false;
    bool simulation = false;
    double dwpd = 1.0;
    int guaranteed_years = 5;
    std::string config_file = "2tier.json";
    std::string output_file = "results.txt";
    std::string ssd_failure_trace = "";  // Path to SSD failure trace file
    bool qlc_cache = false;
    int nprocs = 20;
    double box_mttf = 0.0;
    double io_module_mttr = 0.0;
    double rebuild_bw_ratio = 0.2;
    bool no_result = false;
    double target_perf_ratio = 0.8;
    bool single_port_ssd = false;
    bool active_active = false;
};

Arguments parse_arguments(int argc, char* argv[]) {
    Arguments args;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--total_ssds" && i + 1 < argc) args.total_ssds = std::stoi(argv[++i]);
        else if (arg == "--m" && i + 1 < argc) args.m = std::stoi(argv[++i]);
        else if (arg == "--k" && i + 1 < argc) args.k = std::stoi(argv[++i]);
        else if (arg == "--l" && i + 1 < argc) args.l = std::stoi(argv[++i]);
        else if (arg == "--cached_ssds" && i + 1 < argc) args.cached_ssds = std::stoi(argv[++i]);
        else if (arg == "--cached_m" && i + 1 < argc) args.cached_m = std::stoi(argv[++i]);
        else if (arg == "--cached_k" && i + 1 < argc) args.cached_k = std::stoi(argv[++i]);
        else if (arg == "--cached_l" && i + 1 < argc) args.cached_l = std::stoi(argv[++i]);
        else if (arg == "--network_m" && i + 1 < argc) args.network_m = std::stoi(argv[++i]);
        else if (arg == "--network_k" && i + 1 < argc) args.network_k = std::stoi(argv[++i]);
        else if (arg == "--network_l" && i + 1 < argc) args.network_l = std::stoi(argv[++i]);
        else if (arg == "--inter_replicas" && i + 1 < argc) args.inter_replicas = std::stoi(argv[++i]);
        else if (arg == "--intra_replicas" && i + 1 < argc) args.intra_replicas = std::stoi(argv[++i]);
        else if (arg == "--cached_write_ratio" && i + 1 < argc) args.cached_write_ratio = std::stod(argv[++i]);
        else if (arg == "--capacity" && i + 1 < argc) args.capacity = std::stoull(argv[++i]);
        else if (arg == "--qlc") args.qlc = true;
        else if (arg == "--simulation") args.simulation = true;
        else if (arg == "--dwpd" && i + 1 < argc) args.dwpd = std::stod(argv[++i]);
        else if (arg == "--guaranteed_years" && i + 1 < argc) args.guaranteed_years = std::stoi(argv[++i]);
        else if (arg == "--config_file" && i + 1 < argc) args.config_file = argv[++i];
        else if (arg == "--output_file" && i + 1 < argc) args.output_file = argv[++i];
        else if (arg == "--ssd_failure_trace" && i + 1 < argc) args.ssd_failure_trace = argv[++i];
        else if (arg == "--qlc_cache") args.qlc_cache = true;
        else if (arg == "--nprocs" && i + 1 < argc) args.nprocs = std::stoi(argv[++i]);
        else if (arg == "--box_mttf" && i + 1 < argc) args.box_mttf = std::stod(argv[++i]);
        else if (arg == "--io_module_mttr" && i + 1 < argc) args.io_module_mttr = std::stod(argv[++i]);
        else if (arg == "--rebuild_bw_ratio" && i + 1 < argc) args.rebuild_bw_ratio = std::stod(argv[++i]);
        else if (arg == "--no_result") args.no_result = true;
        else if (arg == "--target_perf_ratio" && i + 1 < argc) args.target_perf_ratio = std::stod(argv[++i]);
        else if (arg == "--single_port_ssd") args.single_port_ssd = true;
        else if (arg == "--active_active") args.active_active = true;
    }

    return args;
}

int main(int argc, char* argv[]) {
    Arguments args = parse_arguments(argc, argv);

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

    // Get SSD parameters from options
    double qlc_write_bw = Utils::KMG_to_bytes(options["qlc_write_bw"]);
    double qlc_read_bw = Utils::KMG_to_bytes(options["qlc_read_bw"]);
    double qlc_dwpd = options["qlc_dwpd_limit"];
    double tlc_write_bw = Utils::KMG_to_bytes(options["tlc_write_bw"]);
    double tlc_read_bw = Utils::KMG_to_bytes(options["tlc_read_bw"]);
    double tlc_dwpd = options["tlc_dwpd_limit"];

    double write_bw = args.qlc ? qlc_write_bw : tlc_write_bw;
    double read_bw = args.qlc ? qlc_read_bw : tlc_read_bw;
    double dwpd_limit = args.qlc ? qlc_dwpd : tlc_dwpd;

    // Validation
    int n = args.m + args.k + args.l;
    if (n > args.total_ssds) {
        std::cerr << "Error: The sum of m, k, l should not exceed total_ssds" << std::endl;
        return 1;
    }

    if ((args.total_ssds - args.cached_ssds) % n != 0) {
        std::cerr << "Error: total_ssds should be divisible by the sum of m, k, l" << std::endl;
        return 1;
    }

    // Setup cached network parameters
    int cached_network_m = args.cached_network_m;
    if (cached_network_m == 0) {
        cached_network_m = args.network_m;
    }

    // Handle replica configuration
    if (args.inter_replicas > 1) {
        args.network_m = 1;
        args.network_k = args.inter_replicas - 1;
        args.network_l = 0;
        std::cout << "Network m and k are calculated as 1 and " << args.inter_replicas - 1 << std::endl;
    }

    if (args.cached_ssds > 0 && args.intra_replicas > 1) {
        args.cached_m = 1;
        args.cached_k = args.intra_replicas - 1;
        args.cached_l = 0;
        std::cout << "Cached m and k are calculated as 1 and " << args.intra_replicas - 1 << std::endl;
    }

    // Print configuration
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Total SSDs: " << args.total_ssds << std::endl;
    std::cout << "  m=" << args.m << ", k=" << args.k << ", l=" << args.l << std::endl;
    std::cout << "  Network: m=" << args.network_m << ", k=" << args.network_k << std::endl;
    std::cout << "  Simulation: " << (args.simulation ? "enabled" : "disabled") << std::endl;

    if (args.simulation) {
        // Prepare simulation parameters
        std::map<std::string, nlohmann::json> params_and_results;
        params_and_results["total_ssds"] = args.total_ssds;
        params_and_results["m"] = args.m;
        params_and_results["k"] = args.k;
        params_and_results["l"] = args.l;
        params_and_results["cached_ssds"] = args.cached_ssds;
        params_and_results["cached_m"] = args.cached_m;
        params_and_results["cached_k"] = args.cached_k;
        params_and_results["cached_l"] = args.cached_l;
        params_and_results["cached_network_m"] = cached_network_m;
        params_and_results["cached_network_k"] = args.cached_network_k;
        params_and_results["cached_network_l"] = args.cached_network_l;
        params_and_results["inter_replicas"] = args.inter_replicas;
        params_and_results["intra_replicas"] = args.intra_replicas;
        params_and_results["cached_write_ratio"] = args.cached_write_ratio;
        params_and_results["network_m"] = args.network_m;
        params_and_results["network_k"] = args.network_k;
        params_and_results["network_l"] = args.network_l;
        params_and_results["capacity"] = args.capacity;
        params_and_results["qlc"] = args.qlc;
        params_and_results["simulation"] = args.simulation;
        params_and_results["dwpd"] = args.dwpd;
        params_and_results["guaranteed_years"] = args.guaranteed_years;
        params_and_results["dwpd_limit"] = dwpd_limit;
        params_and_results["ssd_read_bw"] = read_bw;
        params_and_results["ssd_write_bw"] = write_bw;

        if (args.qlc_cache) {
            params_and_results["qlc_cache"] = true;
            params_and_results["cached_dwpd_limit"] = qlc_dwpd;
            params_and_results["cached_ssd_read_bw"] = qlc_read_bw;
            params_and_results["cached_ssd_write_bw"] = qlc_write_bw;
        } else {
            params_and_results["qlc_cache"] = false;
            params_and_results["cached_dwpd_limit"] = tlc_dwpd;
            params_and_results["cached_ssd_read_bw"] = tlc_read_bw;
            params_and_results["cached_ssd_write_bw"] = tlc_write_bw;
        }

        params_and_results["config_file"] = args.config_file;
        params_and_results["nprocs"] = args.nprocs;
        params_and_results["box_mttf"] = args.box_mttf;
        params_and_results["io_module_mttr"] = args.io_module_mttr;
        params_and_results["active_active"] = args.active_active;
        params_and_results["single_port_ssd"] = args.single_port_ssd;
        params_and_results["rebuild_bw_ratio"] = args.rebuild_bw_ratio;
        params_and_results["target_perf_ratio"] = args.target_perf_ratio;

        // Set SSD failure trace file (CLI overrides JSON config)
        if (!args.ssd_failure_trace.empty()) {
            options["ssd_failure_trace"] = args.ssd_failure_trace;
        }

        // Run simulation
        int num_simulations = 80000;
        monte_carlo_simulation(params_and_results, hardware_graph, num_simulations, options, costs);

        // Output results
        std::cout << "\nSimulation Results:" << std::endl;
        std::cout << "  Availability: " << params_and_results["availability"] << std::endl;
        std::cout << "  Availability (nines): " << params_and_results["avail_nines"] << std::endl;

        // Write to output file if needed
        if (!args.no_result) {
            Utils::write_results_to_csv(args.output_file, params_and_results);
            std::cout << "  Results written to: " << args.output_file << std::endl;
        }
    }

    return 0;
}
