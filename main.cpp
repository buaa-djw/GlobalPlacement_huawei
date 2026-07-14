#include "db/PlacementDB.h"
#include "evaluator/HPWLEvaluator.h"
#include "evaluator/DensityEvaluator.h"
#include "evaluator/ObjectiveEvaluator.h"
#include "optimizer/GlobalPlacer.h"
#include "grid/BinGrid.h"
#include "parser/LimboBookshelfAdapter.h"
#include "utils/Logger.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {
struct Options {
    std::string aux_path;
    int bins_x = 64;
    int bins_y = 64;
    double target_density = 0.9;
    std::string log_path = "../result/global_placer.log";
    LogLevel log_level = LogLevel::Info;
    GlobalPlacerConfig placer_config;
};

void printUsage(const char* argv0) {
    std::cerr << "Usage:\n  " << argv0
              << " --aux <path> [--bins <nx> <ny>] [--target-density <value>]"
              << " [--log <path>] [--log-level <debug|info|warn|error>]"
              << " [--iterations <int>] [--density-weight <value>]"
              << " [--initial-step-fraction <value>] [--minimum-step-fraction <value>]"
              << " [--backtrack-factor <value>] [--line-search-trials <int>]"
              << " [--stall-iterations <int>] [--report-interval <int>]\n";
}

Options parseArgs(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--aux" && i + 1 < argc) opt.aux_path = argv[++i];
        else if (arg == "--bins" && i + 2 < argc) { opt.bins_x = std::stoi(argv[++i]); opt.bins_y = std::stoi(argv[++i]); }
        else if (arg == "--target-density" && i + 1 < argc) opt.target_density = std::stod(argv[++i]);
        else if (arg == "--log" && i + 1 < argc) opt.log_path = argv[++i];
        else if (arg == "--log-level" && i + 1 < argc) opt.log_level = parseLogLevel(argv[++i]);
        else if (arg == "--iterations" && i + 1 < argc) opt.placer_config.max_iterations = std::stoi(argv[++i]);
        else if (arg == "--density-weight" && i + 1 < argc) opt.placer_config.density_weight = std::stod(argv[++i]);
        else if (arg == "--initial-step-fraction" && i + 1 < argc) opt.placer_config.initial_step_fraction = std::stod(argv[++i]);
        else if (arg == "--minimum-step-fraction" && i + 1 < argc) opt.placer_config.minimum_step_fraction = std::stod(argv[++i]);
        else if (arg == "--backtrack-factor" && i + 1 < argc) opt.placer_config.backtrack_factor = std::stod(argv[++i]);
        else if (arg == "--line-search-trials" && i + 1 < argc) opt.placer_config.max_line_search_trials = std::stoi(argv[++i]);
        else if (arg == "--stall-iterations" && i + 1 < argc) opt.placer_config.max_stall_iterations = std::stoi(argv[++i]);
        else if (arg == "--report-interval" && i + 1 < argc) opt.placer_config.report_interval = std::stoi(argv[++i]);
        else if (arg == "--help" || arg == "-h") { printUsage(argv[0]); std::exit(EXIT_SUCCESS); }
        else throw std::invalid_argument("unknown or incomplete argument: " + arg);
    }
    if (opt.aux_path.empty()) throw std::invalid_argument("missing required --aux <path>");
    opt.placer_config.bins_x = opt.bins_x; opt.placer_config.bins_y = opt.bins_y; opt.placer_config.target_density = opt.target_density;
    return opt;
}

std::string binGridSummary(const BinGrid& grid) {
    const Box& c = grid.coreBounds();
    std::ostringstream os;
    os << std::fixed << std::setprecision(3)
       << "========== BinGrid Summary ==========" << '\n'
       << "Core             : (" << c.lx << ", " << c.ly << ") - (" << c.ux << ", " << c.uy << ")\n"
       << "Grid             : " << grid.numBinsX() << " x " << grid.numBinsY() << '\n'
       << "Bin size         : " << grid.binWidth() << " x " << grid.binHeight() << '\n'
       << "Target density   : " << grid.targetDensity() << '\n'
       << "Movable area     : " << grid.totalMovableArea() << '\n'
       << "Fixed area       : " << grid.totalFixedArea() << '\n'
       << "Total overflow   : " << grid.totalOverflow() << '\n'
       << "Overflow bins    : " << grid.overflowBinCount() << '\n'
       << "Max utilization  : " << grid.maxUtilization() << '\n'
       << "=====================================";
    return os.str();
}
std::string densitySummary(const DensityMetrics& metrics) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6)
       << "========== Density Evaluation ==========" << '\n'
       << "Density penalty        : " << metrics.penalty << '\n'
       << "Total overflow         : " << metrics.total_overflow << '\n'
       << "Overflow ratio         : " << metrics.overflow_ratio << '\n'
       << "Overflow bins          : " << metrics.overflow_bin_count << '\n'
       << "Movable area           : " << metrics.total_movable_area << '\n'
       << "Fixed area             : " << metrics.total_fixed_area << '\n'
       << "Movable capacity       : " << metrics.total_movable_capacity << '\n'
       << "Maximum utilization    : " << metrics.max_utilization << '\n'
       << "Average utilization    : " << metrics.average_utilization << '\n'
       << "Zero-capacity bins     : " << metrics.zero_capacity_bin_count << '\n'
       << "Occupied zero-cap bins : " << metrics.zero_capacity_occupied_bin_count << '\n'
       << "========================================";
    return os.str();
}

void writeDensitySummary(std::ostream& os, const DensityMetrics& metrics) {
    os << std::fixed << std::setprecision(6)
       << "\n[Density Evaluation]\n"
       << "Density penalty: " << metrics.penalty << '\n'
       << "Total overflow: " << metrics.total_overflow << '\n'
       << "Overflow ratio: " << metrics.overflow_ratio << '\n'
       << "Overflow bins: " << metrics.overflow_bin_count << '\n'
       << "Total movable area: " << metrics.total_movable_area << '\n'
       << "Total fixed area: " << metrics.total_fixed_area << '\n'
       << "Total movable capacity: " << metrics.total_movable_capacity << '\n'
       << "Maximum utilization: " << metrics.max_utilization << '\n'
       << "Average utilization: " << metrics.average_utilization << '\n'
       << "Zero-capacity bins: " << metrics.zero_capacity_bin_count << '\n'
       << "Zero-capacity occupied bins: " << metrics.zero_capacity_occupied_bin_count << '\n';
}

}

int main(int argc, char** argv) {
    bool logger_ready = false;
    try {
        const Options opt = parseArgs(argc, argv);
        Logger::instance().initialize(opt.log_path, opt.log_level, true);
        logger_ready = true;
        LOG_INFO("program started");
        LOG_INFO("resolved command-line arguments: aux=" << opt.aux_path << ", bins=" << opt.bins_x << "x" << opt.bins_y << ", target_density=" << opt.target_density);
        LOG_INFO("log file path: " << opt.log_path);
        LOG_INFO("benchmark path: " << opt.aux_path);
        if (opt.bins_x <= 0 || opt.bins_y <= 0 || opt.target_density <= 0.0 || opt.target_density > 1.0) {
            printUsage(argv[0]);
            throw std::invalid_argument("illegal BinGrid arguments: bins=" +
                                        std::to_string(opt.bins_x) + "x" + std::to_string(opt.bins_y) +
                                        ", target_density=" + std::to_string(opt.target_density));
        }

        //读取文件
        PlacementDB db;
        LimboBookshelfAdapter adapter(db);
        LOG_INFO("Bookshelf parsing started");
        if (!adapter.read(opt.aux_path)) LOG_FATAL("Failed to read Bookshelf design: " << opt.aux_path);
        LOG_INFO("Bookshelf parsing completed");
        LOG_INFO("cell/net/pin/row counts: cells=" << db.cells().size() << ", nets=" << db.nets().size() << ", pins=" << db.pins().size() << ", rows=" << db.rows().size());

        //初始化布局
        LOG_INFO("Generating deterministic initial placement");
        db.initializeSimplePlacement();
        LOG_INFO(
            "Deterministic initial placement completed"
        );

        //线长计算
        LOG_INFO("HPWL evaluation started");
        HPWLEvaluator hpwl_evaluator;
        double total_hpwl = 0.0;
        for (std::size_t i = 0; i < db.nets().size(); ++i) {
            const int net_id = static_cast<int>(i);
            const double net_hpwl = hpwl_evaluator.netHPWL(db, net_id);
            db.addNetHPWL(net_id, net_hpwl);
            total_hpwl += net_hpwl;
        }
        LOG_INFO("total HPWL: " << std::fixed << std::setprecision(3) << total_hpwl);

        //开始搭建bin，并计算容量
        LOG_INFO("BinGrid construction started");
        BinGrid grid(db, opt.bins_x, opt.bins_y, opt.target_density);
        const Box& core = grid.coreBounds();
        LOG_INFO("core bounds: (" << core.lx << ", " << core.ly << ") - (" << core.ux << ", " << core.uy << ")");
        LOG_INFO("bin count: " << grid.numBins() << ", bin size=" << grid.binWidth() << "x" << grid.binHeight() << ", target density=" << grid.targetDensity());
        LOG_INFO("BinGrid totals: movable_area=" << grid.totalMovableArea() << ", fixed_area=" << grid.totalFixedArea() << ", total_overflow=" << grid.totalOverflow() << ", overflow_bins=" << grid.overflowBinCount() << ", max_utilization=" << grid.maxUtilization());

        //计算每一个bin的密度
        LOG_INFO("Density evaluation started");
        DensityEvaluator density_evaluator;
        const DensityMetrics density_metrics = density_evaluator.evaluate(grid);
        LOG_INFO("Density evaluation completed");
        LOG_INFO("Density evaluation summary: penalty=" << density_metrics.penalty
                 << ", total_overflow=" << density_metrics.total_overflow
                 << ", overflow_ratio=" << density_metrics.overflow_ratio
                 << ", overflow_bins=" << density_metrics.overflow_bin_count
                 << ", total_movable_area=" << density_metrics.total_movable_area
                 << ", total_fixed_area=" << density_metrics.total_fixed_area
                 << ", total_movable_capacity=" << density_metrics.total_movable_capacity
                 << ", max_utilization=" << density_metrics.max_utilization
                 << ", average_utilization=" << density_metrics.average_utilization
                 << ", zero_capacity_bins=" << density_metrics.zero_capacity_bin_count
                 << ", zero_capacity_occupied_bins=" << density_metrics.zero_capacity_occupied_bin_count);

        
        //开始进行优化
        ObjectiveEvaluator objective_evaluator;
        ObjectiveMetrics initial_objective = objective_evaluator.evaluate(db, grid, opt.placer_config.density_weight);

        //开始全局布局
        GlobalPlacer placer(opt.placer_config);
        GlobalPlacerResult placement_result = placer.optimize(db);

        //记录最终结果
        BinGrid final_grid(db, opt.bins_x, opt.bins_y, opt.target_density);
        const DensityMetrics final_density_metrics = density_evaluator.evaluate(final_grid);
        const ObjectiveMetrics final_objective = objective_evaluator.evaluate(db, final_grid, opt.placer_config.density_weight);
        const double improvement = (initial_objective.total_cost - final_objective.total_cost) / std::max(std::abs(initial_objective.total_cost), 1e-12);

        //打印summary
        std::filesystem::create_directories("../result");
        const std::string placement_summary = "../result/placementdb_summary.txt";
        std::ofstream fout(placement_summary);
        if (!fout) LOG_FATAL("cannot open summary output: " << placement_summary);
        db.printSummary(fout);
        fout << "\n[Initial Placement Evaluation]\n" << std::fixed << std::setprecision(3) << "Total HPWL: " << total_hpwl << "\n";
        writeDensitySummary(fout, density_metrics);
        fout << "\n[Global Placement Result]\n"
             << "Initial total cost: " << initial_objective.total_cost << "\n"
             << "Final total cost: " << final_objective.total_cost << "\n"
             << "Accepted iterations: " << placement_result.accepted_iterations << "\n"
             << "Termination reason: " << placement_result.termination_reason << "\n";
        fout.close();

        const std::string grid_summary_path = "../result/bin_grid_summary.txt";
        const std::string final_grid_summary_path = "../result/final_bin_grid_summary.txt";
        std::ofstream gout(grid_summary_path);
        std::ofstream fgout(final_grid_summary_path);

        if (!gout) LOG_FATAL("cannot open BinGrid summary output: " << grid_summary_path);
        gout << binGridSummary(grid) << '\n' << densitySummary(density_metrics) << '\n';
        LOG_INFO("summary output path: " << grid_summary_path);

        if (!fgout) LOG_FATAL("cannot open BinGrid summary output: " << final_grid_summary_path);
        fgout << binGridSummary(final_grid) << '\n' << densitySummary(final_density_metrics) << '\n';
        LOG_INFO("summary output path: " << final_grid_summary_path);

        db.printSummary(std::cout);
        std::cout << "========== Initial Placement Evaluation ==========" << '\n'
                  << std::fixed << std::setprecision(3) << "Total HPWL: " << total_hpwl << '\n'
                  << "==================================================" << '\n'
                  << binGridSummary(grid) << '\n'
                  << densitySummary(density_metrics) << '\n'
                  << binGridSummary(final_grid) << '\n'
                  << densitySummary(final_density_metrics) << '\n'
                  << "========== Global Placement Result ==========\n"
                  << "Initial HPWL: " << initial_objective.hpwl << "\n"
                  << "Final HPWL: " << final_objective.hpwl << "\n"
                  << "HPWL change: " << (final_objective.hpwl - initial_objective.hpwl) << "\n"
                  << "Initial density penalty: " << initial_objective.density_penalty << "\n"
                  << "Final density penalty: " << final_objective.density_penalty << "\n"
                  << "Initial total cost: " << initial_objective.total_cost << "\n"
                  << "Final total cost: " << final_objective.total_cost << "\n"
                  << "Total cost improvement ratio: " << improvement << "\n"
                  << "Initial overflow ratio: " << initial_objective.overflow_ratio << "\n"
                  << "Final overflow ratio: " << final_objective.overflow_ratio << "\n"
                  << "Attempted iterations: " << placement_result.attempted_iterations << "\n"
                  << "Accepted iterations: " << placement_result.accepted_iterations << "\n"
                  << "Termination reason: " << placement_result.termination_reason << "\n"
                  << "=============================================\n";
        LOG_INFO("program completed successfully");
        Logger::instance().flush();
        return EXIT_SUCCESS;
    } catch (const FatalLogError&) {
        if (logger_ready) Logger::instance().flush();
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        if (logger_ready) { LOG_ERROR("Unhandled exception: " << e.what()); Logger::instance().flush(); }
        else { std::cerr << "Error before logger initialization: " << e.what() << '\n'; printUsage(argv[0]); }
        return EXIT_FAILURE;
    } catch (...) {
        if (logger_ready) { LOG_ERROR("Unhandled unknown exception"); Logger::instance().flush(); }
        else std::cerr << "Unknown error before logger initialization\n";
        return EXIT_FAILURE;
    }
}
