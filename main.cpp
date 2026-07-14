#include "db/PlacementDB.h"
#include "evaluator/HPWLEvaluator.h"
#include "evaluator/DensityEvaluator.h"
#include "evaluator/ObjectiveEvaluator.h"
#include "optimizer/GlobalPlacer.h"
#include "grid/BinGrid.h"
#include "initializer/QuadraticInitialPlacer.h"
#include "initializer/SimplePackingInitialPlacer.h"
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
#include <memory>

namespace
{
    struct Options
    {
        std::string aux_path;
        int bins_x = 64;
        int bins_y = 64;
        double target_density = 0.9;
        std::string log_path = "../result/global_placer.log";
        LogLevel log_level = LogLevel::Info;
        GlobalPlacerConfig placer_config;
        InitialPlacementMethod initial_method = InitialPlacementMethod::SimplePacking;
        SimplePackingInitialPlacerConfig simple_initial_config;
        QuadraticInitialPlacerConfig quadratic_initial_config;
    };

    void printUsage(const char *argv0)
    {
        std::cerr << "Usage:\n  " << argv0
                  << " --aux <path> [--bins <nx> <ny>] [--target-density <value>]\n"
                  << " --log <path> [--log-level <debug|info|warn|error>]\n"
                  << " --optimizer <direct|moreau> [--iterations <int>] [--density-weight <value>]\n"
                  << " --initial-placement <simple|quadratic> [--initial-max-net-degree <int>]\n"
                  << " --initial-anchor-weight-ratio <value> [--initial-blend-factor <value>]\n"
                  << " --initial-cg-iterations <int> [--initial-cg-tolerance <value>]\n"
                  << " --initial-require-convergence [--initial-ignore-pin-offsets]\n"
                  << " --zero-capacity-repulsion <value>\n"
                  << " --initial-step-fraction <value> [--minimum-step-fraction <value>]\n"
                  << " --backtrack-factor <value> [--line-search-trials <int>]\n"
                  << " --stall-iterations <int> [--report-interval <int>]\n"
                  << " --moreau-mu <value> [--moreau-inner-iterations <int>]\n"
                  << " --moreau-initial-step-fraction <value> [--moreau-minimum-step-fraction <value>]\n"
                  << " --moreau-backtrack-factor <value> [--moreau-line-search-trials <int>]\n"
                  << " --moreau-stall-iterations <int> [--moreau-residual-tolerance <value>]\n"
                  << " --moreau-displacement-tolerance <value>\n";
    }

    Options parseArgs(int argc, char **argv)
    {
        Options opt;
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            if (arg == "--aux" && i + 1 < argc)
                opt.aux_path = argv[++i];
            else if (arg == "--bins" && i + 2 < argc)
            {
                opt.bins_x = std::stoi(argv[++i]);
                opt.bins_y = std::stoi(argv[++i]);
            }
            else if (arg == "--target-density" && i + 1 < argc)
                opt.target_density = std::stod(argv[++i]);
            else if (arg == "--log" && i + 1 < argc)
                opt.log_path = argv[++i];
            else if (arg == "--log-level" && i + 1 < argc)
                opt.log_level = parseLogLevel(argv[++i]);
            else if (arg == "--optimizer" && i + 1 < argc)
            {
                const std::string v = argv[++i];
                if (v == "direct")
                    opt.placer_config.method = GlobalPlacementMethod::DirectSubgradient;
                else if (v == "moreau")
                    opt.placer_config.method = GlobalPlacementMethod::MoreauProximal;
                else
                    throw std::invalid_argument("invalid optimizer: " + v + " (expected direct or moreau)");
            }
            else if (arg == "--initial-placement" && i + 1 < argc)
            {
                const std::string v = argv[++i];
                if (v == "simple") opt.initial_method = InitialPlacementMethod::SimplePacking;
                else if (v == "quadratic") opt.initial_method = InitialPlacementMethod::QuadraticConnectivity;
                else throw std::invalid_argument("invalid initial placement: " + v + " (expected simple or quadratic)");
            }
            else if (arg == "--initial-max-net-degree" && i + 1 < argc)
                opt.quadratic_initial_config.maximum_net_degree = std::stoi(argv[++i]);
            else if (arg == "--initial-anchor-weight-ratio" && i + 1 < argc)
                opt.quadratic_initial_config.anchor_weight_ratio = std::stod(argv[++i]);
            else if (arg == "--initial-blend-factor" && i + 1 < argc)
                opt.quadratic_initial_config.blend_factor = std::stod(argv[++i]);
            else if (arg == "--initial-cg-iterations" && i + 1 < argc)
                opt.quadratic_initial_config.solver.maximum_iterations = std::stoi(argv[++i]);
            else if (arg == "--initial-cg-tolerance" && i + 1 < argc)
                opt.quadratic_initial_config.solver.relative_tolerance = std::stod(argv[++i]);
            else if (arg == "--initial-require-convergence")
                opt.quadratic_initial_config.require_solver_convergence = true;
            else if (arg == "--initial-ignore-pin-offsets")
                opt.quadratic_initial_config.use_pin_offsets = false;
            else if (arg == "--iterations" && i + 1 < argc)
                opt.placer_config.max_iterations = std::stoi(argv[++i]);
            else if (arg == "--density-weight" && i + 1 < argc)
                opt.placer_config.density_weight = std::stod(argv[++i]);
            else if (arg == "--zero-capacity-repulsion" && i + 1 < argc)
                opt.placer_config.zero_capacity_repulsion = std::stod(argv[++i]);
            else if (arg == "--initial-step-fraction" && i + 1 < argc)
                opt.placer_config.initial_step_fraction = std::stod(argv[++i]);
            else if (arg == "--minimum-step-fraction" && i + 1 < argc)
                opt.placer_config.minimum_step_fraction = std::stod(argv[++i]);
            else if (arg == "--backtrack-factor" && i + 1 < argc)
                opt.placer_config.backtrack_factor = std::stod(argv[++i]);
            else if (arg == "--line-search-trials" && i + 1 < argc)
                opt.placer_config.max_line_search_trials = std::stoi(argv[++i]);
            else if (arg == "--stall-iterations" && i + 1 < argc)
                opt.placer_config.max_stall_iterations = std::stoi(argv[++i]);
            else if (arg == "--report-interval" && i + 1 < argc)
                opt.placer_config.report_interval = std::stoi(argv[++i]);
            else if (arg == "--moreau-mu" && i + 1 < argc)
                opt.placer_config.moreau.mu = std::stod(argv[++i]);
            else if (arg == "--moreau-inner-iterations" && i + 1 < argc)
                opt.placer_config.moreau.max_inner_iterations = std::stoi(argv[++i]);
            else if (arg == "--moreau-initial-step-fraction" && i + 1 < argc)
                opt.placer_config.moreau.initial_step_fraction = std::stod(argv[++i]);
            else if (arg == "--moreau-minimum-step-fraction" && i + 1 < argc)
                opt.placer_config.moreau.minimum_step_fraction = std::stod(argv[++i]);
            else if (arg == "--moreau-backtrack-factor" && i + 1 < argc)
                opt.placer_config.moreau.backtrack_factor = std::stod(argv[++i]);
            else if (arg == "--moreau-line-search-trials" && i + 1 < argc)
                opt.placer_config.moreau.max_line_search_trials = std::stoi(argv[++i]);
            else if (arg == "--moreau-stall-iterations" && i + 1 < argc)
                opt.placer_config.moreau.max_stall_iterations = std::stoi(argv[++i]);
            else if (arg == "--moreau-residual-tolerance" && i + 1 < argc)
                opt.placer_config.moreau.residual_rms_tolerance = std::stod(argv[++i]);
            else if (arg == "--moreau-displacement-tolerance" && i + 1 < argc)
                opt.placer_config.moreau.displacement_rms_tolerance = std::stod(argv[++i]);
            else if (arg == "--help" || arg == "-h")
            {
                printUsage(argv[0]);
                std::exit(EXIT_SUCCESS);
            }
            else
                throw std::invalid_argument("unknown or incomplete argument: " + arg);
        }
        if (opt.aux_path.empty())
            throw std::invalid_argument("missing required --aux <path>");
        opt.placer_config.bins_x = opt.bins_x;
        opt.placer_config.bins_y = opt.bins_y;
        opt.placer_config.target_density = opt.target_density;
        opt.quadratic_initial_config.simple_reference = opt.simple_initial_config;
        return opt;
    }


    std::unique_ptr<InitialPlacer> createInitialPlacer(const Options& options)
    {
        if (options.initial_method == InitialPlacementMethod::QuadraticConnectivity)
        {
            return std::make_unique<QuadraticInitialPlacer>(options.quadratic_initial_config);
        }
        return std::make_unique<SimplePackingInitialPlacer>(options.simple_initial_config);
    }

    void logInitialPlacementResult(const InitialPlacementResult& r)
    {
        LOG_INFO("initial placement method: " << r.method);
        LOG_INFO("initial movable cell count: " << r.movable_cell_count);
        LOG_INFO("initial active net count: " << r.active_net_count);
        LOG_INFO("initial skipped high-degree net count: " << r.skipped_high_degree_net_count);
        LOG_INFO("initial x solver iterations: " << r.x_solver_iterations);
        LOG_INFO("initial y solver iterations: " << r.y_solver_iterations);
        LOG_INFO("initial x relative residual: " << r.x_relative_residual);
        LOG_INFO("initial y relative residual: " << r.y_relative_residual);
        LOG_INFO("initial x converged: " << (r.x_converged ? "true" : "false"));
        LOG_INFO("initial y converged: " << (r.y_converged ? "true" : "false"));
        LOG_INFO("initial used fallback: " << (r.used_fallback ? "true" : "false"));
        LOG_INFO("initial message: " << r.message);
    }

    std::string binGridSummary(const BinGrid &grid)
    {
        const Box &c = grid.coreBounds();
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
    std::string densitySummary(const DensityMetrics &metrics)
    {
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

    void writeDensitySummary(std::ostream &os, const DensityMetrics &metrics)
    {
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

int main(int argc, char **argv)
{
    bool logger_ready = false;
    try
    {
        const Options opt = parseArgs(argc, argv);
        Logger::instance().initialize(opt.log_path, opt.log_level, true);
        logger_ready = true;
        LOG_INFO("program started");
        LOG_INFO("resolved command-line arguments: aux=" << opt.aux_path << ", bins=" << opt.bins_x << "x" << opt.bins_y << ", target_density=" << opt.target_density << ", zero_capacity_repulsion=" << opt.placer_config.zero_capacity_repulsion << ", optimizer=" << (opt.placer_config.method == GlobalPlacementMethod::MoreauProximal ? "moreau" : "direct") << ", moreau_mu=" << opt.placer_config.moreau.mu);
        LOG_INFO("log file path: " << opt.log_path);
        LOG_INFO("benchmark path: " << opt.aux_path);
        if (opt.bins_x <= 0 || opt.bins_y <= 0 || opt.target_density <= 0.0 || opt.target_density > 1.0)
        {
            printUsage(argv[0]);
            throw std::invalid_argument("illegal BinGrid arguments: bins=" +
                                        std::to_string(opt.bins_x) + "x" + std::to_string(opt.bins_y) +
                                        ", target_density=" + std::to_string(opt.target_density));
        }

        // 读取文件
        PlacementDB db;
        LimboBookshelfAdapter adapter(db);
        LOG_INFO("Bookshelf parsing started");
        if (!adapter.read(opt.aux_path))
            LOG_FATAL("Failed to read Bookshelf design: " << opt.aux_path);
        LOG_INFO("Bookshelf parsing completed");
        LOG_INFO("cell/net/pin/row counts: cells=" << db.cells().size() << ", nets=" << db.nets().size() << ", pins=" << db.pins().size() << ", rows=" << db.rows().size());

        // 初始化布局
        LOG_INFO("Initial placement started");
        std::unique_ptr<InitialPlacer> initial_placer = createInitialPlacer(opt);
        const InitialPlacementResult initial_result = initial_placer->place(db);
        logInitialPlacementResult(initial_result);
        LOG_INFO("Initial placement completed");

        // 线长计算
        LOG_INFO("HPWL evaluation started");
        HPWLEvaluator hpwl_evaluator;
        double total_hpwl = 0.0;
        for (std::size_t i = 0; i < db.nets().size(); ++i)
        {
            const int net_id = static_cast<int>(i);
            const double net_hpwl = hpwl_evaluator.netHPWL(db, net_id);
            db.addNetHPWL(net_id, net_hpwl);
            total_hpwl += net_hpwl;
        }
        LOG_INFO("total HPWL: " << std::fixed << std::setprecision(3) << total_hpwl);

        // 开始搭建bin，并计算容量
        LOG_INFO("BinGrid construction started");
        BinGrid grid(db, opt.bins_x, opt.bins_y, opt.target_density);
        const Box &core = grid.coreBounds();
        LOG_INFO("core bounds: (" << core.lx << ", " << core.ly << ") - (" << core.ux << ", " << core.uy << ")");
        LOG_INFO("bin count: " << grid.numBins() << ", bin size=" << grid.binWidth() << "x" << grid.binHeight() << ", target density=" << grid.targetDensity());
        LOG_INFO("BinGrid totals: movable_area=" << grid.totalMovableArea() << ", fixed_area=" << grid.totalFixedArea() << ", total_overflow=" << grid.totalOverflow() << ", overflow_bins=" << grid.overflowBinCount() << ", max_utilization=" << grid.maxUtilization());

        // 计算每一个bin的密度
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

        // 开始进行优化
        ObjectiveEvaluator objective_evaluator;
        ObjectiveMetrics initial_objective = objective_evaluator.evaluate(db, grid, opt.placer_config.density_weight);

        // 开始全局布局
        GlobalPlacer placer(opt.placer_config);
        GlobalPlacerResult placement_result = placer.optimize(db);

        // 记录最终结果
        BinGrid final_grid(db, opt.bins_x, opt.bins_y, opt.target_density);
        const DensityMetrics final_density_metrics = density_evaluator.evaluate(final_grid);
        const ObjectiveMetrics final_objective = objective_evaluator.evaluate(db, final_grid, opt.placer_config.density_weight);
        const double improvement = (initial_objective.total_cost - final_objective.total_cost) / std::max(std::abs(initial_objective.total_cost), 1e-12);

        LOG_INFO("Final HPWL evaluation started");
        HPWLEvaluator final_hpwl;
        double final_total_hpwl = 0.0;
        final_total_hpwl = final_hpwl.updateAllNetHPWL(db);
        LOG_INFO("final total HPWL: " << std::fixed << std::setprecision(3) << final_total_hpwl);

        // 打印summary
        std::filesystem::create_directories("../result");
        const std::string placement_summary = "../result/placementdb_summary.txt";
        std::ofstream fout(placement_summary);
        if (!fout)
            LOG_FATAL("cannot open summary output: " << placement_summary);
        db.printSummary(fout);
        fout << "\n[HPWL Results]\n"
             << std::fixed << std::setprecision(3) << "Total HPWL: " << total_hpwl << "\n"
             << std::fixed << std::setprecision(3) << "Final Total HPWL: " << final_total_hpwl << "\n";

        writeDensitySummary(fout, density_metrics);
        fout << "\n[Global Placement Result]\n"
             << "Initial total cost: " << initial_objective.total_cost << "\n"
             << "Final total cost: " << final_objective.total_cost << "\n"
             << "Optimizer method: " << (opt.placer_config.method == GlobalPlacementMethod::MoreauProximal ? "moreau" : "direct") << "\n"
             << "Density weight: " << opt.placer_config.density_weight << "\n"
             << "Zero-capacity repulsion: " << opt.placer_config.zero_capacity_repulsion << "\n"
             << "Moreau mu: " << opt.placer_config.moreau.mu << "\n"
             << "Accepted iterations: " << placement_result.accepted_iterations << "\n"
             << "Termination reason: " << placement_result.termination_reason << "\n";
        fout.close();



        const std::string grid_summary_path = "../result/bin_grid_summary.txt";
        const std::string final_grid_summary_path = "../result/final_bin_grid_summary.txt";
        std::ofstream gout(grid_summary_path);
        std::ofstream fgout(final_grid_summary_path);

        if (!gout)
            LOG_FATAL("cannot open BinGrid summary output: " << grid_summary_path);
        gout << binGridSummary(grid) << '\n'
             << densitySummary(density_metrics) << '\n';
        LOG_INFO("summary output path: " << grid_summary_path);

        if (!fgout)
            LOG_FATAL("cannot open BinGrid summary output: " << final_grid_summary_path);
        fgout << binGridSummary(final_grid) << '\n'
              << densitySummary(final_density_metrics) << '\n';
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
                  << "Optimizer method: " << (opt.placer_config.method == GlobalPlacementMethod::MoreauProximal ? "moreau" : "direct") << "\n"
                  << "Density weight: " << opt.placer_config.density_weight << "\n"
                  << "Zero-capacity repulsion: " << opt.placer_config.zero_capacity_repulsion << "\n"
                  << "Moreau mu: " << opt.placer_config.moreau.mu << "\n"
                  << "Initial HPWL: " << initial_objective.hpwl << "\n"
                  << "Final HPWL: " << final_objective.hpwl << "\n"
                  << "HPWL change: " << (final_objective.hpwl - initial_objective.hpwl) << "\n"
                  << "Initial density penalty: " << initial_objective.density_penalty << "\n"
                  << "Final density penalty: " << final_objective.density_penalty << "\n"
                  << "Initial normalized HPWL: " << initial_objective.normalized_hpwl << "\n"
                  << "Final normalized HPWL: " << final_objective.normalized_hpwl << "\n"
                  << "Initial normalized density penalty: " << initial_objective.normalized_density_penalty << "\n"
                  << "Final normalized density penalty: " << final_objective.normalized_density_penalty << "\n"
                  << "Initial total cost: " << initial_objective.total_cost << "\n"
                  << "Final total cost: " << final_objective.total_cost << "\n"
                  << "Total cost improvement ratio: " << improvement << "\n"
                  << "Initial overflow ratio: " << initial_objective.overflow_ratio << "\n"
                  << "Final overflow ratio: " << final_objective.overflow_ratio << "\n"
                  << "Proximal term: " << (placement_result.history.empty() ? 0.0 : placement_result.history.back().proximal_term) << "\n"
                  << "Proximal objective: " << (placement_result.history.empty() ? final_objective.total_cost : placement_result.history.back().proximal_objective) << "\n"
                  << "Displacement RMS: " << (placement_result.history.empty() ? 0.0 : placement_result.history.back().proximal_displacement_rms) << "\n"
                  << "Proximal residual RMS: " << (placement_result.history.empty() ? 0.0 : placement_result.history.back().proximal_residual_rms) << "\n"
                  << "Attempted iterations: " << placement_result.attempted_iterations << "\n"
                  << "Accepted iterations: " << placement_result.accepted_iterations << "\n"
                  << "Termination reason: " << placement_result.termination_reason << "\n"
                  << "=============================================\n";
        LOG_INFO("program completed successfully");
        Logger::instance().flush();
        return EXIT_SUCCESS;
    }
    catch (const FatalLogError &)
    {
        if (logger_ready)
            Logger::instance().flush();
        return EXIT_FAILURE;
    }
    catch (const std::exception &e)
    {
        if (logger_ready)
        {
            LOG_ERROR("Unhandled exception: " << e.what());
            Logger::instance().flush();
        }
        else
        {
            std::cerr << "Error before logger initialization: " << e.what() << '\n';
            printUsage(argv[0]);
        }
        return EXIT_FAILURE;
    }
    catch (...)
    {
        if (logger_ready)
        {
            LOG_ERROR("Unhandled unknown exception");
            Logger::instance().flush();
        }
        else
            std::cerr << "Unknown error before logger initialization\n";
        return EXIT_FAILURE;
    }
}
