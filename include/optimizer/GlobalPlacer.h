#pragma once

#include "evaluator/ObjectiveEvaluator.h"
#include "optimizer/MoreauProximalSolver.h"

#include <string>
#include <vector>

class BinGrid;
class PlacementDB;

/** @brief Available global placement optimization methods. */
enum class GlobalPlacementMethod {
    DirectSubgradient,
    MoreauProximal
};

/** @brief User-tunable controls for the simple nonsmooth global placer. */
struct GlobalPlacerConfig {
    int bins_x = 64;
    int bins_y = 64;
    double target_density = 0.9;
    double density_weight = 1.0;
    double zero_capacity_repulsion = 2.0;
    GlobalPlacementMethod method = GlobalPlacementMethod::DirectSubgradient;
    MoreauProximalConfig moreau;
    int max_iterations = 20;
    double initial_step_fraction = 0.10;
    double minimum_step_fraction = 1e-5;
    double backtrack_factor = 0.5;
    int max_line_search_trials = 8;
    int max_stall_iterations = 5;
    double relative_improvement_tolerance = 1e-10;
    int report_interval = 1;
};

/** @brief Per-iteration optimization diagnostics. */
struct GlobalPlacerIteration {
    int iteration = 0;
    bool accepted = false;
    int line_search_trials = 0;
    double hpwl = 0.0;
    double density_penalty = 0.0;
    double normalized_hpwl = 0.0;
    double normalized_density_penalty = 0.0;
    double weighted_density_penalty = 0.0;
    double total_cost = 0.0;
    double total_overflow = 0.0;
    double overflow_ratio = 0.0;
    double step_fraction = 0.0;
    double accepted_step = 0.0;
    double hpwl_gradient_rms = 0.0;
    double density_direction_rms = 0.0;
    double combined_direction_rms = 0.0;
    std::size_t active_net_count = 0;
    std::size_t active_overflow_bin_count = 0;
    std::size_t zero_capacity_overflow_bin_count = 0;
    double maximum_direction_norm = 0.0;
    double relative_improvement = 0.0;
    bool used_moreau = false;
    double moreau_mu = 0.0;
    int proximal_inner_iterations = 0;
    int proximal_accepted_inner_iterations = 0;
    double proximal_term = 0.0;
    double proximal_objective = 0.0;
    double proximal_displacement_rms = 0.0;
    double proximal_residual_rms = 0.0;
    std::string proximal_termination_reason;
};

/** @brief Final result and history from GlobalPlacer::optimize(). */
struct GlobalPlacerResult {
    ObjectiveMetrics initial_metrics;
    ObjectiveMetrics final_metrics;
    ObjectiveMetrics best_metrics;
    std::vector<GlobalPlacerIteration> history;
    int attempted_iterations = 0;
    int accepted_iterations = 0;
    bool converged = false;
    std::string termination_reason;
};

/** @brief Coordinate-only snapshot indexed by Cell::id. */
struct CellPosition { double x = 0.0; double y = 0.0; };

/**
 * @brief Backtracking nonsmooth global placer for HPWL plus density overflow.
 *
 * The placer supports DirectSubgradient and MoreauProximal optimization modes.
 * It intentionally excludes overlap, legalizing, PR+, and file output.
 * Fixed/terminal cells are never updated and are checked after optimize.
 */
class GlobalPlacer {
public:
    /** @brief Validate and store the configuration. */
    explicit GlobalPlacer(GlobalPlacerConfig config);

    /** @brief Optimize movable cell coordinates in-place and restore the best state. */
    GlobalPlacerResult optimize(PlacementDB& db) const;

private:
    GlobalPlacerConfig config_;
};
