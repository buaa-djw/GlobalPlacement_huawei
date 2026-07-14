#pragma once

#include "evaluator/ObjectiveEvaluator.h"
#include "optimizer/SubgradientEvaluator.h"

#include <string>
#include <vector>

class BinGrid;
class PlacementDB;

/** @brief Configuration for the Moreau proximal inner solver. */
struct MoreauProximalConfig {
    /** Moreau proximal parameter μ. Smaller μ produces stronger attraction to the outer anchor. */
    double mu = 1.0;
    /** Maximum number of inner proximal descent iterations. */
    int max_inner_iterations = 10;
    /** Initial line-search step as a fraction of the smaller bin dimension. */
    double initial_step_fraction = 0.02;
    /** Minimum admissible step fraction before termination. */
    double minimum_step_fraction = 1e-6;
    /** Multiplicative backtracking factor in (0, 1). */
    double backtrack_factor = 0.5;
    /** Maximum line-search trials per inner iteration. */
    int max_line_search_trials = 8;
    /** Maximum consecutive failed inner iterations. */
    int max_stall_iterations = 3;
    /** Required relative decrease in proximal objective. */
    double relative_improvement_tolerance = 1e-10;
    /** Residual RMS convergence tolerance. */
    double residual_rms_tolerance = 1e-12;
    /** Displacement RMS convergence tolerance. */
    double displacement_rms_tolerance = 1e-8;
    /** Logging interval for inner iterations. */
    int report_interval = 1;
};

/** @brief Objective and first-order diagnostics for a proximal state. */
struct MoreauProximalMetrics {
    ObjectiveMetrics objective;
    double proximal_term = 0.0;
    double proximal_objective = 0.0;
    double squared_displacement = 0.0;
    double displacement_rms = 0.0;
    double objective_direction_rms = 0.0;
    double proximal_direction_rms = 0.0;
    double combined_residual_rms = 0.0;
    double normalized_search_direction_rms = 0.0;
    double maximum_search_direction_norm = 0.0;
};

/** @brief Per-inner-iteration Moreau proximal diagnostics. */
struct MoreauProximalIteration {
    int inner_iteration = 0;
    bool accepted = false;
    int line_search_trials = 0;
    double step_fraction = 0.0;
    double accepted_step = 0.0;
    double objective_cost = 0.0;
    double proximal_term = 0.0;
    double proximal_objective = 0.0;
    double squared_displacement = 0.0;
    double displacement_rms = 0.0;
    double objective_direction_rms = 0.0;
    double proximal_direction_rms = 0.0;
    double combined_residual_rms = 0.0;
    double normalized_search_direction_rms = 0.0;
    double maximum_search_direction_norm = 0.0;
    double relative_proximal_improvement = 0.0;
};

/** @brief Result of one approximate Moreau proximal subproblem solve. */
struct MoreauProximalResult {
    ObjectiveMetrics anchor_metrics;
    ObjectiveMetrics final_metrics;
    ObjectiveMetrics best_metrics;
    MoreauProximalMetrics anchor_proximal_metrics;
    MoreauProximalMetrics final_proximal_metrics;
    MoreauProximalMetrics best_proximal_metrics;
    std::vector<MoreauProximalIteration> history;
    int attempted_inner_iterations = 0;
    int accepted_inner_iterations = 0;
    bool accepted = false;
    bool converged = false;
    std::string termination_reason;
};

/** @brief Approximate solver for Q_mu(z;x)=F(z)+||z-x||^2/(2 mu L_core^2). */
class MoreauProximalSolver {
public:
    /** @brief Validate and store proximal and subgradient configuration. */
    MoreauProximalSolver(MoreauProximalConfig config, SubgradientConfig subgradient_config);

    /**
     * @brief Approximately solve the Moreau proximal subproblem.
     *
     * Entry state: db and grid represent the outer anchor x_k.
     * Successful exit: db and grid represent the best approximate proximal point z_k.
     * Failed/no-progress exit: db and grid are restored to the anchor x_k.
     */
    MoreauProximalResult solve(PlacementDB& db, BinGrid& grid) const;

private:
    MoreauProximalConfig config_;
    SubgradientConfig subgradient_config_;
};
