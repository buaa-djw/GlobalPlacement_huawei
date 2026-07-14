#pragma once

#include <cstddef>
#include <vector>

class BinGrid;
class PlacementDB;

/** @brief Configuration for normalized HPWL and density surrogate directions. */
struct SubgradientConfig {
    double density_weight = 1.0;
    double zero_capacity_repulsion = 2.0;
    double numerical_epsilon = 1e-12;
    double extreme_tolerance = 1e-9;
};

/** @brief Diagnostics produced while building a placement search direction. */
struct SubgradientMetrics {
    std::size_t movable_cell_count = 0;
    std::size_t active_net_count = 0;
    std::size_t active_overflow_bin_count = 0;
    std::size_t zero_capacity_overflow_bin_count = 0;
    double core_reference_length = 0.0;
    double hpwl_subgradient_rms = 0.0;
    double hpwl_direction_rms = 0.0;
    double density_direction_rms = 0.0;
    double combined_direction_rms_before_normalization = 0.0;
    double combined_direction_rms = 0.0;
    double maximum_direction_norm = 0.0;
};

/** @brief Full vector result of normalized subgradient direction evaluation. */
struct SubgradientResult {
    /// Normalized HPWL subgradient of ObjectiveMetrics::normalized_hpwl.
    std::vector<double> hpwl_gradient_x;
    std::vector<double> hpwl_gradient_y;
    /// HPWL descent direction: negative normalized HPWL subgradient.
    std::vector<double> hpwl_direction_x;
    std::vector<double> hpwl_direction_y;
    /// Normalized-density-scale surrogate descent direction.
    std::vector<double> density_direction_x;
    std::vector<double> density_direction_y;
    /// Final RMS-normalized direction consumed by GlobalPlacer.
    std::vector<double> combined_direction_x;
    std::vector<double> combined_direction_y;
    SubgradientMetrics metrics;
};

/**
 * @brief Read-only evaluator for first-order global-placement search data.
 *
 * SubgradientEvaluator does not mutate PlacementDB, rebuild BinGrid, call the
 * objective evaluator, log, or terminate the process. It computes the
 * normalized HPWL subgradient, a normalized-density-scale surrogate descent
 * direction, and their density-weighted final direction.
 */
class SubgradientEvaluator {
public:
    /** @brief Validate and store the evaluator configuration. */
    explicit SubgradientEvaluator(SubgradientConfig config);

    /** @brief Compute normalized HPWL, density surrogate, and combined directions. */
    SubgradientResult evaluate(const PlacementDB& db, const BinGrid& grid) const;

private:
    SubgradientConfig config_;
};
