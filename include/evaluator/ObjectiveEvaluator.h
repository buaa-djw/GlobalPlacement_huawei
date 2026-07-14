#pragma once

#include <cstddef>

class PlacementDB;
class BinGrid;

/** @brief Objective values for HPWL + density_weight * density penalty. */
struct ObjectiveMetrics
{
    double hpwl = 0.0;
    double density_penalty = 0.0;

    // Dimensionless metrics.
    double normalized_hpwl = 0.0;
    double normalized_density_penalty = 0.0;

    double density_weight = 0.0;
    double weighted_density_penalty = 0.0;

    double total_cost = 0.0;
    double total_overflow = 0.0;
    double overflow_ratio = 0.0;
    
    std::size_t overflow_bin_count = 0;
    std::size_t zero_capacity_bin_count = 0;
    std::size_t zero_capacity_occupied_bin_count = 0;
};

/**
 * @brief Read-only evaluator for the current nonsmooth placement objective.
 *
 * The caller must rebuild BinGrid after changing cell coordinates. This class
 * does not move cells, rebuild grids, mutate nets, log, write files, or exit.
 */
class ObjectiveEvaluator
{
public:
    /** @brief Evaluate HPWL + density_weight * sum overflow^2. */
    ObjectiveMetrics evaluate(const PlacementDB &db, const BinGrid &grid, double density_weight) const;
};
