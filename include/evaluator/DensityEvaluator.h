#pragma once

#include "grid/BinGrid.h"

#include <cstddef>

/**
 * @brief Aggregated density metrics for one placement state.
 *
 * The density penalty is the mathematical objective \f$\sum_b O_b^2\f$,
 * where \f$O_b\f$ is the overflow already stored by BinGrid for bin \f$b\f$.
 * Overflow ratio is normalized by the total movable cell area intersecting the
 * placement core, not by bin area, so it reports the fraction of movable area
 * that is currently over capacity.
 */
struct DensityMetrics {
    double penalty = 0.0;

    double total_overflow = 0.0;
    double overflow_ratio = 0.0;

    double total_movable_area = 0.0;
    double total_fixed_area = 0.0;
    double total_movable_capacity = 0.0;

    double max_utilization = 0.0;
    double average_utilization = 0.0;

    std::size_t bin_count = 0;
    std::size_t overflow_bin_count = 0;
    std::size_t valid_utilization_bin_count = 0;
    std::size_t zero_capacity_bin_count = 0;
    std::size_t zero_capacity_occupied_bin_count = 0;
};

/**
 * @brief Read-only evaluator for bin-grid density quality metrics.
 *
 * DensityEvaluator intentionally does not rebuild the grid, move cells, or
 * write reports.  The caller owns the PlacementDB -> BinGrid::rebuild() step,
 * then this evaluator performs a single O(B) pass over the existing bins and
 * returns DensityMetrics.  Invalid numeric state is reported with exceptions
 * instead of process termination so library users and tests can decide how to
 * log, flush, and exit.
 */
class DensityEvaluator {
public:
    /**
     * @brief Evaluates density-related metrics for the current BinGrid.
     *
     * This method does not rebuild or modify the grid. The caller must ensure
     * that BinGrid::rebuild() has already been called after any placement
     * change.
     *
     * @param grid Current read-only bin grid.
     * @return DensityMetrics Density penalty and related statistics.
     *
     * @throws std::invalid_argument if grid data contains invalid values.
     */
    DensityMetrics evaluate(const BinGrid& grid) const;
};
