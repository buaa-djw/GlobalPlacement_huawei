#include "density/ExactDensityEvaluator.h"

#include "density/ExactBinGrid.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{

constexpr double kRelativeTolerance = 1e-9;

bool near(double lhs, double rhs)
{
    return std::abs(lhs - rhs) <=
           kRelativeTolerance *
               std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

void requireFiniteNonnegative(
    double value,
    const std::string& name
)
{
    if (!std::isfinite(value) ||
        value < -kRelativeTolerance) {
        throw std::runtime_error(
            "ExactDensityEvaluator: invalid " + name
        );
    }
}

} // namespace

ExactDensityMetrics ExactDensityEvaluator::evaluate(
    const ExactBinGrid& grid
) const
{
    ExactDensityMetrics metrics;

    metrics.bin_count = grid.bins().size();

    metrics.total_cell_area =
        grid.totalCellArea();

    metrics.total_movable_area =
        grid.totalMovableArea();

    metrics.total_non_movable_area =
        grid.totalNonMovableArea();

    double area_sum = 0.0;
    double overflow_sum = 0.0;
    double penalty_sum = 0.0;

    for (const Bin& bin : grid.bins()) {
        const double bin_area =
            bin.bounds.area();

        if (!(bin_area > 0.0) ||
            !std::isfinite(bin_area)) {
            throw std::runtime_error(
                "ExactDensityEvaluator: "
                "bin area must be positive"
            );
        }

        requireFiniteNonnegative(
            bin.total_area,
            "bin total area"
        );

        requireFiniteNonnegative(
            bin.movable_area,
            "bin movable area"
        );

        requireFiniteNonnegative(
            bin.fixed_area,
            "bin non-movable area"
        );

        if (!near(
                bin.total_area,
                bin.movable_area + bin.fixed_area)) {
            throw std::runtime_error(
                "ExactDensityEvaluator: "
                "bin area split mismatch"
            );
        }

        const double overflow =
            std::max(
                0.0,
                bin.total_area -
                    bin.target_capacity
            );

        if (!near(overflow, bin.overflow)) {
            throw std::runtime_error(
                "ExactDensityEvaluator: "
                "bin overflow mismatch"
            );
        }

        metrics.total_capacity +=
            bin.target_capacity;

        overflow_sum += overflow;
        penalty_sum += overflow * overflow;

        if (overflow > 0.0) {
            ++metrics.overflow_bin_count;
        }

        metrics.maximum_bin_density =
            std::max(
                metrics.maximum_bin_density,
                bin.total_area / bin_area
            );

        area_sum += bin.total_area;
    }

    metrics.total_overflow =
        overflow_sum;

    metrics.quadratic_penalty =
        penalty_sum;

    if (metrics.total_movable_area > 0.0) {
        metrics.overflow_ratio =
            metrics.total_overflow /
            metrics.total_movable_area;

        metrics.overflow_ratio_defined = true;
    } else {
        metrics.overflow_ratio = 0.0;
        metrics.overflow_ratio_defined = false;
    }

    const double core_area =
        grid.coreBounds().area();

    if (!(core_area > 0.0) ||
        !std::isfinite(core_area)) {
        throw std::runtime_error(
            "ExactDensityEvaluator: "
            "invalid placement core area"
        );
    }

    metrics.average_bin_density =
        area_sum / core_area;

    requireFiniteNonnegative(
        metrics.total_overflow,
        "total overflow"
    );

    requireFiniteNonnegative(
        metrics.quadratic_penalty,
        "quadratic penalty"
    );

    return metrics;
}