#include "density/ExactDensityEvaluator.h"

#include "density/ExactBinGrid.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
constexpr double kRelativeTolerance = 1e-9;

void requireFiniteNonnegative(double value, const std::string& name)
{
    if (!std::isfinite(value) || value < -kRelativeTolerance) {
        throw std::runtime_error("ExactDensityEvaluator: invalid " + name);
    }
}
} // namespace

ExactDensityResult ExactDensityEvaluator::evaluate(const ExactBinGrid& grid) const
{
    const double total_cell_area = grid.rawTotalCellArea();
    if (!(total_cell_area > 0.0) || !std::isfinite(total_cell_area)) {
        throw std::runtime_error(
            "ExactDensityEvaluator: equation (18) total cell area must be positive");
    }

    ExactDensityResult result;
    double overflow_sum = 0.0;

    for (const Bin& bin : grid.bins()) {
        const double bin_area = bin.bounds.area();
        if (!(bin_area > 0.0) || !std::isfinite(bin_area)) {
            throw std::runtime_error("ExactDensityEvaluator: bin area must be positive");
        }

        requireFiniteNonnegative(bin.total_area, "bin density area");
        const double capacity = grid.targetDensity() * bin_area;
        const double overflow = std::max(0.0, bin.total_area - capacity);
        overflow_sum += overflow;
        result.penalty += overflow * overflow;
    }

    result.ofr = overflow_sum / total_cell_area;
    requireFiniteNonnegative(result.penalty, "penalty");
    requireFiniteNonnegative(result.ofr, "OFR");
    return result;
}
