#include "density/ExactDensitySubgradient.h"

#include "db/Cell.h"
#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/PaperOverlapFunction.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
constexpr double kRelativeTolerance = 1e-9;

void requireFiniteNonnegative(double value, const char* name)
{
    if (!std::isfinite(value) || value < -kRelativeTolerance) {
        throw std::runtime_error(name);
    }
}

int clampIndex(int value, int lo, int hi)
{
    return std::max(lo, std::min(value, hi));
}
} // namespace

ExactDensitySubgradientResult ExactDensitySubgradient::evaluate(
    const PlacementDB& db,
    const ExactBinGrid& grid,
    std::mt19937_64& random_engine
) const
{
    const double total_cell_area = grid.rawTotalCellArea();
    if (!(total_cell_area > 0.0) || !std::isfinite(total_cell_area)) {
        throw std::runtime_error(
            "ExactDensitySubgradient: equation (18) total cell area must be positive");
    }

    ExactDensitySubgradientResult result;
    result.grad_x.assign(db.cells().size(), 0.0);
    result.grad_y.assign(db.cells().size(), 0.0);

    std::vector<double> overflow(static_cast<std::size_t>(grid.numBins()), 0.0);
    double overflow_sum = 0.0;
    for (const Bin& bin : grid.bins()) {
        const double capacity = grid.targetDensity() * bin.bounds.area();
        const double value = std::max(0.0, bin.total_area - capacity);
        overflow[static_cast<std::size_t>(bin.id)] = value;
        overflow_sum += value;
        result.penalty += value * value;
    }
    result.ofr = overflow_sum / total_cell_area;

    const Box& core = grid.coreBounds();
    const double bw = grid.binWidth();
    const double bh = grid.binHeight();

    for (const Cell& cell : db.cells()) {
        if (!cell.isMovable()) {
            continue;
        }

        int ix_begin = static_cast<int>(std::floor((cell.x - core.lx) / bw)) - 1;
        int iy_begin = static_cast<int>(std::floor((cell.y - core.ly) / bh)) - 1;
        int ix_end = static_cast<int>(std::ceil((cell.x + cell.width - core.lx) / bw));
        int iy_end = static_cast<int>(std::ceil((cell.y + cell.height - core.ly) / bh));

        ix_begin = clampIndex(ix_begin, 0, grid.numBinsX() - 1);
        iy_begin = clampIndex(iy_begin, 0, grid.numBinsY() - 1);
        ix_end = clampIndex(ix_end, 0, grid.numBinsX() - 1);
        iy_end = clampIndex(iy_end, 0, grid.numBinsY() - 1);

        const double cell_center_x = cell.x + 0.5 * cell.width;
        const double cell_center_y = cell.y + 0.5 * cell.height;

        for (int iy = iy_begin; iy <= iy_end; ++iy) {
            for (int ix = ix_begin; ix <= ix_end; ++ix) {
                const Bin& bin = grid.bin(iy * grid.numBinsX() + ix);
                const double bin_center_x = 0.5 * (bin.bounds.lx + bin.bounds.ux);
                const double bin_center_y = 0.5 * (bin.bounds.ly + bin.bounds.uy);

                AxisOverlapSubgradient ox = paperAxisOverlapSubgradient(
                    cell_center_x, cell.width, bin_center_x, bin.bounds.width(), random_engine);
                AxisOverlapSubgradient oy = paperAxisOverlapSubgradient(
                    cell_center_y, cell.height, bin_center_y, bin.bounds.height(), random_engine);

                const double bin_overflow = overflow[static_cast<std::size_t>(bin.id)];
                if (bin_overflow <= 0.0) {
                    continue;
                }

                result.grad_x[static_cast<std::size_t>(cell.id)] +=
                    2.0 * bin_overflow * ox.subgradient * oy.overlap;
                result.grad_y[static_cast<std::size_t>(cell.id)] +=
                    2.0 * bin_overflow * ox.overlap * oy.subgradient;
            }
        }
    }

    requireFiniteNonnegative(result.penalty, "ExactDensitySubgradient: invalid penalty");
    requireFiniteNonnegative(result.ofr, "ExactDensitySubgradient: invalid OFR");
    for (double value : result.grad_x) {
        if (!std::isfinite(value)) throw std::runtime_error("ExactDensitySubgradient: non-finite grad_x");
    }
    for (double value : result.grad_y) {
        if (!std::isfinite(value)) throw std::runtime_error("ExactDensitySubgradient: non-finite grad_y");
    }
    return result;
}
