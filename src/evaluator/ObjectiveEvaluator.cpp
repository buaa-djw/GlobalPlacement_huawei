#include "evaluator/ObjectiveEvaluator.h"

#include "db/PlacementDB.h"
#include "evaluator/DensityEvaluator.h"
#include "evaluator/HPWLEvaluator.h"
#include "grid/BinGrid.h"
#include "db/Net.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
    constexpr double kEps = 1e-12;
    void requireFiniteNonNegative(double value, const char *name)
    {
        if (!std::isfinite(value) || value < -kEps)
        {
            throw std::runtime_error(std::string("ObjectiveEvaluator: invalid ") + name);
        }
    }
}

namespace
{

    constexpr double kObjectiveEpsilon = 1e-30;

    std::size_t countActiveNets(
        const PlacementDB &db)
    {
        std::size_t count = 0;

        for (const Net &net : db.nets())
        {
            if (net.pin_ids.size() > 1)
            {
                ++count;
            }
        }

        return count;
    }

    double computeNormalizedHpwl(
        const PlacementDB &db,
        const BinGrid &grid,
        double hpwl)
    {
        const Box &core = grid.coreBounds();

        const double reference_length =
            std::hypot(
                core.width(),
                core.height());

        if (!std::isfinite(reference_length) ||
            reference_length <= 0.0)
        {
            throw std::runtime_error(
                "ObjectiveEvaluator: invalid core reference length");
        }

        const std::size_t active_net_count =
            countActiveNets(db);

        const double denominator = static_cast<double>(std::max<std::size_t>(active_net_count, 1)) * reference_length;

        return hpwl / denominator;
    }

    double computeNormalizedDensityPenalty(const BinGrid &grid)
    {
        if (grid.numBins() <= 0)
        {
            throw std::runtime_error("ObjectiveEvaluator: BinGrid contains no bins");
        }

        double normalized_sum = 0.0;

        for (const Bin &bin : grid.bins())
        {
            const double bin_area = bin.bounds.area();

            if (!std::isfinite(bin_area) || bin_area <= 0.0)
            {
                throw std::runtime_error("ObjectiveEvaluator: invalid bin area");
            }

            if (!std::isfinite(bin.overflow) || bin.overflow < 0.0)
            {
                throw std::runtime_error("ObjectiveEvaluator: invalid bin overflow");
            }

            const double normalized_overflow = bin.overflow / std::max(bin_area, kObjectiveEpsilon);

            normalized_sum += normalized_overflow * normalized_overflow;
        }

        return normalized_sum / static_cast<double>(grid.numBins());
    }

} // namespace

ObjectiveMetrics ObjectiveEvaluator::evaluate(
    const PlacementDB &db,
    const BinGrid &grid,
    double density_weight, bool emit_density_warnings) const
{
    if (!std::isfinite(density_weight) || density_weight < 0.0)
    {
        throw std::invalid_argument("ObjectiveEvaluator: density weight must be finite and non-negative");
    }

    HPWLEvaluator hpwl_evaluator;
    DensityEvaluator density_evaluator;

    ObjectiveMetrics metrics;

    metrics.hpwl = hpwl_evaluator.totalHPWL(db);

    const DensityMetrics density = density_evaluator.evaluate(grid, emit_density_warnings);

    metrics.density_penalty = density.penalty;

    metrics.normalized_hpwl = computeNormalizedHpwl(db, grid, metrics.hpwl);

    metrics.normalized_density_penalty = computeNormalizedDensityPenalty(grid);

    metrics.density_weight = density_weight;

    metrics.weighted_density_penalty = density_weight * metrics.normalized_density_penalty;

    metrics.total_cost = metrics.normalized_hpwl + metrics.weighted_density_penalty;

    metrics.total_overflow = density.total_overflow;

    metrics.overflow_ratio = density.overflow_ratio;

    metrics.overflow_bin_count = density.overflow_bin_count;

    metrics.zero_capacity_bin_count = density.zero_capacity_bin_count;

    metrics.zero_capacity_occupied_bin_count = density.zero_capacity_occupied_bin_count;

    if (!std::isfinite(metrics.normalized_hpwl) || metrics.normalized_hpwl < 0.0)
    {
        throw std::runtime_error("ObjectiveEvaluator: invalid normalized HPWL");
    }

    if (!std::isfinite(metrics.normalized_density_penalty) || metrics.normalized_density_penalty < 0.0)
    {
        throw std::runtime_error("ObjectiveEvaluator: invalid normalized density penalty");
    }

    if (!std::isfinite(metrics.total_cost) || metrics.total_cost < 0.0)
    {
        throw std::runtime_error("ObjectiveEvaluator: invalid total cost");
    }

    return metrics;
}