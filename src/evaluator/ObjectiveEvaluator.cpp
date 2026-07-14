#include "evaluator/ObjectiveEvaluator.h"

#include "db/PlacementDB.h"
#include "evaluator/DensityEvaluator.h"
#include "evaluator/HPWLEvaluator.h"
#include "grid/BinGrid.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace {
constexpr double kEps = 1e-12;
void requireFiniteNonNegative(double value, const char* name) {
    if (!std::isfinite(value) || value < -kEps) {
        throw std::runtime_error(std::string("ObjectiveEvaluator: invalid ") + name);
    }
}
}

ObjectiveMetrics ObjectiveEvaluator::evaluate(const PlacementDB& db, const BinGrid& grid, double density_weight) const {
    if (!std::isfinite(density_weight) || density_weight < 0.0) {
        throw std::invalid_argument("ObjectiveEvaluator: density_weight must be finite and non-negative");
    }

    ObjectiveMetrics metrics;
    metrics.hpwl = HPWLEvaluator{}.totalHPWL(db);
    const DensityMetrics density = DensityEvaluator{}.evaluate(grid);

    metrics.density_penalty = density.penalty;
    metrics.density_weight = density_weight;
    metrics.weighted_density_penalty = density_weight * density.penalty;
    metrics.total_cost = metrics.hpwl + metrics.weighted_density_penalty;
    metrics.total_overflow = density.total_overflow;
    metrics.overflow_ratio = density.overflow_ratio;
    metrics.overflow_bin_count = density.overflow_bin_count;
    metrics.zero_capacity_bin_count = density.zero_capacity_bin_count;
    metrics.zero_capacity_occupied_bin_count = density.zero_capacity_occupied_bin_count;

    requireFiniteNonNegative(metrics.hpwl, "HPWL");
    requireFiniteNonNegative(metrics.density_penalty, "density penalty");
    requireFiniteNonNegative(metrics.weighted_density_penalty, "weighted density penalty");
    requireFiniteNonNegative(metrics.total_cost, "total cost");
    requireFiniteNonNegative(metrics.total_overflow, "total overflow");
    requireFiniteNonNegative(metrics.overflow_ratio, "overflow ratio");
    return metrics;
}
