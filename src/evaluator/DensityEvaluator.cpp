#include "evaluator/DensityEvaluator.h"

#include "utils/Logger.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
/**
 * @brief Shared tolerance for tiny floating-point residue around zero.
 *
 * BinGrid area arithmetic may create values such as -1e-15 after subtraction in
 * future implementations.  The evaluator treats only such tiny residue as zero;
 * clear negative values remain invalid and are surfaced to the caller.
 */
constexpr double kDensityEpsilon = 1e-12;

std::string valueToString(double value) {
    std::ostringstream os;
    os << value;
    return os.str();
}

[[noreturn]] void throwInvalidBin(int bin_id, const char* field, double value) {
    throw std::invalid_argument("DensityEvaluator: bin " + std::to_string(bin_id) +
                                " has invalid " + field + ": " + valueToString(value));
}

[[noreturn]] void throwInvalidGrid(const std::string& message) {
    throw std::invalid_argument("DensityEvaluator: " + message);
}

double checkedNonNegative(double value, int bin_id, const char* field) {
    if (!std::isfinite(value)) throwInvalidBin(bin_id, field, value);
    if (value < -kDensityEpsilon) throwInvalidBin(bin_id, field, value);
    return (value < 0.0) ? 0.0 : value;
}

void checkTotal(double value, const char* field) {
    if (!std::isfinite(value) || value < -kDensityEpsilon) {
        throwInvalidGrid(std::string("invalid ") + field + ": " + valueToString(value));
    }
}
}

DensityMetrics DensityEvaluator::evaluate(const BinGrid& grid, bool emit_warnings) const {
    const int nx = grid.numBinsX();
    const int ny = grid.numBinsY();
    const int num_bins = grid.numBins();
    if (nx <= 0 || ny <= 0 || num_bins <= 0) {
        throwInvalidGrid("grid must contain at least one bin");
    }
    if (num_bins != nx * ny || static_cast<std::size_t>(num_bins) != grid.bins().size()) {
        throwInvalidGrid("bin count does not match grid dimensions");
    }
    if (!std::isfinite(grid.targetDensity()) || grid.targetDensity() <= 0.0 || grid.targetDensity() > 1.0) {
        throwInvalidGrid("invalid target density: " + valueToString(grid.targetDensity()));
    }

    LOG_DEBUG("Density evaluation started for " << num_bins << " bins");

    DensityMetrics metrics;
    metrics.bin_count = static_cast<std::size_t>(num_bins);

    for (const Bin& bin : grid.bins()) {
        const double movable_area = checkedNonNegative(bin.movable_area, bin.id, "movable area");
        const double fixed_area = checkedNonNegative(bin.fixed_area, bin.id, "fixed area");
        const double movable_capacity = checkedNonNegative(bin.movable_capacity, bin.id, "movable capacity");
        const double overflow = checkedNonNegative(bin.overflow, bin.id, "overflow");

        if (std::isnan(bin.utilization) || bin.utilization < -kDensityEpsilon ||
            (std::isinf(bin.utilization) && bin.utilization < 0.0)) {
            throwInvalidBin(bin.id, "utilization", bin.utilization);
        }
        if (std::isinf(bin.utilization) && !(movable_capacity <= kDensityEpsilon && movable_area > kDensityEpsilon)) {
            throwInvalidBin(bin.id, "utilization", bin.utilization);
        }

        metrics.total_movable_area += movable_area;
        metrics.total_fixed_area += fixed_area;
        metrics.total_movable_capacity += movable_capacity;
        metrics.total_overflow += overflow;
        metrics.penalty += overflow * overflow;
        if (overflow > kDensityEpsilon) ++metrics.overflow_bin_count;

        if (movable_capacity <= kDensityEpsilon) {
            ++metrics.zero_capacity_bin_count;
            // Occupied zero-capacity bins mean no movable capacity remains, yet
            // movable area is present; they may legitimately create +infinity
            // utilization and must still contribute overflow/penalty.
            if (movable_area > kDensityEpsilon) ++metrics.zero_capacity_occupied_bin_count;
        } else {
            // Average utilization excludes zero-capacity bins so infinity or
            // undefined ratios cannot pollute the mean of usable bins.
            metrics.average_utilization += bin.utilization;
            ++metrics.valid_utilization_bin_count;
        }
        metrics.max_utilization = std::max(metrics.max_utilization, bin.utilization);
    }

    checkTotal(metrics.total_movable_area, "total movable area");
    checkTotal(metrics.total_fixed_area, "total fixed area");
    checkTotal(metrics.total_movable_capacity, "total movable capacity");
    checkTotal(metrics.total_overflow, "total overflow");
    checkTotal(metrics.penalty, "density penalty");

    const double grid_total_movable_area = grid.totalMovableArea();
    const double grid_total_fixed_area = grid.totalFixedArea();
    const double grid_total_overflow = grid.totalOverflow();
    checkTotal(grid_total_movable_area, "grid total movable area");
    checkTotal(grid_total_fixed_area, "grid total fixed area");
    checkTotal(grid_total_overflow, "grid total overflow");
    if (std::abs(grid_total_movable_area - metrics.total_movable_area) > kDensityEpsilon) {
        throwInvalidGrid("grid total movable area disagrees with bins");
    }
    if (std::abs(grid_total_fixed_area - metrics.total_fixed_area) > kDensityEpsilon) {
        throwInvalidGrid("grid total fixed area disagrees with bins");
    }
    if (std::abs(grid_total_overflow - metrics.total_overflow) > kDensityEpsilon) {
        throwInvalidGrid("grid total overflow disagrees with bins");
    }

    if (metrics.valid_utilization_bin_count > 0) {
        metrics.average_utilization /= static_cast<double>(metrics.valid_utilization_bin_count);
    }
    if (metrics.total_movable_area <= kDensityEpsilon) {
        metrics.overflow_ratio = 0.0;
        if (emit_warnings) LOG_WARN("Density evaluation found no movable cell area; overflow ratio is set to 0");
    } else {
        metrics.overflow_ratio = metrics.total_overflow / metrics.total_movable_area;
    }
    if (emit_warnings && metrics.zero_capacity_bin_count > 0) {
        LOG_WARN("Density evaluation found " << metrics.zero_capacity_bin_count << " zero-capacity bins");
    }
    if (emit_warnings && metrics.zero_capacity_occupied_bin_count > 0) {
        LOG_WARN("Density evaluation found " << metrics.zero_capacity_occupied_bin_count << " occupied zero-capacity bins");
    }
    if (emit_warnings && std::isinf(metrics.max_utilization)) {
        LOG_WARN("Density evaluation maximum utilization is infinity");
    }
    if (!std::isfinite(metrics.overflow_ratio) || metrics.overflow_ratio < -kDensityEpsilon) {
        throwInvalidGrid("invalid overflow ratio: " + valueToString(metrics.overflow_ratio));
    }
    if (!(std::isfinite(metrics.max_utilization) || metrics.max_utilization == std::numeric_limits<double>::infinity())) {
        throwInvalidGrid("invalid maximum utilization: " + valueToString(metrics.max_utilization));
    }
    if (!std::isfinite(metrics.average_utilization) || metrics.average_utilization < -kDensityEpsilon) {
        throwInvalidGrid("invalid average utilization: " + valueToString(metrics.average_utilization));
    }

    LOG_DEBUG("Density evaluation completed: penalty=" << metrics.penalty
              << ", total_overflow=" << metrics.total_overflow);
    return metrics;
}
