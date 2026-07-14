#include "optimizer/SubgradientEvaluator.h"

#include "db/Cell.h"
#include "db/Net.h"
#include "db/Pin.h"
#include "db/PlacementDB.h"
#include "geometry/Box.h"
#include "geometry/Point.h"
#include "grid/Bin.h"
#include "grid/BinGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double kSqrtHalf = 0.7071067811865475244;

void validateConfig(const SubgradientConfig& c) {
    if (!std::isfinite(c.density_weight) || c.density_weight < 0.0) throw std::invalid_argument("SubgradientEvaluator: invalid density_weight");
    if (!std::isfinite(c.zero_capacity_repulsion) || c.zero_capacity_repulsion < 1.0) throw std::invalid_argument("SubgradientEvaluator: invalid zero_capacity_repulsion");
    if (!std::isfinite(c.numerical_epsilon) || c.numerical_epsilon <= 0.0) throw std::invalid_argument("SubgradientEvaluator: invalid numerical_epsilon");
    if (!std::isfinite(c.extreme_tolerance) || c.extreme_tolerance <= 0.0) throw std::invalid_argument("SubgradientEvaluator: invalid extreme_tolerance");
}

bool isExtreme(double value, double extreme, double tolerance) {
    return std::abs(value - extreme) <= tolerance * std::max(1.0, std::abs(extreme));
}

void addUnique(std::vector<int>& ids, int id) {
    if (std::find(ids.begin(), ids.end(), id) == ids.end()) ids.push_back(id);
}

std::pair<double, double> deterministicTieDirection(int cell_id, int bin_id) {
    static const std::array<std::pair<double, double>, 8> dirs{{
        {1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}, {0.0, -1.0},
        {kSqrtHalf, kSqrtHalf}, {kSqrtHalf, -kSqrtHalf}, {-kSqrtHalf, kSqrtHalf}, {-kSqrtHalf, -kSqrtHalf}}};
    std::uint32_t h = static_cast<std::uint32_t>(cell_id) * 2654435761u;
    h ^= static_cast<std::uint32_t>(bin_id) * 2246822519u;
    h ^= h >> 16;
    return dirs[h & 7u];
}

double movableRms(const PlacementDB& db, const std::vector<double>& x, const std::vector<double>& y) {
    double sum = 0.0;
    std::size_t n = 0;
    for (const Cell& c : db.cells()) if (c.isMovable()) { sum += x.at(c.id) * x.at(c.id) + y.at(c.id) * y.at(c.id); ++n; }
    return n == 0 ? 0.0 : std::sqrt(sum / static_cast<double>(n));
}

void ensureFinite(double value, const std::string& what) {
    if (!std::isfinite(value)) throw std::runtime_error("SubgradientEvaluator: non-finite " + what);
}

void validateOutput(const PlacementDB& db, const SubgradientResult& r, double eps) {
    const std::size_t n = db.cells().size();
    const std::array<const std::vector<double>*, 10> vectors{{&r.hpwl_gradient_x, &r.hpwl_gradient_y, &r.hpwl_direction_x, &r.hpwl_direction_y, &r.density_direction_x, &r.density_direction_y, &r.objective_direction_x, &r.objective_direction_y, &r.combined_direction_x, &r.combined_direction_y}};
    for (const auto* v : vectors) if (v->size() != n) throw std::runtime_error("SubgradientEvaluator: vector size mismatch");
    for (const Cell& c : db.cells()) {
        for (const auto* v : vectors) {
            const double value = v->at(c.id);
            if (!std::isfinite(value)) throw std::runtime_error("SubgradientEvaluator: non-finite direction for cell " + c.name);
            if (!c.isMovable() && std::abs(value) > eps) throw std::runtime_error("SubgradientEvaluator: non-movable cell has non-zero direction: " + c.name);
        }
    }
    const SubgradientMetrics& m = r.metrics;
    ensureFinite(m.core_reference_length, "core_reference_length");
    ensureFinite(m.hpwl_subgradient_rms, "hpwl_subgradient_rms");
    ensureFinite(m.hpwl_direction_rms, "hpwl_direction_rms");
    ensureFinite(m.density_direction_rms, "density_direction_rms");
    ensureFinite(m.combined_direction_rms_before_normalization, "combined_direction_rms_before_normalization");
    ensureFinite(m.combined_direction_rms, "combined_direction_rms");
    ensureFinite(m.maximum_direction_norm, "maximum_direction_norm");
    if (m.combined_direction_rms > eps && std::abs(m.combined_direction_rms - 1.0) > 1e-8) throw std::runtime_error("SubgradientEvaluator: combined RMS is not normalized");
}
} // namespace

SubgradientEvaluator::SubgradientEvaluator(SubgradientConfig config) : config_(config) { validateConfig(config_); }

SubgradientResult SubgradientEvaluator::evaluate(const PlacementDB& db, const BinGrid& grid) const {
    const std::size_t n = db.cells().size();
    SubgradientResult r;
    r.hpwl_gradient_x.assign(n, 0.0); r.hpwl_gradient_y.assign(n, 0.0);
    r.hpwl_direction_x.assign(n, 0.0); r.hpwl_direction_y.assign(n, 0.0);
    r.density_direction_x.assign(n, 0.0); r.density_direction_y.assign(n, 0.0);
    r.objective_direction_x.assign(n, 0.0); r.objective_direction_y.assign(n, 0.0);
    r.combined_direction_x.assign(n, 0.0); r.combined_direction_y.assign(n, 0.0);

    for (const Cell& c : db.cells()) if (c.isMovable()) ++r.metrics.movable_cell_count;
    const Box& core = grid.coreBounds();
    r.metrics.core_reference_length = std::hypot(core.width(), core.height());
    if (!std::isfinite(r.metrics.core_reference_length) || r.metrics.core_reference_length <= 0.0) throw std::runtime_error("SubgradientEvaluator: invalid core reference length");

    for (const Net& net : db.nets()) {
        if (net.pin_ids.size() <= 1) continue;
        ++r.metrics.active_net_count;
        double minx = std::numeric_limits<double>::infinity(), maxx = -minx, miny = minx, maxy = -minx;
        for (int pid : net.pin_ids) {
            const Point p = db.pinPosition(pid);
            minx = std::min(minx, p.x); maxx = std::max(maxx, p.x); miny = std::min(miny, p.y); maxy = std::max(maxy, p.y);
        }
        std::vector<int> minxs, maxxs, minys, maxys;
        for (int pid : net.pin_ids) {
            const Point p = db.pinPosition(pid);
            const Cell& c = db.cell(db.pin(pid).cell_id);
            if (!c.isMovable()) continue;
            if (isExtreme(p.x, minx, config_.extreme_tolerance)) addUnique(minxs, c.id);
            if (isExtreme(p.x, maxx, config_.extreme_tolerance)) addUnique(maxxs, c.id);
            if (isExtreme(p.y, miny, config_.extreme_tolerance)) addUnique(minys, c.id);
            if (isExtreme(p.y, maxy, config_.extreme_tolerance)) addUnique(maxys, c.id);
        }
        for (int id : minxs) r.hpwl_gradient_x.at(id) -= 1.0 / static_cast<double>(minxs.size());
        for (int id : maxxs) r.hpwl_gradient_x.at(id) += 1.0 / static_cast<double>(maxxs.size());
        for (int id : minys) r.hpwl_gradient_y.at(id) -= 1.0 / static_cast<double>(minys.size());
        for (int id : maxys) r.hpwl_gradient_y.at(id) += 1.0 / static_cast<double>(maxys.size());
    }
    if (r.metrics.active_net_count > 0) {
        const double scale = 1.0 / (static_cast<double>(r.metrics.active_net_count) * r.metrics.core_reference_length);
        for (const Cell& c : db.cells()) if (c.isMovable()) {
            r.hpwl_gradient_x.at(c.id) *= scale; r.hpwl_gradient_y.at(c.id) *= scale;
            r.hpwl_direction_x.at(c.id) = -r.hpwl_gradient_x.at(c.id); r.hpwl_direction_y.at(c.id) = -r.hpwl_gradient_y.at(c.id);
        }
    }
    r.metrics.hpwl_subgradient_rms = movableRms(db, r.hpwl_gradient_x, r.hpwl_gradient_y);
    r.metrics.hpwl_direction_rms = movableRms(db, r.hpwl_direction_x, r.hpwl_direction_y);

    if (grid.numBins() <= 0) throw std::runtime_error("SubgradientEvaluator: BinGrid contains no bins");
    /**
     * The density direction below is a surrogate descent direction whose scale
     * matches the normalized density objective. It is not the exact analytic
     * gradient of cell-rectangle overlap. GlobalPlacer decides whether to
     * accept it through ObjectiveEvaluator plus backtracking line search.
     */
    for (const Bin& b : grid.bins()) {
        if (!std::isfinite(b.overflow) || b.overflow < 0.0) throw std::runtime_error("SubgradientEvaluator: invalid bin overflow");
        if (b.overflow <= config_.numerical_epsilon) continue;
        ++r.metrics.active_overflow_bin_count;
        if (b.movable_capacity <= config_.numerical_epsilon && b.movable_area > config_.numerical_epsilon) ++r.metrics.zero_capacity_overflow_bin_count;
        const double bin_area = b.bounds.area();
        if (!std::isfinite(bin_area) || bin_area <= 0.0) throw std::runtime_error("SubgradientEvaluator: invalid bin area");
        const double coeff = 2.0 * b.overflow / (static_cast<double>(grid.numBins()) * bin_area * bin_area);
        if (!std::isfinite(coeff)) throw std::runtime_error("SubgradientEvaluator: non-finite density coefficient");
        const double bc_x = 0.5 * (b.bounds.lx + b.bounds.ux), bc_y = 0.5 * (b.bounds.ly + b.bounds.uy);
        const double repulsion = (b.movable_capacity <= config_.numerical_epsilon && b.movable_area > config_.numerical_epsilon) ? config_.zero_capacity_repulsion : 1.0;
        for (int cid : b.cell_ids) {
            const Cell& c = db.cell(cid);
            if (!c.isMovable()) continue;
            const Box cb{c.x, c.y, c.x + c.width, c.y + c.height};
            const double overlap = rectangleOverlapArea(cb, b.bounds);
            if (overlap <= 0.0) continue;
            const double cell_area = c.width * c.height;
            const double fraction = overlap / std::max(cell_area, config_.numerical_epsilon);
            double vx = 0.5 * (cb.lx + cb.ux) - bc_x, vy = 0.5 * (cb.ly + cb.uy) - bc_y;
            double norm = std::hypot(vx, vy);
            if (norm <= config_.numerical_epsilon) { auto d = deterministicTieDirection(c.id, b.id); vx = d.first; vy = d.second; norm = 1.0; }
            const double magnitude = coeff * fraction * repulsion;
            r.density_direction_x.at(c.id) += magnitude * vx / norm;
            r.density_direction_y.at(c.id) += magnitude * vy / norm;
        }
    }
    r.metrics.density_direction_rms = movableRms(db, r.density_direction_x, r.density_direction_y);

    /**
     * HPWL and density components have already been scaled to the normalized
     * objective. Forcing each component to RMS=1 separately would erase their
     * relative normalized-objective scales and make density_weight lose its
     * ObjectiveEvaluator meaning. We therefore combine with the same
     * density_weight first, then apply one final movable RMS normalization so
     * step size remains geometrically stable.
     */
    for (const Cell& c : db.cells()) if (c.isMovable()) {
        r.objective_direction_x.at(c.id) = r.hpwl_direction_x.at(c.id) + config_.density_weight * r.density_direction_x.at(c.id);
        r.objective_direction_y.at(c.id) = r.hpwl_direction_y.at(c.id) + config_.density_weight * r.density_direction_y.at(c.id);
        r.combined_direction_x.at(c.id) = r.objective_direction_x.at(c.id);
        r.combined_direction_y.at(c.id) = r.objective_direction_y.at(c.id);
    }
    r.metrics.combined_direction_rms_before_normalization = movableRms(db, r.combined_direction_x, r.combined_direction_y);
    if (r.metrics.combined_direction_rms_before_normalization > config_.numerical_epsilon) {
        for (const Cell& c : db.cells()) if (c.isMovable()) {
            r.combined_direction_x.at(c.id) /= r.metrics.combined_direction_rms_before_normalization;
            r.combined_direction_y.at(c.id) /= r.metrics.combined_direction_rms_before_normalization;
            r.metrics.maximum_direction_norm = std::max(r.metrics.maximum_direction_norm, std::hypot(r.combined_direction_x.at(c.id), r.combined_direction_y.at(c.id)));
        }
    } else {
        std::fill(r.combined_direction_x.begin(), r.combined_direction_x.end(), 0.0);
        std::fill(r.combined_direction_y.begin(), r.combined_direction_y.end(), 0.0);
    }
    r.metrics.combined_direction_rms = movableRms(db, r.combined_direction_x, r.combined_direction_y);
    validateOutput(db, r, config_.numerical_epsilon);
    return r;
}
