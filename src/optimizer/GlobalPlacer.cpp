#include "optimizer/GlobalPlacer.h"

#include "db/Cell.h"
#include "db/Net.h"
#include "db/Pin.h"
#include "db/PlacementDB.h"
#include "geometry/Box.h"
#include "geometry/Point.h"
#include "grid/BinGrid.h"
#include "utils/Logger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace
{
    constexpr double kEps = 1e-12;
    constexpr double kSqrtHalf = 0.7071067811865475244;

    void validateConfig(const GlobalPlacerConfig &c)
    {
        if (c.bins_x <= 0 || c.bins_y <= 0)
            throw std::invalid_argument("GlobalPlacer: bins must be positive");
        if (!std::isfinite(c.target_density) || c.target_density <= 0.0 || c.target_density > 1.0)
            throw std::invalid_argument("GlobalPlacer: invalid target_density");
        if (!std::isfinite(c.density_weight) || c.density_weight < 0.0)
            throw std::invalid_argument("GlobalPlacer: invalid density_weight");
        if (c.max_iterations < 0)
            throw std::invalid_argument("GlobalPlacer: max_iterations must be non-negative");
        if (!std::isfinite(c.initial_step_fraction) || c.initial_step_fraction <= 0.0)
            throw std::invalid_argument("GlobalPlacer: invalid initial_step_fraction");
        if (!std::isfinite(c.minimum_step_fraction) || c.minimum_step_fraction <= 0.0 || c.minimum_step_fraction > c.initial_step_fraction)
            throw std::invalid_argument("GlobalPlacer: invalid minimum_step_fraction");
        if (!std::isfinite(c.backtrack_factor) || c.backtrack_factor <= 0.0 || c.backtrack_factor >= 1.0)
            throw std::invalid_argument("GlobalPlacer: invalid backtrack_factor");
        if (c.max_line_search_trials <= 0 || c.max_stall_iterations <= 0)
            throw std::invalid_argument("GlobalPlacer: invalid line search limits");
        if (!std::isfinite(c.relative_improvement_tolerance) || c.relative_improvement_tolerance < 0.0)
            throw std::invalid_argument("GlobalPlacer: invalid tolerance");
        if (c.report_interval <= 0)
            throw std::invalid_argument("GlobalPlacer: report_interval must be positive");
    }

    std::vector<CellPosition> savePositions(const PlacementDB &db)
    {
        std::vector<CellPosition> out(db.cells().size());
        for (const Cell &c : db.cells())
            out[static_cast<std::size_t>(c.id)] = {c.x, c.y};
        return out;
    }

    void restorePositions(PlacementDB &db, const std::vector<CellPosition> &p)
    {
        if (p.size() != db.cells().size())
            throw std::invalid_argument("GlobalPlacer: snapshot size mismatch");
        for (Cell &c : const_cast<std::vector<Cell> &>(db.cells()))
        {
            c.x = p[static_cast<std::size_t>(c.id)].x;
            c.y = p[static_cast<std::size_t>(c.id)].y;
        }
    }

    double movableRms(const PlacementDB &db, const std::vector<double> &x, const std::vector<double> &y)
    {
        double sum = 0.0;
        std::size_t n = 0;
        for (const Cell &c : db.cells())
            if (c.isMovable())
            {
                sum += x[c.id] * x[c.id] + y[c.id] * y[c.id];
                ++n;
            }
        return n ? std::sqrt(sum / static_cast<double>(n)) : 0.0;
    }

    bool isExtreme(double v, double e) { return std::abs(v - e) <= 1e-9 * std::max(1.0, std::abs(e)); }

    void addUnique(std::vector<int> &v, int id)
    {
        if (std::find(v.begin(), v.end(), id) == v.end())
            v.push_back(id);
    }

    std::pair<double, double> tieDirection(int cell_id, int bin_id)
    {
        static const std::array<std::pair<double, double>, 8> dirs{{{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {kSqrtHalf, kSqrtHalf}, {kSqrtHalf, -kSqrtHalf}, {-kSqrtHalf, kSqrtHalf}, {-kSqrtHalf, -kSqrtHalf}}};
        const unsigned h = static_cast<unsigned>(cell_id) * 2654435761u ^ static_cast<unsigned>(bin_id) * 2246822519u;
        return dirs[h & 7u];
    }

    void projectMovables(PlacementDB &db, const Box &core)
    {
        for (Cell &c : const_cast<std::vector<Cell> &>(db.cells()))
            if (c.isMovable())
            {
                if (c.width > core.width() + kEps || c.height > core.height() + kEps)
                    throw std::runtime_error("GlobalPlacer: movable cell larger than core: " + c.name);
                c.x = std::clamp(c.x, core.lx, core.ux - c.width);
                c.y = std::clamp(c.y, core.ly, core.uy - c.height);
                if (!std::isfinite(c.x) || !std::isfinite(c.y))
                    throw std::runtime_error("GlobalPlacer: non-finite movable coordinate");
            }
    }

    bool acceptCost(double trial, double cur, double tol)
    {
        if (!std::isfinite(trial))
            return false;
        return tol == 0.0 ? trial < cur : trial <= cur - tol * std::max(1.0, std::abs(cur));
    }
}

GlobalPlacer::GlobalPlacer(GlobalPlacerConfig config) : config_(config) { validateConfig(config_); }

void GlobalPlacer::computeHpwlSubgradient(const PlacementDB &db, std::vector<double> &gx, std::vector<double> &gy) const
{
    gx.assign(db.cells().size(), 0.0);
    gy.assign(db.cells().size(), 0.0);
    for (const Net &net : db.nets())
    {
        if (net.pin_ids.size() <= 1)
            continue;
        double minx = std::numeric_limits<double>::infinity(), maxx = -minx, miny = minx, maxy = -minx;
        for (int pid : net.pin_ids)
        {
            Point p = db.pinPosition(pid);
            minx = std::min(minx, p.x);
            maxx = std::max(maxx, p.x);
            miny = std::min(miny, p.y);
            maxy = std::max(maxy, p.y);
        }
        std::vector<int> minxs, maxxs, minys, maxys;
        for (int pid : net.pin_ids)
        {
            Point p = db.pinPosition(pid);
            const Cell &c = db.cell(db.pin(pid).cell_id);
            if (!c.isMovable())
                continue;
            if (isExtreme(p.x, minx))
                addUnique(minxs, c.id);
            if (isExtreme(p.x, maxx))
                addUnique(maxxs, c.id);
            if (isExtreme(p.y, miny))
                addUnique(minys, c.id);
            if (isExtreme(p.y, maxy))
                addUnique(maxys, c.id);
        }
        for (int id : minxs)
            gx[id] -= 1.0 / static_cast<double>(minxs.size());
        for (int id : maxxs)
            gx[id] += 1.0 / static_cast<double>(maxxs.size());
        for (int id : minys)
            gy[id] -= 1.0 / static_cast<double>(minys.size());
        for (int id : maxys)
            gy[id] += 1.0 / static_cast<double>(maxys.size());
    }
}

void GlobalPlacer::computeDensityDirection(const PlacementDB &db, const BinGrid &grid, std::vector<double> &dx, std::vector<double> &dy) const
{
    dx.assign(db.cells().size(), 0.0);
    dy.assign(db.cells().size(), 0.0);
    // Surrogate direction: pushes movable cells out of overflowing bins. This is
    // not the analytic gradient of sum(overflow^2); line search validates it.
    // It avoids utilization and 1/capacity because zero-capacity bins can occur.
    for (const Bin &b : grid.bins())
        if (b.overflow > 0.0)
        {
            const double bc_x = 0.5 * (b.bounds.lx + b.bounds.ux), bc_y = 0.5 * (b.bounds.ly + b.bounds.uy);
            for (int cid : b.cell_ids)
            {
                const Cell &c = db.cell(cid);
                if (!c.isMovable())
                    continue;
                const Box cb{c.x, c.y, c.x + c.width, c.y + c.height};
                const double ov = rectangleOverlapArea(cb, b.bounds);
                if (ov <= 0.0)
                    continue;
                double vx = 0.5 * (cb.lx + cb.ux) - bc_x, vy = 0.5 * (cb.ly + cb.uy) - bc_y;
                double norm = std::sqrt(vx * vx + vy * vy);
                if (norm <= kEps)
                {
                    auto d = tieDirection(cid, b.id);
                    vx = d.first;
                    vy = d.second;
                    norm = 1.0;
                }
                const double mag = (b.overflow / std::max(b.bounds.area(), kEps)) * (ov / std::max(c.width * c.height, kEps));
                dx[cid] += mag * vx / norm;
                dy[cid] += mag * vy / norm;
            }
        }
}

GlobalPlacerResult GlobalPlacer::optimize(PlacementDB &db) const
{
    const std::size_t movable_count = std::count_if(db.cells().begin(), db.cells().end(), [](const Cell &c)
                                                    { return c.isMovable(); });
    LOG_INFO("Global placement started: movable_cells=" << movable_count << " bins=" << config_.bins_x << "x" << config_.bins_y << " target_density=" << config_.target_density << " density_weight=" << config_.density_weight << " max_iterations=" << config_.max_iterations);
    BinGrid grid(db, config_.bins_x, config_.bins_y, config_.target_density);
    ObjectiveEvaluator eval;
    GlobalPlacerResult r;
    const auto fixed_initial = savePositions(db);
    r.initial_metrics = eval.evaluate(db, grid, config_.density_weight);
    r.best_metrics = r.initial_metrics;
    auto best_pos = savePositions(db);
    ObjectiveMetrics current = r.initial_metrics;
    if (config_.max_iterations == 0)
    {
        r.final_metrics = r.initial_metrics;
        r.best_metrics = r.initial_metrics;
        r.termination_reason = "zero iterations requested";
        restorePositions(db, best_pos);
        return r;
    }
    double step_fraction = config_.initial_step_fraction;
    int stalls = 0;
    for (int it = 1; it <= config_.max_iterations; ++it)
    {
        r.attempted_iterations = it;
        const auto start_pos = savePositions(db);
        const ObjectiveMetrics start = current;
        std::vector<double> gx, gy, wx, wy, dx, dy, cx(db.cells().size(), 0.0), cy(db.cells().size(), 0.0);
        computeHpwlSubgradient(db, gx, gy);
        computeDensityDirection(db, grid, dx, dy);
        wx = gx;
        wy = gy;
        for (std::size_t i = 0; i < wx.size(); ++i)
        {
            wx[i] = -wx[i];
            wy[i] = -wy[i];
        }
        const double wr = movableRms(db, wx, wy), dr = movableRms(db, dx, dy);
        for (const Cell &c : db.cells())
            if (c.isMovable())
            {
                if (wr > kEps)
                {
                    cx[c.id] += wx[c.id] / wr;
                    cy[c.id] += wy[c.id] / wr;
                }
                if (dr > kEps)
                {
                    cx[c.id] += config_.density_weight * dx[c.id] / dr;
                    cy[c.id] += config_.density_weight * dy[c.id] / dr;
                }
            }
        double cr = movableRms(db, cx, cy);
        if (cr <= kEps)
        {
            r.converged = true;
            r.termination_reason = "zero search direction";
            break;
        }
        for (const Cell &c : db.cells())
            if (c.isMovable())
            {
                cx[c.id] /= cr;
                cy[c.id] /= cr;
            }
        cr = movableRms(db, cx, cy);
        GlobalPlacerIteration rec;
        rec.iteration = it;
        rec.step_fraction = step_fraction;
        rec.hpwl_gradient_rms = wr;
        rec.density_direction_rms = dr;
        rec.combined_direction_rms = cr;
        double trial_frac = step_fraction;
        for (int tr = 1; tr <= config_.max_line_search_trials; ++tr)
        {
            restorePositions(db, start_pos);
            const double step = trial_frac * std::min(grid.binWidth(), grid.binHeight());
            for (Cell &c : const_cast<std::vector<Cell> &>(db.cells()))
                if (c.isMovable())
                {
                    c.x += step * cx[c.id];
                    c.y += step * cy[c.id];
                }
            projectMovables(db, grid.coreBounds());
            grid.rebuild(db);
            const ObjectiveMetrics tm = eval.evaluate(db, grid, config_.density_weight);
            rec.line_search_trials = tr;
            const bool ok = acceptCost(tm.total_cost, start.total_cost, config_.relative_improvement_tolerance);
            LOG_DEBUG("[GlobalPlacer] trial=" << tr << " step=" << step << " total=" << tm.total_cost << " delta=" << (tm.total_cost - start.total_cost) << " accepted=" << ok);
            if (ok)
            {
                current = tm;
                rec.accepted = true;
                rec.accepted_step = step;
                rec.hpwl = tm.hpwl;
                rec.density_penalty = tm.density_penalty;
                rec.weighted_density_penalty = tm.weighted_density_penalty;
                rec.total_cost = tm.total_cost;
                rec.total_overflow = tm.total_overflow;
                rec.overflow_ratio = tm.overflow_ratio;
                rec.relative_improvement = (start.total_cost - tm.total_cost) / std::max(std::abs(start.total_cost), kEps);
                ++r.accepted_iterations;
                stalls = 0;
                if (tm.total_cost < r.best_metrics.total_cost)
                {
                    r.best_metrics = tm;
                    best_pos = savePositions(db);
                }
                break;
            }
            trial_frac *= config_.backtrack_factor;
        }
        if (!rec.accepted)
        {
            restorePositions(db, start_pos);
            grid.rebuild(db);
            current = eval.evaluate(db, grid, config_.density_weight);
            rec.hpwl = current.hpwl;
            rec.density_penalty = current.density_penalty;
            rec.weighted_density_penalty = current.weighted_density_penalty;
            rec.total_cost = current.total_cost;
            rec.total_overflow = current.total_overflow;
            rec.overflow_ratio = current.overflow_ratio;
            ++stalls;
            step_fraction *= config_.backtrack_factor;
        }
        r.history.push_back(rec);
        if (it % config_.report_interval == 0)
            LOG_INFO("[GlobalPlacer] iter=" << it << " accepted=" << rec.accepted << " hpwl=" << rec.hpwl << " density=" << rec.density_penalty << " weighted_density=" << rec.weighted_density_penalty << " total=" << rec.total_cost << " overflow=" << rec.total_overflow << " overflow_ratio=" << rec.overflow_ratio << " step_fraction=" << rec.step_fraction << " accepted_step=" << rec.accepted_step << " wire_rms=" << rec.hpwl_gradient_rms << " density_rms=" << rec.density_direction_rms << " direction_rms=" << rec.combined_direction_rms << " line_search_trials=" << rec.line_search_trials << " relative_improvement=" << rec.relative_improvement);
        if (!rec.accepted && step_fraction < config_.minimum_step_fraction)
        {
            r.termination_reason = "step size below minimum";
            break;
        }
        if (!rec.accepted && stalls >= config_.max_stall_iterations)
        {
            r.termination_reason = "line search stalled";
            break;
        }
        if (it == config_.max_iterations)
            r.termination_reason = "maximum iterations reached";
    }

    restorePositions(db, best_pos);
    grid.rebuild(db);
    r.final_metrics = eval.evaluate(db, grid, config_.density_weight);
    if (!std::isfinite(r.final_metrics.total_cost) || r.final_metrics.total_cost > r.initial_metrics.total_cost + 1e-8)
        throw std::runtime_error("GlobalPlacer: invalid final objective");
    for (const Cell &c : db.cells())
    {
        if (!std::isfinite(c.x) || !std::isfinite(c.y))
            throw std::runtime_error("GlobalPlacer: non-finite coordinate at final state");
        if (!c.isMovable() && (c.x != fixed_initial[c.id].x || c.y != fixed_initial[c.id].y))
            throw std::runtime_error("GlobalPlacer: fixed cell moved: " + c.name);
        if (c.isMovable() && (c.x < grid.coreBounds().lx - kEps || c.y < grid.coreBounds().ly - kEps || c.x + c.width > grid.coreBounds().ux + kEps || c.y + c.height > grid.coreBounds().uy + kEps))
            throw std::runtime_error("GlobalPlacer: movable outside core");
    }
    if (r.termination_reason.empty())
        r.termination_reason = "maximum iterations reached";
    LOG_INFO("Global placement completed: attempted_iterations=" << r.attempted_iterations << " accepted_iterations=" << r.accepted_iterations << " initial_total=" << r.initial_metrics.total_cost << " final_total=" << r.final_metrics.total_cost << " relative_improvement=" << ((r.initial_metrics.total_cost - r.final_metrics.total_cost) / std::max(std::abs(r.initial_metrics.total_cost), kEps)) << " termination_reason=" << r.termination_reason);
    return r;
}
