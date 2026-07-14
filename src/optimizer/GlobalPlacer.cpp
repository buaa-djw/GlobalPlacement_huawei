#include "optimizer/GlobalPlacer.h"


#include "db/Cell.h"
#include "db/PlacementDB.h"
#include "geometry/Box.h"
#include "grid/BinGrid.h"
#include "utils/Logger.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
    constexpr double kEps = 1e-12;

    void validateConfig(const GlobalPlacerConfig &c)
    {
        if (c.bins_x <= 0 || c.bins_y <= 0)
            throw std::invalid_argument("GlobalPlacer: bins must be positive");
        if (!std::isfinite(c.target_density) || c.target_density <= 0.0 || c.target_density > 1.0)
            throw std::invalid_argument("GlobalPlacer: invalid target_density");
        if (!std::isfinite(c.density_weight) || c.density_weight < 0.0)
            throw std::invalid_argument("GlobalPlacer: invalid density_weight");
        if (!std::isfinite(c.zero_capacity_repulsion) || c.zero_capacity_repulsion < 1.0)
            throw std::invalid_argument("GlobalPlacer: invalid zero_capacity_repulsion");
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

GlobalPlacerResult GlobalPlacer::optimize(PlacementDB &db) const
{
    const std::size_t movable_count = std::count_if(db.cells().begin(), db.cells().end(), [](const Cell &c)
                                                    { return c.isMovable(); });
    LOG_INFO("Global placement started: movable_cells=" << movable_count << " bins=" << config_.bins_x << "x" << config_.bins_y << " target_density=" << config_.target_density << " density_weight=" << config_.density_weight << " zero_capacity_repulsion=" << config_.zero_capacity_repulsion << " max_iterations=" << config_.max_iterations << " optimizer=" << (config_.method == GlobalPlacementMethod::MoreauProximal ? "moreau" : "direct"));
    BinGrid grid(db, config_.bins_x, config_.bins_y, config_.target_density);
    ObjectiveEvaluator eval;
    SubgradientConfig subgradient_config;
    subgradient_config.density_weight = config_.density_weight;
    subgradient_config.zero_capacity_repulsion = config_.zero_capacity_repulsion;
    SubgradientEvaluator subgradient_evaluator(subgradient_config);
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
    if (config_.method == GlobalPlacementMethod::MoreauProximal)
    {
        int stalls = 0;
        for (int it = 1; it <= config_.max_iterations; ++it)
        {
            r.attempted_iterations = it;
            const auto start_pos = savePositions(db);
            const ObjectiveMetrics start = current;
            MoreauProximalSolver solver(config_.moreau, subgradient_config);
            const MoreauProximalResult pr = solver.solve(db, grid);
            GlobalPlacerIteration rec;
            rec.iteration = it; rec.used_moreau = true; rec.moreau_mu = config_.moreau.mu;
            rec.proximal_inner_iterations = pr.attempted_inner_iterations;
            rec.proximal_accepted_inner_iterations = pr.accepted_inner_iterations;
            rec.proximal_term = pr.final_proximal_metrics.proximal_term;
            rec.proximal_objective = pr.final_proximal_metrics.proximal_objective;
            rec.proximal_displacement_rms = pr.final_proximal_metrics.displacement_rms;
            rec.proximal_residual_rms = pr.final_proximal_metrics.combined_residual_rms;
            rec.proximal_termination_reason = pr.termination_reason;
            rec.accepted = pr.accepted; rec.line_search_trials = pr.attempted_inner_iterations;
            current = pr.accepted ? pr.final_metrics : start;
            if (!pr.accepted) { restorePositions(db, start_pos); grid.rebuild(db); current = eval.evaluate(db, grid, config_.density_weight); ++stalls; }
            else { ++r.accepted_iterations; stalls = 0; if (current.total_cost < r.best_metrics.total_cost) { r.best_metrics = current; best_pos = savePositions(db); } }
            rec.hpwl=current.hpwl; rec.density_penalty=current.density_penalty; rec.normalized_hpwl=current.normalized_hpwl; rec.normalized_density_penalty=current.normalized_density_penalty; rec.weighted_density_penalty=current.weighted_density_penalty; rec.total_cost=current.total_cost; rec.total_overflow=current.total_overflow; rec.overflow_ratio=current.overflow_ratio;
            rec.relative_improvement=(start.total_cost-current.total_cost)/std::max(std::abs(start.total_cost), kEps);
            r.history.push_back(rec);
            if (it % config_.report_interval == 0) LOG_INFO("[GlobalPlacer] iter=" << it << " optimizer=moreau accepted=" << rec.accepted << " raw_hpwl=" << rec.hpwl << " normalized_hpwl=" << rec.normalized_hpwl << " raw_density=" << rec.density_penalty << " normalized_density=" << rec.normalized_density_penalty << " total_normalized=" << rec.total_cost << " density_weight=" << config_.density_weight << " zero_capacity_repulsion=" << config_.zero_capacity_repulsion << " overflow_ratio=" << rec.overflow_ratio << " moreau_mu=" << rec.moreau_mu << " proximal_inner_iterations=" << rec.proximal_inner_iterations << " proximal_term=" << rec.proximal_term << " proximal_objective=" << rec.proximal_objective << " proximal_residual_rms=" << rec.proximal_residual_rms << " proximal_termination_reason=" << rec.proximal_termination_reason);
            if (!rec.accepted && stalls >= config_.max_stall_iterations) { r.termination_reason = "line search stalled"; break; }
            if (it == config_.max_iterations) r.termination_reason = "maximum iterations reached";
        }
        restorePositions(db, best_pos); grid.rebuild(db);
        r.final_metrics = eval.evaluate(db, grid, config_.density_weight);
        for (const Cell &c : db.cells()) { if (!std::isfinite(c.x) || !std::isfinite(c.y)) throw std::runtime_error("GlobalPlacer: non-finite coordinate at final state"); if (!c.isMovable() && (c.x != fixed_initial[c.id].x || c.y != fixed_initial[c.id].y)) throw std::runtime_error("GlobalPlacer: fixed cell moved: " + c.name); if (c.isMovable() && (c.x < grid.coreBounds().lx - kEps || c.y < grid.coreBounds().ly - kEps || c.x + c.width > grid.coreBounds().ux + kEps || c.y + c.height > grid.coreBounds().uy + kEps)) throw std::runtime_error("GlobalPlacer: movable outside core"); }
        if (r.termination_reason.empty()) r.termination_reason = "maximum iterations reached";
        LOG_INFO("Global placement completed: optimizer=moreau attempted_iterations=" << r.attempted_iterations << " accepted_iterations=" << r.accepted_iterations << " initial_total=" << r.initial_metrics.total_cost << " final_total=" << r.final_metrics.total_cost << " termination_reason=" << r.termination_reason);
        return r;
    }

    double step_fraction = config_.initial_step_fraction;
    int stalls = 0;
    for (int it = 1; it <= config_.max_iterations; ++it)
    {
        r.attempted_iterations = it;
        const auto start_pos = savePositions(db);
        const ObjectiveMetrics start = current;
        const SubgradientResult gradient = subgradient_evaluator.evaluate(db, grid);
        if (gradient.metrics.combined_direction_rms <= kEps)
        {
            r.converged = true;
            r.termination_reason = "zero search direction";
            break;
        }
        const std::vector<double>& cx = gradient.combined_direction_x;
        const std::vector<double>& cy = gradient.combined_direction_y;
        GlobalPlacerIteration rec;
        rec.iteration = it;
        rec.step_fraction = step_fraction;
        rec.hpwl_gradient_rms = gradient.metrics.hpwl_subgradient_rms;
        rec.density_direction_rms = gradient.metrics.density_direction_rms;
        rec.combined_direction_rms = gradient.metrics.combined_direction_rms;
        rec.active_net_count = gradient.metrics.active_net_count;
        rec.active_overflow_bin_count = gradient.metrics.active_overflow_bin_count;
        rec.zero_capacity_overflow_bin_count = gradient.metrics.zero_capacity_overflow_bin_count;
        rec.maximum_direction_norm = gradient.metrics.maximum_direction_norm;
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
                rec.normalized_hpwl = tm.normalized_hpwl;
                rec.normalized_density_penalty = tm.normalized_density_penalty;
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
            rec.normalized_hpwl = current.normalized_hpwl;
            rec.normalized_density_penalty = current.normalized_density_penalty;
            rec.weighted_density_penalty = current.weighted_density_penalty;
            rec.total_cost = current.total_cost;
            rec.total_overflow = current.total_overflow;
            rec.overflow_ratio = current.overflow_ratio;
            ++stalls;
            step_fraction *= config_.backtrack_factor;
        }
        r.history.push_back(rec);
        if (it % config_.report_interval == 0)
            LOG_INFO("[GlobalPlacer] iter=" << it << " optimizer=direct accepted=" << rec.accepted << " raw_hpwl=" << rec.hpwl << " normalized_hpwl=" << rec.normalized_hpwl << " raw_density=" << rec.density_penalty << " normalized_density=" << rec.normalized_density_penalty << " weighted_density=" << rec.weighted_density_penalty << " total_normalized=" << rec.total_cost << " density_weight=" << config_.density_weight << " zero_capacity_repulsion=" << config_.zero_capacity_repulsion << " overflow=" << rec.total_overflow << " overflow_ratio=" << rec.overflow_ratio << " active_overflow_bins=" << rec.active_overflow_bin_count << " zero_capacity_overflow_bins=" << rec.zero_capacity_overflow_bin_count << " step_fraction=" << rec.step_fraction << " accepted_step=" << rec.accepted_step << " hpwl_rms=" << rec.hpwl_gradient_rms << " density_rms=" << rec.density_direction_rms << " direction_rms=" << rec.combined_direction_rms << " max_direction_norm=" << rec.maximum_direction_norm << " line_search_trials=" << rec.line_search_trials << " relative_improvement=" << rec.relative_improvement);
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
