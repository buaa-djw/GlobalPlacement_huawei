#include "optimizer/MoreauProximalSolver.h"

#include "db/Cell.h"
#include "db/PlacementDB.h"
#include "geometry/Box.h"
#include "grid/BinGrid.h"
#include "utils/Logger.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double kEps = 1e-12;
struct ProximalCellPosition { double x = 0.0; double y = 0.0; };

void validateConfig(const MoreauProximalConfig& c) {
    if (!std::isfinite(c.mu) || c.mu <= 0.0) throw std::invalid_argument("MoreauProximalSolver: invalid mu");
    if (c.max_inner_iterations < 0) throw std::invalid_argument("MoreauProximalSolver: invalid max_inner_iterations");
    if (!std::isfinite(c.initial_step_fraction) || c.initial_step_fraction <= 0.0) throw std::invalid_argument("MoreauProximalSolver: invalid initial_step_fraction");
    if (!std::isfinite(c.minimum_step_fraction) || c.minimum_step_fraction <= 0.0 || c.minimum_step_fraction > c.initial_step_fraction) throw std::invalid_argument("MoreauProximalSolver: invalid minimum_step_fraction");
    if (!std::isfinite(c.backtrack_factor) || c.backtrack_factor <= 0.0 || c.backtrack_factor >= 1.0) throw std::invalid_argument("MoreauProximalSolver: invalid backtrack_factor");
    if (c.max_line_search_trials <= 0 || c.max_stall_iterations <= 0) throw std::invalid_argument("MoreauProximalSolver: invalid line search limits");
    if (!std::isfinite(c.relative_improvement_tolerance) || c.relative_improvement_tolerance < 0.0) throw std::invalid_argument("MoreauProximalSolver: invalid relative_improvement_tolerance");
    if (!std::isfinite(c.residual_rms_tolerance) || c.residual_rms_tolerance < 0.0) throw std::invalid_argument("MoreauProximalSolver: invalid residual_rms_tolerance");
    if (!std::isfinite(c.displacement_rms_tolerance) || c.displacement_rms_tolerance < 0.0) throw std::invalid_argument("MoreauProximalSolver: invalid displacement_rms_tolerance");
    if (c.report_interval <= 0) throw std::invalid_argument("MoreauProximalSolver: invalid report_interval");
}

std::vector<ProximalCellPosition> savePositions(const PlacementDB& db) { std::vector<ProximalCellPosition> p(db.cells().size()); for (const Cell& c: db.cells()) p.at(c.id)={c.x,c.y}; return p; }
void restorePositions(PlacementDB& db, const std::vector<ProximalCellPosition>& p) { if (p.size()!=db.cells().size()) throw std::invalid_argument("MoreauProximalSolver: snapshot size mismatch"); for (Cell& c: const_cast<std::vector<Cell>&>(db.cells())) { c.x=p.at(c.id).x; c.y=p.at(c.id).y; } }
std::size_t movableCount(const PlacementDB& db) { return std::count_if(db.cells().begin(), db.cells().end(), [](const Cell& c){return c.isMovable();}); }
double computeSquaredDisplacement(const PlacementDB& db, const std::vector<ProximalCellPosition>& a) { double s=0; for (const Cell& c: db.cells()) if (c.isMovable()) { double dx=c.x-a.at(c.id).x, dy=c.y-a.at(c.id).y; s += dx*dx+dy*dy; } return s; }
double movableRms(const PlacementDB& db, const std::vector<double>& x, const std::vector<double>& y) { double s=0; std::size_t n=0; for (const Cell& c: db.cells()) if (c.isMovable()) { s += x.at(c.id)*x.at(c.id)+y.at(c.id)*y.at(c.id); ++n; } return n?std::sqrt(s/static_cast<double>(n)):0.0; }
void ensureFinite(double v, const char* n) { if (!std::isfinite(v)) throw std::runtime_error(std::string("MoreauProximalSolver: non-finite ")+n); }
bool accept(double trial, double cur, double tol) { if (!std::isfinite(trial)) return false; return tol==0.0 ? trial < cur : trial <= cur - tol*std::max(1.0,std::abs(cur)); }
void project(PlacementDB& db, const Box& core) { for (Cell& c: const_cast<std::vector<Cell>&>(db.cells())) if (c.isMovable()) { c.x=std::clamp(c.x, core.lx, core.ux-c.width); c.y=std::clamp(c.y, core.ly, core.uy-c.height); } }

MoreauProximalMetrics metricsFor(const PlacementDB& db, const BinGrid& grid, const std::vector<ProximalCellPosition>& anchor, const MoreauProximalConfig& cfg, double density_weight) {
    MoreauProximalMetrics m; m.objective = ObjectiveEvaluator{}.evaluate(db, grid, density_weight, false);
    const double l = std::hypot(grid.coreBounds().width(), grid.coreBounds().height()); if (!std::isfinite(l)||l<=0) throw std::runtime_error("MoreauProximalSolver: invalid core reference length");
    m.squared_displacement = computeSquaredDisplacement(db, anchor);
    m.proximal_term = m.squared_displacement/(2.0*cfg.mu*l*l);
    m.proximal_objective = m.objective.total_cost + m.proximal_term;
    m.displacement_rms = std::sqrt(m.squared_displacement/static_cast<double>(std::max<std::size_t>(movableCount(db),1)));
    ensureFinite(m.proximal_term,"proximal_term"); ensureFinite(m.proximal_objective,"proximal_objective"); ensureFinite(m.displacement_rms,"displacement_rms"); return m;
}

MoreauProximalIteration recFrom(int it, bool accepted, int trials, double sf, double step, const MoreauProximalMetrics& m, double rel) {
    MoreauProximalIteration r; r.inner_iteration=it; r.accepted=accepted; r.line_search_trials=trials; r.step_fraction=sf; r.accepted_step=step; r.objective_cost=m.objective.total_cost; r.proximal_term=m.proximal_term; r.proximal_objective=m.proximal_objective; r.squared_displacement=m.squared_displacement; r.displacement_rms=m.displacement_rms; r.objective_direction_rms=m.objective_direction_rms; r.proximal_direction_rms=m.proximal_direction_rms; r.combined_residual_rms=m.combined_residual_rms; r.normalized_search_direction_rms=m.normalized_search_direction_rms; r.maximum_search_direction_norm=m.maximum_search_direction_norm; r.relative_proximal_improvement=rel; return r;
}
}

MoreauProximalSolver::MoreauProximalSolver(MoreauProximalConfig config, SubgradientConfig subgradient_config) : config_(config), subgradient_config_(subgradient_config) { validateConfig(config_); (void)SubgradientEvaluator(subgradient_config_); }

MoreauProximalResult MoreauProximalSolver::solve(PlacementDB& db, BinGrid& grid) const {
    const auto anchor = savePositions(db); const std::size_t movables = movableCount(db);
    MoreauProximalResult r; r.anchor_metrics = ObjectiveEvaluator{}.evaluate(db, grid, subgradient_config_.density_weight, true); r.anchor_proximal_metrics = metricsFor(db, grid, anchor, config_, subgradient_config_.density_weight); r.best_metrics=r.anchor_metrics; r.best_proximal_metrics=r.anchor_proximal_metrics; auto best_pos=anchor;
    LOG_INFO("[Moreau] solve started: mu="<<config_.mu<<" inner_iterations="<<config_.max_inner_iterations<<" anchor_total="<<r.anchor_metrics.total_cost<<" anchor_proximal="<<r.anchor_proximal_metrics.proximal_objective);
    if (config_.max_inner_iterations == 0) { r.final_metrics=r.anchor_metrics; r.final_proximal_metrics=r.anchor_proximal_metrics; r.termination_reason="zero inner iterations requested"; LOG_INFO("[Moreau] solve completed: accepted=0 attempted_inner_iterations=0 accepted_inner_iterations=0 anchor_objective="<<r.anchor_metrics.total_cost<<" final_objective="<<r.final_metrics.total_cost<<" anchor_proximal="<<r.anchor_proximal_metrics.proximal_objective<<" final_proximal="<<r.final_proximal_metrics.proximal_objective<<" displacement_rms=0 termination_reason="<<r.termination_reason); return r; }
    double step_fraction=config_.initial_step_fraction; int stalls=0; const double l2=std::pow(std::hypot(grid.coreBounds().width(), grid.coreBounds().height()),2);
    for (int it=1; it<=config_.max_inner_iterations; ++it) {
        r.attempted_inner_iterations=it; const auto start_pos=savePositions(db); MoreauProximalMetrics cur=metricsFor(db, grid, anchor, config_, subgradient_config_.density_weight);
        const SubgradientResult sg=SubgradientEvaluator(subgradient_config_).evaluate(db, grid); std::vector<double> px(db.cells().size(),0),py(px.size(),0),rx(px.size(),0),ry(px.size(),0),sx(px.size(),0),sy(px.size(),0);
        for (const Cell& c: db.cells()) if (c.isMovable()) { px.at(c.id)=-(c.x-anchor.at(c.id).x)/(config_.mu*l2); py.at(c.id)=-(c.y-anchor.at(c.id).y)/(config_.mu*l2); rx.at(c.id)=sg.objective_direction_x.at(c.id)+px.at(c.id); ry.at(c.id)=sg.objective_direction_y.at(c.id)+py.at(c.id); }
        cur.objective_direction_rms=movableRms(db, sg.objective_direction_x, sg.objective_direction_y); cur.proximal_direction_rms=movableRms(db, px, py); cur.combined_residual_rms=movableRms(db, rx, ry);
        if (cur.combined_residual_rms <= config_.residual_rms_tolerance) { r.converged=true; r.termination_reason="proximal residual below tolerance"; break; }
        for (const Cell& c: db.cells()) if (c.isMovable()) { sx.at(c.id)=rx.at(c.id)/cur.combined_residual_rms; sy.at(c.id)=ry.at(c.id)/cur.combined_residual_rms; cur.maximum_search_direction_norm=std::max(cur.maximum_search_direction_norm,std::hypot(sx.at(c.id),sy.at(c.id))); }
        cur.normalized_search_direction_rms=movableRms(db,sx,sy);
        bool ok=false; MoreauProximalMetrics accepted_m; double accepted_step=0, trial_frac=step_fraction; int trials=0;
        for (int tr=1; tr<=config_.max_line_search_trials; ++tr) { trials=tr; restorePositions(db,start_pos); double step=trial_frac*std::min(grid.binWidth(),grid.binHeight()); for (Cell& c: const_cast<std::vector<Cell>&>(db.cells())) if (c.isMovable()) { c.x += step*sx.at(c.id); c.y += step*sy.at(c.id); } project(db, grid.coreBounds()); grid.rebuild(db); MoreauProximalMetrics tm=metricsFor(db,grid,anchor,config_,subgradient_config_.density_weight); if (accept(tm.proximal_objective, cur.proximal_objective, config_.relative_improvement_tolerance)) { ok=true; accepted_m=tm; accepted_m.objective_direction_rms=cur.objective_direction_rms; accepted_m.proximal_direction_rms=cur.proximal_direction_rms; accepted_m.combined_residual_rms=cur.combined_residual_rms; accepted_m.normalized_search_direction_rms=cur.normalized_search_direction_rms; accepted_m.maximum_search_direction_norm=cur.maximum_search_direction_norm; accepted_step=step; break; } trial_frac *= config_.backtrack_factor; }
        if (ok) { ++r.accepted_inner_iterations; stalls=0; double rel=(cur.proximal_objective-accepted_m.proximal_objective)/std::max(std::abs(cur.proximal_objective),kEps); r.history.push_back(recFrom(it,true,trials,step_fraction,accepted_step,accepted_m,rel)); if (accepted_m.proximal_objective < r.best_proximal_metrics.proximal_objective) { r.best_proximal_metrics=accepted_m; r.best_metrics=accepted_m.objective; best_pos=savePositions(db); } }
        else { restorePositions(db,start_pos); grid.rebuild(db); ++stalls; step_fraction *= config_.backtrack_factor; r.history.push_back(recFrom(it,false,trials,step_fraction,0,cur,0)); }
        const auto& last=r.history.back(); if (it % config_.report_interval==0) LOG_INFO("[Moreau] inner_iter="<<it<<" accepted="<<last.accepted<<" objective="<<last.objective_cost<<" proximal_term="<<last.proximal_term<<" proximal_objective="<<last.proximal_objective<<" displacement_rms="<<last.displacement_rms<<" objective_direction_rms="<<last.objective_direction_rms<<" proximal_direction_rms="<<last.proximal_direction_rms<<" residual_rms="<<last.combined_residual_rms<<" step_fraction="<<last.step_fraction<<" accepted_step="<<last.accepted_step<<" line_search_trials="<<last.line_search_trials<<" relative_improvement="<<last.relative_proximal_improvement);
        if (!ok && step_fraction < config_.minimum_step_fraction) { r.termination_reason="proximal step size below minimum"; break; }
        if (!ok && stalls >= config_.max_stall_iterations) { r.termination_reason="proximal line search stalled"; break; }
        if (ok && accepted_m.displacement_rms <= config_.displacement_rms_tolerance) { r.converged=true; r.termination_reason="proximal displacement below tolerance"; break; }
    }
    if (r.termination_reason.empty()) r.termination_reason="maximum inner iterations reached";
    const double tol=1e-10*std::max(1.0,std::abs(r.anchor_metrics.total_cost));
    r.accepted = r.best_proximal_metrics.proximal_objective < r.anchor_proximal_metrics.proximal_objective && r.best_metrics.total_cost <= r.anchor_metrics.total_cost + tol;
    restorePositions(db, r.accepted ? best_pos : anchor); grid.rebuild(db); if (!r.accepted) r.termination_reason="no proximal improvement";
    r.final_metrics=ObjectiveEvaluator{}.evaluate(db,grid,subgradient_config_.density_weight,true); r.final_proximal_metrics=metricsFor(db,grid,anchor,config_,subgradient_config_.density_weight);
    for (const Cell& c: db.cells()) if (!c.isMovable() && (c.x != anchor.at(c.id).x || c.y != anchor.at(c.id).y)) throw std::runtime_error("MoreauProximalSolver: fixed or terminal cell moved: "+c.name);
    LOG_INFO("[Moreau] solve completed: accepted="<<r.accepted<<" attempted_inner_iterations="<<r.attempted_inner_iterations<<" accepted_inner_iterations="<<r.accepted_inner_iterations<<" anchor_objective="<<r.anchor_metrics.total_cost<<" final_objective="<<r.final_metrics.total_cost<<" anchor_proximal="<<r.anchor_proximal_metrics.proximal_objective<<" final_proximal="<<r.final_proximal_metrics.proximal_objective<<" displacement_rms="<<r.final_proximal_metrics.displacement_rms<<" termination_reason="<<r.termination_reason);
    (void)movables; return r;
}
