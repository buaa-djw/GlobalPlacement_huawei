#include "optimizer/SubgradientEvaluator.h"

#include "db/PlacementDB.h"
#include "grid/BinGrid.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double EPS = 1e-9;

void near(double a, double b, double eps = EPS) { assert(std::abs(a - b) <= eps); }
void row(PlacementDB& db, double w = 100.0, double h = 100.0) { Row r; r.x_start = 0; r.x_end = w; r.y = 0; r.height = h; r.site_width = 1; r.site_spacing = 1; r.num_sites = static_cast<int>(w); db.addRow(r); }
int cell(PlacementDB& db, const std::string& n, double x, double y, double w = 2.0, double h = 2.0, CellType t = CellType::Standard, bool fixed = false) { int id = db.addCell(n, w, h, t); db.setCellLocation(n, x, y, fixed); return id; }
void net2(PlacementDB& db, int a, int b) { int n = db.addNet("n" + std::to_string(db.nets().size())); db.addPin(a, n, 0, 0, "B"); db.addPin(b, n, 0, 0, "B"); }
void net3(PlacementDB& db, int a, int b, int c) { int n = db.addNet("n" + std::to_string(db.nets().size())); db.addPin(a, n, 0, 0, "B"); db.addPin(b, n, 0, 0, "B"); db.addPin(c, n, 0, 0, "B"); }
double norm(double x, double y) { return std::hypot(x, y); }

void assertFiniteResult(const SubgradientResult& r) {
    std::vector<const std::vector<double>*> vs{&r.hpwl_gradient_x, &r.hpwl_gradient_y, &r.hpwl_direction_x, &r.hpwl_direction_y, &r.density_direction_x, &r.density_direction_y, &r.objective_direction_x, &r.objective_direction_y, &r.combined_direction_x, &r.combined_direction_y};
    for (auto* v : vs) for (double x : *v) assert(std::isfinite(x));
    assert(std::isfinite(r.metrics.core_reference_length));
    assert(std::isfinite(r.metrics.hpwl_subgradient_rms));
    assert(std::isfinite(r.metrics.density_direction_rms));
    assert(std::isfinite(r.metrics.combined_direction_rms));
    assert(std::isfinite(r.metrics.maximum_direction_norm));
}

void assertCellZero(const SubgradientResult& r, int id) {
    near(r.hpwl_gradient_x[id], 0); near(r.hpwl_gradient_y[id], 0); near(r.hpwl_direction_x[id], 0); near(r.hpwl_direction_y[id], 0);
    near(r.density_direction_x[id], 0); near(r.density_direction_y[id], 0); near(r.objective_direction_x[id], 0); near(r.objective_direction_y[id], 0); near(r.combined_direction_x[id], 0); near(r.combined_direction_y[id], 0);
}
}

int main() {
    { // two-pin normalized HPWL
        PlacementDB db; row(db, 100, 100); int a = cell(db, "a", 0, 0), b = cell(db, "b", 10, 0); net2(db, a, b); BinGrid g(db, 1, 1, 1); auto r = SubgradientEvaluator({}).evaluate(db, g);
        const double scale = 1.0 / std::hypot(100.0, 100.0); near(r.hpwl_gradient_x[a], -scale); near(r.hpwl_gradient_x[b], scale); assert(r.hpwl_direction_x[a] > 0); assert(r.hpwl_direction_x[b] < 0); assert(r.metrics.active_net_count == 1);
    }
    { // fixed cell zero
        PlacementDB db; row(db); int a = cell(db, "a", 0, 0), f = cell(db, "f", 10, 0, 2, 2, CellType::Standard, true); net2(db, a, f); BinGrid g(db, 1, 1, 1); auto r = SubgradientEvaluator({}).evaluate(db, g); assertCellZero(r, f);
    }
    { // max extreme split
        PlacementDB db; row(db); int l = cell(db, "l", 0, 0), r0 = cell(db, "r0", 10, 0), r1 = cell(db, "r1", 10, 5); net3(db, l, r0, r1); BinGrid g(db, 1, 1, 1); auto r = SubgradientEvaluator({}).evaluate(db, g); double s = 1.0 / std::hypot(100.0, 100.0); near(r.hpwl_gradient_x[r0], 0.5 * s); near(r.hpwl_gradient_x[r1], 0.5 * s); near(r.hpwl_gradient_x[r0] + r.hpwl_gradient_x[r1], s);
    }
    { // no nets
        PlacementDB db; row(db); cell(db, "a", 0, 0); BinGrid g(db, 1, 1, 1); auto r = SubgradientEvaluator({}).evaluate(db, g); assert(r.metrics.active_net_count == 0); near(r.metrics.hpwl_subgradient_rms, 0); near(r.metrics.hpwl_direction_rms, 0); assertFiniteResult(r);
    }
    { // no overflow
        PlacementDB db; row(db); cell(db, "a", 0, 0, 1, 1); BinGrid g(db, 2, 2, 1); auto r = SubgradientEvaluator({}).evaluate(db, g); assert(r.metrics.active_overflow_bin_count == 0); for (double x : r.density_direction_x) near(x, 0); for (double y : r.density_direction_y) near(y, 0);
    }
    { // overflow pushes outward
        PlacementDB db; row(db, 10, 10); int a = cell(db, "a", 0, 0, 8, 8), b = cell(db, "b", 1, 1, 8, 8); BinGrid g(db, 1, 1, 0.5); auto r = SubgradientEvaluator({}).evaluate(db, g); assert(r.metrics.active_overflow_bin_count == 1); assert(norm(r.density_direction_x[a], r.density_direction_y[a]) > 0); assert((r.density_direction_x[a] * ((0 + 4) - 5) + r.density_direction_y[a] * ((0 + 4) - 5)) > 0);
    }
    { // tie deterministic
        PlacementDB db; row(db, 10, 10); int a = cell(db, "a", 4, 4, 2, 2); cell(db, "b", 0, 0, 8, 8); BinGrid g(db, 1, 1, 0.5); SubgradientEvaluator e({}); auto r1 = e.evaluate(db, g); auto r2 = e.evaluate(db, g); assert(norm(r1.density_direction_x[a], r1.density_direction_y[a]) > 0); near(r1.density_direction_x[a], r2.density_direction_x[a], 0); near(r1.density_direction_y[a], r2.density_direction_y[a], 0);
    }
    { // zero-cap repulsion larger
        PlacementDB normal; row(normal, 10, 10); int n = cell(normal, "m", 4, 4, 4, 4); BinGrid gn(normal, 1, 1, 0.1);
        PlacementDB zero; row(zero, 10, 10); cell(zero, "fix", 0, 0, 10, 10, CellType::Terminal, true); int z = cell(zero, "m", 4, 4, 4, 4); BinGrid gz(zero, 1, 1, 1.0);
        SubgradientConfig c; c.zero_capacity_repulsion = 2.0; auto rn = SubgradientEvaluator(c).evaluate(normal, gn); auto rz = SubgradientEvaluator(c).evaluate(zero, gz); assert(norm(rz.density_direction_x[z], rz.density_direction_y[z]) > norm(rn.density_direction_x[n], rn.density_direction_y[n]));
    }
    { // density_weight zero combines normalized HPWL only
        PlacementDB db; row(db); int a = cell(db, "a", 0, 0), b = cell(db, "b", 10, 0); net2(db, a, b); BinGrid g(db, 1, 1, 1); SubgradientConfig c; c.density_weight = 0; auto r = SubgradientEvaluator(c).evaluate(db, g); double rms = r.metrics.hpwl_direction_rms; near(r.combined_direction_x[a], r.hpwl_direction_x[a] / rms); near(r.combined_direction_x[b], r.hpwl_direction_x[b] / rms);
    }
    { // objective direction composition
        PlacementDB db; row(db, 10, 10); int a = cell(db, "a", 0, 0, 8, 8), b = cell(db, "b", 1, 1, 8, 8); net2(db, a, b); BinGrid g(db, 1, 1, 0.5); SubgradientConfig c; c.density_weight = 3.0; auto r = SubgradientEvaluator(c).evaluate(db, g); for (const Cell& ce : db.cells()) if (ce.isMovable()) { near(r.objective_direction_x[ce.id], r.hpwl_direction_x[ce.id] + c.density_weight * r.density_direction_x[ce.id]); near(r.objective_direction_y[ce.id], r.hpwl_direction_y[ce.id] + c.density_weight * r.density_direction_y[ce.id]); } double rms = r.metrics.combined_direction_rms_before_normalization; if (rms > 0) { near(r.combined_direction_x[a], r.objective_direction_x[a] / rms); near(r.combined_direction_y[a], r.objective_direction_y[a] / rms); }
    }
    { // non-movable Standard fixed, Terminal, TerminalNI zero and finite
        PlacementDB db; row(db); int m = cell(db, "m", 0, 0), f = cell(db, "f", 2, 0, 2, 2, CellType::Standard, true), t = cell(db, "t", 4, 0, 2, 2, CellType::Terminal, true), ni = cell(db, "ni", 6, 0, 2, 2, CellType::TerminalNI, true); net2(db, m, f); BinGrid g(db, 1, 1, 1); auto r = SubgradientEvaluator({}).evaluate(db, g); assertFiniteResult(r); assertCellZero(r, f); assertCellZero(r, t); assertCellZero(r, ni);
    }
    { // invalid config
        auto throws = [](SubgradientConfig c){ bool ok=false; try { SubgradientEvaluator e(c); } catch (const std::invalid_argument&) { ok=true; } assert(ok); };
        SubgradientConfig c; c.density_weight = -1; throws(c); c = {}; c.density_weight = std::numeric_limits<double>::quiet_NaN(); throws(c); c = {}; c.zero_capacity_repulsion = 0.5; throws(c); c = {}; c.zero_capacity_repulsion = std::numeric_limits<double>::infinity(); throws(c); c = {}; c.numerical_epsilon = 0; throws(c); c = {}; c.extreme_tolerance = 0; throws(c);
    }
    { // read-only and deterministic
        PlacementDB db; row(db); int a = cell(db, "a", 0, 0, 2, 3); db.cell(a).orientation = "FN"; int b = cell(db, "b", 10, 0); net2(db, a, b); auto before = db.cells(); BinGrid g(db, 1, 1, 1); SubgradientEvaluator e({}); auto r1 = e.evaluate(db, g); auto r2 = e.evaluate(db, g); for (std::size_t i = 0; i < db.cells().size(); ++i) { const Cell& x = db.cells()[i]; const Cell& y = before[i]; near(x.x, y.x); near(x.y, y.y); near(x.width, y.width); near(x.height, y.height); assert(x.is_fixed == y.is_fixed); assert(x.type == y.type); assert(x.orientation == y.orientation); }
        near(r1.metrics.combined_direction_rms, r2.metrics.combined_direction_rms, 0); for (std::size_t i=0;i<r1.combined_direction_x.size();++i) { near(r1.combined_direction_x[i], r2.combined_direction_x[i], 0); near(r1.combined_direction_y[i], r2.combined_direction_y[i], 0); }
    }
    return 0;
}
