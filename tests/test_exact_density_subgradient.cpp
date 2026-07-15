#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "density/ExactDensitySubgradient.h"
#include "density/PaperOverlapFunction.h"

#include <cassert>
#include <cmath>
#include <random>

static bool near(double a, double b, double eps = 1e-9) { return std::abs(a - b) <= eps * std::max(1.0, std::abs(b)); }
static void row(PlacementDB& db, double w = 10, double h = 10) { Row r; r.y = 0; r.height = h; r.site_width = 1; r.site_spacing = 1; r.x_start = 0; r.num_sites = static_cast<int>(w); db.addRow(r); }
static int cell(PlacementDB& db, const char* n, double x, double y, double w, double h, bool f = false, CellType t = CellType::Standard) { int id = db.addCell(n, w, h, t); db.setCellLocation(n, x, y, f); return id; }
static double rel(double a, double b) { return std::abs(a - b) / std::max(1.0, std::abs(b)); }

int main()
{
    {
        std::mt19937_64 rng(1);
        assert(near(paperAxisOverlapLength(0, 2, 0, 10), 2));
        assert(near(paperAxisOverlapLength(-5, 2, 0, 10), 1));
        assert(near(paperAxisOverlapLength(5, 2, 0, 10), 1));
        assert(near(paperAxisOverlapLength(7, 2, 0, 10), 0));
        assert(near(paperAxisOverlapLength(0, 10, 0, 2), 2));
        assert(near(paperAxisOverlapLength(-5, 10, 0, 2), 1));
        assert(near(paperAxisOverlapLength(5, 10, 0, 2), 1));
        assert(near(paperAxisOverlapLength(7, 10, 0, 2), 0));
        assert(near(paperAxisOverlapSubgradient(-5, 2, 0, 10, rng).subgradient, 1));
        assert(near(paperAxisOverlapSubgradient(5, 2, 0, 10, rng).subgradient, -1));
        assert(near(paperAxisOverlapSubgradient(0, 2, 0, 10, rng).subgradient, 0));
        assert(near(paperAxisOverlapSubgradient(7, 2, 0, 10, rng).subgradient, 0));
    }
    {
        for (double z : {4.0, 6.0}) { std::mt19937_64 a(7), b(7); auto ga = paperAxisOverlapSubgradient(z, 2, 0, 10, a).subgradient; auto gb = paperAxisOverlapSubgradient(z, 2, 0, 10, b).subgradient; assert(ga >= -1 && ga <= 0); assert(ga == gb); }
        for (double z : {-4.0, -6.0}) { std::mt19937_64 a(9), b(9); auto ga = paperAxisOverlapSubgradient(z, 2, 0, 10, a).subgradient; auto gb = paperAxisOverlapSubgradient(z, 2, 0, 10, b).subgradient; assert(ga >= 0 && ga <= 1); assert(ga == gb); }
        std::mt19937_64 a(11), b(11); double g = paperAxisOverlapSubgradient(0, 10, 0, 10, a).subgradient; assert(g >= -1 && g <= 1); assert(g == paperAxisOverlapSubgradient(0, 10, 0, 10, b).subgradient);
    }
    {
        PlacementDB db; row(db, 20, 10); int m = cell(db, "touch", 10, 2, 2, 2); cell(db, "fill", 0, 0, 10, 10, true);
        ExactBinGrid grid(db, 2, 1, 0.5); assert(near(grid.bin(0).total_area, 100));
        std::mt19937_64 rng(3); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        assert(r.grad_x[m] <= 0); assert(std::abs(r.grad_x[m]) > 0); assert(near(r.grad_y[m], 0));
    }
    {
        PlacementDB db; row(db, 10, 10); int m = cell(db, "m", 0, 0, 1, 1);
        ExactBinGrid grid(db, 1, 1, 0.9); std::mt19937_64 rng(1); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        assert(near(r.penalty, 0)); assert(near(r.grad_x[m], 0)); assert(near(r.grad_y[m], 0));
    }
    {
        PlacementDB db; row(db, 10, 10); int m = cell(db, "m", 2, 3, 6, 4);
        ExactBinGrid grid(db, 1, 1, 0.1); std::mt19937_64 rng(1); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        double D = 24, cap = 10, ov = 14; assert(near(r.penalty, ov * ov)); assert(near(r.ofr, ov / D)); assert(near(r.grad_x[m], 0)); assert(near(r.grad_y[m], 0));
    }
    {
        PlacementDB db; row(db, 20, 10); int m = cell(db, "m", 7, 2, 6, 4);
        ExactBinGrid grid(db, 2, 1, 0.1); std::mt19937_64 rng(1); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        // left bin: overflow 2, dD/dx=+4; right bin: overflow 2, dD/dx=-4.
        assert(near(r.grad_x[m], 0)); assert(near(r.grad_y[m], 0));
    }
    {
        PlacementDB db; row(db, 10, 10); int f = cell(db, "f", 0, 0, 10, 10, true); int m = cell(db, "m", 2, 2, 2, 2);
        ExactBinGrid grid(db, 1, 1, 0.5); std::mt19937_64 rng(1); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        assert(near(r.grad_x[f], 0)); assert(near(r.grad_y[f], 0)); assert(near(r.ofr, 54.0 / 104.0));
    }
    {
        PlacementDB db; row(db, 10, 10); cell(db, "m", 0, 0, 2, 2); cell(db, "f", 0, 0, 2, 2, true); cell(db, "t", 100, 100, 2, 2, true, CellType::Terminal); cell(db, "o", 200, 200, 2, 2);
        ExactBinGrid grid(db, 1, 1, 0.01); std::mt19937_64 rng(1); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        assert(near(r.ofr, 7.0 / 16.0));
    }
    {
        PlacementDB db; row(db, 10, 10); int m = cell(db, "m", 3.3, 2.7, 4.2, 3.1);
        ExactBinGrid grid(db, 1, 1, 0.05); std::mt19937_64 rng(1); auto r = ExactDensitySubgradient().evaluate(db, grid, rng);
        double eps = 1e-5;
        db.setCellPosition(m, 3.3 + eps, 2.7); double pp = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.05)).penalty;
        db.setCellPosition(m, 3.3 - eps, 2.7); double pm = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.05)).penalty;
        db.setCellPosition(m, 3.3, 2.7 + eps); double qp = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.05)).penalty;
        db.setCellPosition(m, 3.3, 2.7 - eps); double qm = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.05)).penalty;
        assert(rel((pp - pm) / (2 * eps), r.grad_x[m]) <= 1e-5); assert(rel((qp - qm) / (2 * eps), r.grad_y[m]) <= 1e-5);
    }
    {
        PlacementDB db; row(db, 10, 10); cell(db, "m", 3.3, 2.7, 4.2, 3.1);
        ExactBinGrid grid(db, 1, 1, 0.05); std::mt19937_64 rng(1); auto s = ExactDensitySubgradient().evaluate(db, grid, rng); auto e = ExactDensityEvaluator().evaluate(grid);
        assert(near(s.penalty, e.penalty)); assert(near(s.ofr, e.ofr));
        // Python-style differentiable toy reference: overlap area is 13.02, capacity is 5.
        assert(near(s.penalty, (13.02 - 5.0) * (13.02 - 5.0), 1e-10)); assert(near(s.ofr, (13.02 - 5.0) / 13.02, 1e-10));
        assert(near(s.grad_x[0], 0)); assert(near(s.grad_y[0], 0));
    }
}
