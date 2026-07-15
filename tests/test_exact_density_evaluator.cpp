#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "db/PlacementDB.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

static bool near(double a, double b, double eps = 1e-9) { return std::abs(a - b) <= eps * std::max(1.0, std::abs(b)); }
static void row(PlacementDB& db, double w = 20, double h = 10)
{
    Row r; r.y = 0; r.height = h; r.site_width = 1; r.site_spacing = 1; r.x_start = 0; r.num_sites = static_cast<int>(w); db.addRow(r);
}
static void cell(PlacementDB& db, const char* n, double x, double y, double w, double h, bool f = false, CellType t = CellType::Standard)
{
    db.addCell(n, w, h, t); db.setCellLocation(n, x, y, f);
}
int main()
{
    {
        PlacementDB db; row(db); cell(db, "c", 0, 0, 2, 2);
        ExactDensityResult r = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.9));
        assert(near(r.penalty, 0)); assert(near(r.ofr, 0));
    }
    {
        PlacementDB db; row(db); cell(db, "c", 0, 0, 20, 10);
        ExactDensityResult r = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.5));
        assert(near(r.penalty, 10000)); assert(near(r.ofr, 0.5));
    }
    {
        PlacementDB db; row(db); cell(db, "a", 0, 0, 10, 10); cell(db, "b", 10, 0, 10, 10);
        ExactDensityResult r = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 2, 1, 0.5));
        assert(near(r.penalty, 5000)); assert(near(r.ofr, 0.5));
    }
    {
        PlacementDB db; row(db);
        cell(db, "m", 0, 0, 10, 10);
        cell(db, "f", 0, 0, 10, 10, true);
        cell(db, "t", 100, 100, 10, 10, true, CellType::Terminal);
        cell(db, "o", 200, 200, 10, 10);
        ExactDensityResult r = ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.5));
        assert(near(r.penalty, 10000));
        assert(near(r.ofr, 100.0 / 400.0));
    }
    {
        PlacementDB db; row(db); bool threw = false;
        try { (void)ExactDensityEvaluator().evaluate(ExactBinGrid(db, 1, 1, 0.5)); } catch (const std::runtime_error&) { threw = true; }
        assert(threw);
    }
}
