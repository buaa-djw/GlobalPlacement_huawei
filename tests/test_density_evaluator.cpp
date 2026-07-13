#include "db/PlacementDB.h"
#include "evaluator/DensityEvaluator.h"
#include "grid/BinGrid.h"
#include "utils/Logger.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
constexpr double EPS = 1e-9;

void expectNearImpl(double actual, double expected, int line) { if (std::abs(actual - expected) >= EPS) { std::cerr << "expectNear failed at line " << line << ": actual=" << actual << " expected=" << expected << "\n"; assert(false); } }
#define expectNear(actual, expected) expectNearImpl((actual), (expected), __LINE__)

void addCoreRow(PlacementDB& db, double lx = 0.0, double ly = 0.0, double ux = 100.0, double h = 100.0) {
    Row row;
    row.x_start = lx;
    row.x_end = ux;
    row.y = ly;
    row.height = h;
    row.site_width = 1.0;
    row.site_spacing = 1.0;
    row.num_sites = static_cast<int>(ux - lx);
    db.addRow(row);
}

void addCell(
    PlacementDB& db,
    const std::string& name,
    double x,
    double y,
    double width,
    double height,
    CellType type = CellType::Standard,
    bool fixed = false
) {
    db.addCell(name, width, height, type);
    db.setCellLocation(name, x, y, fixed);
}

DensityMetrics eval(const BinGrid& grid) { return DensityEvaluator{}.evaluate(grid); }

void expectThrowsInvalid(const BinGrid& grid) {
    bool thrown = false;
    try {
        (void)eval(grid);
    } catch (const std::invalid_argument&) {
        thrown = true;
    }
    assert(thrown);
}

bool fileContains(const std::string& path, const std::string& needle) {
    std::ifstream in(path);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return text.find(needle) != std::string::npos;
}

int runLoggedDensityFailure(const std::string& log_path) {
    try {
        Logger::instance().initialize(log_path, LogLevel::Debug, false);
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "m", 0, 0, 1, 1);
        BinGrid grid(db, 1, 1, 1.0);
        auto& bins = const_cast<std::vector<Bin>&>(grid.bins());
        bins[0].movable_area = -1.0;
        (void)eval(grid);
        Logger::instance().flush();
        return 0;
    } catch (const std::exception& e) {
        LOG_ERROR("Caught density evaluator exception: " << e.what());
        Logger::instance().flush();
        Logger::instance().shutdown();
        return 2;
    }
}
}

int main() {
    // 1: no overflow.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "m", 0, 0, 5, 5);
        BinGrid grid(db, 1, 1, 1.0);
        const DensityMetrics m = eval(grid);
        expectNear(m.penalty, 0.0);
        expectNear(m.total_overflow, 0.0);
        expectNear(m.overflow_ratio, 0.0);
        assert(m.overflow_bin_count == 0);
    }
    // 2: one overflowing bin: 80 - 50 = 30, penalty 900.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "m", 0, 0, 8, 10);
        BinGrid grid(db, 1, 1, 0.5);
        const DensityMetrics m = eval(grid);
        expectNear(m.total_overflow, 30.0);
        expectNear(m.penalty, 900.0);
        assert(m.overflow_bin_count == 1);
    }
    // 3 and 4: two bins overflow by 10 and 20; total movable area is 120.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 20, 10);
        addCell(db, "m0", 0, 0, 6, 10);
        addCell(db, "m1", 10, 0, 6, 10);
        BinGrid grid(db, 2, 1, 0.5);
        const DensityMetrics m = eval(grid);
        expectNear(m.total_overflow, 20.0);
        expectNear(m.penalty, 200.0);
        expectNear(m.overflow_ratio, 20.0 / 120.0);
        assert(m.overflow_bin_count == 2);
    }
    // Explicit multi-bin 10 and 20 overflow case.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 20, 10);
        addCell(db, "m0", 0, 0, 6, 10);
        addCell(db, "m1", 10, 0, 7, 10);
        BinGrid grid(db, 2, 1, 0.5);
        const DensityMetrics m = eval(grid);
        expectNear(m.total_overflow, 30.0);
        expectNear(m.penalty, 500.0);
        expectNear(m.overflow_ratio, 30.0 / 130.0);
    }
    // 5: no movable cell area returns a zero overflow ratio and does not crash.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        BinGrid grid(db, 1, 1, 1.0);
        const DensityMetrics m = eval(grid);
        expectNear(m.total_movable_area, 0.0);
        expectNear(m.overflow_ratio, 0.0);
    }
    // 6: average utilization only over positive-capacity bins.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 20, 10);
        addCell(db, "m0", 0, 0, 5, 10);   // 50 / 100 = 0.5
        addCell(db, "m1", 10, 0, 10, 10); // 100 / 100 = 1.0
        BinGrid grid(db, 2, 1, 1.0);
        const DensityMetrics m = eval(grid);
        expectNear(m.average_utilization, 0.75);
        assert(m.valid_utilization_bin_count == 2);
    }
    // 7: zero-capacity empty bin is excluded from utilization average.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "f", 0, 0, 10, 10, CellType::Standard, true);
        BinGrid grid(db, 1, 1, 1.0);
        const DensityMetrics m = eval(grid);
        assert(m.zero_capacity_bin_count == 1);
        assert(m.zero_capacity_occupied_bin_count == 0);
        assert(m.valid_utilization_bin_count == 0);
        expectNear(m.average_utilization, 0.0);
    }
    // 8: zero-capacity occupied bin overflows and may have infinite utilization.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "f", 0, 0, 10, 10, CellType::Standard, true);
        addCell(db, "m", 0, 0, 1, 10);
        BinGrid grid(db, 1, 1, 1.0);
        const DensityMetrics m = eval(grid);
        expectNear(m.total_overflow, 10.0);
        expectNear(m.penalty, 100.0);
        assert(m.zero_capacity_occupied_bin_count == 1);
        assert(std::isinf(m.max_utilization));
    }
    // 9: fixed cells reduce movable capacity and the evaluator uses that capacity.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "f", 0, 0, 5, 10, CellType::Standard, true);
        addCell(db, "m", 0, 0, 8, 10);
        BinGrid grid(db, 1, 1, 1.0);
        const DensityMetrics m = eval(grid);
        expectNear(m.total_fixed_area, 50.0);
        expectNear(m.total_movable_capacity, 50.0);
        expectNear(m.total_overflow, 30.0);
        expectNear(m.penalty, 900.0);
    }
    // 10: rebuild after movement changes the density metrics.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 20, 10);
        addCell(db, "m", 0, 0, 5, 10);
        BinGrid grid(db, 2, 1, 1.0);
        expectNear(eval(grid).total_overflow, 0.0);
        db.setCellLocation("m", 0, 0, false);
        db.cell(db.getCellId("m")).width = 15.0;
        grid.rebuild(db);
        expectNear(eval(grid).total_overflow, 0.0);
        BinGrid tighter(db, 2, 1, 0.5);
        expectNear(eval(tighter).total_overflow, 50.0);
    }
    // 11 and 12: invalid negative/NaN states are rejected; +inf utilization is only allowed for zero-cap occupied bins.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "m", 0, 0, 1, 1);
        BinGrid grid(db, 1, 1, 1.0);
        auto& b = const_cast<std::vector<Bin>&>(grid.bins())[0];
        b.movable_area = -1.0;
        expectThrowsInvalid(grid);
        b.movable_area = 1.0;
        b.movable_capacity = std::numeric_limits<double>::quiet_NaN();
        expectThrowsInvalid(grid);
        b.movable_capacity = 1.0;
        b.utilization = std::numeric_limits<double>::infinity();
        expectThrowsInvalid(grid);
    }
    // 13: repeat calls are stable.
    {
        PlacementDB db;
        addCoreRow(db, 0, 0, 10, 10);
        addCell(db, "m", 0, 0, 8, 10);
        BinGrid grid(db, 1, 1, 0.5);
        const DensityMetrics a = eval(grid);
        const DensityMetrics b = eval(grid);
        expectNear(a.penalty, b.penalty);
        expectNear(a.total_overflow, b.total_overflow);
        expectNear(a.overflow_ratio, b.overflow_ratio);
    }
    // 14: BinGrid rejects empty/invalid grids before DensityEvaluator can accept one.
    {
        PlacementDB db;
        bool thrown = false;
        try {
            BinGrid grid(db, 1, 1, 1.0);
            (void)grid;
        } catch (const std::exception&) {
            thrown = true;
        }
        assert(thrown);
    }
    // Logging integration for caught density exceptions.
    {
        const std::string log_path = "result/test_density_evaluator_error.log";
        const int rc = runLoggedDensityFailure(log_path);
        assert(rc != 0);
        assert(fileContains(log_path, "ERROR"));
        assert(fileContains(log_path, "DensityEvaluator: bin 0 has invalid movable area"));
    }
    return 0;
}
