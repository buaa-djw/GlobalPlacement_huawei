#include "db/PlacementDB.h"
#include "evaluator/HPWLEvaluator.h"
#include "geometry/Point.h"

#include <cassert>
#include <cmath>
#include <string>

namespace
{
    bool near(double actual, double expected)
    {
        return std::abs(actual - expected) < 1e-9;
    }

    void expectNear(double actual, double expected)
    {
        assert(near(actual, expected));
    }

    int addPlacedCell(
        PlacementDB &db,
        const std::string &name,
        double width,
        double height,
        double x,
        double y,
        CellType type = CellType::Standard,
        bool fixed = false)
    {
        const int cell_id = db.addCell(name, width, height, type);
        db.setCellLocation(name, x, y, fixed);
        return cell_id;
    }
}

int main()
{
    HPWLEvaluator evaluator;

    {
        PlacementDB db;
        const int cell = addPlacedCell(db, "single", 10.0, 10.0, 0.0, 0.0);
        const int net = db.addNet("single_net");
        db.addPin(cell, net, 0.0, 0.0, "I");
        expectNear(evaluator.netHPWL(db, net), 0.0);
    }

    {
        PlacementDB db;
        const int a = addPlacedCell(db, "a", 10.0, 10.0, 0.0, 0.0);
        const int b = addPlacedCell(db, "b", 10.0, 10.0, 20.0, 30.0);
        const int net = db.addNet("two_pin_net");
        db.addPin(a, net, 0.0, 0.0, "I");
        db.addPin(b, net, 0.0, 0.0, "O");
        expectNear(evaluator.netHPWL(db, net), 50.0);
    }

    {
        PlacementDB db;
        const int a = addPlacedCell(db, "a", 10.0, 10.0, 0.0, 0.0);
        const int b = addPlacedCell(db, "b", 8.0, 4.0, 10.0, 20.0);
        const int net = db.addNet("offset_net");
        db.addPin(a, net, -2.0, 3.0, "I"); // (3, 8)
        db.addPin(b, net, 6.0, -5.0, "O"); // (20, 17)
        expectNear(evaluator.netHPWL(db, net), 26.0);
    }

    {
        PlacementDB db;
        const int c0 = addPlacedCell(db, "c0", 2.0, 2.0, 0.0, 0.0);
        const int c1 = addPlacedCell(db, "c1", 2.0, 2.0, 10.0, 5.0);
        const int c2 = addPlacedCell(db, "c2", 2.0, 2.0, -4.0, 20.0);
        const int c3 = addPlacedCell(db, "c3", 2.0, 2.0, 30.0, -10.0);
        const int net = db.addNet("multi_pin_net");
        db.addPin(c0, net, 0.0, 0.0, "I");
        db.addPin(c1, net, 0.0, 0.0, "I");
        db.addPin(c2, net, 0.0, 0.0, "I");
        db.addPin(c3, net, 0.0, 0.0, "O");
        expectNear(evaluator.netHPWL(db, net), 64.0);
    }

    {
        PlacementDB db;
        const int a = addPlacedCell(db, "a", 10.0, 10.0, 0.0, 0.0);
        const int b = addPlacedCell(db, "b", 10.0, 10.0, 20.0, 30.0);
        const int c = addPlacedCell(db, "c", 4.0, 4.0, -10.0, 0.0);
        const int d = addPlacedCell(db, "d", 4.0, 4.0, -5.0, 5.0);
        const int n0 = db.addNet("n0");
        db.addPin(a, n0, 0.0, 0.0, "I");
        db.addPin(b, n0, 0.0, 0.0, "O");
        const int n1 = db.addNet("n1");
        db.addPin(c, n1, 0.0, 0.0, "I");
        db.addPin(d, n1, 0.0, 0.0, "O");
        const int n2 = db.addNet("n2");
        db.addPin(a, n2, 1.0, 1.0, "I");
        expectNear(evaluator.netHPWL(db, n0), 50.0);
        expectNear(evaluator.netHPWL(db, n1), 10.0);
        expectNear(evaluator.netHPWL(db, n2), 0.0);
        expectNear(evaluator.totalHPWL(db), 60.0);
    }

    {
        PlacementDB db;
        const int movable = addPlacedCell(db, "movable", 10.0, 10.0, 0.0, 0.0);
        const int fixed = addPlacedCell(db, "fixed", 10.0, 10.0, 100.0, 50.0, CellType::Standard, true);
        const int net = db.addNet("fixed_net");
        db.addPin(movable, net, 0.0, 0.0, "I");
        db.addPin(fixed, net, 0.0, 0.0, "O");
        expectNear(evaluator.netHPWL(db, net), 150.0);
    }

    {
        PlacementDB db;
        const int cell = addPlacedCell(db, "pin_position", 10.0, 20.0, 3.0, 4.0);
        const int net = db.addNet("pin_position_net");
        const int pin = db.addPin(cell, net, -1.5, 2.25, "I");
        const Point pos = db.pinPosition(pin);
        expectNear(pos.x, 6.5);
        expectNear(pos.y, 16.25);
    }

    return 0;
}
