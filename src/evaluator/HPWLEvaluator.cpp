#include "evaluator/HPWLEvaluator.h"

#include "db/Net.h"
#include "db/PlacementDB.h"
#include "geometry/Point.h"

#include <limits>

namespace {
double computeNetHPWL(const PlacementDB& db, const Net& net) {
    if (net.pin_ids.size() <= 1) {
        return 0.0;
    }

    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (const int pin_id : net.pin_ids) {
        const Point pos = db.pinPosition(pin_id);
        if (pos.x < min_x) min_x = pos.x;
        if (pos.x > max_x) max_x = pos.x;
        if (pos.y < min_y) min_y = pos.y;
        if (pos.y > max_y) max_y = pos.y;
    }

    return (max_x - min_x) + (max_y - min_y);
}
}

double HPWLEvaluator::netHPWL(const PlacementDB& db, int net_id) const {
    return computeNetHPWL(db, db.net(net_id));
}

double HPWLEvaluator::totalHPWL(const PlacementDB& db) const {
    double total = 0.0;
    for (const Net& net : db.nets()) {
        total += computeNetHPWL(db, net);
    }
    return total;
}

double HPWLEvaluator::updateAllNetHPWL(PlacementDB& db) const {
    double total = 0.0;

    for (const Net& net : db.nets()) {
        const double hpwl = computeNetHPWL(db, net);

        db.net(net.id).hpwl = hpwl;
        total += hpwl;
    }

    return total;
}