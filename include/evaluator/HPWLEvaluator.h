#pragma once

class PlacementDB;

class HPWLEvaluator {
public:
    double netHPWL(const PlacementDB& db, int net_id) const;
    double totalHPWL(const PlacementDB& db) const;
};
