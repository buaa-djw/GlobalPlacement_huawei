#pragma once

#include <random>
#include <vector>

class PlacementDB;

struct WirelengthCoordinateSnapshot
{
    std::vector<double> pin_x;
    std::vector<double> pin_y;
};

struct ExactWirelengthSubgradientResult
{
    // Equations (3)-(8), must equal exact HPWL.
    double hpwl = 0.0;

    // Indexed directly by cell_id.
    std::vector<double> grad_x;
    std::vector<double> grad_y;
};

class ExactWirelengthSubgradient
{
public:
    WirelengthCoordinateSnapshot capture(const PlacementDB& db) const;

    ExactWirelengthSubgradientResult evaluate(
        const PlacementDB& db,
        const WirelengthCoordinateSnapshot* previous,
        std::mt19937_64& random_engine
    ) const;
};
