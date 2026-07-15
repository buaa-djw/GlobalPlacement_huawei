#pragma once

#include "wirelength/ExactWirelengthSubgradient.h"

#include <random>
#include <vector>

class ExactBinGrid;
class PlacementDB;

struct NonsmoothObjectiveResult
{
    // Equation (12).
    double objective = 0.0;

    // Equation (3), also exact HPWL.
    double hpwl = 0.0;

    // Density penalty in equation (12), before lambda.
    double density_penalty = 0.0;

    // Equation (18).
    double ofr = 0.0;

    // Subgradient of equation (12), indexed by cell_id.
    std::vector<double> grad_x;
    std::vector<double> grad_y;
};

class NonsmoothObjectiveEvaluator
{
public:
    NonsmoothObjectiveResult evaluate(
        const PlacementDB& db,
        const ExactBinGrid& grid,
        double lambda,
        const WirelengthCoordinateSnapshot* previous_wire_positions,
        std::mt19937_64& random_engine
    ) const;
};
