#pragma once

#include <random>
#include <vector>

class ExactBinGrid;
class PlacementDB;

struct ExactDensitySubgradientResult
{
    // Equation (12), before multiplication by lambda.
    double penalty = 0.0;

    // Equation (18).
    double ofr = 0.0;

    // Indexed directly by cell_id.
    std::vector<double> grad_x;
    std::vector<double> grad_y;
};

class ExactDensitySubgradient
{
public:
    ExactDensitySubgradientResult evaluate(
        const PlacementDB& db,
        const ExactBinGrid& grid,
        std::mt19937_64& random_engine
    ) const;
};
