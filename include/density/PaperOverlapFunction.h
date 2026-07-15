#pragma once

#include <random>

struct AxisOverlapSubgradient
{
    double overlap = 0.0;
    double subgradient = 0.0;
};

double paperAxisOverlapLength(
    double cell_center,
    double cell_size,
    double bin_center,
    double bin_size
);

AxisOverlapSubgradient paperAxisOverlapSubgradient(
    double cell_center,
    double cell_size,
    double bin_center,
    double bin_size,
    std::mt19937_64& random_engine
);
