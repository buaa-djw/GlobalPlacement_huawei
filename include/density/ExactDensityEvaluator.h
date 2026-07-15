#pragma once

#include <cstddef>

class ExactBinGrid;

struct ExactDensityMetrics
{
    std::size_t bin_count = 0;
    std::size_t overflow_bin_count = 0;

    double total_cell_area = 0.0;
    double total_movable_area = 0.0;
    double total_non_movable_area = 0.0;

    double total_capacity = 0.0;
    double total_overflow = 0.0;
    double overflow_ratio = 0.0;
    bool overflow_ratio_defined = false;

    double quadratic_penalty = 0.0;
    double maximum_bin_density = 0.0;
    double average_bin_density = 0.0;
};

class ExactDensityEvaluator
{
public:
    ExactDensityMetrics evaluate(
        const ExactBinGrid& grid
    ) const;
};