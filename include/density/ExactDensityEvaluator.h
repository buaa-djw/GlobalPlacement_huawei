#pragma once

class ExactBinGrid;

struct ExactDensityResult
{
    // Equation (12), before multiplication by lambda.
    double penalty = 0.0;

    // Equation (18).
    double ofr = 0.0;
};

class ExactDensityEvaluator
{
public:
    ExactDensityResult evaluate(
        const ExactBinGrid& grid
    ) const;
};
