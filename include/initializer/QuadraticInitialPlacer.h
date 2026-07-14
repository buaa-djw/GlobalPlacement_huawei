#pragma once

#include "initializer/SimplePackingInitialPlacer.h"
#include "solver/PreconditionedConjugateGradient.h"

/** @brief Configuration for connection-aware quadratic initial placement. */
struct QuadraticInitialPlacerConfig {
    int maximum_net_degree = 1000;
    double anchor_weight_ratio = 1e-3;
    double minimum_anchor_weight = 1e-12;
    double blend_factor = 0.5;
    bool use_pin_offsets = true;
    bool project_to_core = true;
    bool require_solver_convergence = false;
    double numerical_epsilon = 1e-12;
    SimplePackingInitialPlacerConfig simple_reference;
    ConjugateGradientConfig solver;
};

/** @brief Builds and solves a star-model quadratic system for initial coordinates. */
class QuadraticInitialPlacer final : public InitialPlacer {
public:
    explicit QuadraticInitialPlacer(QuadraticInitialPlacerConfig config);
    InitialPlacementResult place(PlacementDB& db) const override;
private:
    QuadraticInitialPlacerConfig config_;
};
