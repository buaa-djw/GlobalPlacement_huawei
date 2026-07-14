#pragma once

#include "initializer/InitialPlacer.h"

/** @brief Configuration for deterministic row-interval packing initialization. */
struct SimplePackingInitialPlacerConfig {
    double numerical_epsilon = 1e-12;
};

/** @brief Places movable cells into row free intervals without evaluating objectives. */
class SimplePackingInitialPlacer final : public InitialPlacer {
public:
    explicit SimplePackingInitialPlacer(SimplePackingInitialPlacerConfig config = {});
    InitialPlacementResult place(PlacementDB& db) const override;
private:
    SimplePackingInitialPlacerConfig config_;
};
