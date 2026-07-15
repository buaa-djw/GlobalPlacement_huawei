#include "objective/NonsmoothObjectiveEvaluator.h"

#include "db/Cell.h"
#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensitySubgradient.h"
#include "wirelength/ExactWirelengthSubgradient.h"

#include <cmath>
#include <stdexcept>

NonsmoothObjectiveResult NonsmoothObjectiveEvaluator::evaluate(
    const PlacementDB& db,
    const ExactBinGrid& grid,
    double lambda,
    const WirelengthCoordinateSnapshot* previous_wire_positions,
    std::mt19937_64& random_engine) const
{
    if (!std::isfinite(lambda) || lambda < 0.0) {
        throw std::invalid_argument("NonsmoothObjectiveEvaluator: lambda must be finite and nonnegative");
    }

    ExactWirelengthSubgradient wirelength;
    ExactWirelengthSubgradientResult wire =
        wirelength.evaluate(db, previous_wire_positions, random_engine);

    ExactBinGrid current_grid(db, grid.numBinsX(), grid.numBinsY(), grid.targetDensity());
    ExactDensitySubgradient density_subgradient;
    ExactDensitySubgradientResult density =
        density_subgradient.evaluate(db, current_grid, random_engine);

    if (wire.grad_x.size() != density.grad_x.size() ||
        wire.grad_y.size() != density.grad_y.size() ||
        wire.grad_x.size() != db.cells().size() ||
        wire.grad_y.size() != db.cells().size()) {
        throw std::runtime_error("NonsmoothObjectiveEvaluator: gradient size mismatch");
    }

    NonsmoothObjectiveResult result;
    result.hpwl = wire.hpwl;
    result.density_penalty = density.penalty;
    result.ofr = density.ofr;
    result.objective = result.hpwl + lambda * result.density_penalty;
    result.grad_x.assign(db.cells().size(), 0.0);
    result.grad_y.assign(db.cells().size(), 0.0);

    if (!std::isfinite(result.objective) || !std::isfinite(result.hpwl) ||
        !std::isfinite(result.density_penalty) || !std::isfinite(result.ofr)) {
        throw std::runtime_error("NonsmoothObjectiveEvaluator: non-finite objective component");
    }

    for (const Cell& cell : db.cells()) {
        const std::size_t id = static_cast<std::size_t>(cell.id);
        if (cell.isMovable()) {
            result.grad_x[id] = wire.grad_x[id] + lambda * density.grad_x[id];
            result.grad_y[id] = wire.grad_y[id] + lambda * density.grad_y[id];
        }
        if (!std::isfinite(result.grad_x[id]) || !std::isfinite(result.grad_y[id])) {
            throw std::runtime_error("NonsmoothObjectiveEvaluator: non-finite gradient");
        }
    }
    return result;
}
