#pragma once

#include <cstddef>
#include <string>

class PlacementDB;

/** @brief Available deterministic initial placement algorithms. */
enum class InitialPlacementMethod {
    SimplePacking,
    QuadraticConnectivity
};

/** @brief Diagnostics returned by an initial placer after updating movable cells. */
struct InitialPlacementResult {
    std::string method;
    std::size_t movable_cell_count = 0;
    std::size_t active_net_count = 0;
    std::size_t skipped_high_degree_net_count = 0;
    int x_solver_iterations = 0;
    int y_solver_iterations = 0;
    double x_relative_residual = 0.0;
    double y_relative_residual = 0.0;
    bool x_converged = false;
    bool y_converged = false;
    bool used_fallback = false;
    std::string message;
};

/** @brief Base interface for modules that create movable-cell initial coordinates. */
class InitialPlacer {
public:
    virtual ~InitialPlacer() = default;
    virtual InitialPlacementResult place(PlacementDB& db) const = 0;
};
