#pragma once

#include <string>
#include <vector>
class CsrMatrix;

/** @brief Configuration for Jacobi-preconditioned conjugate gradient. */
struct ConjugateGradientConfig { int maximum_iterations = 500; double relative_tolerance = 1e-6; double absolute_tolerance = 1e-10; double numerical_epsilon = 1e-18; };

/** @brief Result and diagnostics from a PCG solve. */
struct ConjugateGradientResult { std::vector<double> solution; int iterations = 0; double initial_residual_norm = 0.0; double final_residual_norm = 0.0; double relative_residual = 0.0; bool converged = false; std::string termination_reason; };

/** @brief Sparse SPD linear solver using Jacobi preconditioning. */
class PreconditionedConjugateGradient {
public:
    explicit PreconditionedConjugateGradient(ConjugateGradientConfig config);
    ConjugateGradientResult solve(const CsrMatrix& matrix, const std::vector<double>& rhs, const std::vector<double>& initial_solution) const;
private:
    ConjugateGradientConfig config_;
};
