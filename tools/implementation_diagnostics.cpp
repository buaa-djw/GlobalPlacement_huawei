#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "density/ExactDensitySubgradient.h"
#include "evaluator/HPWLEvaluator.h"
#include "objective/NonsmoothObjectiveEvaluator.h"
#include "parser/LimboBookshelfAdapter.h"
#include "wirelength/ExactWirelengthSubgradient.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double kDiagnosticZeroTolerance = 1e-12;
constexpr double kRelativeTolerance = 1e-9;

struct Options
{
    std::string aux;

    int bins_x = 64;
    int bins_y = 64;

    double target_density = 0.9;

    // This is only a diagnostic lambda.
    // It is not the paper's automatically initialized lambda.
    double lambda = 1.0;

    std::uint64_t seed = 1;

    std::size_t sample_gradients = 10;
};

struct PositionSnapshot
{
    std::vector<double> x;
    std::vector<double> y;
};

struct GradientSummary
{
    std::size_t movable_cell_count = 0;
    std::size_t nonzero_movable_cell_count = 0;

    double movable_l1_norm = 0.0;
    double movable_l2_norm = 0.0;
    double movable_linf_norm = 0.0;

    double non_movable_linf_norm = 0.0;
};

void printUsage(const char* executable)
{
    std::cerr
        << "Usage: " << executable
        << " --aux <benchmark.aux>"
        << " [--bins <nx> <ny>]"
        << " [--target-density <value>]"
        << " [--lambda <value>]"
        << " [--seed <integer>]"
        << " [--samples <count>]\n";
}

Options parseOptions(int argc, char** argv)
{
    Options options;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];

        if (argument == "--aux" && index + 1 < argc) {
            options.aux = argv[++index];
        } else if (
            argument == "--bins" &&
            index + 2 < argc) {
            options.bins_x = std::stoi(argv[++index]);
            options.bins_y = std::stoi(argv[++index]);
        } else if (
            argument == "--target-density" &&
            index + 1 < argc) {
            options.target_density =
                std::stod(argv[++index]);
        } else if (
            argument == "--lambda" &&
            index + 1 < argc) {
            options.lambda =
                std::stod(argv[++index]);
        } else if (
            argument == "--seed" &&
            index + 1 < argc) {
            options.seed =
                static_cast<std::uint64_t>(
                    std::stoull(argv[++index])
                );
        } else if (
            argument == "--samples" &&
            index + 1 < argc) {
            options.sample_gradients =
                static_cast<std::size_t>(
                    std::stoull(argv[++index])
                );
        } else if (
            argument == "-h" ||
            argument == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument(
                "unknown or incomplete argument: " +
                argument
            );
        }
    }

    if (options.aux.empty()) {
        throw std::invalid_argument("missing --aux");
    }

    if (options.bins_x <= 0 ||
        options.bins_y <= 0) {
        throw std::invalid_argument(
            "bin counts must be positive"
        );
    }

    if (!std::isfinite(options.target_density) ||
        options.target_density <= 0.0 ||
        options.target_density > 1.0) {
        throw std::invalid_argument(
            "target density must satisfy 0 < value <= 1"
        );
    }

    if (!std::isfinite(options.lambda) ||
        options.lambda < 0.0) {
        throw std::invalid_argument(
            "lambda must be finite and nonnegative"
        );
    }

    return options;
}

PositionSnapshot capturePositions(
    const PlacementDB& database
)
{
    PositionSnapshot snapshot;

    snapshot.x.reserve(database.cells().size());
    snapshot.y.reserve(database.cells().size());

    for (const Cell& cell : database.cells()) {
        snapshot.x.push_back(cell.x);
        snapshot.y.push_back(cell.y);
    }

    return snapshot;
}

bool positionsUnchanged(
    const PlacementDB& database,
    const PositionSnapshot& snapshot
)
{
    if (snapshot.x.size() != database.cells().size() ||
        snapshot.y.size() != database.cells().size()) {
        return false;
    }

    for (const Cell& cell : database.cells()) {
        const std::size_t id =
            static_cast<std::size_t>(cell.id);

        if (cell.x != snapshot.x[id] ||
            cell.y != snapshot.y[id]) {
            return false;
        }
    }

    return true;
}

double maxAbsDifference(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs
)
{
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error(
            "diagnostics: vector size mismatch"
        );
    }

    double maximum = 0.0;

    for (std::size_t index = 0;
         index < lhs.size();
         ++index) {
        maximum = std::max(
            maximum,
            std::abs(lhs[index] - rhs[index])
        );
    }

    return maximum;
}

double relativeTolerance(double reference)
{
    return kRelativeTolerance *
           std::max(1.0, std::abs(reference));
}

GradientSummary summarizeGradient(
    const PlacementDB& database,
    const std::vector<double>& grad_x,
    const std::vector<double>& grad_y
)
{
    if (grad_x.size() != database.cells().size() ||
        grad_y.size() != database.cells().size()) {
        throw std::runtime_error(
            "diagnostics: gradient size mismatch"
        );
    }

    GradientSummary summary;
    double l2_squared = 0.0;

    for (const Cell& cell : database.cells()) {
        const std::size_t id =
            static_cast<std::size_t>(cell.id);

        const double gx = grad_x[id];
        const double gy = grad_y[id];

        if (!std::isfinite(gx) ||
            !std::isfinite(gy)) {
            throw std::runtime_error(
                "diagnostics: non-finite gradient"
            );
        }

        const double local_linf =
            std::max(std::abs(gx), std::abs(gy));

        if (cell.isMovable()) {
            ++summary.movable_cell_count;

            summary.movable_l1_norm +=
                std::abs(gx) + std::abs(gy);

            l2_squared += gx * gx + gy * gy;

            summary.movable_linf_norm =
                std::max(
                    summary.movable_linf_norm,
                    local_linf
                );

            if (local_linf >
                kDiagnosticZeroTolerance) {
                ++summary.nonzero_movable_cell_count;
            }
        } else {
            summary.non_movable_linf_norm =
                std::max(
                    summary.non_movable_linf_norm,
                    local_linf
                );
        }
    }

    summary.movable_l2_norm =
        std::sqrt(l2_squared);

    return summary;
}

void printGradientSummary(
    const std::string& name,
    const GradientSummary& summary
)
{
    std::cout
        << '[' << name << " Gradient]\n"
        << "movable_cell_count "
        << summary.movable_cell_count << '\n'
        << "nonzero_movable_cell_count "
        << summary.nonzero_movable_cell_count << '\n'
        << "movable_l1_norm "
        << summary.movable_l1_norm << '\n'
        << "movable_l2_norm "
        << summary.movable_l2_norm << '\n'
        << "movable_linf_norm "
        << summary.movable_linf_norm << '\n'
        << "non_movable_linf_norm "
        << summary.non_movable_linf_norm
        << "\n\n";
}

double gradientRecompositionError(
    const ExactWirelengthSubgradientResult& wire,
    const ExactDensitySubgradientResult& density,
    const NonsmoothObjectiveResult& objective,
    double lambda
)
{
    if (wire.grad_x.size() != density.grad_x.size() ||
        wire.grad_y.size() != density.grad_y.size() ||
        objective.grad_x.size() != wire.grad_x.size() ||
        objective.grad_y.size() != wire.grad_y.size()) {
        throw std::runtime_error(
            "diagnostics: gradient recomposition size mismatch"
        );
    }

    double maximum_error = 0.0;

    for (std::size_t id = 0;
         id < wire.grad_x.size();
         ++id) {
        const double expected_x =
            wire.grad_x[id] +
            lambda * density.grad_x[id];

        const double expected_y =
            wire.grad_y[id] +
            lambda * density.grad_y[id];

        maximum_error = std::max(
            maximum_error,
            std::abs(
                objective.grad_x[id] -
                expected_x
            )
        );

        maximum_error = std::max(
            maximum_error,
            std::abs(
                objective.grad_y[id] -
                expected_y
            )
        );
    }

    return maximum_error;
}

void printSampleGradients(
    const PlacementDB& database,
    const ExactWirelengthSubgradientResult& wire,
    const ExactDensitySubgradientResult& density,
    const NonsmoothObjectiveResult& objective,
    std::size_t maximum_samples
)
{
    std::cout << "[Sample Movable Gradients]\n";

    std::size_t printed = 0;

    for (const Cell& cell : database.cells()) {
        if (!cell.isMovable()) {
            continue;
        }

        const std::size_t id =
            static_cast<std::size_t>(cell.id);

        const double largest_component = std::max({
            std::abs(wire.grad_x[id]),
            std::abs(wire.grad_y[id]),
            std::abs(density.grad_x[id]),
            std::abs(density.grad_y[id]),
            std::abs(objective.grad_x[id]),
            std::abs(objective.grad_y[id])
        });

        if (largest_component <=
            kDiagnosticZeroTolerance) {
            continue;
        }

        std::cout
            << "cell_id " << cell.id
            << " name " << cell.name
            << " wire_g=("
            << wire.grad_x[id] << ','
            << wire.grad_y[id] << ')'
            << " density_g=("
            << density.grad_x[id] << ','
            << density.grad_y[id] << ')'
            << " objective_g=("
            << objective.grad_x[id] << ','
            << objective.grad_y[id] << ')'
            << '\n';

        ++printed;

        if (printed >= maximum_samples) {
            break;
        }
    }

    if (printed == 0) {
        std::cout << "none\n";
    }

    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Options options =
            parseOptions(argc, argv);

        PlacementDB database;
        LimboBookshelfAdapter parser(database);

        if (!parser.read(options.aux)) {
            throw std::runtime_error(
                "failed to read Bookshelf benchmark"
            );
        }

        const PositionSnapshot positions_before =
            capturePositions(database);

        const HPWLEvaluator hpwl_evaluator;
        const double reference_hpwl =
            hpwl_evaluator.totalHPWL(database);

        ExactBinGrid grid(
            database,
            options.bins_x,
            options.bins_y,
            options.target_density
        );

        const ExactDensityResult density_value =
            ExactDensityEvaluator().evaluate(grid);

        /*
         * The component RNG and objective RNG use the same seed.
         *
         * Components are called in the same order used by
         * NonsmoothObjectiveEvaluator:
         *
         *   1. wirelength
         *   2. density
         *
         * This makes gradient recomposition directly comparable.
         */
        std::mt19937_64 component_rng(options.seed);

        const ExactWirelengthSubgradientResult wire =
            ExactWirelengthSubgradient().evaluate(
                database,
                nullptr,
                component_rng
            );

        const ExactDensitySubgradientResult density =
            ExactDensitySubgradient().evaluate(
                database,
                grid,
                component_rng
            );

        std::mt19937_64 objective_rng(options.seed);

        const NonsmoothObjectiveResult objective =
            NonsmoothObjectiveEvaluator().evaluate(
                database,
                grid,
                options.lambda,
                nullptr,
                objective_rng
            );

        /*
         * Repeat with the same seed to check reproducibility.
         */
        std::mt19937_64 repeat_rng(options.seed);

        const NonsmoothObjectiveResult repeated_objective =
            NonsmoothObjectiveEvaluator().evaluate(
                database,
                grid,
                options.lambda,
                nullptr,
                repeat_rng
            );

        const double hpwl_error =
            std::abs(wire.hpwl - reference_hpwl);

        const double density_penalty_error =
            std::abs(
                density.penalty -
                density_value.penalty
            );

        const double density_ofr_error =
            std::abs(
                density.ofr -
                density_value.ofr
            );

        const double expected_objective =
            wire.hpwl +
            options.lambda * density.penalty;

        const double objective_value_error =
            std::abs(
                objective.objective -
                expected_objective
            );

        const double objective_hpwl_error =
            std::abs(
                objective.hpwl -
                wire.hpwl
            );

        const double objective_penalty_error =
            std::abs(
                objective.density_penalty -
                density.penalty
            );

        const double objective_ofr_error =
            std::abs(
                objective.ofr -
                density.ofr
            );

        const double gradient_recomposition_error =
            gradientRecompositionError(
                wire,
                density,
                objective,
                options.lambda
            );

        double reproducibility_error = 0.0;

        reproducibility_error = std::max(
            reproducibility_error,
            std::abs(
                repeated_objective.objective -
                objective.objective
            )
        );

        reproducibility_error = std::max(
            reproducibility_error,
            std::abs(
                repeated_objective.hpwl -
                objective.hpwl
            )
        );

        reproducibility_error = std::max(
            reproducibility_error,
            std::abs(
                repeated_objective.density_penalty -
                objective.density_penalty
            )
        );

        reproducibility_error = std::max(
            reproducibility_error,
            std::abs(
                repeated_objective.ofr -
                objective.ofr
            )
        );

        reproducibility_error = std::max(
            reproducibility_error,
            maxAbsDifference(
                repeated_objective.grad_x,
                objective.grad_x
            )
        );

        reproducibility_error = std::max(
            reproducibility_error,
            maxAbsDifference(
                repeated_objective.grad_y,
                objective.grad_y
            )
        );

        const GradientSummary wire_summary =
            summarizeGradient(
                database,
                wire.grad_x,
                wire.grad_y
            );

        const GradientSummary density_summary =
            summarizeGradient(
                database,
                density.grad_x,
                density.grad_y
            );

        const GradientSummary objective_summary =
            summarizeGradient(
                database,
                objective.grad_x,
                objective.grad_y
            );

        const bool coordinates_unchanged =
            positionsUnchanged(
                database,
                positions_before
            );

        const double objective_gradient_tolerance =
            kRelativeTolerance *
            std::max(
                1.0,
                objective_summary.movable_linf_norm
            );

        bool passed = true;

        passed =
            passed &&
            hpwl_error <=
                relativeTolerance(reference_hpwl);

        passed =
            passed &&
            density_penalty_error <=
                relativeTolerance(density_value.penalty);

        passed =
            passed &&
            density_ofr_error <=
                relativeTolerance(density_value.ofr);

        passed =
            passed &&
            objective_value_error <=
                relativeTolerance(expected_objective);

        passed =
            passed &&
            objective_hpwl_error <=
                relativeTolerance(wire.hpwl);

        passed =
            passed &&
            objective_penalty_error <=
                relativeTolerance(density.penalty);

        passed =
            passed &&
            objective_ofr_error <=
                relativeTolerance(density.ofr);

        passed =
            passed &&
            gradient_recomposition_error <=
                objective_gradient_tolerance;

        passed =
            passed &&
            wire_summary.non_movable_linf_norm <=
                kDiagnosticZeroTolerance;

        passed =
            passed &&
            density_summary.non_movable_linf_norm <=
                kDiagnosticZeroTolerance;

        passed =
            passed &&
            objective_summary.non_movable_linf_norm <=
                kDiagnosticZeroTolerance;

        passed =
            passed &&
            reproducibility_error <=
                objective_gradient_tolerance;

        passed =
            passed &&
            coordinates_unchanged;

        std::cout
            << std::setprecision(12)
            << std::scientific;

        std::cout
            << "==================[Implementation Diagnostics]===================\n"
            << "[Configuration]\n"
            << "bins_x " << options.bins_x << '\n'
            << "bins_y " << options.bins_y << '\n'
            << "target_density "
            << options.target_density << '\n'
            << "diagnostic_lambda "
            << options.lambda << '\n'
            << "random_seed "
            << options.seed << "\n\n";

        std::cout
            << "[Value Consistency]\n"
            << "reference_hpwl "
            << reference_hpwl << '\n'
            << "weighted_l1_hpwl "
            << wire.hpwl << '\n'
            << "hpwl_abs_error "
            << hpwl_error << '\n'
            << "density_evaluator_penalty "
            << density_value.penalty << '\n'
            << "density_subgradient_penalty "
            << density.penalty << '\n'
            << "density_penalty_abs_error "
            << density_penalty_error << '\n'
            << "density_evaluator_ofr "
            << density_value.ofr << '\n'
            << "density_subgradient_ofr "
            << density.ofr << '\n'
            << "density_ofr_abs_error "
            << density_ofr_error << '\n'
            << "internal_objective "
            << objective.objective << '\n'
            << "recomposed_objective "
            << expected_objective << '\n'
            << "objective_abs_error "
            << objective_value_error << "\n\n";

        std::cout
            << "[Objective Component Consistency]\n"
            << "objective_hpwl_abs_error "
            << objective_hpwl_error << '\n'
            << "objective_penalty_abs_error "
            << objective_penalty_error << '\n'
            << "objective_ofr_abs_error "
            << objective_ofr_error << '\n'
            << "gradient_recomposition_linf_error "
            << gradient_recomposition_error
            << "\n\n";

        printGradientSummary(
            "Wirelength",
            wire_summary
        );

        printGradientSummary(
            "Density",
            density_summary
        );

        printGradientSummary(
            "Objective",
            objective_summary
        );

        std::cout
            << "[Invariants]\n"
            << "coordinates_unchanged "
            << (coordinates_unchanged ? "true" : "false")
            << '\n'
            << "same_seed_reproducibility_error "
            << reproducibility_error << '\n'
            << "overall_status "
            << (passed ? "PASS" : "FAIL")
            << "\n\n";

        if (density_summary.nonzero_movable_cell_count == 0) {
            std::cout
                << "[Diagnostic Note]\n"
                << "The density subgradient is zero for every movable cell "
                << "at the current placement.\n"
                << "For the original adaptec2 placement this is expected "
                << "when movable cells have no effective overlap with the "
                << "placement bins.\n\n";
        }

        printSampleGradients(
            database,
            wire,
            density,
            objective,
            options.sample_gradients
        );

        return passed ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        printUsage(argv[0]);
        return 1;
    }
}