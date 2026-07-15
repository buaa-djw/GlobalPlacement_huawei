#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "evaluator/HPWLEvaluator.h"
#include "parser/LimboBookshelfAdapter.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace
{

struct Options
{
    std::string aux;
    int bins_x = 64;
    int bins_y = 64;
    double target_density = 0.9;
};

void printUsage(const char* executable)
{
    std::cerr
        << "Usage: "
        << executable
        << " --aux <benchmark.aux>"
        << " [--bins <nx> <ny>]"
        << " [--target-density <value>]\n";
}

Options parseOptions(int argc, char** argv)
{
    Options options;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];

        if (argument == "--aux" &&
            index + 1 < argc) {
            options.aux = argv[++index];
        } else if (argument == "--bins" &&
                   index + 2 < argc) {
            options.bins_x =
                std::stoi(argv[++index]);

            options.bins_y =
                std::stoi(argv[++index]);
        } else if (
            argument == "--target-density" &&
            index + 1 < argc) {
            options.target_density =
                std::stod(argv[++index]);
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
        throw std::invalid_argument(
            "missing --aux"
        );
    }

    return options;
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

        std::size_t standard_count = 0;
        std::size_t terminal_count = 0;
        std::size_t terminal_ni_count = 0;

        std::size_t movable_count = 0;
        std::size_t non_movable_count = 0;
        std::size_t fixed_flag_count = 0;

        std::map<std::string, std::size_t>
            orientation_counts;

        for (const Cell& cell : database.cells()) {
            switch (cell.type) {
            case CellType::Standard:
                ++standard_count;
                break;

            case CellType::Terminal:
                ++terminal_count;
                break;

            case CellType::TerminalNI:
                ++terminal_ni_count;
                break;
            }

            if (cell.isMovable()) {
                ++movable_count;
            } else {
                ++non_movable_count;
            }

            if (cell.is_fixed) {
                ++fixed_flag_count;
            }

            ++orientation_counts[cell.orientation];
        }

        const Box core =
            database.coreBounds();

        const HPWLEvaluator hpwl_evaluator;

        const double input_hpwl =
            hpwl_evaluator.totalHPWL(database);

        std::cout
            << std::fixed
            << std::setprecision(6);

        std::cout
            << "==================[Database]===================\n"
            << "cells "
            << database.cells().size() << '\n'
            << "standard_cells "
            << standard_count << '\n'
            << "terminal_cells "
            << terminal_count << '\n'
            << "terminal_NI_cells "
            << terminal_ni_count << '\n'
            << "movable_cells "
            << movable_count << '\n'
            << "non_movable_cells "
            << non_movable_count << '\n'
            << "fixed_flag_cells "
            << fixed_flag_count << '\n'
            << "nets "
            << database.nets().size() << '\n'
            << "pins "
            << database.pins().size() << '\n'
            << "rows "
            << database.rows().size() << '\n'
            << "core_bounds "
            << core.lx << ' '
            << core.ly << ' '
            << core.ux << ' '
            << core.uy << "\n\n";

        std::cout
            << "==================[Input Placement]===================\n"
            << "input_hpwl "
            << input_hpwl
            << "\n\n";

        std::cout
            << "==================[Orientation Statistics]===================\n";

        bool has_non_n_orientation = false;

        for (const auto& [orientation, count] :
             orientation_counts) {
            std::cout
                << orientation << ' '
                << count << '\n';

            if (orientation != "N") {
                has_non_n_orientation = true;
            }
        }

        if (has_non_n_orientation) {
            std::cerr
                << "Warning: non-N orientations exist; "
                << "pin offsets are still interpreted using "
                << "N-orientation semantics.\n";
        }

        ExactBinGrid grid(
            database,
            options.bins_x,
            options.bins_y,
            options.target_density
        );

        std::cout
            << "\n==================[Core Diagnostics]===================\n"
            << "movable_cell_count "
            << grid.movableCellCount() << '\n'
            << "movable_inside_core_count "
            << grid.movableInsideCoreCount() << '\n'
            << "movable_outside_core_count "
            << grid.movableOutsideCoreCount() << '\n'
            << "non_movable_cell_count "
            << grid.nonMovableCellCount() << '\n'
            << "non_movable_inside_core_count "
            << grid.nonMovableInsideCoreCount() << '\n'
            << "non_movable_outside_core_count "
            << grid.nonMovableOutsideCoreCount() << '\n'
            << "raw_movable_area "
            << grid.rawMovableArea() << '\n'
            << "raw_non_movable_area "
            << grid.rawNonMovableArea() << '\n'
            << "clipped_movable_area "
            << grid.totalMovableArea() << '\n'
            << "clipped_non_movable_area "
            << grid.totalNonMovableArea() << '\n';

        if (grid.movableCellCount() > 0 &&
            grid.totalMovableArea() <= 0.0) {
            std::cerr
                << "Warning: the database contains movable cells, "
                << "but none has positive overlap with the "
                << "placement core. Density OFR is undefined for "
                << "this input placement. An initial placement is "
                << "required before density optimization.\n";
        }

        const ExactDensityMetrics metrics =
            ExactDensityEvaluator().evaluate(grid);

        std::cout
            << "\n==================[Exact Density]===================\n"
            << "bins_x "
            << grid.numBinsX() << '\n'
            << "bins_y "
            << grid.numBinsY() << '\n'
            << "bin_width "
            << grid.binWidth() << '\n'
            << "bin_height "
            << grid.binHeight() << '\n'
            << "target_density "
            << grid.targetDensity() << '\n'
            << "total_cell_area "
            << metrics.total_cell_area << '\n'
            << "total_movable_area "
            << metrics.total_movable_area << '\n'
            << "total_non_movable_area "
            << metrics.total_non_movable_area << '\n'
            << "total_capacity "
            << metrics.total_capacity << '\n'
            << "total_overflow "
            << metrics.total_overflow << '\n';

        std::cout << "overflow_ratio ";

        if (metrics.overflow_ratio_defined) {
            std::cout
                << metrics.overflow_ratio
                << '\n';
        } else {
            std::cout << "N/A\n";
        }

        std::cout
            << "quadratic_penalty "
            << metrics.quadratic_penalty << '\n'
            << "overflow_bin_count "
            << metrics.overflow_bin_count << '\n'
            << "maximum_bin_density "
            << metrics.maximum_bin_density << '\n'
            << "average_bin_density "
            << metrics.average_bin_density << '\n';

        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        printUsage(argv[0]);
        return 1;
    }
}