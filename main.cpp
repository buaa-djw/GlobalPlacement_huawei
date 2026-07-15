#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "evaluator/HPWLEvaluator.h"
#include "parser/LimboBookshelfAdapter.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
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

        const HPWLEvaluator hpwl_evaluator;
        const double hpwl = hpwl_evaluator.totalHPWL(database);

        ExactBinGrid grid(
            database,
            options.bins_x,
            options.bins_y,
            options.target_density
        );

        const ExactDensityResult density =
            ExactDensityEvaluator().evaluate(grid);

        std::cout
            << std::fixed
            << std::setprecision(6)
            << "==================[Placement Result]===================\n"
            << "HPWL "
            << hpwl
            << '\n'
            << std::setprecision(12)
            << "OFR "
            << density.ofr
            << '\n';

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