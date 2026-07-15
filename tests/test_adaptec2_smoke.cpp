#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "evaluator/HPWLEvaluator.h"
#include "parser/LimboBookshelfAdapter.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

int main()
{
    const char* root_env = std::getenv("BENCHMARK_ROOT");
    const std::filesystem::path root = root_env ? root_env : "testbench/ispd2005";
    const std::filesystem::path aux = root / "adaptec2" / "adaptec2.aux";
    if (!std::filesystem::exists(aux)) {
        std::cerr << "Skipping adaptec2 smoke check: benchmark not found at " << aux << '\n';
        return 0;
    }

    PlacementDB db;
    LimboBookshelfAdapter parser(db);
    assert(parser.read(aux.string()));
    ExactBinGrid grid(db, 64, 64, 0.9);
    const ExactDensityResult density = ExactDensityEvaluator().evaluate(grid);
    const double hpwl = HPWLEvaluator().totalHPWL(db);
    assert(std::abs(hpwl - 147771030.0) <= 1.0);
    assert(std::abs(density.ofr - 0.04781179245485895) <= 1e-9);
    return 0;
}
