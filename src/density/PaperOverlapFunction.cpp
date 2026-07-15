#include "density/PaperOverlapFunction.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRelativeBoundaryTolerance = 1e-12;

bool near(double a, double b)
{
    return std::abs(a - b) <=
           kRelativeBoundaryTolerance *
               std::max({1.0, std::abs(a), std::abs(b)});
}

void requireSize(double value, const char* name)
{
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(name);
    }
}

double randomCos(double lo, double hi, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> dist(lo, hi);
    return std::cos(dist(rng));
}
} // namespace

double paperAxisOverlapLength(
    double cell_center,
    double cell_size,
    double bin_center,
    double bin_size
)
{
    requireSize(cell_size, "paperAxisOverlapLength: invalid cell size");
    requireSize(bin_size, "paperAxisOverlapLength: invalid bin size");
    if (!std::isfinite(cell_center) || !std::isfinite(bin_center)) {
        throw std::invalid_argument("paperAxisOverlapLength: invalid center");
    }

    const double z = cell_center - bin_center;
    const double d = std::abs(z);
    const double lower = std::abs(bin_size - cell_size) * 0.5;
    const double upper = (bin_size + cell_size) * 0.5;
    const double plateau = std::min(bin_size, cell_size);

    if (d <= lower || near(d, lower)) {
        return plateau;
    }
    if (d < upper && !near(d, upper)) {
        return upper - d;
    }
    return 0.0;
}

AxisOverlapSubgradient paperAxisOverlapSubgradient(
    double cell_center,
    double cell_size,
    double bin_center,
    double bin_size,
    std::mt19937_64& random_engine
)
{
    const double overlap = paperAxisOverlapLength(
        cell_center, cell_size, bin_center, bin_size);

    const double z = cell_center - bin_center;
    const double d = std::abs(z);
    const double lower = std::abs(bin_size - cell_size) * 0.5;
    const double upper = (bin_size + cell_size) * 0.5;

    double subgradient = 0.0;
    if (near(lower, 0.0) && near(z, 0.0)) {
        subgradient = randomCos(0.0, kPi, random_engine);
    } else if (near(d, lower) || near(d, upper)) {
        if (z >= 0.0) {
            subgradient = randomCos(0.5 * kPi, kPi, random_engine);
        } else {
            subgradient = randomCos(0.0, 0.5 * kPi, random_engine);
        }
    } else if (d < lower) {
        subgradient = 0.0;
    } else if (d < upper) {
        subgradient = (z < 0.0) ? 1.0 : -1.0;
    } else {
        subgradient = 0.0;
    }

    return AxisOverlapSubgradient{overlap, subgradient};
}
