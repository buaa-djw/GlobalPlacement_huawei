#include "wirelength/ExactWirelengthSubgradient.h"

#include "db/Cell.h"
#include "db/Net.h"
#include "db/Pin.h"
#include "db/PlacementDB.h"
#include "geometry/Point.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTieRelativeTolerance = 1e-12;

struct AxisPair
{
    int first_pin = -1;
    int second_pin = -1;
    double weight = 0.0;
};

struct AxisPin
{
    int pin_id = -1;
    double coordinate = 0.0;
};

bool nearTie(double diff, double a, double b)
{
    return std::abs(diff) <=
           kTieRelativeTolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

void requireFinite(double value, const char* name)
{
    if (!std::isfinite(value)) {
        throw std::runtime_error(name);
    }
}

void addUniqueSortedPair(std::vector<AxisPair>& pairs, int a, int b, double weight)
{
    if (a == b) {
        throw std::runtime_error("ExactWirelengthSubgradient: self pin pair");
    }
    if (b < a) std::swap(a, b);
    pairs.push_back(AxisPair{a, b, weight});
}

std::vector<AxisPair> buildAxisPairs(
    const PlacementDB& db,
    const Net& net,
    bool x_axis)
{
    const std::size_t p = net.pin_ids.size();
    std::vector<AxisPair> pairs;
    if (p <= 1) return pairs;

    if (p == 2) {
        addUniqueSortedPair(pairs, net.pin_ids[0], net.pin_ids[1], 1.0);
    } else if (p == 3) {
        for (std::size_t i = 0; i < p; ++i) {
            for (std::size_t j = i + 1; j < p; ++j) {
                addUniqueSortedPair(pairs, net.pin_ids[i], net.pin_ids[j], 0.5);
            }
        }
    } else {
        std::vector<AxisPin> axis_pins;
        axis_pins.reserve(p);
        for (int pin_id : net.pin_ids) {
            const Point pos = db.pinPosition(pin_id);
            const double coord = x_axis ? pos.x : pos.y;
            requireFinite(coord, "ExactWirelengthSubgradient: non-finite pin coordinate");
            axis_pins.push_back(AxisPin{pin_id, coord});
        }
        std::sort(axis_pins.begin(), axis_pins.end(), [](const AxisPin& lhs, const AxisPin& rhs) {
            if (lhs.coordinate != rhs.coordinate) return lhs.coordinate < rhs.coordinate;
            return lhs.pin_id < rhs.pin_id;
        });
        const int min_pin = axis_pins.front().pin_id;
        const int max_pin = axis_pins.back().pin_id;
        const double weight = 1.0 / static_cast<double>(p - 1);
        for (std::size_t i = 1; i < p; ++i) {
            addUniqueSortedPair(pairs, min_pin, axis_pins[i].pin_id, weight);
        }
        for (std::size_t i = 1; i + 1 < p; ++i) {
            addUniqueSortedPair(pairs, max_pin, axis_pins[i].pin_id, weight);
        }
    }

    std::sort(pairs.begin(), pairs.end(), [](const AxisPair& lhs, const AxisPair& rhs) {
        return std::pair<int, int>{lhs.first_pin, lhs.second_pin} <
               std::pair<int, int>{rhs.first_pin, rhs.second_pin};
    });
    pairs.erase(std::unique(pairs.begin(), pairs.end(), [](const AxisPair& lhs, const AxisPair& rhs) {
        return lhs.first_pin == rhs.first_pin && lhs.second_pin == rhs.second_pin;
    }), pairs.end());
    return pairs;
}

double randomCos(double lo, double hi, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> dist(lo, hi);
    return std::cos(dist(rng));
}

double pairSubgradient(
    double first,
    double second,
    int first_pin,
    int second_pin,
    const WirelengthCoordinateSnapshot* previous,
    bool x_axis,
    std::mt19937_64& rng)
{
    const double diff = first - second;
    if (!nearTie(diff, first, second)) return diff > 0.0 ? 1.0 : -1.0;
    if (previous == nullptr) return 1.0;
    const std::vector<double>& prev = x_axis ? previous->pin_x : previous->pin_y;
    const double prev_first = prev[static_cast<std::size_t>(first_pin)];
    const double prev_second = prev[static_cast<std::size_t>(second_pin)];
    const double prev_diff = prev_first - prev_second;
    if (!nearTie(prev_diff, prev_first, prev_second)) {
        if (prev_diff > 0.0) return randomCos(0.0, kPi / 3.0, rng);
        return randomCos(2.0 * kPi / 3.0, kPi, rng);
    }
    return 1.0;
}

void evaluateAxis(
    const PlacementDB& db,
    const Net& net,
    const WirelengthCoordinateSnapshot* previous,
    bool x_axis,
    std::mt19937_64& rng,
    ExactWirelengthSubgradientResult& result)
{
    const std::vector<AxisPair> pairs = buildAxisPairs(db, net, x_axis);
    for (const AxisPair& pair : pairs) {
        const Point first_pos = db.pinPosition(pair.first_pin);
        const Point second_pos = db.pinPosition(pair.second_pin);
        const double first = x_axis ? first_pos.x : first_pos.y;
        const double second = x_axis ? second_pos.x : second_pos.y;
        requireFinite(first, "ExactWirelengthSubgradient: non-finite first coordinate");
        requireFinite(second, "ExactWirelengthSubgradient: non-finite second coordinate");
        const double diff = first - second;
        result.hpwl += pair.weight * std::abs(diff);
        const double g = pair.weight * pairSubgradient(
            first, second, pair.first_pin, pair.second_pin, previous, x_axis, rng);
        const int first_cell = db.pin(pair.first_pin).cell_id;
        const int second_cell = db.pin(pair.second_pin).cell_id;
        if (db.cell(first_cell).isMovable()) {
            (x_axis ? result.grad_x : result.grad_y)[static_cast<std::size_t>(first_cell)] += g;
        }
        if (db.cell(second_cell).isMovable()) {
            (x_axis ? result.grad_x : result.grad_y)[static_cast<std::size_t>(second_cell)] -= g;
        }
    }
}
} // namespace

WirelengthCoordinateSnapshot ExactWirelengthSubgradient::capture(const PlacementDB& db) const
{
    WirelengthCoordinateSnapshot snapshot;
    snapshot.pin_x.reserve(db.pins().size());
    snapshot.pin_y.reserve(db.pins().size());
    for (const Pin& pin : db.pins()) {
        const Point pos = db.pinPosition(pin.id);
        requireFinite(pos.x, "ExactWirelengthSubgradient::capture: non-finite pin x");
        requireFinite(pos.y, "ExactWirelengthSubgradient::capture: non-finite pin y");
        snapshot.pin_x.push_back(pos.x);
        snapshot.pin_y.push_back(pos.y);
    }
    return snapshot;
}

ExactWirelengthSubgradientResult ExactWirelengthSubgradient::evaluate(
    const PlacementDB& db,
    const WirelengthCoordinateSnapshot* previous,
    std::mt19937_64& random_engine) const
{
    if (previous != nullptr &&
        (previous->pin_x.size() != db.pins().size() || previous->pin_y.size() != db.pins().size())) {
        throw std::invalid_argument("ExactWirelengthSubgradient: previous snapshot pin count mismatch");
    }
    if (previous != nullptr) {
        for (double v : previous->pin_x) requireFinite(v, "ExactWirelengthSubgradient: non-finite previous pin x");
        for (double v : previous->pin_y) requireFinite(v, "ExactWirelengthSubgradient: non-finite previous pin y");
    }

    ExactWirelengthSubgradientResult result;
    result.grad_x.assign(db.cells().size(), 0.0);
    result.grad_y.assign(db.cells().size(), 0.0);

    for (const Net& net : db.nets()) {
        evaluateAxis(db, net, previous, true, random_engine, result);
        evaluateAxis(db, net, previous, false, random_engine, result);
    }

    if (!std::isfinite(result.hpwl) || result.hpwl < 0.0) {
        throw std::runtime_error("ExactWirelengthSubgradient: invalid HPWL");
    }
    for (const Cell& cell : db.cells()) {
        if (!cell.isMovable()) {
            result.grad_x[static_cast<std::size_t>(cell.id)] = 0.0;
            result.grad_y[static_cast<std::size_t>(cell.id)] = 0.0;
        }
    }
    for (double v : result.grad_x) requireFinite(v, "ExactWirelengthSubgradient: non-finite grad_x");
    for (double v : result.grad_y) requireFinite(v, "ExactWirelengthSubgradient: non-finite grad_y");
    return result;
}
