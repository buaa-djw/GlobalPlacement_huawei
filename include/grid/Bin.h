#pragma once

#include "geometry/Box.h"

#include <vector>

/**
 * @brief Per-bin density accounting record.
 *
 * Bins store only stable integer cell ids, never Cell objects or pointers, so a
 * rebuild remains valid after the PlacementDB vector storage changes.
 */
struct Bin {
    int id = -1;
    int ix = -1;
    int iy = -1;
    Box bounds;
    double movable_area = 0.0;
    double fixed_area = 0.0;
    double physical_capacity = 0.0;
    double movable_capacity = 0.0;
    double overflow = 0.0;
    double utilization = 0.0;
    std::vector<int> cell_ids;
};
