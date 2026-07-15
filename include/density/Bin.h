#pragma once
#include "geometry/Box.h"
#include <vector>

struct Bin {
    int id = -1;
    int ix = -1;
    int iy = -1;
    Box bounds;
    double total_area = 0.0;
    double movable_area = 0.0;
    double fixed_area = 0.0;
    double target_capacity = 0.0;
    double overflow = 0.0;
    std::vector<int> overlapping_cell_ids;
};
