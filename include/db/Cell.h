#pragma once

#include <string>
#include <vector>

struct Cell {
    int id = -1;
    std::string name;
    double width = 0.0;
    double height = 0.0;
    double x = 0.0;
    double y = 0.0;
    bool is_terminal = false;
    bool is_fixed = false;
    std::vector<int> pin_ids;
};
