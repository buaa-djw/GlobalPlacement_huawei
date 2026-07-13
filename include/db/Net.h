#pragma once

#include <string>
#include <vector>

struct Net {
    int id = -1;
    double hpwl = 0.0;
    std::string name;
    std::vector<int> pin_ids;
};
