#pragma once

#include <string>

struct Pin {
    int id = -1;
    int cell_id = -1;
    int net_id = -1;
    std::string cell_name;
    std::string net_name;
    double offset_x = 0.0;
    double offset_y = 0.0;
    std::string direction;
};
