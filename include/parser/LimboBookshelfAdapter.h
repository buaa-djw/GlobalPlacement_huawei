#pragma once

#include <string>

class PlacementDB;

class LimboBookshelfAdapter {
public:
    explicit LimboBookshelfAdapter(PlacementDB& db);
    bool read(const std::string& aux_path);
private:
    PlacementDB& db_;
};
