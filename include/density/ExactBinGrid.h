#pragma once
#include "density/Bin.h"
#include "geometry/Box.h"
#include <vector>
class PlacementDB;
class ExactBinGrid {
public:
    ExactBinGrid(const PlacementDB& db, int bins_x, int bins_y, double target_density);
    void rebuild(const PlacementDB& db);
    int numBinsX() const { return bins_x_; }
    int numBinsY() const { return bins_y_; }
    int numBins() const { return static_cast<int>(bins_.size()); }
    double binWidth() const { return bin_width_; }
    double binHeight() const { return bin_height_; }
    double targetDensity() const { return target_density_; }
    const Box& coreBounds() const { return core_; }
    const std::vector<Bin>& bins() const { return bins_; }
    const Bin& bin(int id) const;
    double totalCellArea() const { return total_cell_area_; }
    double totalMovableArea() const { return total_movable_area_; }
    double totalFixedArea() const { return total_fixed_area_; }
private:
    int bins_x_ = 0, bins_y_ = 0;
    double target_density_ = 0.0, bin_width_ = 0.0, bin_height_ = 0.0;
    Box core_;
    std::vector<Bin> bins_;
    double total_cell_area_ = 0.0, total_movable_area_ = 0.0, total_fixed_area_ = 0.0;
};
