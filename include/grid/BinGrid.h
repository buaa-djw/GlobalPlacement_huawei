#pragma once

#include "geometry/Box.h"
#include "grid/Bin.h"

#include <vector>

class PlacementDB;

/**
 * @brief Regular 2-D bin grid for geometric density projection.
 *
 * BinGrid is deliberately optimization-free: it projects cell rectangles onto a
 * fixed core grid, subtracts fixed/terminal blockage from movable capacity, and
 * reports overflow statistics.  Invalid inputs are reported with standard
 * exceptions; the executable top level owns logging and process termination so
 * library code remains testable and reusable.
 */
class BinGrid {
public:
    /** @brief Build a grid from the database core and immediately project cells. */
    BinGrid(const PlacementDB& db, int num_bins_x, int num_bins_y, double target_density);

    /** @brief Recompute all per-bin areas and statistics using current cell positions. */
    void rebuild(const PlacementDB& db);

    int numBinsX() const { return num_bins_x_; }
    int numBinsY() const { return num_bins_y_; }
    int numBins() const { return static_cast<int>(bins_.size()); }
    double binWidth() const { return bin_width_; }
    double binHeight() const { return bin_height_; }
    double targetDensity() const { return target_density_; }
    const Box& coreBounds() const { return core_; }
    int binId(int ix, int iy) const;
    const Bin& bin(int bin_id) const;
    const Bin& bin(int ix, int iy) const { return bin(binId(ix, iy)); }
    const std::vector<Bin>& bins() const { return bins_; }
    double totalMovableArea() const { return total_movable_area_; }
    double totalFixedArea() const { return total_fixed_area_; }
    double totalOverflow() const { return total_overflow_; }
    int overflowBinCount() const { return overflow_bin_count_; }
    double maxUtilization() const { return max_utilization_; }

private:
    void validateParameters() const;
    void initializeBins();
    void projectCell(const PlacementDB& db, int cell_id);
    void computeCapacitiesAndStats();

    int num_bins_x_ = 0;
    int num_bins_y_ = 0;
    double target_density_ = 0.0;
    Box core_;
    double bin_width_ = 0.0;
    double bin_height_ = 0.0;
    std::vector<Bin> bins_;
    double total_movable_area_ = 0.0;
    double total_fixed_area_ = 0.0;
    double total_overflow_ = 0.0;
    int overflow_bin_count_ = 0;
    double max_utilization_ = 0.0;
};
