#pragma once

#include "density/Bin.h"
#include "geometry/Box.h"

#include <cstddef>
#include <vector>

class PlacementDB;

/**
 * @brief Exact cell-bin overlap grid.
 *
 * All cells contribute to total bin density.
 *
 * A cell is classified as movable only when Cell::isMovable() is true.
 * The existing fixed_area field in Bin therefore represents the area of all
 * non-movable objects, including terminals and fixed standard cells.
 */
class ExactBinGrid
{
public:
    ExactBinGrid(
        const PlacementDB& db,
        int bins_x,
        int bins_y,
        double target_density
    );

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

    // Area accumulated inside the placement core.
    double totalCellArea() const { return total_cell_area_; }
    double totalMovableArea() const { return total_movable_area_; }

    // Kept for compatibility. This means all non-movable object area.
    double totalFixedArea() const { return total_non_movable_area_; }
    double totalNonMovableArea() const {
        return total_non_movable_area_;
    }

    // Raw object area before clipping against the placement core.
    double rawMovableArea() const { return raw_movable_area_; }
    double rawNonMovableArea() const {
        return raw_non_movable_area_;
    }

    std::size_t movableCellCount() const {
        return movable_cell_count_;
    }

    std::size_t nonMovableCellCount() const {
        return non_movable_cell_count_;
    }

    std::size_t movableInsideCoreCount() const {
        return movable_inside_core_count_;
    }

    std::size_t movableOutsideCoreCount() const {
        return movable_outside_core_count_;
    }

    std::size_t nonMovableInsideCoreCount() const {
        return non_movable_inside_core_count_;
    }

    std::size_t nonMovableOutsideCoreCount() const {
        return non_movable_outside_core_count_;
    }

private:
    int bins_x_ = 0;
    int bins_y_ = 0;

    double target_density_ = 0.0;
    double bin_width_ = 0.0;
    double bin_height_ = 0.0;

    Box core_;
    std::vector<Bin> bins_;

    // Areas after clipping objects to the core.
    double total_cell_area_ = 0.0;
    double total_movable_area_ = 0.0;
    double total_non_movable_area_ = 0.0;

    // Areas before clipping.
    double raw_movable_area_ = 0.0;
    double raw_non_movable_area_ = 0.0;

    std::size_t movable_cell_count_ = 0;
    std::size_t non_movable_cell_count_ = 0;

    std::size_t movable_inside_core_count_ = 0;
    std::size_t movable_outside_core_count_ = 0;

    std::size_t non_movable_inside_core_count_ = 0;
    std::size_t non_movable_outside_core_count_ = 0;
};