#include "grid/BinGrid.h"

#include "db/Cell.h"
#include "db/PlacementDB.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
constexpr double kUtilEpsilon = 1e-30;

void requireFinite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("BinGrid: non-finite ") + name);
    }
}

int clampIndex(int value, int lo, int hi) {
    return std::max(lo, std::min(value, hi));
}
}

BinGrid::BinGrid(const PlacementDB& db, int num_bins_x, int num_bins_y, double target_density)
    : num_bins_x_(num_bins_x), num_bins_y_(num_bins_y), target_density_(target_density), core_(db.coreBounds()) {
    validateParameters();
    bin_width_ = core_.width() / static_cast<double>(num_bins_x_);
    bin_height_ = core_.height() / static_cast<double>(num_bins_y_);
    requireFinite(bin_width_, "bin width");
    requireFinite(bin_height_, "bin height");
    if (bin_width_ <= 0.0 || bin_height_ <= 0.0) {
        throw std::invalid_argument("BinGrid: bin dimensions must be positive");
    }
    initializeBins();
    rebuild(db);
}

void BinGrid::validateParameters() const {
    if (num_bins_x_ <= 0 || num_bins_y_ <= 0) {
        throw std::invalid_argument("BinGrid: num_bins_x and num_bins_y must be positive");
    }
    requireFinite(target_density_, "target density");
    if (target_density_ <= 0.0 || target_density_ > 1.0) {
        throw std::invalid_argument("BinGrid: target density must satisfy 0 < density <= 1");
    }
    if (!core_.valid()) {
        throw std::invalid_argument("BinGrid: core bounds must have positive width and height");
    }
}

void BinGrid::initializeBins() {
    bins_.clear();
    bins_.reserve(static_cast<std::size_t>(num_bins_x_) * static_cast<std::size_t>(num_bins_y_));
    for (int iy = 0; iy < num_bins_y_; ++iy) {
        for (int ix = 0; ix < num_bins_x_; ++ix) {
            Bin bin;
            bin.id = binId(ix, iy);
            bin.ix = ix;
            bin.iy = iy;
            // Last row/column are snapped to the core upper bounds to avoid a
            // floating-point accumulation gap at the placement boundary.
            bin.bounds.lx = core_.lx + static_cast<double>(ix) * bin_width_;
            bin.bounds.ux = (ix == num_bins_x_ - 1) ? core_.ux : core_.lx + static_cast<double>(ix + 1) * bin_width_;
            bin.bounds.ly = core_.ly + static_cast<double>(iy) * bin_height_;
            bin.bounds.uy = (iy == num_bins_y_ - 1) ? core_.uy : core_.ly + static_cast<double>(iy + 1) * bin_height_;
            bins_.push_back(std::move(bin));
        }
    }
}

int BinGrid::binId(int ix, int iy) const {
    if (ix < 0 || ix >= num_bins_x_ || iy < 0 || iy >= num_bins_y_) {
        throw std::out_of_range("BinGrid::binId: bin indices out of range");
    }
    return iy * num_bins_x_ + ix;
}

const Bin& BinGrid::bin(int bin_id) const {
    if (bin_id < 0 || bin_id >= static_cast<int>(bins_.size())) {
        throw std::out_of_range("BinGrid::bin: bin id out of range");
    }
    return bins_[static_cast<std::size_t>(bin_id)];
}

void BinGrid::rebuild(const PlacementDB& db) {
    for (Bin& bin : bins_) {
        bin.movable_area = 0.0;
        bin.fixed_area = 0.0;
        bin.physical_capacity = 0.0;
        bin.movable_capacity = 0.0;
        bin.overflow = 0.0;
        bin.utilization = 0.0;
        bin.cell_ids.clear();
    }
    total_movable_area_ = total_fixed_area_ = total_overflow_ = max_utilization_ = 0.0;
    overflow_bin_count_ = 0;
    for (std::size_t i = 0; i < db.cells().size(); ++i) {
        projectCell(db, static_cast<int>(i));
    }
    computeCapacitiesAndStats();
}

void BinGrid::projectCell(const PlacementDB& db, int cell_id) {
    const Cell& cell = db.cell(cell_id);
    requireFinite(cell.x, "cell x");
    requireFinite(cell.y, "cell y");
    requireFinite(cell.width, "cell width");
    requireFinite(cell.height, "cell height");
    if (cell.width < 0.0 || cell.height < 0.0) {
        throw std::invalid_argument("BinGrid: cell has negative width or height: " + cell.name);
    }

    const Box cell_box{cell.x, cell.y, cell.x + cell.width, cell.y + cell.height};
    const double clipped_lx = std::max(cell_box.lx, core_.lx);
    const double clipped_ly = std::max(cell_box.ly, core_.ly);
    const double clipped_ux = std::min(cell_box.ux, core_.ux);
    const double clipped_uy = std::min(cell_box.uy, core_.uy);
    if (clipped_ux <= clipped_lx || clipped_uy <= clipped_ly) {
        return; // Completely outside core or zero area: no positive bin overlap.
    }

    // Use ceil(upper)-1 for the ending index: if a cell's right/top edge lies
    // exactly on a bin boundary, the bin on the far side has zero-area overlap
    // and must not receive the cell id.  Clamp handles cells crossing the core.
    int ix_begin = static_cast<int>(std::floor((clipped_lx - core_.lx) / bin_width_));
    int iy_begin = static_cast<int>(std::floor((clipped_ly - core_.ly) / bin_height_));
    int ix_end = static_cast<int>(std::ceil((clipped_ux - core_.lx) / bin_width_)) - 1;
    int iy_end = static_cast<int>(std::ceil((clipped_uy - core_.ly) / bin_height_)) - 1;
    ix_begin = clampIndex(ix_begin, 0, num_bins_x_ - 1);
    ix_end = clampIndex(ix_end, 0, num_bins_x_ - 1);
    iy_begin = clampIndex(iy_begin, 0, num_bins_y_ - 1);
    iy_end = clampIndex(iy_end, 0, num_bins_y_ - 1);

    for (int iy = iy_begin; iy <= iy_end; ++iy) {
        for (int ix = ix_begin; ix <= ix_end; ++ix) {
            Bin& bin_ref = bins_[static_cast<std::size_t>(binId(ix, iy))];
            // The exact rectangle intersection protects against floating-point
            // edge cases in the candidate index range and defines membership.
            const double overlap = rectangleOverlapArea(cell_box, bin_ref.bounds);
            if (overlap <= 0.0) continue;
            if (cell.isMovable()) bin_ref.movable_area += overlap;
            else bin_ref.fixed_area += overlap;
            bin_ref.cell_ids.push_back(cell_id);
        }
    }
}

void BinGrid::computeCapacitiesAndStats() {
    for (Bin& bin : bins_) {
        bin.physical_capacity = target_density_ * bin.bounds.area();
        // Fixed and terminal cells consume physical target capacity first; only
        // the remaining non-negative capacity is available to movable cells.
        bin.movable_capacity = std::max(0.0, bin.physical_capacity - bin.fixed_area);
        bin.overflow = std::max(0.0, bin.movable_area - bin.movable_capacity);
        if (bin.movable_capacity <= 0.0) {
            bin.utilization = (bin.movable_area > 0.0) ? std::numeric_limits<double>::infinity() : 0.0;
        } else {
            bin.utilization = bin.movable_area / std::max(bin.movable_capacity, kUtilEpsilon);
        }
        total_movable_area_ += bin.movable_area;
        total_fixed_area_ += bin.fixed_area;
        total_overflow_ += bin.overflow;
        if (bin.overflow > 0.0) ++overflow_bin_count_;
        max_utilization_ = std::max(max_utilization_, bin.utilization);
    }
}
