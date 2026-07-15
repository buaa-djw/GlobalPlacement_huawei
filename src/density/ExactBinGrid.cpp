#include "density/ExactBinGrid.h"

#include "db/Cell.h"
#include "db/PlacementDB.h"
#include "density/PaperOverlapFunction.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{

constexpr double kRelativeTolerance = 1e-9;

bool near(double lhs, double rhs)
{
    return std::abs(lhs - rhs) <=
           kRelativeTolerance *
               std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

void requireFiniteNonnegative(
    double value,
    const std::string& name
)
{
    if (!std::isfinite(value) ||
        value < -kRelativeTolerance) {
        throw std::runtime_error(
            "ExactBinGrid: invalid " + name
        );
    }
}

double paperOverlapArea(
    const Cell& cell,
    const Bin& bin
)
{
    const double cell_center_x = cell.x + 0.5 * cell.width;
    const double cell_center_y = cell.y + 0.5 * cell.height;
    const double bin_center_x = 0.5 * (bin.bounds.lx + bin.bounds.ux);
    const double bin_center_y = 0.5 * (bin.bounds.ly + bin.bounds.uy);

    const double overlap_x = paperAxisOverlapLength(
        cell_center_x, cell.width, bin_center_x, bin.bounds.width());
    const double overlap_y = paperAxisOverlapLength(
        cell_center_y, cell.height, bin_center_y, bin.bounds.height());

    return overlap_x * overlap_y;
}


Box clipToCore(
    const Box& object_box,
    const Box& core
)
{
    return Box{
        std::max(object_box.lx, core.lx),
        std::max(object_box.ly, core.ly),
        std::min(object_box.ux, core.ux),
        std::min(object_box.uy, core.uy)
    };
}

} // namespace

ExactBinGrid::ExactBinGrid(
    const PlacementDB& db,
    int bins_x,
    int bins_y,
    double target_density
)
    : bins_x_(bins_x),
      bins_y_(bins_y),
      target_density_(target_density)
{
    if (bins_x_ <= 0 || bins_y_ <= 0) {
        throw std::invalid_argument(
            "ExactBinGrid: bins_x and bins_y must be positive"
        );
    }

    if (!std::isfinite(target_density_) ||
        target_density_ <= 0.0 ||
        target_density_ > 1.0) {
        throw std::invalid_argument(
            "ExactBinGrid: target_density must satisfy "
            "0 < value <= 1"
        );
    }

    rebuild(db);
}

const Bin& ExactBinGrid::bin(int id) const
{
    if (id < 0 ||
        id >= static_cast<int>(bins_.size())) {
        throw std::out_of_range(
            "ExactBinGrid::bin: invalid bin id"
        );
    }

    return bins_[static_cast<std::size_t>(id)];
}

void ExactBinGrid::rebuild(const PlacementDB& db)
{
    core_ = db.coreBounds();

    if (!core_.valid() ||
        !std::isfinite(core_.area())) {
        throw std::runtime_error(
            "ExactBinGrid: invalid core bounds"
        );
    }

    bin_width_ =
        core_.width() / static_cast<double>(bins_x_);

    bin_height_ =
        core_.height() / static_cast<double>(bins_y_);

    if (!(bin_width_ > 0.0) ||
        !(bin_height_ > 0.0) ||
        !std::isfinite(bin_width_) ||
        !std::isfinite(bin_height_)) {
        throw std::runtime_error(
            "ExactBinGrid: invalid bin dimensions"
        );
    }

    bins_.clear();
    bins_.reserve(
        static_cast<std::size_t>(bins_x_) *
        static_cast<std::size_t>(bins_y_)
    );

    for (int iy = 0; iy < bins_y_; ++iy) {
        for (int ix = 0; ix < bins_x_; ++ix) {
            Bin bin;

            bin.id = iy * bins_x_ + ix;
            bin.ix = ix;
            bin.iy = iy;

            bin.bounds.lx =
                core_.lx +
                static_cast<double>(ix) * bin_width_;

            bin.bounds.ly =
                core_.ly +
                static_cast<double>(iy) * bin_height_;

            bin.bounds.ux =
                (ix == bins_x_ - 1)
                    ? core_.ux
                    : core_.lx +
                          static_cast<double>(ix + 1) *
                              bin_width_;

            bin.bounds.uy =
                (iy == bins_y_ - 1)
                    ? core_.uy
                    : core_.ly +
                          static_cast<double>(iy + 1) *
                              bin_height_;

            const double bin_area = bin.bounds.area();

            if (!(bin_area > 0.0) ||
                !std::isfinite(bin_area)) {
                throw std::runtime_error(
                    "ExactBinGrid: bin area must be positive"
                );
            }

            bin.target_capacity =
                target_density_ * bin_area;

            bins_.push_back(std::move(bin));
        }
    }

    total_cell_area_ = 0.0;
    total_movable_area_ = 0.0;
    total_non_movable_area_ = 0.0;

    raw_movable_area_ = 0.0;
    raw_non_movable_area_ = 0.0;

    movable_cell_count_ = 0;
    non_movable_cell_count_ = 0;

    movable_inside_core_count_ = 0;
    movable_outside_core_count_ = 0;

    non_movable_inside_core_count_ = 0;
    non_movable_outside_core_count_ = 0;

    double expected_total_area = 0.0;
    double expected_movable_area = 0.0;
    double expected_non_movable_area = 0.0;

    for (const Cell& cell : db.cells()) {
        if (!std::isfinite(cell.x) ||
            !std::isfinite(cell.y) ||
            !std::isfinite(cell.width) ||
            !std::isfinite(cell.height) ||
            cell.width < 0.0 ||
            cell.height < 0.0) {
            throw std::runtime_error(
                "ExactBinGrid: invalid cell geometry for " +
                cell.name
            );
        }

        const double raw_area =
            cell.width * cell.height;

        requireFiniteNonnegative(
            raw_area,
            "raw cell area"
        );

        if (cell.isMovable()) {
            ++movable_cell_count_;
            raw_movable_area_ += raw_area;
        } else {
            ++non_movable_cell_count_;
            raw_non_movable_area_ += raw_area;
        }

        const Box cell_box{
            cell.x,
            cell.y,
            cell.x + cell.width,
            cell.y + cell.height
        };

        const Box clipped =
            clipToCore(cell_box, core_);

        if (!clipped.valid()) {
            if (cell.isMovable()) {
                ++movable_outside_core_count_;
            } else {
                ++non_movable_outside_core_count_;
            }

            continue;
        }

        const double clipped_area =
            clipped.area();

        requireFiniteNonnegative(
            clipped_area,
            "clipped cell area"
        );

        if (cell.isMovable()) {
            ++movable_inside_core_count_;
            expected_movable_area += clipped_area;
        } else {
            ++non_movable_inside_core_count_;
            expected_non_movable_area += clipped_area;
        }

        expected_total_area += clipped_area;

        int ix_begin = static_cast<int>(
            std::floor(
                (clipped.lx - core_.lx) /
                bin_width_
            )
        );

        int iy_begin = static_cast<int>(
            std::floor(
                (clipped.ly - core_.ly) /
                bin_height_
            )
        );

        int ix_end = static_cast<int>(
            std::ceil(
                (clipped.ux - core_.lx) /
                bin_width_
            )
        ) - 1;

        int iy_end = static_cast<int>(
            std::ceil(
                (clipped.uy - core_.ly) /
                bin_height_
            )
        ) - 1;

        ix_begin = std::clamp(
            ix_begin,
            0,
            bins_x_ - 1
        );

        iy_begin = std::clamp(
            iy_begin,
            0,
            bins_y_ - 1
        );

        ix_end = std::clamp(
            ix_end,
            0,
            bins_x_ - 1
        );

        iy_end = std::clamp(
            iy_end,
            0,
            bins_y_ - 1
        );

        for (int iy = iy_begin; iy <= iy_end; ++iy) {
            for (int ix = ix_begin; ix <= ix_end; ++ix) {
                Bin& bin =
                    bins_[static_cast<std::size_t>(
                        iy * bins_x_ + ix
                    )];

                const double overlap_area =
                    paperOverlapArea(cell, bin);

                if (overlap_area <= 0.0) {
                    continue;
                }

                requireFiniteNonnegative(
                    overlap_area,
                    "cell-bin overlap area"
                );

                bin.total_area += overlap_area;

                if (cell.isMovable()) {
                    bin.movable_area += overlap_area;
                } else {
                    // In the current data model this means all
                    // non-movable objects, including terminals.
                    bin.fixed_area += overlap_area;
                }

                bin.overlapping_cell_ids.push_back(
                    cell.id
                );
            }
        }
    }

    for (Bin& bin : bins_) {
        if (!near(
                bin.total_area,
                bin.movable_area + bin.fixed_area)) {
            throw std::runtime_error(
                "ExactBinGrid: bin area split mismatch "
                "for bin " +
                std::to_string(bin.id)
            );
        }

        bin.overflow = std::max(
            0.0,
            bin.total_area - bin.target_capacity
        );

        total_cell_area_ += bin.total_area;
        total_movable_area_ += bin.movable_area;
        total_non_movable_area_ += bin.fixed_area;
    }

    if (!near(
            total_cell_area_,
            expected_total_area)) {
        throw std::runtime_error(
            "ExactBinGrid: total clipped area "
            "conservation mismatch"
        );
    }

    if (!near(
            total_movable_area_,
            expected_movable_area)) {
        throw std::runtime_error(
            "ExactBinGrid: movable clipped area "
            "conservation mismatch"
        );
    }

    if (!near(
            total_non_movable_area_,
            expected_non_movable_area)) {
        throw std::runtime_error(
            "ExactBinGrid: non-movable clipped area "
            "conservation mismatch"
        );
    }

    if (!near(
            total_cell_area_,
            total_movable_area_ +
                total_non_movable_area_)) {
        throw std::runtime_error(
            "ExactBinGrid: total area does not equal "
            "movable plus non-movable area"
        );
    }
}