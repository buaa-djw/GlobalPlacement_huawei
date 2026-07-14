#include "initializer/SimplePackingInitialPlacer.h"
#include "db/PlacementDB.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

SimplePackingInitialPlacer::SimplePackingInitialPlacer(SimplePackingInitialPlacerConfig config) : config_(config) {
    if (!std::isfinite(config_.numerical_epsilon) || config_.numerical_epsilon <= 0.0) throw std::invalid_argument("SimplePackingInitialPlacer: invalid numerical epsilon");
}

InitialPlacementResult SimplePackingInitialPlacer::place(PlacementDB& db) const {
    const double kEpsilon = config_.numerical_epsilon;

    /**
     * @brief One-dimensional interval on a placement row.
     */
    struct Interval {
        double lx = 0.0;
        double ux = 0.0;
    };

    /**
     * @brief Remaining usable segment on one placement row.
     *
     * next_x records the first position that has not yet been occupied by
     * previously initialized movable cells.
     */
    struct FreeSlot {
        std::size_t row_index = 0;
        double lx = 0.0;
        double ux = 0.0;
        double next_x = 0.0;
    };

    InitialPlacementResult result;
    result.method = "simple";
    for (const Cell& cell : db.cells()) if (cell.isMovable()) ++result.movable_cell_count;
    if (db.cells().empty()) { result.message = "empty database"; return result; }

    if (db.rows().empty()) {
        throw std::runtime_error(
            "SimplePackingInitialPlacer: "
            "no placement rows are available"
        );
    }

    /*
     * Align a coordinate to a row's placement-site grid.
     *
     * Site alignment is relative to row.x_start rather than coordinate zero,
     * because different Bookshelf subrows may have different origins.
     */
    const auto alignToSite =
        [kEpsilon](
            double value,
            double row_origin,
            double site_spacing
        ) -> double {
            if (!std::isfinite(value) ||
                !std::isfinite(row_origin) ||
                !std::isfinite(site_spacing)) {
                throw std::runtime_error(
                    "SimplePackingInitialPlacer: "
                    "non-finite row alignment value"
                );
            }

            if (site_spacing <= kEpsilon) {
                return value;
            }

            const double relative =
                (value - row_origin) / site_spacing;

            const double site_index =
                std::ceil(relative - kEpsilon);

            return row_origin +
                   std::max(0.0, site_index) *
                       site_spacing;
        };

    /*
     * Sort row references by y coordinate and then by x origin.
     *
     * The original db.rows() vector remains unchanged.
     */
    std::vector<const Row*> sorted_rows;
    sorted_rows.reserve(db.rows().size());

    for (const Row& row : db.rows()) {
        if (!std::isfinite(row.x_start) ||
            !std::isfinite(row.x_end) ||
            !std::isfinite(row.y) ||
            !std::isfinite(row.height)) {
            throw std::runtime_error(
                "SimplePackingInitialPlacer: "
                "placement row contains non-finite geometry"
            );
        }

        if (row.x_end <= row.x_start + kEpsilon ||
            row.height <= kEpsilon) {
            throw std::runtime_error(
                "SimplePackingInitialPlacer: "
                "placement row has non-positive dimensions"
            );
        }

        sorted_rows.push_back(&row);
    }

    std::sort(
        sorted_rows.begin(),
        sorted_rows.end(),
        [](const Row* lhs, const Row* rhs) {
            if (lhs->y != rhs->y) {
                return lhs->y < rhs->y;
            }

            if (lhs->x_start != rhs->x_start) {
                return lhs->x_start < rhs->x_start;
            }

            return lhs->id < rhs->id;
        }
    );

    /*
     * Store free intervals for each row after subtracting fixed obstacles.
     */
    std::vector<std::vector<Interval>> row_free_intervals(
        sorted_rows.size()
    );

    for (std::size_t row_index = 0;
         row_index < sorted_rows.size();
         ++row_index) {
        const Row& row = *sorted_rows[row_index];

        const double row_ly = row.y;
        const double row_uy = row.y + row.height;

        std::vector<Interval> blocked_intervals;

        /*
         * Only cells explicitly marked fixed by the .pl file are physical
         * obstacles here. Their coordinates must remain unchanged.
         */
        for (const Cell& cell : db.cells()) {
            if (!cell.is_fixed) {
                continue;
            }

            if (!std::isfinite(cell.x) ||
                !std::isfinite(cell.y) ||
                !std::isfinite(cell.width) ||
                !std::isfinite(cell.height)) {
                throw std::runtime_error(
                    "SimplePackingInitialPlacer: "
                    "fixed cell has invalid geometry: " +
                    cell.name
                );
            }

            const double cell_lx = cell.x;
            const double cell_ux = cell.x + cell.width;
            const double cell_ly = cell.y;
            const double cell_uy = cell.y + cell.height;

            /*
             * A fixed cell blocks this row only when their vertical spans
             * overlap by a positive amount.
             */
            const double overlap_y =
                std::min(row_uy, cell_uy) -
                std::max(row_ly, cell_ly);

            if (overlap_y <= kEpsilon) {
                continue;
            }

            const double blocked_lx =
                std::max(row.x_start, cell_lx);

            const double blocked_ux =
                std::min(row.x_end, cell_ux);

            if (blocked_ux <= blocked_lx + kEpsilon) {
                continue;
            }

            blocked_intervals.push_back(
                Interval{blocked_lx, blocked_ux}
            );
        }

        /*
         * Sort and merge overlapping fixed-cell intervals.
         */
        std::sort(
            blocked_intervals.begin(),
            blocked_intervals.end(),
            [](const Interval& lhs, const Interval& rhs) {
                if (lhs.lx != rhs.lx) {
                    return lhs.lx < rhs.lx;
                }

                return lhs.ux < rhs.ux;
            }
        );

        std::vector<Interval> merged_blocked;

        for (const Interval& blocked : blocked_intervals) {
            if (merged_blocked.empty() ||
                blocked.lx >
                    merged_blocked.back().ux + kEpsilon) {
                merged_blocked.push_back(blocked);
            } else {
                merged_blocked.back().ux =
                    std::max(
                        merged_blocked.back().ux,
                        blocked.ux
                    );
            }
        }

        /*
         * Construct free row intervals as the complement of the merged
         * fixed-cell intervals.
         */
        double cursor = row.x_start;

        for (const Interval& blocked : merged_blocked) {
            if (blocked.lx > cursor + kEpsilon) {
                row_free_intervals[row_index].push_back(
                    Interval{cursor, blocked.lx}
                );
            }

            cursor = std::max(cursor, blocked.ux);
        }

        if (cursor < row.x_end - kEpsilon) {
            row_free_intervals[row_index].push_back(
                Interval{cursor, row.x_end}
            );
        }
    }

    /*
     * Create a pool of all currently free row slots.
     *
     * The key is the usable width after site alignment. std::multimap allows
     * selection of the smallest slot that is still large enough, which
     * reduces fragmentation compared with always selecting the first row.
     */
    std::multimap<double, FreeSlot> free_slots;

    const auto insertFreeSlot =
        [&](const FreeSlot& slot) {
            const Row& row =
                *sorted_rows[slot.row_index];

            const double spacing =
                row.site_spacing > kEpsilon
                    ? row.site_spacing
                    : (
                        row.site_width > kEpsilon
                            ? row.site_width
                            : 1.0
                    );

            const double aligned_x =
                alignToSite(
                    std::max(slot.next_x, slot.lx),
                    row.x_start,
                    spacing
                );

            const double available_width =
                slot.ux - aligned_x;

            if (available_width > kEpsilon) {
                free_slots.emplace(
                    available_width,
                    slot
                );
            }
        };

    for (std::size_t row_index = 0;
         row_index < row_free_intervals.size();
         ++row_index) {
        for (const Interval& interval :
             row_free_intervals[row_index]) {
            insertFreeSlot(
                FreeSlot{
                    row_index,
                    interval.lx,
                    interval.ux,
                    interval.lx
                }
            );
        }
    }

    /*
     * Collect only true movable cells.
     *
     * Terminal and TerminalNI cells are not initialized here, even when their
     * is_fixed flag is false. Their type already excludes them from global
     * placement through Cell::isMovable().
     */
    std::vector<int> movable_cell_ids;
    movable_cell_ids.reserve(db.cells().size());

    for (const Cell& cell : db.cells()) {
        if (cell.isMovable()) {
            movable_cell_ids.push_back(cell.id);
        }
    }

    /*
     * Place larger cells first. This deterministic ordering reduces the risk
     * that large cells are left with only small fragmented slots.
     */
    std::stable_sort(
        movable_cell_ids.begin(),
        movable_cell_ids.end(),
        [&](int lhs_id, int rhs_id) {
            const Cell& lhs =
                db.cells()[static_cast<std::size_t>(lhs_id)];

            const Cell& rhs =
                db.cells()[static_cast<std::size_t>(rhs_id)];

            if (lhs.height != rhs.height) {
                return lhs.height > rhs.height;
            }

            if (lhs.width != rhs.width) {
                return lhs.width > rhs.width;
            }

            return lhs.id < rhs.id;
        }
    );

    std::size_t placed_count = 0;

    for (const int cell_id : movable_cell_ids) {
        const Cell& cell =
            db.cells()[static_cast<std::size_t>(cell_id)];

        if (!std::isfinite(cell.width) ||
            !std::isfinite(cell.height) ||
            cell.width < 0.0 ||
            cell.height < 0.0) {
            throw std::runtime_error(
                "SimplePackingInitialPlacer: "
                "movable cell has invalid dimensions: " +
                cell.name
            );
        }

        /*
         * Slots belonging to rows that are too short are removed temporarily.
         * They are restored after searching for the current cell.
         */
        std::vector<std::pair<double, FreeSlot>>
            deferred_slots;

        auto slot_it =
            free_slots.lower_bound(
                std::max(0.0, cell.width - kEpsilon)
            );

        bool placed = false;
        FreeSlot selected_slot;
        double selected_x = 0.0;

        while (slot_it != free_slots.end()) {
            const FreeSlot candidate =
                slot_it->second;

            const Row& row =
                *sorted_rows[candidate.row_index];

            /*
             * This simple initializer supports cells that fit within one row.
             * Multi-row movable macros should be handled by a dedicated macro
             * placement strategy.
             */
            if (cell.height >
                row.height + kEpsilon) {
                deferred_slots.push_back(*slot_it);
                slot_it =
                    free_slots.erase(slot_it);
                continue;
            }

            const double spacing =
                row.site_spacing > kEpsilon
                    ? row.site_spacing
                    : (
                        row.site_width > kEpsilon
                            ? row.site_width
                            : 1.0
                    );

            const double candidate_x =
                alignToSite(
                    std::max(
                        candidate.next_x,
                        candidate.lx
                    ),
                    row.x_start,
                    spacing
                );

            if (candidate_x + cell.width <=
                candidate.ux + kEpsilon) {
                selected_slot = candidate;
                selected_x = candidate_x;
                placed = true;
                break;
            }

            /*
             * The multimap key should normally prevent this case, but remove
             * stale or unusable slots defensively.
             */
            slot_it = free_slots.erase(slot_it);
        }

        for (const auto& deferred :
             deferred_slots) {
            free_slots.insert(deferred);
        }

        if (!placed) {
            throw std::runtime_error(
                "SimplePackingInitialPlacer: "
                "cannot place movable cell '" +
                cell.name +
                "' with size " +
                std::to_string(cell.width) +
                " x " +
                std::to_string(cell.height)
            );
        }

        const Row& selected_row =
            *sorted_rows[selected_slot.row_index];

        /*
         * Remove the selected slot before inserting its remaining portion.
         */
        free_slots.erase(slot_it);

        db.setCellPosition(cell.id, selected_x, selected_row.y);

        selected_slot.next_x =
            selected_x + cell.width;

        insertFreeSlot(selected_slot);

        ++placed_count;
    }

    if (placed_count != movable_cell_ids.size()) {
        throw std::runtime_error(
            "SimplePackingInitialPlacer: "
            "not all movable cells were initialized"
        );
    }
    result.message = "simple packing completed";
    return result;
}
