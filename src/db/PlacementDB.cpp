#include "db/PlacementDB.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <map>
#include <utility>

int PlacementDB::addCell(const std::string& name, double width, double height, CellType type) {

    if (name.empty()) throw std::invalid_argument("PlacementDB::addCell: empty cell name");
    if (hasCell(name)) throw std::runtime_error("PlacementDB::addCell: duplicate cell '" + name + "'");
    if (!std::isfinite(width) ||!std::isfinite(height) ||width < 0.0 ||height < 0.0) {
        throw std::invalid_argument("PlacementDB::addCell: invalid dimensions for cell '" +name + "'");
    }

    const int id = static_cast<int>(cells_.size());

    Cell cell;
    cell.id = id;
    cell.name = name;
    cell.width = width;
    cell.height = height;

    // Coordinates will be initialized after all Bookshelf files are parsed.
    cell.x = 0.0;
    cell.y = 0.0;

    cell.type = type;

    // Terminal type does not automatically imply a valid fixed placement.
    cell.is_fixed = false;

    cell.orientation = "N";
    
    cells_.push_back(std::move(cell));
    cell_name_to_id_[name] = id;

    return id;
}

int PlacementDB::addNet(const std::string& name) {
    if (name.empty()) throw std::invalid_argument("PlacementDB::addNet: empty net name");
    if (hasNet(name)) throw std::runtime_error("PlacementDB::addNet: duplicate net '" + name + "'");
    const int id = static_cast<int>(nets_.size());
    nets_.push_back(Net{id, 0.0 ,name, {}});
    net_name_to_id_[name] = id;
    return id;
}

void PlacementDB::addNetHPWL(int id, double hpwl) {
    checkId(id, static_cast<int>(nets_.size()), "net");
    if(hpwl < 0.0) throw std::invalid_argument("PlacementDB::addNetHPWL: negative HPWL value");
    nets_[id].hpwl = hpwl;
}

int PlacementDB::addPin(int cell_id, int net_id, double offset_x, double offset_y, const std::string& direction) {
    checkId(cell_id, static_cast<int>(cells_.size()), "cell");
    checkId(net_id, static_cast<int>(nets_.size()), "net");
    const int id = static_cast<int>(pins_.size());
    pins_.push_back(Pin{id, cell_id, net_id, cells_[cell_id].name, nets_[net_id].name, offset_x, offset_y, direction});
    cells_[cell_id].pin_ids.push_back(id);
    nets_[net_id].pin_ids.push_back(id);
    return id;
}

int PlacementDB::addRow(const Row& row) {
    Row copy = row;
    copy.id = static_cast<int>(rows_.size());
    if (copy.x_end == 0.0 && copy.num_sites > 0) copy.x_end = copy.x_start + copy.num_sites * copy.site_spacing;
    rows_.push_back(copy);
    return copy.id;
}

bool PlacementDB::hasCell(const std::string& name) const { return cell_name_to_id_.count(name) != 0; }
bool PlacementDB::hasNet(const std::string& name) const { return net_name_to_id_.count(name) != 0; }

int PlacementDB::getCellId(const std::string& name) const {
    auto it = cell_name_to_id_.find(name);
    if (it == cell_name_to_id_.end()) throw std::runtime_error("PlacementDB::getCellId: unknown cell '" + name + "'");
    return it->second;
}

int PlacementDB::getNetId(const std::string& name) const {
    auto it = net_name_to_id_.find(name);
    if (it == net_name_to_id_.end()) throw std::runtime_error("PlacementDB::getNetId: unknown net '" + name + "'");
    return it->second;
}

double PlacementDB::getNetHPWL(int id) {
    checkId(id, static_cast<int>(nets_.size()), "net");
    return nets_[id].hpwl;
}

Cell& PlacementDB::cell(int id) { checkId(id, static_cast<int>(cells_.size()), "cell"); return cells_[id]; }
const Cell& PlacementDB::cell(int id) const { checkId(id, static_cast<int>(cells_.size()), "cell"); return cells_[id]; }
Net& PlacementDB::net(int id) { checkId(id, static_cast<int>(nets_.size()), "net"); return nets_[id]; }
const Net& PlacementDB::net(int id) const { checkId(id, static_cast<int>(nets_.size()), "net"); return nets_[id]; }
Pin& PlacementDB::pin(int id) { checkId(id, static_cast<int>(pins_.size()), "pin"); return pins_[id]; }
const Pin& PlacementDB::pin(int id) const { checkId(id, static_cast<int>(pins_.size()), "pin"); return pins_[id]; }

Point PlacementDB::pinPosition(int pin_id) const {
    checkId(pin_id, static_cast<int>(pins_.size()), "pin");
    const Pin& p = pins_[pin_id];
    checkId(p.cell_id, static_cast<int>(cells_.size()), "pin cell");
    const Cell& c = cells_[p.cell_id];

    // Bookshelf pin offsets are interpreted relative to the cell center.
    // Orientation is not fully represented in the current database yet, so all
    // cells are evaluated with N-orientation semantics in this stage.
    return Point{c.x + 0.5 * c.width + p.offset_x,
                 c.y + 0.5 * c.height + p.offset_y};
}


Box PlacementDB::coreBounds() const {
    double lx = std::numeric_limits<double>::infinity();
    double ly = std::numeric_limits<double>::infinity();
    double ux = -std::numeric_limits<double>::infinity();
    double uy = -std::numeric_limits<double>::infinity();

    if (!rows_.empty()) {
        for (const Row& row : rows_) {
            lx = std::min(lx, row.x_start);
            ux = std::max(ux, row.x_end);
            ly = std::min(ly, row.y);
            uy = std::max(uy, row.y + row.height);
        }
    } else {
        if (cells_.empty()) {
            throw std::runtime_error("PlacementDB::coreBounds: database has neither rows nor cells");
        }
        for (const Cell& cell : cells_) {
            lx = std::min(lx, cell.x);
            ux = std::max(ux, cell.x + cell.width);
            ly = std::min(ly, cell.y);
            uy = std::max(uy, cell.y + cell.height);
        }
    }

    const Box core{lx, ly, ux, uy};
    if (!core.valid()) {
        throw std::runtime_error("PlacementDB::coreBounds: computed core has non-positive area");
    }
    return core;
}

void PlacementDB::setCellLocation(const std::string& name, double x, double y, bool fixed, const std::string& orientation ) {
    if (!hasCell(name)) {
        std::cerr << "Warning: .pl references unknown cell '" << name << "'\n";
        return;
    }
    if (!std::isfinite(x) ||!std::isfinite(y)) {
        throw std::invalid_argument("PlacementDB::setCellLocation: non-finite coordinate for '" + name + "'");
    }

    Cell& cell =cells_[static_cast<std::size_t>(getCellId(name))];
    cell.x = x;
    cell.y = y;
    cell.is_fixed = fixed;
    cell.orientation = orientation.empty() ? "N" : orientation;
}
/*
void PlacementDB::initializeSimplePlacement() {
    if (cells_.empty()) {
        return;
    }

    const Box core = coreBounds();

    if (!core.valid()) {
        throw std::runtime_error(
            "PlacementDB::initializeSimplePlacement: invalid core"
        );
    }

    std::size_t cells_to_initialize = 0;

    for (const Cell& cell : cells_) {
        if (cell.isMovable()) {
            ++cells_to_initialize;
        }
    }

    if (cells_to_initialize == 0) {
        return;
    }

    const double aspect_ratio =
        core.width() / core.height();

    const std::size_t num_columns =
        std::max<std::size_t>(
            1,
            static_cast<std::size_t>(
                std::ceil(
                    std::sqrt(
                        static_cast<double>(
                            cells_to_initialize
                        ) * aspect_ratio
                    )
                )
            )
        );

    const std::size_t num_rows =
        static_cast<std::size_t>(
            std::ceil(
                static_cast<double>(
                    cells_to_initialize
                ) /
                static_cast<double>(num_columns)
            )
        );

    const double step_x =
        core.width() /
        static_cast<double>(num_columns);

    const double step_y =
        core.height() /
        static_cast<double>(num_rows);

    std::size_t placement_index = 0;

    for (Cell& cell : cells_) {
        // Preserve real fixed positions if a benchmark provides them.
        if (cell.isMovable()) {
            continue;
        }

        if (cell.width > core.width() ||
            cell.height > core.height()) {
            throw std::runtime_error(
                "PlacementDB::initializeSimplePlacement: "
                "cell is larger than the core: " +
                cell.name
            );
        }

        const std::size_t column =
            placement_index % num_columns;

        const std::size_t row =
            placement_index / num_columns;

        const double center_x =
            core.lx +
            (static_cast<double>(column) + 0.5) *
                step_x;

        const double center_y =
            core.ly +
            (static_cast<double>(row) + 0.5) *
                step_y;

        const double maximum_x =
            core.ux - cell.width;

        const double maximum_y =
            core.uy - cell.height;

        cell.x = std::clamp(
            center_x - 0.5 * cell.width,
            core.lx,
            maximum_x
        );

        cell.y = std::clamp(
            center_y - 0.5 * cell.height,
            core.ly,
            maximum_y
        );

        ++placement_index;
    }
}
*/
void PlacementDB::initializeSimplePlacement() {
    constexpr double kEpsilon = 1e-9;

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

    if (cells_.empty()) {
        return;
    }

    if (rows_.empty()) {
        throw std::runtime_error(
            "PlacementDB::initializeSimplePlacement: "
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
                    "PlacementDB::initializeSimplePlacement: "
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
     * The original rows_ vector remains unchanged.
     */
    std::vector<const Row*> sorted_rows;
    sorted_rows.reserve(rows_.size());

    for (const Row& row : rows_) {
        if (!std::isfinite(row.x_start) ||
            !std::isfinite(row.x_end) ||
            !std::isfinite(row.y) ||
            !std::isfinite(row.height)) {
            throw std::runtime_error(
                "PlacementDB::initializeSimplePlacement: "
                "placement row contains non-finite geometry"
            );
        }

        if (row.x_end <= row.x_start + kEpsilon ||
            row.height <= kEpsilon) {
            throw std::runtime_error(
                "PlacementDB::initializeSimplePlacement: "
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
        for (const Cell& cell : cells_) {
            if (!cell.is_fixed) {
                continue;
            }

            if (!std::isfinite(cell.x) ||
                !std::isfinite(cell.y) ||
                !std::isfinite(cell.width) ||
                !std::isfinite(cell.height)) {
                throw std::runtime_error(
                    "PlacementDB::initializeSimplePlacement: "
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
    movable_cell_ids.reserve(cells_.size());

    for (const Cell& cell : cells_) {
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
                cells_[static_cast<std::size_t>(lhs_id)];

            const Cell& rhs =
                cells_[static_cast<std::size_t>(rhs_id)];

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
        Cell& cell =
            cells_[static_cast<std::size_t>(cell_id)];

        if (!std::isfinite(cell.width) ||
            !std::isfinite(cell.height) ||
            cell.width < 0.0 ||
            cell.height < 0.0) {
            throw std::runtime_error(
                "PlacementDB::initializeSimplePlacement: "
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
                "PlacementDB::initializeSimplePlacement: "
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

        cell.x = selected_x;
        cell.y = selected_row.y;

        selected_slot.next_x =
            selected_x + cell.width;

        insertFreeSlot(selected_slot);

        ++placed_count;
    }

    if (placed_count != movable_cell_ids.size()) {
        throw std::runtime_error(
            "PlacementDB::initializeSimplePlacement: "
            "not all movable cells were initialized"
        );
    }
}
/*
void PlacementDB::printSummary(std::ostream& os) const {
    const auto fixed_count = std::count_if(cells_.begin(), cells_.end(), [](const Cell& c){ return c.is_terminal || c.is_fixed; });
    size_t max_degree = 0;
    for (const auto& n : nets_) max_degree = std::max(max_degree, n.pin_ids.size());
    double min_x = std::numeric_limits<double>::infinity(), min_y = min_x;
    double max_x = -min_x, max_y = -min_x;
    for (const auto& c : cells_) { min_x = std::min(min_x, c.x); min_y = std::min(min_y, c.y); max_x = std::max(max_x, c.x + c.width); max_y = std::max(max_y, c.y + c.height); }
    if (cells_.empty()) min_x = min_y = max_x = max_y = 0.0;

    os << "========== PlacementDB Summary ==========\n"
       << "Number of cells:\n  total: " << cells_.size() << "\n  movable: " << cells_.size() - fixed_count
       << "\n  terminal/fixed: " << fixed_count << "\n\n"
       << "Number of nets:\n  total: " << nets_.size() << "\n\n"
       << "Number of pins:\n  total: " << pins_.size() << "\n  average pins per net: " << std::fixed << std::setprecision(2)
       << (nets_.empty() ? 0.0 : static_cast<double>(pins_.size()) / static_cast<double>(nets_.size()))
       << "\n  max pins per net: " << max_degree << "\n\n"
       << "Rows:\n  total: " << rows_.size() << "\n\n"
       << "Placement region:\n  min x: " << min_x << "\n  max x: " << max_x << "\n  min y: " << min_y << "\n  max y: " << max_y << "\n\n"
       << "Sample cells:\n";
    for (size_t i = 0; i < std::min<size_t>(5, cells_.size()); ++i) os << "  [" << cells_[i].id << "] name=" << cells_[i].name << ", width=" << cells_[i].width << ", height=" << cells_[i].height << ", x=" << cells_[i].x << ", y=" << cells_[i].y << ", fixed=" << (cells_[i].is_fixed ? "true" : "false") << "\n";
    os << "\nSample nets:\n";
    for (size_t i = 0; i < std::min<size_t>(5, nets_.size()); ++i) os << "  [" << nets_[i].id << "] name=" << nets_[i].name << ", degree=" << nets_[i].pin_ids.size() << "\n";
    os << "========================================\n";
}
*/

void PlacementDB::printSummary(std::ostream& os) const {
    constexpr size_t kMaxDump = 100;

    const size_t fixed_count = static_cast<size_t>(
        std::count_if(cells_.begin(), cells_.end(), [](const Cell& c) {
            return  c.is_fixed;
        })
    );

    size_t max_degree = 0;
    for (const auto& n : nets_) {
        max_degree = std::max(max_degree, n.pin_ids.size());
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (const auto& c : cells_) {
        min_x = std::min(min_x, c.x);
        min_y = std::min(min_y, c.y);
        max_x = std::max(max_x, c.x + c.width);
        max_y = std::max(max_y, c.y + c.height);
    }

    if (cells_.empty()) {
        min_x = min_y = max_x = max_y = 0.0;
    }

    const double avg_pins_per_net =
        nets_.empty()
            ? 0.0
            : static_cast<double>(pins_.size()) / static_cast<double>(nets_.size());

    os << std::fixed << std::setprecision(3);

    os << "========== PlacementDB Summary ==========\n";

    os << "\n[Basic Statistics]\n";
    os << "Number of cells:\n";
    os << "  total: " << cells_.size() << "\n";
    os << "  movable: " << cells_.size() - fixed_count << "\n";
    os << "  terminal/fixed: " << fixed_count << "\n";

    os << "\nNumber of nets:\n";
    os << "  total: " << nets_.size() << "\n";

    os << "\nNumber of pins:\n";
    os << "  total: " << pins_.size() << "\n";
    os << "  average pins per net: " << avg_pins_per_net << "\n";
    os << "  max pins per net: " << max_degree << "\n";

    os << "\nRows:\n";
    os << "  total: " << rows_.size() << "\n";

    os << "\nPlacement region:\n";
    os << "  min x: " << min_x << "\n";
    os << "  max x: " << max_x << "\n";
    os << "  min y: " << min_y << "\n";
    os << "  max y: " << max_y << "\n";

    os << "\n========================================\n";
    os << "[Cell Dump]\n";
    os << "Dump first " << std::min(kMaxDump, cells_.size())
       << " / " << cells_.size() << " cells\n\n";

    os << "Format:\n";
    os << "  cell_id | name | width | height | x | y | type | orientation | is_fixed | num_pins\n\n";

    for (size_t i = 0; i < std::min(kMaxDump, cells_.size()); ++i) {
        const Cell& c = cells_[i];

        os << "Cell[" << c.id << "]"
           << " | name=" << c.name
           << " | width=" << c.width
           << " | height=" << c.height
           << " | x=" << c.x
           << " | y=" << c.y
           << " | type=" << static_cast<int>(c.type)
           << " | orientation=" << c.orientation
           << " | is_fixed=" << (c.is_fixed ? "true" : "false")
           << " | num_pins=" << c.pin_ids.size()
           << "\n";
    }

    os << "\n========================================\n";
    os << "[Net Dump]\n";
    os << "Dump first " << std::min(kMaxDump, nets_.size())
       << " / " << nets_.size() << " nets\n\n";

    os << "Format:\n";
    os << "  net_id | name | degree | HPWL |pin_ids  \n\n";

    for (size_t i = 0; i < std::min(kMaxDump, nets_.size()); ++i) {
        const Net& n = nets_[i];

        os << "Net[" << n.id << "]"
           << " | name=" << n.name
           << " | degree=" << n.pin_ids.size()
           << " | HPWL=" << n.hpwl
           << " | pin_ids=[";

        for (size_t j = 0; j < n.pin_ids.size(); ++j) {
            os << n.pin_ids[j];
            if (j + 1 < n.pin_ids.size()) {
                os << ", ";
            }
        }

        os << "]\n";
    }

    os << "\n========================================\n";
    os << "[Pin Dump]\n";
    os << "Dump first " << std::min(kMaxDump, pins_.size())
       << " / " << pins_.size() << " pins\n\n";

    os << "Format:\n";
    os << "  pin_id | cell_id | cell_name | net_id | net_name | offset_x | offset_y | direction\n\n";

    for (size_t i = 0; i < std::min(kMaxDump, pins_.size()); ++i) {
        const Pin& p = pins_[i];

        os << "Pin[" << p.id << "]"
           << " | cell_id=" << p.cell_id
           << " | cell_name=" << p.cell_name
           << " | net_id=" << p.net_id
           << " | net_name=" << p.net_name
           << " | offset_x=" << p.offset_x
           << " | offset_y=" << p.offset_y
           << " | direction=" << p.direction;

        if (p.cell_id >= 0 && static_cast<size_t>(p.cell_id) < cells_.size()) {
            const Cell& c = cells_[static_cast<size_t>(p.cell_id)];
            const Point pos = pinPosition(p.id);
            os << " | cell_x=" << c.x
               << " | cell_y=" << c.y
               << " | estimated_pin_x=" << pos.x
               << " | estimated_pin_y=" << pos.y;
        }

        os << "\n";
    }

    os << "\n========================================\n";
    os << "[Row Dump]\n";
    os << "Dump first " << std::min(kMaxDump, rows_.size())
       << " / " << rows_.size() << " rows\n\n";

    os << "Format:\n";
    os << "  row_id | y | height | site_width | site_spacing | x_start | x_end | num_sites\n\n";

    for (size_t i = 0; i < std::min(kMaxDump, rows_.size()); ++i) {
        const Row& r = rows_[i];

        os << "Row[" << r.id << "]"
           << " | y=" << r.y
           << " | height=" << r.height
           << " | site_width=" << r.site_width
           << " | site_spacing=" << r.site_spacing
           << " | x_start=" << r.x_start
           << " | x_end=" << r.x_end
           << " | num_sites=" << r.num_sites
           << "\n";
    }

    os << "\n========== End of PlacementDB Summary ==========\n";
}

void PlacementDB::checkId(int id, int size, const char* kind) {
    if (id < 0 || id >= size) throw std::out_of_range(std::string("PlacementDB: invalid ") + kind + " id " + std::to_string(id));
}
