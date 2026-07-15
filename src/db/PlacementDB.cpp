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
void PlacementDB::setCellPosition(int cell_id, double x, double y) {
    checkId(cell_id, static_cast<int>(cells_.size()), "cell");
    if (!std::isfinite(x) || !std::isfinite(y)) {
        throw std::invalid_argument("PlacementDB::setCellPosition: non-finite coordinate");
    }
    cells_[static_cast<std::size_t>(cell_id)].x = x;
    cells_[static_cast<std::size_t>(cell_id)].y = y;
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

void PlacementDB::printSummary(std::ostream& os) const
{
    const std::size_t movable_count =
        static_cast<std::size_t>(
            std::count_if(
                cells_.begin(),
                cells_.end(),
                [](const Cell& cell) {
                    return cell.isMovable();
                }
            )
        );

    const std::size_t fixed_flag_count =
        static_cast<std::size_t>(
            std::count_if(
                cells_.begin(),
                cells_.end(),
                [](const Cell& cell) {
                    return cell.is_fixed;
                }
            )
        );

    const std::size_t terminal_count =
        static_cast<std::size_t>(
            std::count_if(
                cells_.begin(),
                cells_.end(),
                [](const Cell& cell) {
                    return cell.type ==
                           CellType::Terminal;
                }
            )
        );

    const std::size_t terminal_ni_count =
        static_cast<std::size_t>(
            std::count_if(
                cells_.begin(),
                cells_.end(),
                [](const Cell& cell) {
                    return cell.type ==
                           CellType::TerminalNI;
                }
            )
        );

    std::size_t max_degree = 0;

    for (const Net& net : nets_) {
        max_degree =
            std::max(
                max_degree,
                net.pin_ids.size()
            );
    }

    os << "========== PlacementDB Summary ==========\n"
       << "cells " << cells_.size() << '\n'
       << "movable_cells " << movable_count << '\n'
       << "non_movable_cells "
       << cells_.size() - movable_count << '\n'
       << "fixed_flag_cells "
       << fixed_flag_count << '\n'
       << "terminal_cells "
       << terminal_count << '\n'
       << "terminal_NI_cells "
       << terminal_ni_count << '\n'
       << "nets " << nets_.size() << '\n'
       << "pins " << pins_.size() << '\n'
       << "rows " << rows_.size() << '\n'
       << "max_net_degree " << max_degree << '\n'
       << "========================================\n";
}

void PlacementDB::checkId(int id, int size, const char* kind) {
    if (id < 0 || id >= size) throw std::out_of_range(std::string("PlacementDB: invalid ") + kind + " id " + std::to_string(id));
}
