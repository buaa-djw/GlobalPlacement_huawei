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
