#pragma once

#include "db/Cell.h"
#include "db/Net.h"
#include "db/Pin.h"
#include "db/Row.h"
#include "geometry/Point.h"

#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

class PlacementDB {
public:
    int addCell(const std::string& name, double width, double height, bool is_terminal);
    int addNet(const std::string& name);
    void addNetHPWL(int id, double hpwl);
    int addPin(int cell_id, int net_id, double offset_x, double offset_y, const std::string& direction);
    int addRow(const Row& row);

    bool hasCell(const std::string& name) const;
    bool hasNet(const std::string& name) const;

    int getCellId(const std::string& name) const;
    int getNetId(const std::string& name) const;

    Cell& cell(int id);
    const Cell& cell(int id) const;
    Net& net(int id);
    const Net& net(int id) const;
    Pin& pin(int id);
    const Pin& pin(int id) const;

    Point pinPosition(int pin_id) const;

    const std::vector<Cell>& cells() const { return cells_; }
    const std::vector<Net>& nets() const { return nets_; }
    const std::vector<Pin>& pins() const { return pins_; }
    const std::vector<Row>& rows() const { return rows_; }

    void setCellLocation(const std::string& name, double x, double y, bool fixed);
    void printSummary(std::ostream& os) const;

private:
    static void checkId(int id, int size, const char* kind);
    std::vector<Cell> cells_;
    std::vector<Net> nets_;
    std::vector<Pin> pins_;
    std::vector<Row> rows_;
    std::unordered_map<std::string, int> cell_name_to_id_;
    std::unordered_map<std::string, int> net_name_to_id_;
};
