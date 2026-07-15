#include "db/PlacementDB.h"

#include <cassert>
#include <sstream>

int main()
{
    PlacementDB db;
    int c0 = db.addCell("a", 10, 20, CellType::Standard);
    int c1 = db.addCell("b", 5, 5, CellType::Terminal);
    int n0 = db.addNet("n");
    int p0 = db.addPin(c0, n0, 1.5, -2.0, "I");
    db.addPin(c1, n0, 0.0, 0.0, "O");
    Row row;
    row.y = 0;
    row.height = 10;
    row.site_width = 1;
    row.site_spacing = 1;
    row.x_start = 0;
    row.num_sites = 100;
    db.addRow(row);
    db.setCellLocation("a", 3, 4, false);
    assert(db.getCellId("a") == c0);
    assert(db.cell(c0).pin_ids.size() == 1);
    assert(db.net(n0).pin_ids.size() == 2);
    assert(db.pin(p0).cell_name == "a");
    assert(db.rows().size() == 1);
    std::ostringstream os;
    db.printSummary(os);
    assert(os.str().find("Number of cells") != std::string::npos);
    return 0;
}
