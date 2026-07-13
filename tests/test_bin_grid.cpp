#include "db/PlacementDB.h"
#include "grid/BinGrid.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {
constexpr double EPS = 1e-9;
bool near(double a, double b) { return std::abs(a - b) < EPS; }
void expectNear(double a, double b) { assert(near(a, b)); }

void addCoreRow(PlacementDB& db, double lx=0, double ly=0, double ux=100, double h=80) {
    Row r; r.x_start = lx; r.x_end = ux; r.y = ly; r.height = h; r.site_spacing = 1; r.num_sites = static_cast<int>(ux-lx); db.addRow(r);
}
int addCell(
    PlacementDB& db,
    const std::string& name,
    double x,
    double y,
    double width,
    double height,
    CellType type = CellType::Standard,
    bool fixed = false
) {
    const int id = db.addCell(name, width, height, type);
    db.setCellLocation(name, x, y, fixed);
    return id;
}
bool hasCell(const Bin& b, int id) { for (int x : b.cell_ids) if (x == id) return true; return false; }
}

int main() {
    { PlacementDB db; addCoreRow(db); BinGrid g(db,5,4,1.0); expectNear(g.binWidth(),20); expectNear(g.binHeight(),20); assert(g.numBins()==20); expectNear(g.bin(0,0).bounds.lx,0); expectNear(g.bin(0,0).bounds.uy,20); expectNear(g.bin(4,3).bounds.ux,100); expectNear(g.bin(4,3).bounds.uy,80); }
    { PlacementDB db; addCoreRow(db); int c=addCell(db,"c",1,1,5,5); BinGrid g(db,5,4,1.0); expectNear(g.bin(0,0).movable_area,25); assert(hasCell(g.bin(0,0),c)); for(const Bin& b:g.bins()) if(b.id!=0) expectNear(b.movable_area,0); }
    { PlacementDB db; addCoreRow(db); int c=addCell(db,"c",10,10,20,20); BinGrid g(db,5,4,1.0); expectNear(g.bin(0,0).movable_area,100); expectNear(g.bin(1,0).movable_area,100); expectNear(g.bin(0,1).movable_area,100); expectNear(g.bin(1,1).movable_area,100); expectNear(g.totalMovableArea(),400); assert(hasCell(g.bin(1,1),c)); }
    { PlacementDB db; addCoreRow(db,0,0,20,20); int f=addCell(db,"f",0,0,10,10,CellType::Standard,true); BinGrid g(db,1,1,0.5); (void)f; expectNear(g.bin(0).fixed_area,100); expectNear(g.bin(0).movable_area,0); expectNear(g.bin(0).movable_capacity,100); }
    { PlacementDB db; addCoreRow(db,0,0,10,10); addCell(db,"m",0,0,10,10); BinGrid g(db,1,1,0.5); expectNear(g.bin(0).overflow,50); assert(g.overflowBinCount()==1); expectNear(g.totalOverflow(),50); }
    { PlacementDB db; addCoreRow(db,0,0,10,10); addCell(db,"m",0,0,5,5); BinGrid g(db,1,1,0.5); expectNear(g.totalOverflow(),0); }
    { PlacementDB db; addCoreRow(db,0,0,40,20); int c=addCell(db,"m",0,0,20,10); BinGrid g(db,2,1,1.0); assert(hasCell(g.bin(0,0),c)); assert(!hasCell(g.bin(1,0),c)); }
    { PlacementDB db; addCoreRow(db,0,0,10,10); addCell(db,"m",-5,-5,10,10); BinGrid g(db,1,1,1.0); expectNear(g.totalMovableArea(),25); }
    { PlacementDB db; addCoreRow(db,0,0,10,10); addCell(db,"m",20,20,5,5); BinGrid g(db,1,1,1.0); expectNear(g.totalMovableArea(),0); assert(g.bin(0).cell_ids.empty()); }
    { PlacementDB db; addCoreRow(db,0,0,40,20); addCell(db,"m",0,0,10,10); BinGrid g(db,2,1,1.0); expectNear(g.bin(0,0).movable_area,100); db.setCellLocation("m",25,0,false); g.rebuild(db); expectNear(g.bin(0,0).movable_area,0); expectNear(g.bin(1,0).movable_area,100); }
    { PlacementDB db; addCoreRow(db); bool thrown=false; try{BinGrid g(db,0,1,1);}catch(const std::exception&){thrown=true;} assert(thrown); thrown=false; try{BinGrid g(db,1,0,1);}catch(const std::exception&){thrown=true;} assert(thrown); thrown=false; try{BinGrid g(db,1,1,0);}catch(const std::exception&){thrown=true;} assert(thrown); thrown=false; try{BinGrid g(db,1,1,1.1);}catch(const std::exception&){thrown=true;} assert(thrown); PlacementDB empty; thrown=false; try{BinGrid g(empty,1,1,1);}catch(const std::exception&){thrown=true;} assert(thrown); }
    { PlacementDB db; addCoreRow(db,0,0,20,20); addCell(db,"m1",0,0,10,10); addCell(db,"m2",10,10,20,20); addCell(db,"f",0,10,10,10,CellType::Standard,true); BinGrid g(db,2,2,1.0); expectNear(g.totalMovableArea(),200); expectNear(g.totalFixedArea(),100); }
    return 0;
}
