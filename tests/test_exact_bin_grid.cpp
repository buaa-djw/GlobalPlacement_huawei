#include "density/ExactBinGrid.h"
#include "db/PlacementDB.h"
#include <cassert>
#include <cmath>
#include <stdexcept>
static bool near(double a,double b){return std::abs(a-b)<1e-9;} static void row(PlacementDB&db){Row r; r.y=0; r.height=100; r.site_width=1; r.site_spacing=1; r.x_start=0; r.num_sites=100; db.addRow(r);} static int cell(PlacementDB&db,const char*n,double x,double y,double w,double h,bool f=false){int id=db.addCell(n,w,h,CellType::Standard); db.setCellLocation(n,x,y,f); return id;}
int main(){
 {PlacementDB db; row(db); int c=cell(db,"c",10,10,10,10); ExactBinGrid g(db,2,2,0.9); assert(near(g.bin(0).total_area,100)); assert(g.bin(0).overlapping_cell_ids[0]==c); assert(near(g.totalCellArea(),100));}
 {PlacementDB db; row(db); cell(db,"c",40,10,20,10); ExactBinGrid g(db,2,2,0.9); assert(near(g.bin(0).total_area+g.bin(1).total_area,200));}
 {PlacementDB db; row(db); cell(db,"c",40,40,20,20); ExactBinGrid g(db,2,2,0.9); double s=0; for(auto&b:g.bins()) s+=b.total_area; assert(near(s,400)); for(auto&b:g.bins()) assert(near(b.total_area,100));}
 {PlacementDB db; row(db); cell(db,"c",0,0,50,50); ExactBinGrid g(db,2,2,0.9); assert(g.bin(0).overlapping_cell_ids.size()==1); assert(g.bin(1).overlapping_cell_ids.empty()); assert(g.bin(2).overlapping_cell_ids.empty());}
 {PlacementDB db; row(db); cell(db,"m",0,0,10,10,false); cell(db,"f",10,0,10,10,true); ExactBinGrid g(db,1,1,0.9); assert(near(g.totalMovableArea(),100)); assert(near(g.totalFixedArea(),100)); assert(near(g.totalCellArea(),200));}
 {PlacementDB db; row(db); cell(db,"c",-10,-10,20,20); ExactBinGrid g(db,1,1,0.9); assert(near(g.totalCellArea(),100));}
 {PlacementDB db; row(db); int id=cell(db,"c",200,200,10,10); double x=db.cell(id).x,y=db.cell(id).y; ExactBinGrid g(db,2,2,0.9); assert(near(g.totalCellArea(),0)); assert(near(db.cell(id).x,x)&&near(db.cell(id).y,y));}
 {PlacementDB db; row(db); ExactBinGrid g(db,3,7,0.9); assert(near(g.bin(g.numBins()-1).bounds.ux,100)); assert(near(g.bin(g.numBins()-1).bounds.uy,100));}
 {PlacementDB db; row(db); bool ok=false; try{ExactBinGrid g(db,0,1,0.9);}catch(const std::invalid_argument&){ok=true;} assert(ok); ok=false; try{ExactBinGrid g(db,1,-1,0.9);}catch(const std::invalid_argument&){ok=true;} assert(ok);}
 {PlacementDB db; row(db); bool ok=false; try{ExactBinGrid g(db,1,1,0.0);}catch(const std::invalid_argument&){ok=true;} assert(ok); ok=false; try{ExactBinGrid g(db,1,1,1.1);}catch(const std::invalid_argument&){ok=true;} assert(ok);}
}
