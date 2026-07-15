#include "density/ExactBinGrid.h"
#include "density/ExactDensityEvaluator.h"
#include "db/PlacementDB.h"
#include <cassert>
#include <cmath>
static bool near(double a,double b){return std::abs(a-b)<1e-9;} static void row(PlacementDB&db){Row r; r.y=0; r.height=10; r.site_width=1; r.site_spacing=1; r.x_start=0; r.num_sites=20; db.addRow(r);} static void cell(PlacementDB&db,const char*n,double x,double y,double w,double h,bool f=false){db.addCell(n,w,h,CellType::Standard); db.setCellLocation(n,x,y,f);}
int main(){
 {PlacementDB db; row(db); cell(db,"c",0,0,2,2); ExactBinGrid g(db,1,1,0.9); auto m=ExactDensityEvaluator().evaluate(g); assert(near(m.total_overflow,0)); assert(near(m.quadratic_penalty,0)); assert(near(m.overflow_ratio,0));}
 {PlacementDB db; row(db); cell(db,"c",0,0,20,10); ExactBinGrid g(db,1,1,0.5); auto m=ExactDensityEvaluator().evaluate(g); assert(near(m.total_cell_area,200)); assert(near(m.total_capacity,100)); assert(near(m.total_overflow,100)); assert(near(m.quadratic_penalty,10000)); assert(near(m.overflow_ratio,0.5));}
 {PlacementDB db; row(db); cell(db,"a",0,0,10,10); cell(db,"b",10,0,10,10); ExactBinGrid g(db,2,1,0.5); auto m=ExactDensityEvaluator().evaluate(g); assert(near(m.total_overflow,100)); assert(near(m.quadratic_penalty,5000));}
 {PlacementDB db; row(db); cell(db,"f",0,0,10,10,true); cell(db,"m",0,0,10,10,false); ExactBinGrid g(db,1,1,0.5); auto m=ExactDensityEvaluator().evaluate(g); assert(near(m.total_cell_area,200)); assert(near(m.total_overflow,100));}
 {PlacementDB db; row(db); cell(db,"f",0,0,10,10,true); cell(db,"m",0,0,10,10,false); ExactBinGrid g(db,1,1,0.5); auto m=ExactDensityEvaluator().evaluate(g); assert(near(m.overflow_ratio,1.0));}
 {PlacementDB db; row(db); cell(db,"f",0,0,20,10,true); ExactBinGrid g(db,1,1,0.5); auto m=ExactDensityEvaluator().evaluate(g); assert(near(m.overflow_ratio,0)); assert(near(m.total_overflow,100));}
 {PlacementDB db; row(db); cell(db,"a",-5,0,10,10); cell(db,"b",10,0,10,10); ExactBinGrid g(db,2,1,0.9); double s=0; for(auto&b:g.bins()) s+=b.total_area; assert(near(s,150)); assert(near(s,g.totalCellArea()));}
}
