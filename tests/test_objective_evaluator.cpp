#include "db/PlacementDB.h"
#include "evaluator/HPWLEvaluator.h"
#include "evaluator/ObjectiveEvaluator.h"
#include "grid/BinGrid.h"
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
constexpr double EPS=1e-8;
void row(PlacementDB& db,double ux=100,double uy=100){ Row r; r.x_start=0; r.x_end=ux; r.y=0; r.height=uy; r.site_width=1; r.site_spacing=1; r.num_sites=(int)ux; db.addRow(r);} 
int cell(PlacementDB& db,const char*n,double x,double y,double w,double h,bool fixed=false){int id=db.addCell(n,w,h,CellType::Standard); db.setCellLocation(n,x,y,fixed); return id;}
void pin2(PlacementDB& db,int a,int b){int n=db.addNet("n"+std::to_string(db.nets().size())); db.addPin(a,n,0,0,"B"); db.addPin(b,n,0,0,"B");}
void near(double a,double b){assert(std::abs(a-b)<EPS);} 
}
int main(){
 {PlacementDB db; row(db); int a=cell(db,"a",0,0,1,1), b=cell(db,"b",10,0,1,1); pin2(db,a,b); BinGrid g(db,1,1,1); auto m=ObjectiveEvaluator{}.evaluate(db,g,2); near(m.hpwl,HPWLEvaluator{}.totalHPWL(db)); near(m.density_penalty,0); near(m.weighted_density_penalty,0); near(m.total_cost,m.hpwl);} 
 {PlacementDB db; row(db,10,10); cell(db,"a",0,0,8,10); BinGrid g(db,1,1,0.5); auto m=ObjectiveEvaluator{}.evaluate(db,g,3); assert(m.density_penalty>0); assert(m.total_overflow>0); near(m.total_cost,m.hpwl+3*m.density_penalty);} 
 {PlacementDB db; row(db,10,10); cell(db,"a",0,0,8,10); BinGrid g(db,1,1,0.5); auto m=ObjectiveEvaluator{}.evaluate(db,g,0); near(m.weighted_density_penalty,0); near(m.total_cost,m.hpwl);} 
 {PlacementDB db; row(db,10,10); cell(db,"f",0,0,10,10,true); cell(db,"m",0,0,1,10); BinGrid g(db,1,1,1); auto m=ObjectiveEvaluator{}.evaluate(db,g,1); assert(std::isfinite(m.density_penalty)); assert(std::isfinite(m.total_overflow)); assert(std::isfinite(m.total_cost));}
 {PlacementDB db; row(db); BinGrid g(db,1,1,1); for(double w: {-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}){bool th=false; try{(void)ObjectiveEvaluator{}.evaluate(db,g,w);}catch(const std::invalid_argument&){th=true;} assert(th);} }
 {PlacementDB db; row(db); int a=cell(db,"a",1,2,3,4); db.cell(a).type=CellType::Standard; bool f=db.cell(a).is_fixed; double x=db.cell(a).x,y=db.cell(a).y; CellType t=db.cell(a).type; BinGrid g(db,1,1,1); (void)ObjectiveEvaluator{}.evaluate(db,g,1); near(db.cell(a).x,x); near(db.cell(a).y,y); assert(db.cell(a).is_fixed==f); assert(db.cell(a).type==t);} 
}
