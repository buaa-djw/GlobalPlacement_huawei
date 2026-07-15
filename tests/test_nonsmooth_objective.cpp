#include "db/PlacementDB.h"
#include "density/ExactBinGrid.h"
#include "density/ExactDensitySubgradient.h"
#include "evaluator/HPWLEvaluator.h"
#include "objective/NonsmoothObjectiveEvaluator.h"
#include "wirelength/ExactWirelengthSubgradient.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

static bool near(double a, double b, double eps = 1e-9) { return std::abs(a - b) <= eps * std::max(1.0, std::abs(b)); }
static double rel(double a, double b) { return std::abs(a - b) / std::max(1.0, std::abs(b)); }
static void row(PlacementDB& db, double w = 20, double h = 10) { Row r; r.y = 0; r.height = h; r.site_width = 1; r.site_spacing = 1; r.x_start = 0; r.num_sites = static_cast<int>(w); db.addRow(r); }
static int cell(PlacementDB& db, const char* n, double x, double y, double w, double h, bool f = false, CellType t = CellType::Standard) { int id = db.addCell(n, w, h, t); db.setCellLocation(n, x, y, f); return id; }
static void pin(PlacementDB& db, int c, int n) { db.addPin(c, n, 0, 0, "I"); }

static NonsmoothObjectiveResult eval(PlacementDB& db, double lambda, int bx=2, int by=1, double td=0.1)
{
    ExactBinGrid grid(db, bx, by, td); std::mt19937_64 rng(1); return NonsmoothObjectiveEvaluator().evaluate(db, grid, lambda, nullptr, rng);
}

int main()
{
    {
        PlacementDB db; row(db); int a=cell(db,"a",1,1,2,2), b=cell(db,"b",12,1,2,2); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n);
        ExactBinGrid grid(db,2,1,0.1); std::mt19937_64 r1(1), r2(1); auto obj=NonsmoothObjectiveEvaluator().evaluate(db,grid,0,nullptr,r1); auto w=ExactWirelengthSubgradient().evaluate(db,nullptr,r2);
        assert(near(obj.objective,obj.hpwl)); assert(near(obj.grad_x[a],w.grad_x[a])); assert(near(obj.grad_y[b],w.grad_y[b])); assert(obj.density_penalty >= 0); assert(obj.ofr >= 0);
    }
    {
        PlacementDB db; row(db); int a=cell(db,"a",1,1,6,4), b=cell(db,"b",9,2,6,4); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n); ExactBinGrid grid(db,2,1,0.05); double lambda=3.25;
        std::mt19937_64 rw(5), rd(5), ro(5); auto w=ExactWirelengthSubgradient().evaluate(db,nullptr,rw); auto d=ExactDensitySubgradient().evaluate(db,ExactBinGrid(db,2,1,0.05),rd); auto o=NonsmoothObjectiveEvaluator().evaluate(db,grid,lambda,nullptr,ro);
        assert(near(o.objective,w.hpwl+lambda*d.penalty)); assert(near(o.density_penalty,d.penalty)); assert(near(o.ofr,d.ofr)); for (std::size_t i=0;i<db.cells().size();++i) { assert(near(o.grad_x[i], w.grad_x[i]+lambda*d.grad_x[i])); assert(near(o.grad_y[i], w.grad_y[i]+lambda*d.grad_y[i])); }
    }
    {
        PlacementDB db; row(db); int m=cell(db,"m",1,1,2,2), f=cell(db,"f",10,1,2,2,true); int n=db.addNet("n"); pin(db,m,n); pin(db,f,n); auto o=eval(db,2.0); assert(o.hpwl>0); assert(o.ofr>=0); assert(near(o.grad_x[f],0)); assert(near(o.grad_y[f],0)); assert(!near(o.grad_x[m],0));
    }
    {
        PlacementDB db; row(db,100,100); int a=cell(db,"a",1,1,1,1), b=cell(db,"b",10,10,1,1); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n); auto o=eval(db,4.0,1,1,0.9); std::mt19937_64 rng(1); auto w=ExactWirelengthSubgradient().evaluate(db,nullptr,rng); assert(near(o.density_penalty,0)); assert(near(o.objective,o.hpwl)); assert(near(o.grad_x[a],w.grad_x[a]));
    }
    {
        PlacementDB db; row(db,10,10); int a=cell(db,"a",1,1,8,8); int n=db.addNet("n"); pin(db,a,n); auto o=eval(db,2.0,1,1,0.1); std::mt19937_64 rng(1); auto d=ExactDensitySubgradient().evaluate(db,ExactBinGrid(db,1,1,0.1),rng); assert(near(o.hpwl,0)); assert(near(o.objective,2.0*d.penalty)); assert(near(o.grad_x[a],2.0*d.grad_x[a]));
    }
    {
        PlacementDB db; row(db); int a=cell(db,"a",0,1,2,2), b=cell(db,"b",0,1,2,2); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n); ExactWirelengthSubgradient wl; auto prev=wl.capture(db); prev.pin_x[0]=2; ExactBinGrid grid(db,2,1,0.1); std::mt19937_64 r1(9), r2(9); auto x=NonsmoothObjectiveEvaluator().evaluate(db,grid,1.0,&prev,r1); auto y=NonsmoothObjectiveEvaluator().evaluate(db,grid,1.0,&prev,r2); assert(x.objective==y.objective); assert(x.grad_x[a]==y.grad_x[a]);
    }
    {
        PlacementDB db; row(db); int a=cell(db,"a",1,1,2,2); ExactBinGrid grid(db,2,1,0.1); for(double l:{-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}) { std::mt19937_64 rng(1); bool thrown=false; try { (void)NonsmoothObjectiveEvaluator().evaluate(db,grid,l,nullptr,rng); } catch(const std::invalid_argument&) { thrown=true; } assert(thrown); }
    }
    {
        PlacementDB db; row(db); int a=cell(db,"a",1,1,6,4), b=cell(db,"b",15,1,2,2); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n); auto before=eval(db,1.0); db.setCellPosition(a,3,1); auto after=eval(db,1.0); assert(!near(before.objective,after.objective)); assert(near(db.cell(b).x,15));
    }
    {
        PlacementDB db; row(db,20,20); int a=cell(db,"a",2.3,2.7,3.1,4.2), b=cell(db,"b",14.1,12.2,2.0,2.0); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n); double lambda=0.7; auto r=eval(db,lambda,1,1,0.02); double eps=1e-6, ox=db.cell(a).x, oy=db.cell(a).y; db.setCellPosition(a,ox+eps,oy); double pp=eval(db,lambda,1,1,0.02).objective; db.setCellPosition(a,ox-eps,oy); double pm=eval(db,lambda,1,1,0.02).objective; db.setCellPosition(a,ox,oy); assert(rel((pp-pm)/(2*eps), r.grad_x[a])<=1e-5);
    }
    {
        PlacementDB db; row(db); int a=cell(db,"a",1,1,6,4), b=cell(db,"b",9,2,6,4); int n=db.addNet("n"); pin(db,a,n); pin(db,b,n); ExactBinGrid grid(db,2,1,0.05); std::mt19937_64 r1(3), r2(3), r3(3); auto o=NonsmoothObjectiveEvaluator().evaluate(db,grid,1,nullptr,r1); auto w=ExactWirelengthSubgradient().evaluate(db,nullptr,r2); auto d=ExactDensitySubgradient().evaluate(db,ExactBinGrid(db,2,1,0.05),r3); assert(near(o.hpwl,w.hpwl)); assert(near(o.hpwl,HPWLEvaluator().totalHPWL(db))); assert(near(o.density_penalty,d.penalty)); assert(near(o.ofr,d.ofr));
    }
}
