#include "db/PlacementDB.h"
#include "evaluator/HPWLEvaluator.h"
#include "wirelength/ExactWirelengthSubgradient.h"

#include <cassert>
#include <cmath>
#include <random>
#include <vector>

static bool near(double a, double b, double eps = 1e-9) { return std::abs(a - b) <= eps * std::max(1.0, std::abs(b)); }
static double rel(double a, double b) { return std::abs(a - b) / std::max(1.0, std::abs(b)); }
static int cell(PlacementDB& db, const char* n, double x, double y, double w = 0, double h = 0, bool f = false, CellType t = CellType::Standard) { int id = db.addCell(n, w, h, t); db.setCellLocation(n, x, y, f); return id; }
static void checkHpwl(const PlacementDB& db) { std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.hpwl, HPWLEvaluator().totalHPWL(db))); }

static PlacementDB netFromCoords(const std::vector<double>& xs)
{
    PlacementDB db; int net = db.addNet("n");
    for (std::size_t i = 0; i < xs.size(); ++i) { int c = cell(db, ("c" + std::to_string(i)).c_str(), xs[i], 0); db.addPin(c, net, 0, 0, "I"); }
    return db;
}

int main()
{
    {
        PlacementDB db; int a = cell(db, "a", 0, 0); int n = db.addNet("n"); db.addPin(a, n, 0, 0, "I");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng);
        assert(near(r.hpwl, 0)); assert(near(r.grad_x[a], 0)); assert(near(r.grad_y[a], 0)); checkHpwl(db);
    }
    {
        PlacementDB db; int a = cell(db, "a", 0, 0); int b = cell(db, "b", 10, 20); int n = db.addNet("n"); db.addPin(a, n, 0, 0, "I"); db.addPin(b, n, 0, 0, "O");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng);
        assert(near(r.hpwl, 30)); assert(near(r.grad_x[a], -1)); assert(near(r.grad_x[b], 1)); assert(near(r.grad_y[a], -1)); assert(near(r.grad_y[b], 1)); checkHpwl(db);
    }
    {
        PlacementDB db; int n = db.addNet("n"); for (int i = 0; i < 3; ++i) db.addPin(cell(db, ("c"+std::to_string(i)).c_str(), i*10, i*7), n, 0, 0, "I");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng);
        assert(near(r.hpwl, 34)); assert(near(r.grad_x[0]+r.grad_x[1]+r.grad_x[2], 0)); assert(near(r.grad_y[0]+r.grad_y[1]+r.grad_y[2], 0)); checkHpwl(db);
    }
    {
        PlacementDB db; int n = db.addNet("n"); double xs[] = {0, 10, 20, 30}; double ys[] = {30, 0, 20, 10}; for (int i=0;i<4;++i) db.addPin(cell(db,("c"+std::to_string(i)).c_str(),xs[i],ys[i]), n, 0, 0, "I");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng);
        assert(near(r.hpwl, 60)); assert(near(r.grad_x[0], -1)); assert(near(r.grad_x[3], 1)); assert(near(r.grad_y[1], -1)); assert(near(r.grad_y[0], 1)); checkHpwl(db);
    }
    {
        PlacementDB db; int n = db.addNet("n"); for (int i=0;i<10;++i) db.addPin(cell(db,("c"+std::to_string(i)).c_str(),i*3,i*2), n, 0, 0, "I");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.hpwl, 45)); checkHpwl(db);
    }
    for (auto xs : {std::vector<double>{0,0,1,1}, {0,0,0,1}, {0,1,1,1}, {0,0,0,0}}) { auto db = netFromCoords(xs); std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.hpwl, *std::max_element(xs.begin(),xs.end())-*std::min_element(xs.begin(),xs.end()))); checkHpwl(db); }
    {
        PlacementDB db; int a = cell(db, "a", 10, 0, 10, 0); int b = cell(db, "b", 0, 0, 10, 0); int n = db.addNet("n"); db.addPin(a, n, -10, 0, "I"); db.addPin(b, n, 20, 0, "O");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.hpwl, 20)); assert(near(r.grad_x[a], -1)); assert(near(r.grad_x[b], 1)); checkHpwl(db);
    }
    {
        PlacementDB db; int a = cell(db, "a", 0, 0, 0, 0); int n = db.addNet("n"); db.addPin(a, n, 0, 0, "I"); db.addPin(a, n, 5, 0, "O");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.hpwl, 5)); assert(near(r.grad_x[a], 0)); checkHpwl(db);
    }
    {
        PlacementDB db; int m = cell(db, "m", 0, 0); int f = cell(db, "f", 10, 0, 0, 0, true); int n = db.addNet("n"); db.addPin(m, n, 0, 0, "I"); db.addPin(f, n, 0, 0, "O");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.hpwl, 10)); assert(near(r.grad_x[m], -1)); assert(near(r.grad_x[f], 0)); checkHpwl(db);
    }
    {
        PlacementDB db; int a = cell(db, "a", 0, 0); int b = cell(db, "b", 0, 0); int n = db.addNet("n"); db.addPin(a, n, 0, 0, "I"); db.addPin(b, n, 0, 0, "O");
        std::mt19937_64 rng(1); auto r = ExactWirelengthSubgradient().evaluate(db, nullptr, rng); assert(near(r.grad_x[a], 1)); assert(near(r.grad_x[b], -1));
        WirelengthCoordinateSnapshot prev; prev.pin_x={2,0}; prev.pin_y={0,0}; std::mt19937_64 r1(7), r2(7); auto a1=ExactWirelengthSubgradient().evaluate(db,&prev,r1); auto a2=ExactWirelengthSubgradient().evaluate(db,&prev,r2); assert(a1.grad_x[a]>=0.5 && a1.grad_x[a]<=1.0); assert(a1.grad_x[a]==a2.grad_x[a]);
        prev.pin_x={0,2}; std::mt19937_64 r3(8), r4(8); auto b1=ExactWirelengthSubgradient().evaluate(db,&prev,r3); auto b2=ExactWirelengthSubgradient().evaluate(db,&prev,r4); assert(b1.grad_x[a]>=-1.0 && b1.grad_x[a]<=-0.5); assert(b1.grad_x[a]==b2.grad_x[a]);
        prev.pin_x={0,0}; std::mt19937_64 r5(9); auto e=ExactWirelengthSubgradient().evaluate(db,&prev,r5); assert(near(e.grad_x[a],1));
    }
    {
        PlacementDB db; int a=cell(db,"a",0,0), b=cell(db,"b",0,0); int n=db.addNet("n"); db.addPin(a,n,0,0,"I"); db.addPin(b,n,0,0,"O"); WirelengthCoordinateSnapshot prev; prev.pin_x={2,0}; prev.pin_y={0,0}; std::mt19937_64 r1(1), r2(2); auto x=ExactWirelengthSubgradient().evaluate(db,&prev,r1); auto y=ExactWirelengthSubgradient().evaluate(db,&prev,r2); assert(x.hpwl==y.hpwl);
    }
    {
        PlacementDB db; int n=db.addNet("n"); int c0=cell(db,"c0",1,2), c1=cell(db,"c1",7,11), c2=cell(db,"c2",19,23), c3=cell(db,"c3",31,37); for(int c:{c0,c1,c2,c3}) db.addPin(c,n,0,0,"I"); std::mt19937_64 rng(1); auto r=ExactWirelengthSubgradient().evaluate(db,nullptr,rng); double eps=1e-6;
        for(int c:{c0,c1,c2,c3}) { double ox=db.cell(c).x, oy=db.cell(c).y; db.setCellPosition(c,ox+eps,oy); double pp=HPWLEvaluator().totalHPWL(db); db.setCellPosition(c,ox-eps,oy); double pm=HPWLEvaluator().totalHPWL(db); db.setCellPosition(c,ox,oy); assert(rel((pp-pm)/(2*eps), r.grad_x[c])<=1e-5); }
    }
}
