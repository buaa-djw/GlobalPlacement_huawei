#include "initializer/QuadraticInitialPlacer.h"

#include "db/PlacementDB.h"
#include "solver/CsrMatrix.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
bool isFiniteVector(const std::vector<double>& v){ for(double x:v) if(!std::isfinite(x)) return false; return true; }
void restoreSimple(PlacementDB& db,const std::vector<double>& sx,const std::vector<double>& sy){ for(const Cell& c:db.cells()) if(c.isMovable()) db.setCellPosition(c.id,sx.at(static_cast<size_t>(c.id)),sy.at(static_cast<size_t>(c.id))); }
void checkNonMovable(const PlacementDB& db,const std::vector<std::pair<double,double>>& saved){ for(const Cell& c:db.cells()) if(!c.isMovable() && (c.x!=saved.at(static_cast<size_t>(c.id)).first || c.y!=saved.at(static_cast<size_t>(c.id)).second)) throw std::runtime_error("QuadraticInitialPlacer: fixed or terminal cell moved: "+c.name); }
}

QuadraticInitialPlacer::QuadraticInitialPlacer(QuadraticInitialPlacerConfig c):config_(c){
    if(c.maximum_net_degree<2) throw std::invalid_argument("QuadraticInitialPlacer: maximum_net_degree must be >= 2");
    if(!std::isfinite(c.anchor_weight_ratio)||c.anchor_weight_ratio<0.0) throw std::invalid_argument("QuadraticInitialPlacer: invalid anchor weight ratio");
    if(!std::isfinite(c.minimum_anchor_weight)||c.minimum_anchor_weight<=0.0) throw std::invalid_argument("QuadraticInitialPlacer: invalid minimum anchor weight");
    if(!std::isfinite(c.blend_factor)||c.blend_factor<0.0||c.blend_factor>1.0) throw std::invalid_argument("QuadraticInitialPlacer: invalid blend factor");
    if(!std::isfinite(c.numerical_epsilon)||c.numerical_epsilon<=0.0) throw std::invalid_argument("QuadraticInitialPlacer: invalid epsilon");
    (void)PreconditionedConjugateGradient(c.solver);
}

InitialPlacementResult QuadraticInitialPlacer::place(PlacementDB& db) const{
    InitialPlacementResult result; result.method="quadratic";
    std::vector<std::pair<double,double>> saved(db.cells().size());
    for(const Cell& c:db.cells()){ saved[static_cast<size_t>(c.id)]={c.x,c.y}; if(c.isMovable()) ++result.movable_cell_count; }
    SimplePackingInitialPlacer simple(config_.simple_reference); InitialPlacementResult simple_result=simple.place(db); result.movable_cell_count=simple_result.movable_cell_count;
    std::vector<double> sx(db.cells().size()), sy(db.cells().size()); for(const Cell& c:db.cells()){ sx[static_cast<size_t>(c.id)]=c.x; sy[static_cast<size_t>(c.id)]=c.y; }
    auto fallback=[&](const std::string& msg){ restoreSimple(db,sx,sy); checkNonMovable(db,saved); result.used_fallback=true; result.message=msg; return result; };
    try{
        if(result.movable_cell_count==0) return fallback("no movable cells; retained simple packing reference");
        std::vector<int> cell_to_var(db.cells().size(),-1), net_to_var(db.nets().size(),-1); int var=0;
        for(const Cell& c:db.cells()) if(c.isMovable()) cell_to_var.at(static_cast<size_t>(c.id))=var++;
        const int M=var;
        for(const Net& n:db.nets()){
            const int degree=static_cast<int>(n.pin_ids.size()); if(degree<=1) continue; if(degree>config_.maximum_net_degree){++result.skipped_high_degree_net_count; continue;}
            bool has_movable=false; for(int pid:n.pin_ids) if(db.cell(db.pin(pid).cell_id).isMovable()) {has_movable=true; break;} if(!has_movable) continue;
            net_to_var.at(static_cast<size_t>(n.id))=var++; ++result.active_net_count;
        }
        if(result.active_net_count==0) return fallback("no active nets; retained simple packing reference");
        std::vector<SparseTriplet> trips; std::vector<double> bx(static_cast<size_t>(var),0.0), by(static_cast<size_t>(var),0.0), initx(static_cast<size_t>(var),0.0), inity(static_cast<size_t>(var),0.0), movable_diag(static_cast<size_t>(M),0.0);
        for(const Cell& c:db.cells()) if(c.isMovable()){ int v=cell_to_var[c.id]; initx[v]=sx[c.id]; inity[v]=sy[c.id]; }
        for(const Net& n:db.nets()){
            int nv=net_to_var.at(static_cast<size_t>(n.id)); if(nv<0) continue; int d=static_cast<int>(n.pin_ids.size()); double w=1.0/static_cast<double>(std::max(1,d-1)); double sumx=0.0,sumy=0.0; int cnt=0;
            for(int pid:n.pin_ids){ const Pin& p=db.pin(pid); const Cell& c=db.cell(p.cell_id); double dx=0.5*c.width+(config_.use_pin_offsets?p.offset_x:0.0); double dy=0.5*c.height+(config_.use_pin_offsets?p.offset_y:0.0); if(c.isMovable()){ int cv=cell_to_var[c.id]; trips.push_back({cv,cv,w}); trips.push_back({nv,nv,w}); trips.push_back({cv,nv,-w}); trips.push_back({nv,cv,-w}); bx[cv]-=w*dx; bx[nv]+=w*dx; by[cv]-=w*dy; by[nv]+=w*dy; movable_diag[cv]+=w; sumx+=sx[c.id]+dx; sumy+=sy[c.id]+dy; } else { Point pos=config_.use_pin_offsets?db.pinPosition(pid):Point{c.x+0.5*c.width,c.y+0.5*c.height}; trips.push_back({nv,nv,w}); bx[nv]+=w*pos.x; by[nv]+=w*pos.y; sumx+=pos.x; sumy+=pos.y; } ++cnt; }
            initx[nv]=sumx/std::max(1,cnt); inity[nv]=sumy/std::max(1,cnt);
        }
        double avg=0.0; for(double d:movable_diag) avg+=d; avg/=static_cast<double>(M); double alpha=std::max(config_.minimum_anchor_weight, config_.anchor_weight_ratio*std::max(avg,1.0));
        for(const Cell& c:db.cells()) if(c.isMovable()){ int cv=cell_to_var[c.id]; trips.push_back({cv,cv,alpha}); bx[cv]+=alpha*sx[c.id]; by[cv]+=alpha*sy[c.id]; }
        CsrMatrix A=CsrMatrix::fromTriplets(var,var,std::move(trips),config_.numerical_epsilon*1e-6); auto diag=A.diagonal(); for(double d:diag) if(!std::isfinite(d)||d<=config_.numerical_epsilon) throw std::runtime_error("illegal matrix diagonal");
        PreconditionedConjugateGradient pcg(config_.solver); auto rx=pcg.solve(A,bx,initx); auto ry=pcg.solve(A,by,inity); result.x_solver_iterations=rx.iterations; result.y_solver_iterations=ry.iterations; result.x_relative_residual=rx.relative_residual; result.y_relative_residual=ry.relative_residual; result.x_converged=rx.converged; result.y_converged=ry.converged;
        if(!isFiniteVector(rx.solution)||!isFiniteVector(ry.solution)) return fallback("non-finite quadratic solution; restored simple packing reference");
        if(config_.require_solver_convergence && (!rx.converged || !ry.converged)) return fallback("solver did not converge; restored simple packing reference");
        Box core=db.coreBounds(); double beta=config_.blend_factor;
        for(const Cell& c:db.cells()) if(c.isMovable()){ int cv=cell_to_var[c.id]; double x=(1.0-beta)*sx[c.id]+beta*rx.solution[cv]; double y=(1.0-beta)*sy[c.id]+beta*ry.solution[cv]; if(config_.project_to_core){ if(c.width>core.width()||c.height>core.height()) throw std::runtime_error("movable cell larger than core: "+c.name); x=std::clamp(x,core.lx,core.ux-c.width); y=std::clamp(y,core.ly,core.uy-c.height); } if(!std::isfinite(x)||!std::isfinite(y)) throw std::runtime_error("non-finite projected coordinate"); db.setCellPosition(c.id,x,y); }
        checkNonMovable(db,saved); result.used_fallback=false; result.message="quadratic connectivity initialization completed"; return result;
    } catch(const std::exception& e){ return fallback(std::string("quadratic numerical fallback: ")+e.what()); }
}
