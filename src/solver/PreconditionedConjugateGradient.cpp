#include "solver/PreconditionedConjugateGradient.h"
#include "solver/CsrMatrix.h"

#include <cmath>
#include <stdexcept>

namespace { double dot(const std::vector<double>& a,const std::vector<double>& b){ double s=0.0; for(size_t i=0;i<a.size();++i) s+=a[i]*b[i]; return s;} double norm(const std::vector<double>& a){ return std::sqrt(dot(a,a)); } bool finite(const std::vector<double>& v){ for(double x:v) if(!std::isfinite(x)) return false; return true; } }
PreconditionedConjugateGradient::PreconditionedConjugateGradient(ConjugateGradientConfig c):config_(c){ if(c.maximum_iterations<0||!std::isfinite(c.relative_tolerance)||c.relative_tolerance<0.0||!std::isfinite(c.absolute_tolerance)||c.absolute_tolerance<0.0||!std::isfinite(c.numerical_epsilon)||c.numerical_epsilon<=0.0) throw std::invalid_argument("PreconditionedConjugateGradient: invalid config"); }
ConjugateGradientResult PreconditionedConjugateGradient::solve(const CsrMatrix& A,const std::vector<double>& b,const std::vector<double>& x0) const{
    if(A.rows()!=A.columns()) throw std::invalid_argument("PCG: matrix must be square"); const int n=A.rows(); if(static_cast<int>(b.size())!=n||static_cast<int>(x0.size())!=n) throw std::invalid_argument("PCG: dimension mismatch");
    ConjugateGradientResult res; res.solution=x0; if(n==0){res.converged=true;res.termination_reason="empty system";return res;} if(!finite(b)||!finite(x0)) throw std::invalid_argument("PCG: non-finite input");
    auto diag=A.diagonal(); for(double d:diag) if(!std::isfinite(d)||d<=config_.numerical_epsilon) throw std::runtime_error("PCG: matrix diagonal must be finite positive");
    std::vector<double> Ax,r(n),z(n),p(n),Ap; A.multiply(res.solution,Ax); for(int i=0;i<n;++i){r[i]=b[i]-Ax[i]; z[i]=r[i]/diag[i]; p[i]=z[i];}
    double rho=dot(r,z); if(!std::isfinite(rho)) throw std::runtime_error("PCG: non-finite rho"); res.initial_residual_norm=norm(r); res.final_residual_norm=res.initial_residual_norm; res.relative_residual=res.initial_residual_norm/std::max(res.initial_residual_norm,config_.numerical_epsilon);
    auto converged=[&](){return res.final_residual_norm<=config_.absolute_tolerance||res.final_residual_norm/std::max(res.initial_residual_norm,config_.numerical_epsilon)<=config_.relative_tolerance;}; if(converged()){res.converged=true;res.termination_reason="initial solution converged";return res;}
    for(int it=0; it<config_.maximum_iterations; ++it){ A.multiply(p,Ap); double den=dot(p,Ap); if(!std::isfinite(den)||den<=config_.numerical_epsilon) throw std::runtime_error("PCG: invalid denominator"); double alpha=rho/den; if(!std::isfinite(alpha)) throw std::runtime_error("PCG: invalid alpha"); for(int i=0;i<n;++i){res.solution[i]+=alpha*p[i]; r[i]-=alpha*Ap[i];} if(!finite(res.solution)||!finite(r)) throw std::runtime_error("PCG: non-finite iterate"); res.iterations=it+1; res.final_residual_norm=norm(r); res.relative_residual=res.final_residual_norm/std::max(res.initial_residual_norm,config_.numerical_epsilon); if(converged()){res.converged=true;res.termination_reason="converged";return res;} for(int i=0;i<n;++i) z[i]=r[i]/diag[i]; double rho_new=dot(r,z); if(!std::isfinite(rho_new)) throw std::runtime_error("PCG: non-finite rho"); double beta=rho_new/rho; if(!std::isfinite(beta)) throw std::runtime_error("PCG: non-finite beta"); for(int i=0;i<n;++i) p[i]=z[i]+beta*p[i]; rho=rho_new; }
    res.termination_reason="maximum iterations reached"; return res;
}
