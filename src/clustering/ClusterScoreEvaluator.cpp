#include "clustering/ClusterScoreEvaluator.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
namespace{double term(double area,double e,double d,int id){if(!std::isfinite(area)||area<=0) throw std::invalid_argument("ClusterScoreEvaluator: non-positive area"); if(d<=0) return 0.0; if(e+1e-12<d) throw std::runtime_error("ClusterScoreEvaluator: external degree smaller than internal degree for cluster "+std::to_string(id)); double den=e-d; if(den==0.0) return std::numeric_limits<double>::infinity(); return (1.0/area)*d/den;}}
double ClusterScoreEvaluator::evaluate(const ClusterConnectivity& v) const{ if(v.internal_degree_jk<=0) return 0.0; return term(v.area_j,v.external_degree_j,v.internal_degree_jk,v.cluster_j)+term(v.area_k,v.external_degree_k,v.internal_degree_jk,v.cluster_k);}
