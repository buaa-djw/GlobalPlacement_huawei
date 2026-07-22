#include "clustering/ClusterScoreEvaluator.h"
#include <cassert>
#include <cmath>
#include <limits>
int main(){ClusterScoreEvaluator e; double s=e.evaluate({2,4,3,5,1,0,1}); assert(std::abs(s-(0.5*1/2+0.25*1/4))<1e-12); double a=e.evaluate({2,4,3,5,1}); double b=e.evaluate({20,4,3,5,1}); assert(b<a); assert(e.evaluate({2,4,3,5,0})==0); assert(std::isinf(e.evaluate({2,4,1,5,1}))); bool th=false; try{e.evaluate({0,4,1,5,1});}catch(...){th=true;} assert(th);}
