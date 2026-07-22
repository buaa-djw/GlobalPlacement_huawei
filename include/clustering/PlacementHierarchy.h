#pragma once
#include "clustering/HierarchyLevel.h"
#include "db/PlacementDB.h"
#include <vector>
class PlacementHierarchy{public: void build(const PlacementDB& finest,int current=150); int numLevels()const{return (int)levels_.size();} int coarsestLevelIndex()const{return (int)levels_.size()-1;} const HierarchyLevel& level(int i)const{return levels_.at(i);} const std::vector<HierarchyLevel>& levels()const{return levels_;} private: std::vector<HierarchyLevel> levels_;};
