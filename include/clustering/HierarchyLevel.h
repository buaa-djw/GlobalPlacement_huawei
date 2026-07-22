#pragma once
#include "clustering/ClusterNet.h"
#include "clustering/ClusterNode.h"
#include "clustering/ModifiedBestChoiceClusterer.h"
#include <vector>
struct HierarchyLevel{int level_index=0; std::vector<ClusterNode> clusters; std::vector<ClusterNet> nets; std::vector<int> fine_to_coarse; std::vector<std::vector<int>> coarse_to_fine; ClustererStats stats; double area_conservation_error=0.0; bool mapping_complete=true; bool fixed_preserved=true;};
