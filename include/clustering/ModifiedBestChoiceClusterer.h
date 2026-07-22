#pragma once
#include "clustering/ClusterCandidate.h"
#include "clustering/DynamicClusterHypergraph.h"
#include <vector>
struct ClustererStats{int input_movable_clusters=0,target_movable_clusters=0,output_movable_clusters=0,positive_score_merges=0,stale_candidates=0,invalid_recomputations=0; double top_score_at_stop=0.0; std::string stop_reason; std::vector<std::pair<int,int>> merge_sequence;};
class ModifiedBestChoiceClusterer{public: ClustererStats cluster(DynamicClusterHypergraph& h,int target); std::vector<ClusterCandidate> initialCandidates(DynamicClusterHypergraph& h)const;};
