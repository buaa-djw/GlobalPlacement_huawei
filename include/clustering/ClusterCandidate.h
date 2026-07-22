#pragma once
#include <cstdint>
struct ClusterCandidate{int first_cluster=-1; int second_cluster=-1; int owner_cluster=-1; int partner_cluster=-1; double score=0.0; std::uint64_t first_version=0; std::uint64_t second_version=0;};
struct ClusterCandidateLess{bool operator()(const ClusterCandidate&a,const ClusterCandidate&b)const{if(a.score!=b.score) return a.score<b.score; if(a.first_cluster!=b.first_cluster) return a.first_cluster>b.first_cluster; return a.second_cluster>b.second_cluster;}};
