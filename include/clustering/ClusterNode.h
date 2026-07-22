#pragma once
#include <cstdint>
#include <vector>
struct ClusterNode{int id=-1; bool active=true; bool valid=true; bool movable=true; std::uint64_t version=0; double area=0.0; std::vector<int> original_cell_ids; std::vector<int> incident_net_ids;};
