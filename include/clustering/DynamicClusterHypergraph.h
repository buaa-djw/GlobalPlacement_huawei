#pragma once
#include "clustering/ClusterConnectivityModel.h"
#include "clustering/ClusterNet.h"
#include "clustering/ClusterNode.h"
#include "clustering/ClusterScoreEvaluator.h"
#include "db/PlacementDB.h"
#include <optional>
#include <unordered_map>
#include <vector>
struct BestNeighbor{int cluster=-1; double score=0.0; double internal_degree=0.0;};
class DynamicClusterHypergraph: public ClusterConnectivityModel{public: DynamicClusterHypergraph()=default; explicit DynamicClusterHypergraph(const PlacementDB& db); DynamicClusterHypergraph(std::vector<ClusterNode> c,std::vector<ClusterNet> n); const std::vector<ClusterNode>& clusters()const{return clusters_;} std::vector<ClusterNode>& clusters(){return clusters_;} const std::vector<ClusterNet>& nets()const{return nets_;} std::vector<ClusterNet>& nets(){return nets_;} int activeMovableClusterCount()const; int activeNetCount()const; int activePinCount()const; double activeMovableArea()const; std::vector<int> activeNeighbors(int id)const; int mergeClusters(int j,int k); std::optional<BestNeighbor> findMaximumScoreNeighbor(int id)const; double internalDegree(int j,int k)const override; double externalDegree(int j)const override; bool checkNoDuplicateActiveNetIncidence()const; private: void normalizeAll(); void rebuildIncidentNets(); std::vector<ClusterNode> clusters_; std::vector<ClusterNet> nets_; ClusterScoreEvaluator scorer_;};
class PaperReferencedConnectivityModel: public ClusterConnectivityModel{public: explicit PaperReferencedConnectivityModel(const DynamicClusterHypergraph& h):h_(h){} double internalDegree(int j,int k)const override{return h_.internalDegree(j,k);} double externalDegree(int j)const override{return h_.externalDegree(j);} private: const DynamicClusterHypergraph& h_;};
