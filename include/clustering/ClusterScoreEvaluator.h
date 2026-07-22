#pragma once
struct ClusterConnectivity{double area_j=0, area_k=0, external_degree_j=0, external_degree_k=0, internal_degree_jk=0; int cluster_j=-1, cluster_k=-1;};
class ClusterScoreEvaluator{public: double evaluate(const ClusterConnectivity& value) const;};
