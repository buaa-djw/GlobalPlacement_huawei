#pragma once
class ClusterConnectivityModel{public: virtual ~ClusterConnectivityModel()=default; virtual double internalDegree(int cluster_j,int cluster_k) const=0; virtual double externalDegree(int cluster_j) const=0;};
