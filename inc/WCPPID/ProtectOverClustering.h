#ifndef WIRECELLPID_PROTECTOVERCLUSTERING_H
#define WIRECELLPID_PROTECTOVERCLUSTERING_H

#include "WCPPID/PR3DCluster.h"

namespace WCPPID{
  void Protect_Over_Clustering(WCPPID::PR3DClusterSelection& live_clusters,
			       std::map<WCPPID::PR3DCluster*, int>& map_cluster_parent_id, // cluster to main cluster
			       std::map<int, std::vector<WCPPID::PR3DCluster*> >& map_parentid_clusters // main cluster to clusters
			       );
}

#endif 
