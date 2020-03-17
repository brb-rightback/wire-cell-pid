#include "WCPPID/ProtectOverClustering.h"

using namespace WCP;

void WCPPID::Protect_Over_Clustering(WCPPID::PR3DClusterSelection& live_clusters, std::map<WCPPID::PR3DCluster*, int>& map_cluster_parent_id,  std::map<int, std::vector<WCPPID::PR3DCluster*> >& map_parentid_clusters  ){
  // update live_clusters,  map_parentid_clusters, map_cluster_parent_id;
  // keep the original id for the main cluster ...

  
   WCPPID::PR3DClusterSelection temp_live_clusters;
    std::map<WCPPID::PR3DCluster*, int> temp_map_cluster_parent_id; // cluster to main cluster
    std::map<int, std::vector<WCPPID::PR3DCluster*> > temp_map_parentid_clusters; // main cluster to clusters
}
