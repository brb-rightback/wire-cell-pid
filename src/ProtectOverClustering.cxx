#include "WCPPID/ProtectOverClustering.h"

using namespace WCP;

void WCPPID::Protect_Over_Clustering(std::vector<int>& to_be_checked, WCPPID::PR3DClusterSelection& live_clusters, std::map<WCPPID::PR3DCluster*, int>& map_cluster_parent_id,  std::map<int, std::vector<WCPPID::PR3DCluster*> >& map_parentid_clusters, WCP::ToyCTPointCloud& ct_point_cloud ){
  // update live_clusters,  map_parentid_clusters, map_cluster_parent_id;
  // keep the original id for the main cluster ...
  
   WCPPID::PR3DClusterSelection temp_live_clusters;
   
   std::set<WCPPID::PR3DCluster*> examined_clusters;
   for (size_t i=0;i!=to_be_checked.size(); i++){
     int curr_main_cluster_id = to_be_checked.at(i);
     auto it = map_parentid_clusters.find(curr_main_cluster_id);
     if (it == map_parentid_clusters.end()) continue;
     for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
       WCPPID::PR3DCluster *temp_cluster = *it1;
       examined_clusters.insert(temp_cluster);
       map_cluster_parent_id.erase(temp_cluster);
     }
     map_parentid_clusters.erase(curr_main_cluster_id);
   }
   
   for (auto it = live_clusters.begin(); it != live_clusters.end(); it++){
     if (examined_clusters.find(*it)==examined_clusters.end())
       temp_live_clusters.push_back(*it);
   }
   
   for (auto it = examined_clusters.begin(); it!=examined_clusters.end(); it++){
     delete (*it);
   }
   
   live_clusters = temp_live_clusters;
   
   
}
