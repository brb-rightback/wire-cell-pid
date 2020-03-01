#include "WCPPID/NeutrinoID.h"

using namespace WCP;

void WCPPID::NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, double flash_time){


  // deal with main cluster
  { 
    WCPPID::PR3DCluster *temp_cluster = main_cluster;
    if (temp_cluster->get_point_cloud_steiner()!=0){
      if (temp_cluster->get_point_cloud_steiner()->get_num_points() >= 2){
    	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2); 
    	temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
    	temp_cluster->cal_shortest_path(wcps.second,2);
      }
      if (temp_cluster->get_path_wcps().size()>=2){
    	temp_cluster->collect_charge_trajectory(ct_point_cloud);
    	temp_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
      }
    }
  }
  
  // deal with other clusters ...
  for (auto it1 =other_clusters.begin(); it1!=other_clusters.end();it1++){
    WCPPID::PR3DCluster *temp_cluster = (*it1);
    if (temp_cluster->get_point_cloud_steiner()!=0){
      if (temp_cluster->get_point_cloud_steiner()->get_num_points() >= 2){
    	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2); 
    	temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
    	temp_cluster->cal_shortest_path(wcps.second,2);
      }
      if (temp_cluster->get_path_wcps().size()>=2){
    	temp_cluster->collect_charge_trajectory(ct_point_cloud);
    	temp_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
      }
    }
  }

  
}
