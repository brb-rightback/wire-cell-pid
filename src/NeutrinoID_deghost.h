#include "WCPData/DynamicToyPointCloud.h"

bool sortbysec(const std::pair<WCPPID::PR3DCluster*,double> &a,
	       const std::pair<WCPPID::PR3DCluster*,double> &b){
  return (a.second > b.second);
}


void WCPPID::NeutrinoID::deghost_clusters(){
  std::map<int, WCPPID::ProtoSegmentSelection> map_cluster_id_segments;  // ID segments
  std::map<WCPPID::PR3DCluster*, double> cluster_length_map;
  WCPPID::PR3DClusterSelection ordered_clusters;
  order_clusters(ordered_clusters, map_cluster_id_segments, cluster_length_map);
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  
  DynamicToyPointCloud global_point_cloud(angle_u,angle_v,angle_w);
  DynamicToyPointCloud global_steiner_point_cloud(angle_u,angle_v,angle_w);
  DynamicToyPointCloud global_skeleton_cloud(angle_u,angle_v,angle_w);

  std::vector<WCPPID::PR3DCluster*> to_be_removed_clusters;

  for (size_t i=0;i!= ordered_clusters.size(); i++){
    if (i==0){
      // fill anyway ...
      global_point_cloud.AddPoints(ordered_clusters.at(i)->get_point_cloud());
      global_steiner_point_cloud.AddPoints(ordered_clusters.at(i)->get_point_cloud_steiner());
      auto it = map_cluster_id_segments.find(ordered_clusters.at(i)->get_cluster_id());
      if (it != map_cluster_id_segments.end()){
	for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	  global_skeleton_cloud.AddPoints((*it1)->get_point_vec());
	}
      }
    }else{
      WCPPID::PR3DCluster *cluster = ordered_clusters.at(i);
      int num_dead[3]={0,0,0};
      int num_unique[3]={0,0,0};
      int num_total_points = 0;
      // dist cut ...
      double dis_cut = 1.2*units::cm;
      auto it = map_cluster_id_segments.find(ordered_clusters.at(i)->get_cluster_id());

      if (it != map_cluster_id_segments.end()){
	for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	  PointVector& pts = (*it1)->get_point_vec();
	  num_total_points += pts.size();
	  for (size_t j=0;j!=pts.size();j++){
	    Point test_point = pts.at(j);
	    bool flag_dead = ct_point_cloud->get_closest_dead_chs(test_point,0);
	    if (!flag_dead){
	      std::tuple<double, WCP::PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 0);
	      if (std::get<0>(results)<=dis_cut/2.){
	      }else{
		results = global_steiner_point_cloud.get_closest_2d_point_info(test_point, 0);
		if (std::get<0>(results)<=dis_cut/3.*2.){
		}else{
		  results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 0);
		  if (std::get<0>(results)<=dis_cut*5/4.){
		  }else{
		    num_unique[0]++;
		  }
		}
	      }
	    }else{
	      num_dead[0]++;
	    }

	  
	    // V plane ...
	    flag_dead = ct_point_cloud->get_closest_dead_chs(test_point,1);
	    if (!flag_dead){
	      std::tuple<double, WCP::PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 1);
	      if (std::get<0>(results)<=dis_cut/2.){
	      }else{
		results = global_steiner_point_cloud.get_closest_2d_point_info(test_point, 1);
		if (std::get<0>(results)<=dis_cut/3.*2.){
		}else{
		  results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 1);
		  if (std::get<0>(results)<=dis_cut*5/4.){
		  }else{
		    num_unique[1]++;
		  }
		}
	      }
	    }else{
	      num_dead[1]++;
	    }

	    // W plane ...
	    flag_dead = ct_point_cloud->get_closest_dead_chs(test_point,2);
	    if (!flag_dead){
	      std::tuple<double, WCP::PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 2);
	      if (std::get<0>(results)<=dis_cut/2.){
	      }else{
		results = global_steiner_point_cloud.get_closest_2d_point_info(test_point, 2);
		if (std::get<0>(results)<=dis_cut/3.*2.){
		}else{
		  results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 2);
		  if (std::get<0>(results)<=dis_cut*5/4.){
		  }else{
		    num_unique[2]++;
		  }
		}
	      }
	    }else{
	      num_dead[2]++;
	    }
	    
	    
	  }
	} // loop over segments
      }else{
	// no segments ...
	
      }
      
      bool flag_add = true;
      
      

      double unique_percent_u = num_unique[0]*1.0/num_total_points;
      double unique_percent_v = num_unique[1]*1.0/num_total_points;
      double unique_percent_w = num_unique[2]*1.0/num_total_points;

      double dead_percent_u = num_dead[0]*1.0/num_total_points;
      double dead_percent_v = num_dead[1]*1.0/num_total_points;
      double dead_percent_w = num_dead[2]*1.0/num_total_points;
      
      double max_unique_percent = std::max(unique_percent_u, unique_percent_v);
      max_unique_percent = std::max(max_unique_percent, unique_percent_w);
      double min_unique_percent = std::min(unique_percent_u, unique_percent_v);
      min_unique_percent = std::min(min_unique_percent, unique_percent_w);
      double ave_unique_percent = (unique_percent_u + unique_percent_v + unique_percent_w)/3.;
      double max_dead_percent = std::max(dead_percent_u, dead_percent_v);
      max_dead_percent = std::max(max_dead_percent, dead_percent_w);
      
      
      
      if (max_dead_percent > 0.8 && max_unique_percent < 0.35 && ave_unique_percent < 0.16 && min_unique_percent < 0.08 ||
	  max_unique_percent < 0.1 && ave_unique_percent < 0.05 && min_unique_percent < 0.025)
	flag_add = false;

      //std::cout << flag_add << " " << cluster->get_cluster_id() << " " << num_total_points << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << max_unique_percent << " " << min_unique_percent << " " << ave_unique_percent << " " << max_dead_percent << std::endl;
      
      if (flag_add){
	global_point_cloud.AddPoints(ordered_clusters.at(i)->get_point_cloud());
	global_steiner_point_cloud.AddPoints(ordered_clusters.at(i)->get_point_cloud_steiner());
	auto it = map_cluster_id_segments.find(ordered_clusters.at(i)->get_cluster_id());
	if (it != map_cluster_id_segments.end()){
	  for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	    global_skeleton_cloud.AddPoints((*it1)->get_point_vec());
	  }
	}
      }else{
	to_be_removed_clusters.push_back(cluster);
      }
    }
  } // loop over cluster


  for (auto it = to_be_removed_clusters.begin(); it!= to_be_removed_clusters.end(); it++){
    auto it2 = map_cluster_id_segments.find((*it)->get_cluster_id());
    if (it2 != map_cluster_id_segments.end()){
      for (auto it1 = it2->second.begin(); it1 != it2->second.end(); it1++){
    	del_proto_segment(*it1);
      }
    }
  }
  
  WCPPID::ProtoVertexSelection tmp_vertices;
  for (auto it1 = map_vertex_segments.begin(); it1!=map_vertex_segments.end(); it1++){
    if (it1->second.size()==0) tmp_vertices.push_back(it1->first);
  }
  for (auto it1 = tmp_vertices.begin(); it1!=tmp_vertices.end(); it1++){
    del_proto_vertex(*it1);
  }
  
}


void WCPPID::NeutrinoID::order_clusters(WCPPID::PR3DClusterSelection& ordered_clusters, std::map<int, WCPPID::ProtoSegmentSelection>& map_cluster_id_segments, std::map<WCPPID::PR3DCluster*, double>& map_cluster_total_length){
  
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    double length = sg->get_length();

    //std::cout << sg->get_point_vec().size() << " " << sg->get_point_vec().front() << " " << sg->get_point_vec().back() << std::endl;
    
    if (map_cluster_total_length.find(map_id_cluster[sg->get_cluster_id()])==map_cluster_total_length.end()){
      WCPPID::ProtoSegmentSelection segments;
      segments.push_back(sg);
      map_cluster_id_segments[sg->get_cluster_id()] =segments;
      map_cluster_total_length[map_id_cluster[sg->get_cluster_id()]] = length;
    }else{
      map_cluster_id_segments[sg->get_cluster_id()].push_back(sg);
      map_cluster_total_length[map_id_cluster[sg->get_cluster_id()]] += length;
    }
  }
  
 
  //  std::cout << map_cluster_total_length.size() << " " << map_cluster_id_segments.size() << std::endl;
  std::vector<std::pair<WCPPID::PR3DCluster*, double>> temp_pair_vec;
  for (auto it = map_cluster_total_length.begin(); it!=map_cluster_total_length.end(); it++){
    temp_pair_vec.push_back(std::make_pair(it->first, it->second));
  }
  sort(temp_pair_vec.begin(), temp_pair_vec.end(), sortbysec);
  for (auto it = temp_pair_vec.begin(); it!=temp_pair_vec.end();it++){
    ordered_clusters.push_back(it->first);
    //    std::cout << it->first->get_cluster_id() << " " << map_cluster_total_length[it->first] << " " << map_cluster_id_segments[it->first->get_cluster_id()].size() << std::endl;
  }
}
