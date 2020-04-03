void WCPPID::NeutrinoID::separate_track_shower(){
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    sg->is_shower_trajectory();
    sg->is_shower_topology();
    //  std::cout <<" C: " << sg->get_cluster_id() << " " << sg->get_id() << " " << sg->get_flag_shower_trajectory() << " " << sg->get_flag_shower_topology() << std::endl;
  }
  PR3DClusterSelection all_clusters = other_clusters;
  all_clusters.push_back(main_cluster);
  for (auto it = all_clusters.begin(); it!=all_clusters.end(); it++){
    WCPPID::PR3DCluster* temp_cluster = *it;

    std::map<int, WCPPID::ProtoSegment*> map_id_seg;
    for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
      map_id_seg[sg->get_id()] = sg;
      //std::cout << temp_cluster->get_cluster_id() << " " << sg->get_id() << std::endl;
    }
    
    std::vector<int>& point_sub_cluster_ids = temp_cluster->get_point_sub_cluster_ids();
    std::vector<bool>& point_flag_showers = temp_cluster->get_point_flag_showers();
    point_flag_showers.resize(point_sub_cluster_ids.size(), false);

    
    
    for (size_t i=0;i!=point_sub_cluster_ids.size();i++){ 
      if (point_sub_cluster_ids.at(i) == -1) continue;
      if (map_id_seg.find(point_sub_cluster_ids.at(i)) == map_id_seg.end()) continue;
      //    std::cout << point_sub_cluster_ids.at(i) << " " << map_id_seg[point_sub_cluster_ids.at(i)] << std::endl; 
      if ( map_id_seg[point_sub_cluster_ids.at(i)]->get_flag_shower()) point_flag_showers.at(i) = true; 
    } 
  }
  
}



void WCPPID::NeutrinoID::determine_direction(WCPPID::PR3DCluster* temp_cluster){
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;

    if (sg->get_flag_shower_trajectory()){
      // trajectory shower
    }else if (sg->get_flag_shower_topology()){
      // topology shower
    }else{
      // track ...
    }

    
  }
}
