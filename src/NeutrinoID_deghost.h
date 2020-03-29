bool sortbysec(const std::pair<WCPPID::PR3DCluster*,double> &a,
	       const std::pair<WCPPID::PR3DCluster*,double> &b){
  return (a.second > b.second);
}


void WCPPID::NeutrinoID::deghost_clusters(){
  std::map<int, WCPPID::ProtoSegmentSelection> map_cluster_id_segments;
  std::map<WCPPID::PR3DCluster*, double> map_cluster_max_length;
  WCPPID::PR3DClusterSelection ordered_clusters;
  
  order_clusters(ordered_clusters, map_cluster_id_segments, map_cluster_max_length);
  
 
  
}


void WCPPID::NeutrinoID::order_clusters(WCPPID::PR3DClusterSelection& ordered_clusters, std::map<int, WCPPID::ProtoSegmentSelection>& map_cluster_id_segments, std::map<WCPPID::PR3DCluster*, double>& map_cluster_max_length){


  
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    double length = sg->get_length();
    
    if (map_cluster_max_length.find(map_id_cluster[sg->get_cluster_id()])==map_cluster_max_length.end()){
      WCPPID::ProtoSegmentSelection segments;
      segments.push_back(sg);
      map_cluster_id_segments[sg->get_cluster_id()] =segments;
      map_cluster_max_length[map_id_cluster[sg->get_cluster_id()]] = length;
    }else{
      map_cluster_id_segments[sg->get_cluster_id()].push_back(sg);
      if (length > map_cluster_max_length[map_id_cluster[sg->get_cluster_id()]]) map_cluster_max_length[map_id_cluster[sg->get_cluster_id()]] = length;
    }
  }
  
 
  //  std::cout << map_cluster_max_length.size() << " " << map_cluster_id_segments.size() << std::endl;
  std::vector<std::pair<WCPPID::PR3DCluster*, double>> temp_pair_vec;
  for (auto it = map_cluster_max_length.begin(); it!=map_cluster_max_length.end(); it++){
    temp_pair_vec.push_back(std::make_pair(it->first, it->second));
  }
  sort(temp_pair_vec.begin(), temp_pair_vec.end(), sortbysec);
  for (auto it = temp_pair_vec.begin(); it!=temp_pair_vec.end();it++){
    ordered_clusters.push_back(it->first);
    std::cout << map_cluster_max_length[it->first] << std::endl;
  }
}
