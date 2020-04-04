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

    WCPPID::ProtoVertex *start_v=0, *end_v=0;
    for (auto it = map_segment_vertices[sg].begin(); it!=map_segment_vertices[sg].end(); it++){
      if ((*it)->get_wcpt().index == sg->get_wcpt_vec().front().index) start_v = *it;
      if ((*it)->get_wcpt().index == sg->get_wcpt_vec().back().index) end_v = *it;
    }
    if (start_v==0 || end_v==0){
      std::cout << "Error in finding vertices for a segment" << std::endl; 
    }

    bool flag_print = true;
    
    if (sg->get_flag_shower_trajectory()){
      // trajectory shower
      sg->determine_dir_shower_trajectory(map_vertex_segments[start_v].size(), map_vertex_segments[end_v].size(), flag_print);
    }else if (sg->get_flag_shower_topology()){
      // topology shower
      sg->determine_dir_shower_topology(map_vertex_segments[start_v].size(), map_vertex_segments[end_v].size(), flag_print);
    }else{
      // track ...
      sg->determine_dir_track(map_vertex_segments[start_v].size(), map_vertex_segments[end_v].size(), flag_print);
    }
    
  }
}


void WCPPID::NeutrinoID::determine_main_vertex(WCPPID::PR3DCluster* temp_cluster){
  // update directions ... 
  bool flag_update = true;

  while(flag_update){
    flag_update = false;
    std::set<WCPPID::ProtoVertex* > used_vertices;
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
      WCPPID::ProtoVertex *vtx = it->first;
      if (vtx->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
      if (it->second.size()==1) continue;
      if (used_vertices.find(vtx)!=used_vertices.end()) continue;

      int n_in = 0;
      std::map<WCPPID::ProtoSegment*, bool> map_sg_dir;
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::ProtoSegment *sg = (*it1);
	bool flag_start;
	if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	  flag_start = true;
	else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	  flag_start = false;

	if ((flag_start && sg->get_flag_dir()==-1 || (!flag_start) && sg->get_flag_dir()==1) && (!sg->is_dir_weak()))
	  n_in ++;
	if (sg->get_flag_dir()==0 || sg->is_dir_weak())
	  map_sg_dir[sg] = flag_start;
      }
      
      if (map_sg_dir.size() == 0 ) used_vertices.insert(vtx);
	
      if (n_in>0){
	for (auto it1 = map_sg_dir.begin(); it1!=map_sg_dir.end(); it1++){
	  WCPPID::ProtoSegment *sg = it1->first;
	  bool flag_start = it1->second;
	  if (flag_start){
	    sg->set_flag_dir(1);
	  }else{
	    sg->set_flag_dir(-1);
	  }
	  flag_update = true;
	}
	used_vertices.insert(vtx);
      }
      
      if (flag_update) break; // is this the best?
    }
  } // keep updating

  // examination ...
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    if (it->second.size()==1) continue;

    int n_in = 0;
    int n_in_shower = 0;
    int n_out_tracks = 0;
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      WCPPID::ProtoSegment *sg = (*it1);
      bool flag_start;
      if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	flag_start = true;
      else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	flag_start = false;
      
      if (flag_start && sg->get_flag_dir()==-1 || (!flag_start) && sg->get_flag_dir()==1){
	n_in ++;
	if (sg->get_flag_shower()) n_in_shower ++;
      }
      if (flag_start && sg->get_flag_dir()==1 || (!flag_start) && sg->get_flag_dir()==-1){
	if (!sg->get_flag_shower()) n_out_tracks ++;
      }
    }

    if (n_in > 1 ) std::cout << "Wrong: Multiple particles into a vertex! " << std::endl;
    if (n_in_shower >0 && n_out_tracks > 0) std::cout << "Wrong: one shower in and one track out! " << std::endl;
  }
  

  
  // print ...
  std::cout << "Information: " << std::endl;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    
    if (sg->get_flag_shower_topology()){
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " S_topo "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << std::endl;
    }else if (sg->get_flag_shower_trajectory()){
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " S_traj "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << std::endl;
    }else{
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " Track  "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << std::endl;
    }
  }

  
  // find the main vertex ...
  WCPPID::ProtoVertexSelection main_vertex_candidates;
  for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    bool flag_in = false;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      WCPPID::ProtoSegment *sg = *it1;
      bool flag_start;
      if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	flag_start = true;
      else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	flag_start = false;
      else std::cout << "Something messed up!" << std::endl;

      if (flag_start && sg->get_flag_dir()==-1 && (!sg->is_dir_weak())){
	flag_in = true;
	break;
      }else if ((!flag_start)&& sg->get_flag_dir()==1 && (!sg->is_dir_weak())){
	flag_in = true;
	break;
      }
    }
    if (!flag_in)
      main_vertex_candidates.push_back(vtx);
  }

  std::cout << main_vertex_candidates.size() << std::endl;
  
  
  
}
