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
  {
    bool flag_update = true;
    std::set<WCPPID::ProtoSegment* > used_segments;
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
	  if (used_segments.find(sg) != used_segments.end()) continue;
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
	    sg->set_dir_weak(true);
	    used_segments.insert(sg);
	    flag_update = true;
	  }
	  used_vertices.insert(vtx);
      }
	
	if (flag_update) break; // is this the best?
      }
    } // keep updating
  }

  
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
  std::cout << "Information after initial logic examination: " << std::endl;
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

  //  std::cout << main_vertex_candidates.size() << std::endl;
  
  if (main_vertex_candidates.size()==1){
    main_vertex = main_vertex_candidates.front();
  }else{
    main_vertex = compare_main_vertices(main_vertex_candidates);
  }
  bool flag_check = examine_direction(main_vertex);
  if (!flag_check) std::cout << "Wrong: inconsistency for track directions! " << std::endl;

  std::cout << "Main Vertex " << main_vertex->get_fit_pt() << " connecting to: ";
  for (auto it = map_vertex_segments[main_vertex].begin(); it!=map_vertex_segments[main_vertex].end(); it++){
    std::cout << (*it)->get_id() << ", ";
  }
  std::cout << std::endl;
  
}

WCPPID::ProtoVertex* WCPPID::NeutrinoID::compare_main_vertices(WCPPID::ProtoVertexSelection& vertex_candidates){
  std::map<WCPPID::ProtoVertex*, double> map_vertex_num;
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    map_vertex_num[*it] = 0;
  }

  // find the proton in and out ...
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    
    int n_proton_in = 0;
    int n_proton_out = 0;
    
    for (auto it1 = map_vertex_segments[vtx].begin(); it1 != map_vertex_segments[vtx].end(); it1++){
      WCPPID::ProtoSegment *sg = (*it1);
      if ((sg->is_dir_weak() || sg->get_flag_dir() == 0) && fabs(sg->get_particle_type())==2212){
	WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vtx);

	for (auto it2 = map_vertex_segments[other_vertex].begin(); it2!=map_vertex_segments[other_vertex].end(); it2++){
	  bool flag_start; 
	  if ((*it2)->get_wcpt_vec().front().index == other_vertex->get_wcpt().index) 
	    flag_start = true; 
	  else if ((*it2)->get_wcpt_vec().back().index == other_vertex->get_wcpt().index) 
	    flag_start = false; 
	
	  if ((flag_start && (*it2)->get_flag_dir()==1 || (!flag_start) && (*it2)->get_flag_dir()==-1) && (!(*it2)->is_dir_weak()) && fabs((*it2)->get_particle_type())==2212 ){ 
	    n_proton_out ++; 
	  } 
	
	  if ((flag_start && (*it2)->get_flag_dir()==-1 || (!flag_start) && (*it2)->get_flag_dir()==1) && (!(*it2)->is_dir_weak()) && fabs((*it2)->get_particle_type())==2212){ 
	    n_proton_in ++; 
	  } 
	  if (((*it2)->is_dir_weak() || (*it2)->get_flag_dir() == 0) && fabs((*it2)->get_particle_type())==2212)
	    n_proton_in ++;
	}
      }
    }
    
    // positive is good ...
    map_vertex_num[vtx] -= (n_proton_in - n_proton_out);
    // std::cout << map_vertex_num[vtx] << " " << n_proton_in << " " << n_proton_out << std::endl;
  }

  // whether the vertex is at beginning or not ...
  double min_z = 1e9;
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (vtx->get_fit_pt().z < min_z) min_z = vtx->get_fit_pt().z;
  }
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    // std::cout << vtx->get_fit_pt().z << std::endl;
    map_vertex_num[vtx] -= (vtx->get_fit_pt().z - min_z)/(400*units::cm);  
    //    std::cout << map_vertex_segments[vtx].size() << std::endl;
    // number of tracks, more is good 
    map_vertex_num[vtx] += map_vertex_segments[vtx].size()/4.;

    //    std::cout << map_vertex_num[vtx] << " " << (vtx->get_fit_pt().z - min_z)/(400*units::cm) << " " << map_vertex_segments[vtx].size()/4. << std::endl;
  }
  
  // whether the vetex is at boundary or not ...
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (fid->inside_fiducial_volume(vtx->get_fit_pt(),offset_x))
      map_vertex_num[vtx] +=0.5; // good 
    // std::cout << map_vertex_num[vtx] << " " << fid->inside_fiducial_volume(vtx->get_fit_pt(),offset_x) << std::endl;
  }
  
  double max_val = -1e9; WCPPID::ProtoVertex* max_vertex = 0;
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (map_vertex_num[vtx] > max_val){
      max_val = map_vertex_num[vtx];
      max_vertex = vtx;
    }
  }

  //  std::cout << (max_vertex->get_fit_pt().z-min_z)/(400*units::cm) << std::endl;
  
  return max_vertex;
}



bool WCPPID::NeutrinoID::examine_direction(WCPPID::ProtoVertex* main_vertex){
  bool flag_check = true;

  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;

  std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > segments_to_be_examined;
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    segments_to_be_examined.push_back(std::make_pair(main_vertex, *it));
  }
  used_vertices.insert(main_vertex);

  while(segments_to_be_examined.size()>0){
    std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > temp_segments;
    for (auto it = segments_to_be_examined.begin(); it!= segments_to_be_examined.end(); it++){
      WCPPID::ProtoVertex *prev_vtx = it->first;
      WCPPID::ProtoSegment *current_sg = it->second;
      if (used_segments.find(current_sg)!=used_segments.end()) continue; // looked at it before ...

      if (current_sg->get_flag_dir() ==0 || current_sg->is_dir_weak()){ // change direction ...
	bool flag_start;
	if (current_sg->get_wcpt_vec().front().index == prev_vtx->get_wcpt().index)
	  flag_start = true;
	else if (current_sg->get_wcpt_vec().back().index == prev_vtx->get_wcpt().index)
	  flag_start = false;

	if (flag_start) current_sg->set_flag_dir(1);
	else current_sg->set_flag_dir(-1);
	current_sg->set_dir_weak(true);
      }
      used_segments.insert(current_sg);

      WCPPID::ProtoVertex* curr_vertex = find_other_vertex(current_sg, prev_vtx);
      if (used_vertices.find(curr_vertex) != used_vertices.end()) continue;
      for (auto it1 = map_vertex_segments[curr_vertex].begin(); it1!= map_vertex_segments[curr_vertex].end(); it1++){
	temp_segments.push_back(std::make_pair(curr_vertex, *it1));
      }
      used_vertices.insert(curr_vertex);
    }
    segments_to_be_examined = temp_segments;
  }
  
  // std::cout << used_vertices.size() << " " << used_segments.size() << std::endl;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  // find the muon candidate, make the other ones connecting to main vertex to be pion
  {
    WCPPID::ProtoSegment *muon_sg = 0;
    double muon_length = 0;
    WCPPID::ProtoSegmentSelection pion_sgs;
    for (auto it = map_vertex_segments[main_vertex].begin(); it!=map_vertex_segments[main_vertex].end(); it++){
      WCPPID::ProtoSegment *sg = *it;
      if (abs(sg->get_particle_type()) == 13){
	if (sg->get_length() > muon_length){
	  muon_length = sg->get_length();
	  muon_sg = sg;
	}
	pion_sgs.push_back(sg);
      }
    }
    for (auto it = pion_sgs.begin(); it!= pion_sgs.end(); it++){
      if (*it == muon_sg) continue;
      (*it)->set_particle_type(211);
      (*it)->set_particle_mass(mp.get_mass_pion());
      if ((*it)->get_particle_4mom(3)>0)
	(*it)->cal_4mom_range();
    }
  }
  
  // print ...
  std::cout << "Information after main vertex determination: " << std::endl;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != main_vertex->get_cluster_id()) continue;
    
    if (sg->get_flag_shower_topology()){
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " S_topo "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak() << std::endl;
    }else if (sg->get_flag_shower_trajectory()){
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " S_traj "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak()<< std::endl;
    }else{
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " Track  "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak() << std::endl;
    }
  }
  
  // examination ...
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != main_vertex->get_cluster_id()) continue;
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
    
    if (n_in > 1 ) flag_check = false;
    if (n_in_shower >0 && n_out_tracks > 0) flag_check = false;
  }
  
  return flag_check;
}
