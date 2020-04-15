void WCPPID::NeutrinoID::separate_track_shower(WCPPID::PR3DCluster* temp_cluster){
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    sg->is_shower_trajectory();
    sg->is_shower_topology();
  }
}

void WCPPID::NeutrinoID::separate_track_shower(){
  //  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
  //  WCPPID::ProtoSegment *sg = it->first;
  //  sg->is_shower_trajectory();
  //  sg->is_shower_topology();
    //  std::cout <<" C: " << sg->get_cluster_id() << " " << sg->get_id() << " " << sg->get_flag_shower_trajectory() << " " << sg->get_flag_shower_topology() << std::endl;
  //}
  
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

    bool flag_print = false;
    if (sg->get_cluster_id() == main_cluster->get_cluster_id()) flag_print = true;

    // std::cout << sg << " " << sg->get_id() << " " << sg->get_flag_shower_trajectory() << " " << sg->get_flag_shower_topology() << std::endl;
    
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
    // std::cout << sg->get_id() << " " << sg->get_particle_type() << std::endl;
    
  }
}

void WCPPID::NeutrinoID::improve_maps_one_in(WCPPID::PR3DCluster* temp_cluster, bool flag_strong_check){
  bool flag_update = true;
  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;
  while(flag_update){
    flag_update = false;
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
	
	if ((flag_start && sg->get_flag_dir()==-1 || (!flag_start) && sg->get_flag_dir()==1)){
	  if (flag_strong_check){
	    if (!sg->is_dir_weak()) n_in ++;
	  }else{
	    n_in ++;
	  }
	}
	if (sg->get_flag_dir()==0 || sg->is_dir_weak())
	  map_sg_dir[sg] = flag_start;
      }
      
      if (map_sg_dir.size() == 0 ) used_vertices.insert(vtx); // no outside to change ...
      
      if (n_in>0){
	for (auto it1 = map_sg_dir.begin(); it1!=map_sg_dir.end(); it1++){
	  WCPPID::ProtoSegment *sg = it1->first;
	  bool flag_start = it1->second;
	  if (flag_start){
	    sg->set_flag_dir(1);
	  }else{
	    sg->set_flag_dir(-1);
	  }
	  sg->cal_4mom();
	  sg->set_dir_weak(true);
	  used_segments.insert(sg); // change direction ...
	  flag_update = true;
	}
	used_vertices.insert(vtx); // change vertices ...
      }
      //	if (flag_update) break; // is this the best?
    }
  } // keep updating
}


bool WCPPID::NeutrinoID::examine_maps(WCPPID::ProtoVertex *temp_vertex){
  return examine_maps(temp_vertex->get_cluster_id());
}

bool WCPPID::NeutrinoID::examine_maps(int temp_cluster_id){
  bool flag_return = true;
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster_id) continue;
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

    if (n_in > 1 && n_in != n_in_shower) {
      std::cout << "Wrong: Multiple (" << n_in << ") particles into a vertex! " << std::endl;
      print_segs_info(temp_cluster_id, vtx);
      flag_return = false;
    }
    if (n_in_shower >0 && n_out_tracks > 0) {
      std::cout << "Wrong: " << n_in_shower << " showers in and " << n_out_tracks << " tracks out! " << std::endl;
      print_segs_info(temp_cluster_id, vtx);
      flag_return = false;
    }
  }
  return flag_return;
}

bool WCPPID::NeutrinoID::examine_maps(WCPPID::PR3DCluster* temp_cluster){
  return examine_maps(temp_cluster->get_cluster_id());
}

void WCPPID::NeutrinoID::improve_maps_no_dir_tracks(int temp_cluster_id){
  
  bool flag_update = true;
  while(flag_update){
    flag_update = false;
    for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      if (sg->get_cluster_id() != temp_cluster_id) continue;
      if (sg->get_flag_shower()) continue;
      if (sg->get_flag_dir()!=0) continue;
      
      std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> two_vertices = find_vertices(sg);
      int nshowers[2]={0,0};
      int n_in[2] = {0,0};
      for (auto it1 = map_vertex_segments[two_vertices.first].begin(); it1!=map_vertex_segments[two_vertices.first].end(); it1++){
  	if ((*it1)->get_flag_shower()) nshowers[0] ++;
	bool flag_start;
	if ((*it1)->get_wcpt_vec().front().index == two_vertices.first->get_wcpt().index)
	  flag_start = true;
	else if ((*it1)->get_wcpt_vec().back().index == two_vertices.first->get_wcpt().index)
	  flag_start = false;
	if ((flag_start && (*it1)->get_flag_dir()==-1 || (!flag_start) && (*it1)->get_flag_dir()==1)) n_in[0]++;
      }
      for (auto it1 = map_vertex_segments[two_vertices.second].begin(); it1!=map_vertex_segments[two_vertices.second].end(); it1++){
  	if ((*it1)->get_flag_shower()) nshowers[1] ++;
	bool flag_start;
	if ((*it1)->get_wcpt_vec().front().index == two_vertices.second->get_wcpt().index)
	  flag_start = true;
	else if ((*it1)->get_wcpt_vec().back().index == two_vertices.second->get_wcpt().index)
	  flag_start = false;
	if ((flag_start && (*it1)->get_flag_dir()==-1 || (!flag_start) && (*it1)->get_flag_dir()==1)) n_in[1]++;
      }

     
      if (  (nshowers[0]+1 == map_vertex_segments[two_vertices.first].size() || nshowers[0]>0) && // one shower or nothing
	    (nshowers[1]+1 == map_vertex_segments[two_vertices.second].size() || nshowers[1] >0) && // one shower or nothing
	    (nshowers[0] + nshowers[1] >2) && // 3 showers in total ...
	    (nshowers[0]+1 == map_vertex_segments[two_vertices.first].size() || nshowers[1]+1 == map_vertex_segments[two_vertices.second].size() ) // one side all showers
	    ){
	bool flag_start;
	if (sg->get_wcpt_vec().front().index == two_vertices.first->get_wcpt().index)
	  flag_start = true;
	else if (sg->get_wcpt_vec().back().index == two_vertices.first->get_wcpt().index)
	  flag_start = false;
	if (flag_start){
	  if (n_in[0]>0 && n_in[1]==0){
	    sg->set_flag_dir(1);
	  }else if (n_in[0]==0 && n_in[1]>0){
	    sg->set_flag_dir(-1);
	  }
	}else{
	  if (n_in[0]>0 && n_in[1]==0){
	    sg->set_flag_dir(-1);
	  }else if (n_in[0]==0 && n_in[1]>0){
	    sg->set_flag_dir(1);
	  }
	}
	sg->set_particle_type(11);
      	TPCParams& mp = Singleton<TPCParams>::Instance();
      	sg->set_particle_mass(mp.get_mass_electron());
      	if (sg->get_particle_4mom(3)>0) sg->cal_4mom();
      	flag_update = true;
      }else if ( nshowers[0]+1 == map_vertex_segments[two_vertices.first].size() && nshowers[0]>=2 && sg->get_particle_type() == 2212){
	TVector3 v1 = sg->cal_dir_3vector(two_vertices.first->get_fit_pt(), 5*units::cm);
	double min_angle = 180;
	for (auto it2 = map_vertex_segments[two_vertices.first].begin(); it2 != map_vertex_segments[two_vertices.first].end(); it2++){
	  if (*it2 == sg) continue;
	  TVector3 v2 = (*it2)->cal_dir_3vector(two_vertices.first->get_fit_pt(), 5*units::cm);
	  double angle = fabs(v1.Angle(v2)/3.1415926*180.-180.);
	  if (angle < min_angle) min_angle = angle;
	}
	double dQ_dx_rms = sg->get_rms_dQ_dx();

	if ( dQ_dx_rms > 0.75 * (43e3/units::cm) && min_angle < 30 ||
	     dQ_dx_rms > 0.4 *(43e3/units::cm) && min_angle < 15){
	  bool flag_start;
	  if (sg->get_wcpt_vec().front().index == two_vertices.first->get_wcpt().index)
	    flag_start = true;
	  else if (sg->get_wcpt_vec().back().index == two_vertices.first->get_wcpt().index)
	    flag_start = false;
	  if (flag_start)
	    sg->set_flag_dir(-1);
	  else
	    sg->set_flag_dir(1);
	  sg->set_particle_type(11);
	  TPCParams& mp = Singleton<TPCParams>::Instance();
	  sg->set_particle_mass(mp.get_mass_electron());
	  if (sg->get_particle_4mom(3)>0) sg->cal_4mom();
	  flag_update = true;
	}
      }else if ( nshowers[1]+1 == map_vertex_segments[two_vertices.second].size() && nshowers[1]>=2 && sg->get_particle_type() == 2212){
	TVector3 v1 = sg->cal_dir_3vector(two_vertices.second->get_fit_pt(), 5*units::cm);
	double min_angle = 180;
	for (auto it2 = map_vertex_segments[two_vertices.second].begin(); it2 != map_vertex_segments[two_vertices.second].end(); it2++){
	  if (*it2 == sg) continue;
	  TVector3 v2 = (*it2)->cal_dir_3vector(two_vertices.second->get_fit_pt(), 5*units::cm);
	  double angle = fabs(v1.Angle(v2)/3.1415926*180.-180.);
	  if (angle < min_angle) min_angle = angle;
	}
	double dQ_dx_rms = sg->get_rms_dQ_dx();

	if ( dQ_dx_rms > 0.75 * (43e3/units::cm) && min_angle < 30 ||
	     dQ_dx_rms > 0.4 *(43e3/units::cm) && min_angle < 15){
	  bool flag_start;
	  if (sg->get_wcpt_vec().front().index == two_vertices.second->get_wcpt().index)
	    flag_start = true;
	  else if (sg->get_wcpt_vec().back().index == two_vertices.second->get_wcpt().index)
	    flag_start = false;
	  if (flag_start)
	    sg->set_flag_dir(-1);
	  else
	    sg->set_flag_dir(1);
	  sg->set_particle_type(11);
	  TPCParams& mp = Singleton<TPCParams>::Instance();
	  sg->set_particle_mass(mp.get_mass_electron());
	  if (sg->get_particle_4mom(3)>0) sg->cal_4mom();
	  flag_update = true;
	}
      }

    } // loop over all segments
  } // keep updating
  

  /*
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster_id) continue;
    if (sg->get_flag_shower()) continue;
    if (sg->get_flag_dir()!=0) continue;
    
    std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> two_vertices = find_vertices(sg);
    int nshowers[2]={0,0};
    int n_in[2] = {0,0};
    int n_out[2] = {0,0};
    for (auto it1 = map_vertex_segments[two_vertices.first].begin(); it1!=map_vertex_segments[two_vertices.first].end(); it1++){
      if ((*it1)->get_flag_shower()) nshowers[0] ++;
      bool flag_start;
      if ((*it1)->get_wcpt_vec().front().index == two_vertices.first->get_wcpt().index)
	flag_start = true;
      else if ((*it1)->get_wcpt_vec().back().index == two_vertices.first->get_wcpt().index)
	flag_start = false;
      if ((flag_start && (*it1)->get_flag_dir()==-1 || (!flag_start) && (*it1)->get_flag_dir()==1)) n_in[0]++;
      else if ((flag_start && (*it1)->get_flag_dir()==1 || (!flag_start) && (*it1)->get_flag_dir()==-1)) n_out[0]++;
    }
    for (auto it1 = map_vertex_segments[two_vertices.second].begin(); it1!=map_vertex_segments[two_vertices.second].end(); it1++){
      if ((*it1)->get_flag_shower()) nshowers[1] ++;
      bool flag_start;
      if ((*it1)->get_wcpt_vec().front().index == two_vertices.second->get_wcpt().index)
	flag_start = true;
      else if ((*it1)->get_wcpt_vec().back().index == two_vertices.second->get_wcpt().index)
	flag_start = false;
      if ((flag_start && (*it1)->get_flag_dir()==-1 || (!flag_start) && (*it1)->get_flag_dir()==1)) n_in[1]++;
      else if ((flag_start && (*it1)->get_flag_dir()==1 || (!flag_start) && (*it1)->get_flag_dir()==-1)) n_out[1]++;
    }
    std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " Track  "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak() << " " << nshowers[0] << " " << nshowers[1] << " " << map_vertex_segments[two_vertices.first].size() << " " << map_vertex_segments[two_vertices.second].size() << " " << n_in[0] << " " << n_in[1] << " " << n_out[0] << " " << n_out[1] << " " << sg->get_rms_dQ_dx()/(43e3/units::cm) << std::endl;
  }
  */  
}

void WCPPID::NeutrinoID::fix_maps_multiple_tracks_in(int temp_cluster_id){
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster_id) continue;
    if (it->second.size()==1) continue;

    int n_in = 0;
    int n_in_shower = 0;
    std::vector<WCPPID::ProtoSegment* > in_tracks;
    
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
	else in_tracks.push_back(sg);
      }
    }

    if (n_in > 1 && n_in != n_in_shower) {
      for (auto it1 = in_tracks.begin(); it1 != in_tracks.end(); it1++){
	(*it1)->set_flag_dir(0);
	(*it1)->set_dir_weak(true);
      }
    }
  }
}

void WCPPID::NeutrinoID::fix_maps_shower_in_track_out(int temp_cluster_id){

  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster_id) continue;
    if (it->second.size()==1) continue;
    
    std::vector<WCPPID::ProtoSegment*> in_showers;
    bool flag_turn_shower_dir = false;
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      WCPPID::ProtoSegment *sg = (*it1);
      bool flag_start;
      if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	flag_start = true;
      else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	flag_start = false;
      if (flag_start && sg->get_flag_dir()==-1 || (!flag_start) && sg->get_flag_dir()==1){
	if (sg->get_flag_shower() ) in_showers.push_back(sg);
      }else if (flag_start && sg->get_flag_dir()==1 || (!flag_start) && sg->get_flag_dir()==-1 ){
	if ((!sg->get_flag_shower()) ) {
	  if (!sg->is_dir_weak()) flag_turn_shower_dir = true;
	}
      }
    } // loop segment

    if (flag_turn_shower_dir){
      for (auto it1 = in_showers.begin(); it1!= in_showers.end(); it1++){
	(*it1)->set_flag_dir((*it1)->get_flag_dir() * (-1));
	(*it1)->set_dir_weak(true);
      }
    }
    
  } // loop over all vertices
}

void WCPPID::NeutrinoID::improve_maps_multiple_tracks_in(int temp_cluster_id){
  bool flag_update = true;
  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;
  while(flag_update){
    flag_update = false;
  
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
      WCPPID::ProtoVertex *vtx = it->first;
      if (vtx->get_cluster_id() != temp_cluster_id) continue;
      if (it->second.size()==1) continue;
      if (used_vertices.find(vtx)!=used_vertices.end()) continue;
      
      int n_in = 0;
      int n_in_shower = 0;
      std::vector<ProtoSegment* > in_tracks;
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::ProtoSegment *sg = (*it1);
	bool flag_start;
	if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	  flag_start = true;
	else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	  flag_start = false;
	if (flag_start && sg->get_flag_dir()==-1 || (!flag_start) && sg->get_flag_dir()==1){
	  n_in ++;
	  if (sg->get_flag_shower() ) n_in_shower ++;
	  else in_tracks.push_back(sg);
	}
      }
      
      
      if (n_in > 1 && n_in !=n_in_shower){
	for (auto it1 = in_tracks.begin(); it1!=in_tracks.end(); it1++){
	  WCPPID::ProtoSegment *sg1 = *it1;
	  sg1->set_particle_type(11);
	  TPCParams& mp = Singleton<TPCParams>::Instance();
	  sg1->set_particle_mass(mp.get_mass_electron());
	  if (sg1->get_particle_4mom(3)>0) sg1->cal_4mom();
	  flag_update = true;
	}
      } 
      used_vertices.insert(vtx);
      
    } // loop over everything
  } // while
}

void WCPPID::NeutrinoID::improve_maps_shower_in_track_out(int temp_cluster_id){
  bool flag_update = true;
  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;
  while(flag_update){
    flag_update = false;
  
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
      WCPPID::ProtoVertex *vtx = it->first;
      if (vtx->get_cluster_id() != temp_cluster_id) continue;
      if (it->second.size()==1) continue;
      if (used_vertices.find(vtx)!=used_vertices.end()) continue;
      
      int n_in = 0;
      int n_in_shower = 0;
      std::vector<WCPPID::ProtoSegment*> out_tracks;
      std::map<WCPPID::ProtoSegment*, bool> map_no_dir_segments;
      
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::ProtoSegment *sg = (*it1);
	bool flag_start;
	if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	  flag_start = true;
	else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	  flag_start = false;
	if (flag_start && sg->get_flag_dir()==-1 || (!flag_start) && sg->get_flag_dir()==1){
	  n_in ++;
	  if (sg->get_flag_shower() ) {
	    n_in_shower ++;
	  }
	}else if (flag_start && sg->get_flag_dir()==1 || (!flag_start) && sg->get_flag_dir()==-1 ){
	  if ((!sg->get_flag_shower()) ) {
	    if (sg->is_dir_weak()) out_tracks.push_back(sg);
	  }
	}else if (sg->get_flag_dir()==0){
	  map_no_dir_segments[sg] = flag_start;
	}
      }
      

      
      if (n_in_shower >0 && (out_tracks.size() >0 || map_no_dir_segments.size() >0 )) {
	for (auto it1 = out_tracks.begin(); it1!=out_tracks.end(); it1++){
	  WCPPID::ProtoSegment *sg1 = *it1;
	  sg1->set_particle_type(11);
	  TPCParams& mp = Singleton<TPCParams>::Instance();
	  sg1->set_particle_mass(mp.get_mass_electron());
	  if (sg1->get_particle_4mom(3)>0) sg1->cal_4mom();
	  flag_update = true;
	}
	for (auto it1 = map_no_dir_segments.begin(); it1 != map_no_dir_segments.end(); it1++){
	  WCPPID::ProtoSegment *sg1 = it1->first;
	  if (used_segments.find(sg1)!=used_segments.end()) continue;
	  bool flag_start = it1->second;
	  if (flag_start){
	    sg1->set_flag_dir(1);
	  }else{
	    sg1->set_flag_dir(-1);
	  }
	  if (!sg1->get_flag_shower()){
	    sg1->set_particle_type(11);
	    TPCParams& mp = Singleton<TPCParams>::Instance();
	    sg1->set_particle_mass(mp.get_mass_electron());
	  }
	  if (sg1->get_particle_4mom(3)>0) sg1->cal_4mom();
	  sg1->set_dir_weak(true);
	  used_segments.insert(sg1);
	  flag_update = true;
	}
      } // correct other directions
      used_vertices.insert(vtx);
      
    } // loop over everything
  } // while
}

void WCPPID::NeutrinoID::print_segs_info(int temp_cluster_id, WCPPID::ProtoVertex* spec_vertex){
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster_id) continue;
    if (spec_vertex != 0){
      if (map_vertex_segments[spec_vertex].find(sg)==map_vertex_segments[spec_vertex].end()) continue;
    }
    
    if (sg->get_flag_shower_topology()){
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " S_topo "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak() << std::endl;
    }else if (sg->get_flag_shower_trajectory()){
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " S_traj "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak()<< std::endl;
    }else{
      std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " Track  "  << sg->get_flag_dir() << " " << sg->get_particle_type() << " " << sg->get_particle_mass()/units::MeV << " " << (sg->get_particle_4mom(3)-sg->get_particle_mass())/units::MeV << " " << sg->is_dir_weak() << std::endl;
    }
  }
}

void WCPPID::NeutrinoID::print_segs_info(WCPPID::ProtoVertex* temp_vertex){
  print_segs_info(temp_vertex->get_cluster_id());
}
void WCPPID::NeutrinoID::print_segs_info(WCPPID::PR3DCluster* temp_cluster){
  print_segs_info(temp_cluster->get_cluster_id());
}


void WCPPID::NeutrinoID::determine_main_vertex(WCPPID::PR3DCluster* temp_cluster){
  // update directions ... 
  // improve_maps_one_in(temp_cluster);  
  // examination ...
  // examine_maps(temp_cluster);
  // print ...
  std::cout << "Information after initial logic examination: " << std::endl;
  print_segs_info(temp_cluster);

  
  // find the main vertex ...
  bool flag_save_only_showers = true;
  std::map<ProtoVertex*, std::pair<int, int> > map_vertex_track_shower;
  WCPPID::ProtoVertexSelection main_vertex_candidates;
  for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    bool flag_in = false;
    int ntracks = 0, nshowers = 0;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      WCPPID::ProtoSegment *sg = *it1;
      if (sg->get_flag_shower()) nshowers ++;
      else ntracks ++;
      
      bool flag_start;
      if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
	flag_start = true;
      else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
	flag_start = false;
      else std::cout << "Something messed up!" << std::endl;
      if (flag_start && sg->get_flag_dir()==-1  && (!sg->is_dir_weak()) ){
	flag_in = true;
	break;
      }else if ((!flag_start)&& sg->get_flag_dir()==1 &&(!sg->is_dir_weak()) ){
	flag_in = true;
	break;
      }
    }
    if (!flag_in){
      if (ntracks > 0) 
	flag_save_only_showers = false;
      //	std::cout << vtx << " " << ntracks << " " << nshowers << std::endl;
      map_vertex_track_shower[vtx] = std::make_pair(ntracks, nshowers);
      
    }
  }

  for (auto it = map_vertex_track_shower.begin(); it!= map_vertex_track_shower.end(); it++){
    if (flag_save_only_showers){
      main_vertex_candidates.push_back(it->first);
    }else{
      if (it->second.first >0)
	main_vertex_candidates.push_back(it->first);
    }
  }


  // std::cout << main_vertex_candidates.size() << std::endl;
  
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
    map_vertex_num[vtx] -= (n_proton_in - n_proton_out);   // proton information ...
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
    map_vertex_num[vtx] -= (vtx->get_fit_pt().z - min_z)/(400*units::cm);   // position information
    //    std::cout << map_vertex_segments[vtx].size() << std::endl;
    // number of tracks, more is good 
    map_vertex_num[vtx] += map_vertex_segments[vtx].size()/4.; // number of tracks

    //    std::cout << map_vertex_num[vtx] << " " << (vtx->get_fit_pt().z - min_z)/(400*units::cm) << " " << map_vertex_segments[vtx].size()/4. << std::endl;
  }
  
  // whether the vetex is at boundary or not ...
  for (auto it = vertex_candidates.begin(); it!=vertex_candidates.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (fid->inside_fiducial_volume(vtx->get_fit_pt(),offset_x))
      map_vertex_num[vtx] +=0.5; // good      // fiducial volume ..
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

	//	if (current_sg->get_id()==19) std::cout << current_sg->get_id() << " " << flag_start << " " << map_vertex_segments[prev_vtx].size() << " " << prev_vtx->get_wcpt().index << " " << current_sg->get_wcpt_vec().front().index << " " << current_sg->get_wcpt_vec().back().index << std::endl;
	
	if (flag_start) current_sg->set_flag_dir(1);
	else current_sg->set_flag_dir(-1);
	current_sg->cal_4mom();
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
	(*it)->cal_4mom();
    }
  }

  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != main_vertex->get_cluster_id()) continue;
    if (sg->get_particle_4mom(3)==0 && sg->get_particle_mass() > 0&& (!sg->get_flag_shower_topology())){
      if (!sg->is_dir_weak() ){  // weak direction and not shower
	//      std::cout << sg->get_id() << std::endl;
	sg->cal_4mom();
      } else{  
      	std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> vertices = find_vertices(sg);
      	WCPPID::ProtoVertex *start_v = 0, *end_v = 0;
      	if (vertices.first->get_wcpt().index == sg->get_wcpt_vec().front().index){
      	  start_v = vertices.first;
      	  end_v = vertices.second;
      	}else{
      	  start_v = vertices.second;
      	  end_v = vertices.first;
      	}

      	if (sg->get_flag_dir()==1  && map_vertex_segments[end_v].size()==1  && fid->inside_fiducial_volume(end_v->get_fit_pt(),offset_x) ||
	    sg->get_flag_dir()==-1 && map_vertex_segments[start_v].size()==1 && fid->inside_fiducial_volume(start_v->get_fit_pt(),offset_x))
      	  sg->cal_4mom();
	else if (map_vertex_segments[end_v].size()==2){
	  bool flag_Michel = false;
	  for (auto it1 = map_vertex_segments[end_v].begin(); it1 != map_vertex_segments[end_v].end(); it1++){
	    if ((*it1) == sg) continue;
	    if ( (*it1)->get_flag_shower_trajectory() || (*it1)->get_flag_shower_dQdx() ) flag_Michel = true;
	  }
	  if (flag_Michel) sg->cal_4mom();
	}else if (map_vertex_segments[start_v].size()==2){
	  bool flag_Michel = false;
	  for (auto it1 = map_vertex_segments[start_v].begin(); it1 != map_vertex_segments[start_v].end(); it1++){
	    if ((*it1) == sg) continue;
	    if ( (*it1)->get_flag_shower_trajectory() || (*it1)->get_flag_shower_dQdx() ) flag_Michel = true;
	  }
	  if (flag_Michel) sg->cal_4mom();
	}
      }
    }
  }

  
  
  // print ...
  //  std::cout << "Information after main vertex determination: " << std::endl;
  //print_segs_info(main_vertex);
  // examination ...
  return examine_maps(main_vertex);
  
  
}
