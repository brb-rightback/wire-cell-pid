void WCPPID::NeutrinoID::examine_showers(){

  std::map<WCPPID::ProtoSegment *, WCPPID::WCShower*> map_merge_seg_shower;

  TVector3 drift_dir(1,0,0);
  std::set<WCPPID::WCShower*> del_showers;
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;

    if (map_segment_in_shower.find(sg) != map_segment_in_shower.end() && sg->get_particle_type() !=13) continue;
    
    // 7004_340_17018
    if (sg->get_length() > 45*units::cm && (!sg->is_dir_weak()) || sg->get_length() > 55*units::cm) continue;
    //    std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " " << sg->get_direct_length()/units::cm << std::endl;
    WCPPID::ProtoVertex *vtx = find_other_vertex(sg, main_vertex);

    double daughter_length = calculate_num_daughter_tracks(main_vertex, sg, false).second;
    
    bool flag_checked = false;

    /* { */
    /*   auto tmp_results = calculate_num_daughter_tracks(main_vertex, sg, false,  3*units::cm); */
    /*   std::cout << tmp_results.first << " " << tmp_results.second/units::cm << std::endl; */
    /* } */
    
    /* for (auto it1 = map_vertex_segments[vtx].begin(); it1 != map_vertex_segments[vtx].end(); it1++){ */
    /*   WCPPID::ProtoSegment *sg1 = *it1; */
    /*   if (sg1 == sg) continue; */
    /*   std::cout << sg1->get_length()/units::cm << " " << sg1->get_particle_type() << " " << sg1->is_dir_weak() << std::endl; */
    /* } */
    
    // case I ...
    if (map_vertex_to_shower.find(vtx) != map_vertex_to_shower.end()){
      TVector3 dir1 = sg->cal_dir_3vector(vtx->get_fit_pt(), 15*units::cm);
      TVector3 dir1_1 = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      bool flag_tmp_connected = false;
      
      for (auto it1 = map_vertex_to_shower[vtx].begin(); it1 != map_vertex_to_shower[vtx].end(); it1++){
	WCPPID::WCShower *shower = *it1;
	if (shower->get_start_segment()->get_particle_type()!=11) continue;
	auto pair_result = shower->get_start_vertex();
	double Eshower = 0;
	if (shower->get_kine_best() != 0){ 
	  Eshower = shower->get_kine_best();
	}else{
	  Eshower = shower->get_kine_charge();
	}
	//	std::cout << Eshower << std::endl;
	if (pair_result.second == 1 && Eshower > 60*units::MeV) flag_tmp_connected = true;
      }

      for (auto it1 = map_vertex_to_shower[vtx].begin(); it1 != map_vertex_to_shower[vtx].end(); it1++){
	WCPPID::WCShower *shower = *it1;
	if (shower->get_start_segment()->get_particle_type()!=11) continue;
	auto pair_result = shower->get_start_vertex();
	
	TVector3 dir2 = shower->cal_dir_3vector(shower->get_start_point(), 100*units::cm);
	
	double Eshower = 0;
	if (shower->get_kine_best() != 0){ 
	  Eshower = shower->get_kine_best();
	}else{
	  Eshower = shower->get_kine_charge();
	}
	
	// std::cout << Eshower/units::MeV << " " << dir1.Angle(dir2)/3.1415926*180. << " " << pair_result.second << " " << dir1.Angle(drift_dir)/3.1415926*180. << " " << dir2.Angle(drift_dir)/3.1415926*180. << " " << sg->get_length()/units::cm << std::endl;

	// check shower ..
	Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
	double tmp_total_length = 0;
	double tmp_track_length = 0;
	for (auto it2 = map_seg_vtxs.begin(); it2 != map_seg_vtxs.end(); it2++){
	  WCPPID::ProtoSegment *sg1 = it2->first;
	  if (sg1->get_cluster_id() != shower->get_start_segment()->get_cluster_id()) continue;
	  //	  std::cout << shower->get_start_segment()->get_cluster_id() << " " << sg1->get_particle_type() << " "<< sg1->get_length()/units::cm << " " << sg1->is_dir_weak() << std::endl;
	  tmp_total_length += sg1->get_length();
	  if (!sg1->is_dir_weak()) tmp_track_length += sg1->get_length();
	}
	if (tmp_track_length > 3*units::cm && tmp_track_length > 0.25 * tmp_total_length) continue;
	
	
	if (pair_result.second == 1 && Eshower > 100*units::MeV ) flag_checked = true;

	

	double tmp_angle = std::min(180 - dir1.Angle(dir2)/3.1415926*180., dir1_1.Angle(dir2)/3.1415926*180. );

	//	std::cout << Eshower << " " << pair_result.second << " " << sg->get_particle_type() << " " << sg->is_dir_weak() << " " << " " << sg->get_length()/units::cm << " " << sg->get_closest_point(shower->get_start_point()).first/units::cm << " " << 180 - dir1.Angle(dir2)/3.1415926*180. << " " << dir1_1.Angle(dir2)/3.1415926*180. << " " << daughter_length/units::cm << " " << tmp_angle << std::endl;
	
	if (pair_result.second==2){
	  // 5498_13_700  +  5774_74_3718  
	  if ((!sg->is_dir_weak()) && sg->get_length() > 3*units::cm || flag_tmp_connected) {
	    if (sg->get_closest_point(shower->get_start_point()).first < 8*units::cm && Eshower > 75*units::MeV && tmp_angle < 6){
	    }else
	      continue;
	  }
	  // 7049_1213_60664 
	  if ( sg->get_closest_point(shower->get_start_point()).first > 20*units::cm && Eshower < 150*units::MeV && tmp_angle > 2.5) continue;
	}

	//	std::cout << "haah " << std::endl;
	
	if (Eshower > 800*units::MeV && tmp_angle < 30
	    || Eshower > 150*units::MeV && tmp_angle < 10
	    || Eshower > 150*units::MeV && tmp_angle < 18 && pair_result.second == 1 && sg->is_dir_weak() // 7004_990_49528
	    || Eshower > 100*units::MeV && tmp_angle < 10 && sg->get_length() < 25*units::cm
	    || Eshower > 250*units::MeV && tmp_angle < 15
	    || Eshower > 360*units::MeV && tmp_angle < 25 
	    || Eshower > 100*units::MeV && Eshower <= 150*units::MeV && tmp_angle < 15 && sg->get_length() < 25*units::cm && flag_checked
	    || Eshower > 60*units::MeV &&  pair_result.second == 2 && sg->is_dir_weak() && (tmp_angle < 15 && sg->get_closest_point(shower->get_start_point()).first < 18 * units::cm || tmp_angle < 17.5 && sg->get_closest_point(shower->get_start_point()).first < 6 * units::cm) && sg->get_length() < 15*units::cm  // 7003_1226_61341
	    || Eshower > 60*units::MeV && pair_result.second == 2 && tmp_angle < 7.5 && sg->get_closest_point(shower->get_start_point()).first < 8 * units::cm && sg->get_length() < 20*units::cm // 7004_989_49482
	    ){
	  //std::cout << shower->get_kine_charge()/units::MeV << " " << dir1.Angle(dir2)/3.1415926*180. << " " << pair_result.second << " " << dir1.Angle(drift_dir)/3.1415926*180. << " " << dir2.Angle(drift_dir)/3.1415926*180. << std::endl;
	  map_merge_seg_shower[sg] = shower;
	  continue;
	}
      }
    }

    //    std::cout << daughter_length/units::cm << std::endl;
    
    if (map_merge_seg_shower.find(sg) != map_merge_seg_shower.end()) continue;    
    if (flag_checked) continue;
    if ( (!sg->is_dir_weak()) && sg->get_length()>6*units::cm  && daughter_length < 40*units::cm) continue;

    
    
    // case II?
    if (map_vertex_to_shower.find(main_vertex) != map_vertex_to_shower.end()){
      TVector3 dir1 = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      for (auto it1 = map_vertex_to_shower[main_vertex].begin(); it1 != map_vertex_to_shower[main_vertex].end(); it1++){
	WCPPID::WCShower *shower = *it1;
	if (shower->get_start_segment()->get_particle_type()!=11) continue;
	auto pair_result = shower->get_start_vertex();
	if (pair_result.second == 1 ) continue; //direct connected ...
	TVector3 dir2(shower->get_start_point().x - main_vertex->get_fit_pt().x,
		      shower->get_start_point().y - main_vertex->get_fit_pt().y,
		      shower->get_start_point().z - main_vertex->get_fit_pt().z);
	TVector3 dir3 = shower->cal_dir_3vector(shower->get_start_point(), 100*units::cm);
	double min_dis = sg->get_closest_point(shower->get_start_point()).first;
	  
	//std::cout << "kaka: " << sg->get_id() << " " << sg->get_particle_type() << " " << shower->get_kine_charge() << " " << dir1.Angle(dir2)/3.1415926*180. << " "<< dir1.Angle(dir3)/3.1415926*180. << " " << fabs(90 - dir1.Angle(drift_dir)/3.1415926*180.) << " " << fabs(90 - dir3.Angle(drift_dir)/3.1415926*180.) << " " << sg->get_length()/units::cm << " " << min_dis/units::cm << std::endl;
	
	if ((shower->get_kine_charge() > 80*units::MeV && dir1.Angle(dir2)/3.1415926*180. < 10
	     || shower->get_kine_charge() > 50*units::MeV && dir1.Angle(dir2)/3.1415926*180. < 3
	     || shower->get_kine_charge() > 80*units::MeV && dir1.Angle(dir3)/3.1415926*180. < 6 && dir1.Angle(dir2)/3.1415926*180. < 17.5
	     || shower->get_kine_charge() > 80*units::MeV && dir1.Angle(dir3)/3.1415926*180. < 6 && fabs(90 - dir1.Angle(drift_dir)/3.1415926*180.) < 10 && fabs(90 - dir3.Angle(drift_dir)/3.1415926*180.) < 10 && dir1.Angle(dir2)/3.1415926*180. < 30)
	    && (sg->get_length() > 5*units::cm || sg->get_length() > 3*units::cm && min_dis < 2.0*units::cm)){
	  map_merge_seg_shower[sg] = shower;


		
	  continue;
	}
      }
    } // find shower ...


    // case III
    TVector3 dir1 = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
    for (auto it1 = showers.begin(); it1 != showers.end(); it1++){
      WCPPID::WCShower *shower = *it1;
      if (shower->get_start_vertex().first == vtx || shower->get_start_vertex().first == main_vertex) continue;
      if (shower->get_start_segment()->get_particle_type()!=11) continue;
      if (shower->get_start_vertex().second <=2) {
	TVector3 dir2(shower->get_start_point().x - main_vertex->get_fit_pt().x,
		      shower->get_start_point().y - main_vertex->get_fit_pt().y,
		      shower->get_start_point().z - main_vertex->get_fit_pt().z);
	TVector3 dir3 = shower->cal_dir_3vector(shower->get_start_point(), 100*units::cm);
	//std::cout << dir1.Angle(dir2)/3.1415926*180. << " " << dir2.Angle(dir3)/3.1415926*180. << " " << sg->get_length() << std::endl;
	//std::cout << sg->get_id() << " " << sg->get_particle_type() << " " << shower->get_kine_charge() << " " << dir1.Angle(dir2)/3.1415926*180. << " " << dir2.Angle(dir3)/3.1415926*180. << " "<< sg->get_length()/units::cm << std::endl;
	
	if ( (( shower->get_kine_charge() > 80*units::MeV &&  dir1.Angle(dir2)/3.1415926*180. < 15 ) || (shower->get_kine_charge() > 50*units::MeV && dir1.Angle(dir2)/3.1415926*180. < 5))&&  dir2.Angle(dir3)/3.1415926*180. < 15 && (dir1.Angle(dir2)/3.1415926*180. + dir2.Angle(dir3)/3.1415926*180.) < 25 && sg->get_length() > 5*units::cm ){
	  map_merge_seg_shower[sg] = shower;
	  continue;
	}
      //      std::cout << shower->get_kine_charge() << std::endl;
      }
    }
    
  } // loop over segment

  
  // for the long muons ...
  for (auto it1 = showers.begin(); it1 != showers.end(); it1++){
    WCPPID::ProtoSegment *sg1 = (*it1)->get_start_segment();
    if (map_merge_seg_shower.find(sg1)!=map_merge_seg_shower.end()){
      sg1->set_flag_avoid_muon_check(true);
      del_showers.insert(*it1);
    }
  }
  if (del_showers.size() !=0){
    for (auto it1 = del_showers.begin(); it1!=del_showers.end(); it1 ++){
      WCPPID::WCShower *shower1 = *it1;
      showers.erase(find(showers.begin(), showers.end(), shower1));
      delete shower1;
    }
    del_showers.clear();
    update_shower_maps();
  }
 
  
  // merge ...
  //  std::cout << "Xin: " << map_merge_seg_shower.size() << std::endl;
  std::set<WCPPID::WCShower* > updated_showers;
  for (auto it = map_merge_seg_shower.begin(); it != map_merge_seg_shower.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    WCPPID::WCShower *shower = it->second;

    std::cout << "EM shower modification: " << shower->get_start_segment()->get_id() << " -> " << sg->get_id() << std::endl;
    updated_showers.insert(shower);
    
    auto pair_result = shower->get_start_vertex();
    
    /* for (auto it1 = map_vertex_to_shower[pair_result.first].begin(); it1!= map_vertex_to_shower[pair_result.first].end(); it1++){ */
    /*   WCPPID::WCShower *shower1 = *it1; */
    /*   if (shower == shower1) continue; */
    /*   if (shower1->get_start_vertex().second==1){ */
    /* 	del_showers.insert(shower1); */
    /*   } */
    /* } */
    
    if (pair_result.second != 1){
      shower->add_segment(sg, map_segment_vertices);
      shower->set_start_vertex(main_vertex, 1);
      shower->set_start_segment(sg);
      shower->set_start_point(main_vertex->get_fit_pt());
      std::set<WCPPID::ProtoSegment* > tmp_used_segments;
      shower->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, tmp_used_segments);
      if (sg->get_length()>45*units::cm || sg->is_dir_weak()) sg->set_flag_avoid_muon_check(true);
    }else{
      shower->add_segment(sg, map_segment_vertices);
      shower->set_start_vertex(main_vertex, 1);
      shower->set_start_segment(sg);
      shower->set_start_point(main_vertex->get_fit_pt());
      std::set<WCPPID::ProtoSegment* > tmp_used_segments;
      shower->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, tmp_used_segments);
      //std::cout << kine_charge/units::MeV << std::endl;
    }
    // add the rest of showers;
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
    for (auto it1 = showers.begin(); it1!=showers.end(); it1++){
      WCPPID::WCShower *shower1 = *it1;
      if (shower == shower1) continue;
      if (shower1->get_start_vertex().second==1 && shower1->get_start_vertex().first != main_vertex){
	if (map_vtx_segs.find(shower1->get_start_vertex().first) != map_vtx_segs.end()){
	  shower->add_shower(shower1);
	  del_showers.insert(shower1);
	}
      }
    }
    
    
    sg->set_particle_type(11);
    shower->update_particle_type();
    shower->calculate_kinematics();
    double kine_charge = cal_kine_charge(shower);
    shower->set_kine_charge(kine_charge);
    shower->set_flag_kinematics(true);
  }

 

  
  
  //std::cout << del_showers.size() << std::endl;
  for (auto it1 = del_showers.begin(); it1!=del_showers.end(); it1 ++){
    WCPPID::WCShower *shower1 = *it1;
    if (shower1->get_start_vertex().first != main_vertex){
      showers.erase(find(showers.begin(), showers.end(), shower1));
      delete shower1;
    }
  }
  del_showers.clear();

  // check other showers ...
  for (auto it = updated_showers.begin(); it != updated_showers.end(); it++){
    WCPPID::WCShower *shower = *it;
    
    TVector3 dir1 = shower->cal_dir_3vector(shower->get_start_point(), 25*units::cm);
    TVector3 dir11 = shower->cal_dir_3vector(shower->get_start_point(), 6*units::cm);
	
    for (auto it1 = showers.begin(); it1 != showers.end(); it1++){
      WCPPID::WCShower *shower1 = *it1;
      if (updated_showers.find(shower1) != updated_showers.end()) continue;
      if (shower1->get_start_segment()->get_particle_type()!=11) continue;
      if (shower1->get_start_vertex().second==2) {
	if (del_showers.find(shower1) != del_showers.end()) continue;
	TVector3 dir2(shower1->get_start_point().x - shower->get_start_vertex().first->get_fit_pt().x,
		      shower1->get_start_point().y - shower->get_start_vertex().first->get_fit_pt().y,
		      shower1->get_start_point().z - shower->get_start_vertex().first->get_fit_pt().z);
	TVector3 dir3 = shower1->cal_dir_3vector(shower1->get_start_point(), 25*units::cm);
	
	//
	
	if (dir1.Angle(dir2)/3.1415926*180. < 10 && dir1.Angle(dir3)/3.1415926*180. < 20){
	  //	std::cout << shower1->get_kine_charge() << " " << dir1.Angle(dir2)/3.1415926*180. << " " << dir1.Angle(dir3)/3.1415926*180. << " " << 	dir11.Angle(dir2)/3.1415926*180. << " " << dir11.Angle(dir2)/3.1415926*180. << std::endl;
	  shower->add_shower(shower1);
	  shower->calculate_kinematics();
	  double kine_charge = cal_kine_charge(shower);
	  shower->set_kine_charge(kine_charge);
	  shower->set_flag_kinematics(true);
	  del_showers.insert(shower1);
	}
      }
    }
  }


  
  for (auto it1 = del_showers.begin(); it1!=del_showers.end(); it1 ++){
    WCPPID::WCShower *shower1 = *it1;
    showers.erase(find(showers.begin(), showers.end(), shower1));
    delete shower1;
  }

  /* for (auto it1 = showers.begin(); it1 != showers.end(); it1++){ */
  /*     WCPPID::WCShower *shower1 = *it1; */
  /*     if (shower1->get_start_segment()->get_particle_type()!=11) continue; */
  /*     if (shower1->get_start_vertex().second>2) continue; */
  /*     std::cout << shower1->get_kine_charge() << " " << shower1->get_start_segment()->get_id() << " " << shower1->get_start_vertex().second << std::endl; */
  /* } */
  
  if (map_merge_seg_shower.size()>0)     update_shower_maps();

  
  
  examine_shower_1();
}


void WCPPID::NeutrinoID::examine_shower_1(){
  // multiple gamma, no single gamma is large enough ...

  // if there is a EM shower connecting to it with a large energy, no need to proceed.
  bool flag_skip = false;
  auto it = map_vertex_to_shower.find(main_vertex);
  if (it != map_vertex_to_shower.end()){
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      WCPPID::WCShower  *shower = *it1;
      if (shower->get_start_segment()->get_particle_type() !=11) continue;
      if (shower->get_start_vertex().second != 1) continue;
      double energy =0;
      if (shower->get_kine_best() != 0){ 
	  energy = shower->get_kine_best();
	}else{
	  energy = shower->get_kine_charge();
	}
      if (energy > 80*units::MeV) flag_skip = true;     
    }
  }
  bool flag_added = false;
  
  if (!flag_skip){
    std::set<WCPPID::WCShower* > used_showers; // used showers ...
    std::map<WCPPID::ProtoSegment*, std::set<WCPPID::WCShower*> > map_segment_showers; // association ...
    std::map<WCPPID::ProtoSegment*, WCPPID::WCShower*> map_segment_new_shower;
    std::set<WCPPID::ProtoSegment*> used_segments;
    std::set<WCPPID::WCShower* > del_showers;
      
    for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
      WCPPID::ProtoSegment *sg = *it;
      if ((!sg->is_dir_weak()) || sg->get_particle_type()==2212 && sg->get_medium_dQ_dx()/(43e3/units::cm) > 1.6 || sg->get_particle_type()==11) continue;
    
      // form a new shower ...
      WCPPID::WCShower *shower1 = new WCPPID::WCShower();
      shower1->set_start_vertex(main_vertex, 1);
      shower1->set_start_segment(sg);
      shower1->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, used_segments);
      
      TVector3 dir1 = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      for (auto it1 = showers.begin(); it1 != showers.end(); it1++){
	WCPPID::WCShower *shower = *it1;
	double energy = shower->get_kine_charge();
	double min_dis = shower->get_closest_dis(sg);
	
	
	if (shower->get_start_vertex().second > 2 && min_dis > 3*units::cm || shower->get_start_segment()->get_particle_type() != 11) continue;
	if (shower->get_start_vertex().first == main_vertex && shower->get_start_vertex().second ==1) continue;
	Map_Proto_Vertex_Segments& map_vtx_segs = shower1->get_map_vtx_segs();
	if (shower->get_start_vertex().second==1 && energy > 80*units::MeV){
	  if (map_vtx_segs.find(shower->get_start_vertex().first) == map_vtx_segs.end()) continue;
	  else{
	    TVector3 dir3 = shower->cal_dir_3vector(shower->get_start_point(), 15*units::cm);
	    if (dir3.Angle(dir1)/3.1415926*180. > 30 ) continue;
	    //	    std::cout << dir3.Angle(dir1)/3.1415926*180. << std::endl;
	  }
	}

	
	
	if (used_showers.find(shower) != used_showers.end()) continue;
	WCP::Point shower_point = shower->get_closest_point(main_vertex->get_fit_pt()).second;
	TVector3 dir2(shower_point.x - main_vertex->get_fit_pt().x,
		      shower_point.y - main_vertex->get_fit_pt().y,
		      shower_point.z - main_vertex->get_fit_pt().z);

	double angle = dir1.Angle(dir2)/3.1415926*180.;
	
	
	if (angle < 15 && min_dis < 36*units::cm || angle < 10 && min_dis <46*units::cm || angle < 7.5){
	// consider the WCShower associated to this (also angle)
	// consider WCShower nearby (not necessarily associated)
	  map_segment_showers[sg].insert(shower);
	  used_showers.insert(shower);
	}
      } // loop over shower

      if (map_segment_showers.find(sg) != map_segment_showers.end()){
	map_segment_new_shower[sg] = shower1;
      }else{
	delete shower1;
      }
      
    } // loop over segment 

    
    for (auto it = map_segment_showers.begin(); it != map_segment_showers.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      WCPPID::WCShower *shower1 = map_segment_new_shower[sg];
      std::set<WCPPID::WCShower* >& associated_showers = it->second;
      int num_showers = associated_showers.size();
      double max_energy =0;
      double total_energy = 0;
      for (auto it1 = associated_showers.begin(); it1 != associated_showers.end(); it1++){
	WCPPID::WCShower* shower = *it1;
	double energy = shower->get_kine_charge();
	if (energy > max_energy) max_energy = energy;
	total_energy += energy;
      }

      Map_Proto_Segment_Vertices& map_seg_vtxs = shower1->get_map_seg_vtxs();
      Map_Proto_Vertex_Segments& map_vtx_segs = shower1->get_map_vtx_segs();
      double max_length = 0;
      WCPPID::ProtoSegment *max_sg = 0;
      int n_tracks = 0;
      int n_showers = 0;
      double total_length = 0;
      bool flag_good_track = false;
      for (auto it1 = map_seg_vtxs.begin(); it1 != map_seg_vtxs.end(); it1++){
	WCPPID::ProtoSegment *sg1 = it1->first;
	double length = sg1->get_length();
	double medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
	if (!sg1->is_dir_weak()) {
	  auto pair_vertices = find_vertices(sg1);
	  WCPPID::ProtoVertex *end_vertex = 0;
	  if (sg1->get_flag_dir()==1){
	    if (sg1->get_wcpt_vec().back().index == pair_vertices.first->get_wcpt().index){
	      end_vertex = pair_vertices.first;
	    }else{
	      end_vertex = pair_vertices.second;
	    }
	  }else{
	    if (sg1->get_wcpt_vec().front().index == pair_vertices.first->get_wcpt().index){
	      end_vertex = pair_vertices.first;
	    }else{
	      end_vertex = pair_vertices.second;
	    }
	  }
	  if (map_vertex_segments[end_vertex].size()>1){
	    bool flag_non_ele = false;
	    for (auto it2 = map_vertex_segments[end_vertex].begin(); it2 != map_vertex_segments[end_vertex].end(); it2++){
	      WCPPID::ProtoSegment *sg2 = *it2;
	      if(sg2 == sg1) continue;
	      if (!sg2->get_flag_shower()) flag_non_ele = true;
	    }
	    if ((!flag_non_ele) && map_vertex_segments[end_vertex].size()<=3) flag_good_track == true;
	    //	    std::cout << map_vertex_segments[end_vertex].size() << std::endl;
	  }else{
	    flag_good_track = true;
	  }
	}
	
	if (sg1->get_flag_shower()) n_showers ++;
	n_tracks ++;
	total_length += length;
	if (max_length < length) {
	  max_length = length;
	  max_sg = sg1;
	}
	//	std::cout << sg1->get_length()/units::cm << " " << sg1->get_medium_dQ_dx()/(43e3/units::cm) << " " << sg1->is_dir_weak() << " " << sg1->get_flag_shower() << std::endl;
      }

      bool flag_skip = false;
      for (auto it1 = showers.begin(); it1 != showers.end(); it1++){
	WCPPID::WCShower *shower1 = *it1;
	if (shower1->get_start_vertex().second==1 && associated_showers.find(shower1) == associated_showers.end()){
	  //5774_74_3718
	  if (map_vtx_segs.find(shower1->get_start_vertex().first)!=map_vtx_segs.end() && shower1->get_start_vertex().first != main_vertex && shower1->get_kine_charge() > 60*units::MeV)
	    flag_skip = true;
	}
      }
      

      std::cout << "kaka2: " << num_showers << " " << max_energy << " " << total_energy << " " << shower1->get_total_length()/units::cm << " " << shower1->get_num_segments() << " " << flag_good_track << " " << n_tracks << " " << n_showers << " " << max_length/units::cm << " " << total_length/n_tracks/units::cm << " " << flag_skip << std::endl;
      
      if ( total_length< 70*units::cm && (n_tracks == 1 && total_length< 60*units::cm ||
					  n_tracks == 1 && total_length < 65*units::cm && num_showers > 3 && total_energy > 150*units::MeV || total_length < n_tracks * 36*units::cm) && (total_energy > 50*units::MeV || total_energy/units::MeV > total_length/units::cm * 0.75) && (!flag_skip)){
	// for the new shower,
	shower1->get_start_segment()->set_particle_type(11);
	shower1->get_start_segment()->set_flag_avoid_muon_check(true);
	shower1->update_particle_type();
	for (auto it1 = associated_showers.begin(); it1 != associated_showers.end(); it1++){
	  del_showers.insert(*it1);
	  shower1->add_shower(*it1);
	}
	shower1->calculate_kinematics();
	double kine_charge = cal_kine_charge(shower1);
	shower1->set_kine_charge(kine_charge);
	shower1->set_flag_kinematics(true);
	showers.push_back(shower1);
	std::cout << "Create a new low-energy shower: " << kine_charge << std::endl;
	flag_added = true;
      }else{
	delete shower1;
      }
    } // loop over segment ...

    for (auto it1 = del_showers.begin(); it1!=del_showers.end(); it1 ++){
      WCPPID::WCShower *shower1 = *it1;
      showers.erase(find(showers.begin(), showers.end(), shower1));
      delete shower1;
    }
    update_shower_maps();
  }

  if (!flag_added){
    std::map<WCPPID::WCShower*, std::set<WCPPID::WCShower*> > map_shower_showers;
    auto it = map_vertex_to_shower.find(main_vertex);
    if (it != map_vertex_to_shower.end()){
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::WCShower  *shower = *it1;
	if (shower->get_start_segment()->get_particle_type() !=11) continue;
	if (shower->get_start_vertex().second != 1) continue;

	TVector3 dir1 = shower->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
	for (auto it2 = showers.begin(); it2 != showers.end(); it2++){
	  WCPPID::WCShower *shower1 = *it2;
	  double energy = shower1->get_kine_charge();
	  double min_dis = shower1->get_closest_dis(shower->get_start_segment());
	
	  if (shower1->get_start_vertex().second > 2 && min_dis > 3*units::cm || shower1->get_start_segment()->get_particle_type() != 11) continue;
	  if (shower1->get_start_vertex().first == main_vertex && shower1->get_start_vertex().second ==1) continue;
	  Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
	  if (shower1->get_start_vertex().second==1 && map_vtx_segs.find(shower1->get_start_vertex().first) == map_vtx_segs.end()) continue;

	  if (shower1->get_total_length() < 3*units::cm) continue;
	  
	  WCP::Point shower_point = shower1->get_closest_point(main_vertex->get_fit_pt()).second;
	  TVector3 dir2(shower_point.x - main_vertex->get_fit_pt().x,
			shower_point.y - main_vertex->get_fit_pt().y,
			shower_point.z - main_vertex->get_fit_pt().z);
	  
	  double angle = dir1.Angle(dir2)/3.1415926*180.;
	  TVector3 dir3 = shower1->cal_dir_3vector(shower1->get_start_point(), 30*units::cm);
	  double angle1 = dir2.Angle(dir3)/3.1415926*180.;

	  if (angle < 15 && angle1 < 15 && min_dis < 28*units::cm){
	    //	    std::cout << angle << " " << angle1 << " " << min_dis/units::cm << " " << std::endl;
	    map_shower_showers[shower].insert(shower1);
	  }
	}
      }
    }



    std::set<WCPPID::WCShower*> del_showers;
    WCPPID::WCShower *max_shower = 0;
    double max_energy = 0;
    for (auto it1 = map_shower_showers.begin(); it1 != map_shower_showers.end(); it1++){
      double acc_energy =0;
      WCPPID::WCShower *shower = it1->first;
      if (shower->get_kine_best() != 0){
	acc_energy += shower->get_kine_best();
      }else{
	acc_energy += shower->get_kine_charge();
      }
      for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++){
	WCPPID::WCShower *shower1 = *it2;
	if (shower1->get_kine_best() != 0){
	  acc_energy += shower1->get_kine_best();
	}else{
	  acc_energy += shower1->get_kine_charge();
	}
      }
      if (acc_energy > max_energy){
	max_energy = acc_energy;
	max_shower = shower;
      }
    }

    if (max_shower != 0){
      for (auto it2 = map_shower_showers[max_shower].begin(); it2 != map_shower_showers[max_shower].end(); it2++){
	WCPPID::WCShower *shower1 = *it2;
	max_shower->add_shower(shower1);
	del_showers.insert(shower1);
      }
      max_shower->calculate_kinematics();
      max_shower->get_start_segment()->set_flag_avoid_muon_check(true);
      double kine_charge = cal_kine_charge(max_shower);
      max_shower->set_kine_charge(kine_charge);
      max_shower->set_flag_kinematics(true);
    }
    
    for (auto it1 = del_showers.begin(); it1!=del_showers.end(); it1 ++){
      WCPPID::WCShower *shower1 = *it1;
      showers.erase(find(showers.begin(), showers.end(), shower1));
      delete shower1;
    }
    update_shower_maps();
    
  }

  
  

}
