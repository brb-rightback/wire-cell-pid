bool WCPPID::NeutrinoID::cosmic_tagger(){
  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);

  bool flag_cosmic = false;
  
  bool flag_cosmic_1 = false;
  bool flag_cosmic_2 = false;
  bool flag_cosmic_3 = false;
  bool flag_cosmic_4 = false;
  bool flag_cosmic_5 = false;
  bool flag_cosmic_6 = false;
  bool flag_cosmic_7 = false;
  bool flag_cosmic_8 = false;
  bool flag_cosmic_9 = false;
  bool flag_cosmic_10 = false;
  
  bool flag_print = false;
  
  double tmp_dis = sqrt(pow(main_vertex->get_fit_pt().x  - main_vertex->get_wcpt().x,2) + pow(main_vertex->get_fit_pt().y  - main_vertex->get_wcpt().y,2) + pow(main_vertex->get_fit_pt().z  - main_vertex->get_wcpt().z,2));
  //  std::cout << tmp_dis/units::cm << std::endl;
  // if main vertex is outside the fiducial volume
  Point test_p = main_vertex->get_fit_pt();
  if (tmp_dis > 5*units::cm){
    test_p.x = main_vertex->get_wcpt().x;
    test_p.y = main_vertex->get_wcpt().y;
    test_p.z = main_vertex->get_wcpt().z;
  }
  
  std::vector<double> stm_tol_vec =     {-1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm};
  //  std::cout << fid->inside_fiducial_volume(test_p, offset_x) << " " << fid->inside_fiducial_volume(test_p, offset_x, &stm_tol_vec) << std::endl;
  
  if (! fid->inside_fiducial_volume(test_p, offset_x, &stm_tol_vec)){  // Xiangpan's boundary
    // if (!fid->inside_fiducial_volume(test_p, offset_x)){ // original boundary
    flag_cosmic_1 = true;
    if (flag_print) std::cout << "Xin_A: " << test_p << " " << true << std::endl;
  }

  //  std::cout << "Entering cosmic tagger" << std::endl;
  // single muon, nothing there ...
  {
    double max_length = 0;
    WCPPID::ProtoSegment *muon = 0;
    WCPPID::WCShower *long_muon = 0;
    Int_t valid_tracks = 0;
    for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end() ;it++){
      WCPPID::ProtoSegment *sg = *it;
      double length = sg->get_length();
      double medium_dQ_dx = sg->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
      
      if (sg->get_particle_type() == 13 || (!sg->get_flag_shower()) && medium_dQ_dx < dQ_dx_cut * 1.05 * 43e3/units::cm && sg->get_particle_type()!=211){
	if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()){
	  double long_muon_length = map_segment_in_shower[sg]->get_total_track_length();
	  if (long_muon_length > max_length){
	    long_muon = map_segment_in_shower[sg];
	    max_length = long_muon_length;
	    muon = 0; // reset the muon
	  }
	}else{
	  if (length > max_length){
	    max_length = length;
	    muon = sg;
	    long_muon = 0;
	  }
	}
      }
    }
    for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
      WCPPID::ProtoSegment *sg = *it;
      double length = sg->get_length();
      if (sg == muon ) continue;
      if (long_muon != 0)
	if (sg == long_muon->get_start_segment()) continue;
      
      //  std::cout << "Xin_3: " << length/units::cm << " " << map_vertex_segments[main_vertex].size() << " " << sg->get_particle_type() << " " << sg->is_dir_weak() << std::endl;
      if (length > 2.5*units::cm || length > 0.9*units::cm && sg->get_particle_type()==2212 || (!sg->is_dir_weak())) {
	//	std::cout << "Xin_3: " << sg->get_medium_dQ_dx()/(43e3/units::cm) << " " << length/units::cm << " " << sg->get_particle_type() << " " << sg->is_dir_weak() << std::endl;
	if (length < 15*units::cm && sg->get_medium_dQ_dx()/(43e3/units::cm)  < 0.75 && sg->is_dir_weak()) continue;
	valid_tracks ++;
      }
    }

    //    std::cout << main_vertex << " " << muon << " " << long_muon << " " << valid_tracks << std::endl;

    int connected_showers = 0;
    // high energy showers ...
    for (auto it = map_vertex_to_shower[main_vertex].begin(); it != map_vertex_to_shower[main_vertex].end(); it++){
    //    for (auto it = showers.begin(); it != showers.end(); it++){
      WCPPID::WCShower *shower = *it;
      double Eshower = 0;
      if (shower->get_kine_best() != 0){ 
	Eshower = shower->get_kine_best();
      }else{
	Eshower = shower->get_kine_charge();
      }
      if (shower->get_start_vertex().second >2) continue;
      if (shower == long_muon) continue;
      if (Eshower > 150*units::MeV && (!bad_reconstruction(shower)) )
	valid_tracks ++;
      if (shower->get_start_vertex().second == 1 && Eshower > 70*units::MeV && (!bad_reconstruction(shower)) &&
	  shower->get_start_segment()->get_particle_type()!=13)
	connected_showers ++;
      //      std::cout << Eshower/units::MeV << " " << bad_reconstruction(shower) << " " << shower->get_start_vertex().second << " " << shower->get_start_segment()->get_cluster_id() << " " << valid_tracks << std::endl;
    }

    //    std::cout << " " << map_vertex_to_shower[main_vertex].size() << std::endl;

    if (valid_tracks == 0){
      if ( muon != 0){
	WCPPID::ProtoSegment *sg = muon;
	WCPPID::ProtoVertex *other_vtx = find_other_vertex(sg, main_vertex);
	Point test_p1 = other_vtx->get_fit_pt();
	TVector3 dir = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
	bool flag_inside = fid->inside_fiducial_volume(test_p1, offset_x);
	
	double dQ_dx_front = 0;
	double dQ_dx_end = 0;
	if (sg->get_wcpt_vec().front().index==main_vertex->get_wcpt().index){
	  dQ_dx_front = sg->get_medium_dQ_dx(0,10);
	  dQ_dx_end =sg->get_medium_dQ_dx(int(sg->get_point_vec().size())-10,sg->get_point_vec().size());
	}else{
	  dQ_dx_end = sg->get_medium_dQ_dx(0,10);
	  dQ_dx_front = sg->get_medium_dQ_dx(int(sg->get_point_vec().size())-10,sg->get_point_vec().size());
	}
	int n_muon_tracks = calculate_num_daughter_tracks(main_vertex, sg, false, 3*units::cm).first;
	double total_shower_length = calculate_num_daughter_showers(main_vertex, sg).second;
	//      std::cout << "Xin_1: " << sg->get_length()/units::cm << " " << sg->get_particle_type() << " " << fid->inside_fiducial_volume(test_p1, offset_x) << " " << sg->is_dir_weak() << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << dir.Angle(dir_vertical)/3.1415926*180. << " " << dir.Angle(dir_drift)/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;


	//	std::cout << "Xin_BB: " << sg->get_particle_type() << " " <<  n_muon_tracks << " " << total_shower_length/units::cm << " " <<sg->get_length()/units::cm << " " << main_vertex->get_fit_pt().x/units::cm << " " << main_vertex->get_fit_pt().y/units::cm << " " << main_vertex->get_fit_pt().z/units::cm << " " << fid->inside_fiducial_volume(test_p1, offset_x) << " " << sg->is_dir_weak() << " " << dir.Theta()/3.1415926*180. << " " << dir.Phi()/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;
	
	if (sg->get_particle_type() == 13 && n_muon_tracks <=2 && total_shower_length < 40*units::cm && ( (!flag_inside) && dir.Angle(dir_beam)/3.1415926*180. > 40  || (flag_inside &&  (sg->is_dir_weak() && (!(dQ_dx_end > 1.4 *43e3/units::cm && dQ_dx_end > 1.2 * dQ_dx_front)) || dir.Angle(dir_beam)/3.1415926*180. > 60)) )
	    && (dir.Theta()/3.1415926*180.>=100. || fabs(fabs(dir.Phi()/3.1415926*180.)-90)<=50)
	    ){

	  flag_cosmic_2 = true;
	  if (flag_print) std::cout << "Xin_B: " << sg->get_length()/units::cm << " " << main_vertex->get_fit_pt().x/units::cm << " " << main_vertex->get_fit_pt().y/units::cm << " " << main_vertex->get_fit_pt().z/units::cm << " " << fid->inside_fiducial_volume(test_p1, offset_x) << " " << sg->is_dir_weak() << " " << dir.Theta()/3.1415926*180. << " " << dir.Phi()/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;
	}/* else if (sg->get_particle_type()==13 && flag_inside && (!sg->is_dir_weak()) && */
      /* 		(dir.Angle(dir_drift)/3.1415926*180. > 165     // prolonged towards anode ... */
      /* 		 || dir.Angle(dir_vertical)/3.1415926*180. > 150) // vertical down ... */
      /* 		){ */
      /* 	neutrino_type |= 1UL << 1; */
      /* 	return true; */
      /* } */
      }else if (long_muon != 0){
	
	WCPPID::ProtoSegment *start_sg = long_muon->get_start_segment();
	auto pair_result = long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon);
	WCPPID::ProtoVertex *other_vtx = pair_result.second;
	WCPPID::ProtoSegment *last_sg = pair_result.first;

	/* for (auto it1 = segments_in_long_muon.begin(); it1 != segments_in_long_muon.end(); it1++){ */
	/*   std::cout << *it1 << " " << (*it1)->get_id() << std::endl; */
	/* } */
	
	/* std::cout << "Xin: " << main_vertex << " " << other_vtx <<  " " << last_sg <<" " << segments_in_long_muon.size() <<  std::endl; */
	
	Point test_p1 = other_vtx->get_fit_pt(); 
	TVector3 dir = long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm); 
	bool flag_inside = fid->inside_fiducial_volume(test_p1, offset_x); 
	
	double dQ_dx_front = 0;
	if (start_sg->get_wcpt_vec().front().index == main_vertex->get_wcpt().index){
	  dQ_dx_front = start_sg->get_medium_dQ_dx(0,10); 
	}else{
	  dQ_dx_front = start_sg->get_medium_dQ_dx(int(start_sg->get_point_vec().size())-10,start_sg->get_point_vec().size()); 
	}
	
	double dQ_dx_end = 0; 
	if (last_sg->get_wcpt_vec().front().index==other_vtx->get_wcpt().index){ 
	  dQ_dx_end = last_sg->get_medium_dQ_dx(0,10);
	}else{ 
	  dQ_dx_end = last_sg->get_medium_dQ_dx(int(last_sg->get_point_vec().size())-10,last_sg->get_point_vec().size());
	}

	//	std::cout << "Xin_BB_(long): " << long_muon->get_total_track_length()/units::cm << " " << main_vertex->get_fit_pt().x/units::cm << " " << main_vertex->get_fit_pt().y/units::cm << " " << main_vertex->get_fit_pt().z/units::cm << " " << fid->inside_fiducial_volume(test_p1, offset_x) << " " << last_sg->is_dir_weak() << " " << dir.Theta()/3.1415926*180. << " " << dir.Phi()/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;
	
	if ( ( (!flag_inside) && dir.Angle(dir_beam)/3.1415926*180. > 40
	       || (flag_inside && 
		   (last_sg->is_dir_weak() && (!(dQ_dx_end > 1.4 *43e3/units::cm && dQ_dx_end > 1.2 * dQ_dx_front)) || dir.Angle(dir_beam)/3.1415926*180. > 60)))
	     && (dir.Theta()/3.1415926*180.>=100. || fabs(fabs(dir.Phi()/3.1415926*180.)-90)<=50)
	     ){
	  flag_cosmic_3 = true;
	  if (flag_print) std::cout << "Xin_B_(long): " << long_muon->get_total_track_length()/units::cm << " " << main_vertex->get_fit_pt().x/units::cm << " " << main_vertex->get_fit_pt().y/units::cm << " " << main_vertex->get_fit_pt().z/units::cm << " " << fid->inside_fiducial_volume(test_p1, offset_x) << " " << last_sg->is_dir_weak() << " " << dir.Theta()/3.1415926*180. << " " << dir.Phi()/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;
	}
      }
    }

    //    std::cout << muon << " " << long_muon << std::endl;

    if (muon != 0 ){
      WCPPID::ProtoSegment *sg =muon;
      WCPPID::ProtoVertex *other_vtx = find_other_vertex(sg, main_vertex);
      Point test_p1 = other_vtx->get_fit_pt();
      TVector3 dir = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      bool flag_inside = fid->inside_fiducial_volume(test_p1, offset_x);

      //std::cout << "Xin: " << sg->get_length()/units::cm << " " << flag_inside << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << connected_showers << std::endl;
      if ( (!flag_inside) && dir.Angle(dir_beam)/3.1415926*180. > 100 && connected_showers ==0){
	flag_cosmic_4 = true;
	if (flag_print)  std::cout << "Xin_F: " << sg->get_length()/units::cm << " " << flag_inside << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << connected_showers << std::endl;
      }
    }else if (long_muon != 0){
      WCPPID::ProtoSegment *start_sg = long_muon->get_start_segment();
      auto pair_result = long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon);
      WCPPID::ProtoVertex *other_vtx = pair_result.second;
      Point test_p1 = other_vtx->get_fit_pt(); 
      TVector3 dir = long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm); 
      bool flag_inside = fid->inside_fiducial_volume(test_p1, offset_x); 

      //std::cout << "Xin_(long): " << flag_inside << " " << connected_showers << " " << dir.Angle(dir_beam)/3.1415926*180. << std::endl;
      
      if ( (!flag_inside) && dir.Angle(dir_beam)/3.1415926*180. > 100 && connected_showers ==0){
	flag_cosmic_5 = true;
	if (flag_print)  std::cout << "Xin_F_(long): " << long_muon->get_total_track_length()/units::cm << " " << flag_inside << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << connected_showers << std::endl;
      }
    }
  }


  // stopped muon with a  Michel electron ...
  {
    WCPPID::WCShower *michel_ele = 0;
    double michel_energy = 0;
    WCPPID::ProtoSegment *muon = 0;
    WCPPID::ProtoSegment *muon_2nd = 0;
    WCPPID::WCShower* long_muon = 0;
    double max_length = 0;
    
    Int_t valid_tracks = 0;
    
    for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
      WCPPID::ProtoSegment *sg = *it;
      double length = sg->get_length();
      
      if (sg->get_flag_shower()){
	for (auto it1 = showers.begin(); it1!=showers.end(); it1++){
	  WCPPID::WCShower *shower = *it1;
	  // calculate energy 
	  if (shower->get_start_segment() == sg){
	    double tmp_Eshower = 0;
	    if (shower->get_kine_best() != 0){ 
	      tmp_Eshower = shower->get_kine_best();
	    }else{
	      tmp_Eshower = shower->get_kine_charge();
	    }
	    if (tmp_Eshower > michel_energy){
	      michel_energy = tmp_Eshower;
	      michel_ele = shower;
	    }
	  }
	}
      }else{
	if (sg->get_particle_type()==13 ){
	  if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()){
	    double long_muon_length = map_segment_in_shower[sg]->get_total_track_length();
	    if (long_muon_length > max_length){
	      max_length = long_muon_length;
	      long_muon = map_segment_in_shower[sg];
	      muon = 0;
	    }
	  }else{
	    if (length > max_length){
	      max_length = length;
	      muon = sg;
	      long_muon = 0;
	    }
	  }
	}
      }
    }

    max_length = 0;
    for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
      WCPPID::ProtoSegment *sg = *it;
      double length = sg->get_length();
      if (sg == muon ) continue;
      if (long_muon != 0)
	if (sg == long_muon->get_start_segment()) continue;
      if (michel_ele != 0)
      	if (sg == michel_ele->get_start_segment()) continue;
      if (sg->get_flag_shower()) continue;

      if ((!sg->get_flag_shower()) && (sg->get_particle_type()!=2212)){
	if (length > max_length){
	  max_length = length;
	  muon_2nd = sg;
	}
      }
    }
    
    for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
      WCPPID::ProtoSegment *sg = *it;
      double length = sg->get_length();
      if (sg == muon ) continue;
      if (long_muon != 0)
	if (sg == long_muon->get_start_segment()) continue;
      if (michel_ele != 0)
      	if (sg == michel_ele->get_start_segment()) continue;
      if (sg->get_flag_shower()) continue;
      if (sg == muon_2nd) continue;
      
      if (length > 2.5*units::cm || length > 0.9*units::cm && sg->get_particle_type()==2212 || (!sg->is_dir_weak())) {
      	if (length < 15*units::cm && sg->get_medium_dQ_dx()/(43e3/units::cm)  < 0.75 && sg->is_dir_weak()) continue;
	valid_tracks ++;
      }
    }

  
    
    for (auto it = map_vertex_to_shower[main_vertex].begin(); it != map_vertex_to_shower[main_vertex].end(); it++){
      WCPPID::WCShower *shower = *it;
      if (shower == michel_ele) continue;
      if (shower == long_muon) continue;
      double Eshower = 0;
      if (shower->get_kine_best() != 0){ 
	Eshower = shower->get_kine_best();
      }else{
	Eshower = shower->get_kine_charge();
      }
      if (shower->get_start_vertex().second >2) continue;
      if (Eshower > 150*units::MeV && (!bad_reconstruction(shower))) valid_tracks ++;
      //      std::cout << Eshower/units::MeV << " " << bad_reconstruction(shower) << " " << shower->get_start_vertex().second << " " << shower->get_start_segment()->get_cluster_id() << std::endl;
    }

    //    std::cout << long_muon << " " << muon << " " << michel_ele << " " << muon_2nd << " " << valid_tracks << std::endl;
    if ( (muon!=0 || long_muon !=0) && michel_ele != 0 && muon_2nd !=0 ){
      {
	double length = muon_2nd->get_length();
	if (length > 2.5*units::cm || length > 0.9*units::cm && muon_2nd->get_particle_type()==2212 || (!muon_2nd->is_dir_weak())) {
	  if (length < 15*units::cm && muon_2nd->get_medium_dQ_dx()/(43e3/units::cm)  < 0.75 && muon_2nd->is_dir_weak()){
	  }else	  valid_tracks ++;
	}
      }
      
      TVector3 dir1;
      if (muon != 0){
	dir1 = muon->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      }else if (long_muon !=0){
	dir1 = long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm);
      }
      TVector3 dir2 = michel_ele->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      TVector3 dir3 = muon_2nd->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      double Eshower = 0;
      if (michel_ele->get_kine_best() != 0){ 
	Eshower = michel_ele->get_kine_best();
      }else{
	Eshower = michel_ele->get_kine_charge();
      }
      //      std::cout << dir1.Angle(dir2)/3.1415926*180. << " " << dir1.Angle(dir3)/3.1415926*180. << std::endl;
      if (Eshower < 25*units::cm && dir1.Angle(dir2)/3.1415926*180. > 170){
	valid_tracks -- ;
      }else if (muon_2nd->get_length() < 5*units::cm && dir1.Angle(dir3)/3.1415926*180. > 170){
	valid_tracks -- ;
      } 
      
      if (dir1.Angle(dir3)/3.1415926*180. > 175 && valid_tracks <= 2){
	if (!(muon_2nd->get_length() < 5*units::cm && dir1.Angle(dir3)/3.1415926*180. > 170)) valid_tracks --;
	for (auto it1 = map_vertex_segments[main_vertex].begin(); it1 != map_vertex_segments[main_vertex].end(); it1++){
	  WCPPID::ProtoSegment *sg1 = *it1;
	  if (sg1 == muon || sg1 == muon_2nd) continue;
	  if (sg1 == michel_ele->get_start_segment()) continue;
	  if (long_muon != 0)
	    if (sg1 == long_muon->get_start_segment()) continue;
	  if (sg1->is_dir_weak() && sg1->get_length() < 5*units::cm){
	    valid_tracks --;
	  }
	}
      }
    
      if (valid_tracks <0) valid_tracks = 0;
      //std::cout << Eshower/units::MeV << std::endl;
    }
   
    //    std::cout << muon << " " << long_muon << " " << michel_ele << " " << muon_2nd << " " << valid_tracks << std::endl;

    if ( (muon!=0 || long_muon!=0) && (michel_ele !=0 && valid_tracks == 0 || valid_tracks == 0 && muon_2nd !=0 && michel_ele == 0)){
      // muon related ..
      double dQ_dx_front = 0;
      double dQ_dx_end = 0;
      TVector3 dir;
      bool flag_inside;
      bool flag_weak_dir = false;
      double muon_length;
      int n_muon_tracks = 0;
      double total_shower_length; 
      if (muon != 0){
	WCPPID::ProtoVertex *other_vtx = find_other_vertex(muon, main_vertex);
	Point test_p1 = other_vtx->get_fit_pt();
	dir = muon->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
	flag_inside = fid->inside_fiducial_volume(test_p1, offset_x);
	if (muon->get_wcpt_vec().front().index==main_vertex->get_wcpt().index){
	  dQ_dx_front = muon->get_medium_dQ_dx(0,10);
	  dQ_dx_end = muon->get_medium_dQ_dx(int(muon->get_point_vec().size())-10,muon->get_point_vec().size());
	}else{
	  dQ_dx_end = muon->get_medium_dQ_dx(0,10);
	  dQ_dx_front = muon->get_medium_dQ_dx(int(muon->get_point_vec().size())-10,muon->get_point_vec().size());
	}
	flag_weak_dir = muon->is_dir_weak();
	muon_length = muon->get_length();
	n_muon_tracks = calculate_num_daughter_tracks(main_vertex, muon, false, 3*units::cm).first;
	total_shower_length = calculate_num_daughter_showers(main_vertex, muon).second;
      }else if (long_muon != 0){
	WCPPID::ProtoSegment *start_sg = long_muon->get_start_segment();
        auto pair_result = long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon);
        WCPPID::ProtoVertex *other_vtx = pair_result.second;
        WCPPID::ProtoSegment *last_sg = pair_result.first;
	Point test_p1 = other_vtx->get_fit_pt(); 
	dir = long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm); 
	flag_inside = fid->inside_fiducial_volume(test_p1, offset_x);
	if (start_sg->get_wcpt_vec().front().index == main_vertex->get_wcpt().index){
          dQ_dx_front = start_sg->get_medium_dQ_dx(0,10); 
        }else{
          dQ_dx_front = start_sg->get_medium_dQ_dx(int(start_sg->get_point_vec().size())-10,start_sg->get_point_vec().size()); 
        }
	if (last_sg->get_wcpt_vec().front().index==other_vtx->get_wcpt().index){ 
          dQ_dx_end = last_sg->get_medium_dQ_dx(0,10);
        }else{ 
          dQ_dx_end = last_sg->get_medium_dQ_dx(int(last_sg->get_point_vec().size())-10,last_sg->get_point_vec().size());
        }
	flag_weak_dir = last_sg->is_dir_weak();
	muon_length = long_muon->get_total_track_length();
	n_muon_tracks = 1;
	total_shower_length = 0;
      }

      bool flag_sec = false;
      if (michel_ele !=0){
	double Eshower = 0;
	if (michel_ele->get_kine_best() != 0){ 
	  Eshower = michel_ele->get_kine_best();
	}else{
	  Eshower = michel_ele->get_kine_charge();
	}
	if (Eshower < 70*units::MeV) flag_sec = true;
	else{
	  TVector3 dir2 = michel_ele->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm);
	  double tmp_angle = dir.Angle(dir2)/3.1415926*180.;
	  if (muon_length > 120*units::cm && tmp_angle > 155 || muon_length > 20*units::cm && tmp_angle > 175 || muon_length > 60*units::cm && tmp_angle > 167.5){
	    flag_sec = true;
	    if (flag_print)  std::cout << "Xin_I: " << tmp_angle << " " << muon_length/units::cm << std::endl;
	  }

	}
	
      }else if (muon_2nd !=0){
	if (muon_2nd->is_dir_weak() && muon_2nd->get_length() < 8*units::cm)  flag_sec = true;

	if (muon_2nd->is_dir_weak() && (!fid->inside_fiducial_volume(find_other_vertex(muon_2nd, main_vertex)->get_fit_pt(), offset_x))){
	  TVector3 dir2 = muon_2nd->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
	  if (dir.Angle(dir2)/3.1415926*180. > 170){
	    if (flag_print)       std::cout << "Xin_G: " << std::endl;
	    flag_cosmic_6 = true;
	  }
	}
	//	std::cout << dir.Angle(dir2)/3.1415926*180. << " " << muon_2nd->is_dir_weak() << " " << fid->inside_fiducial_volume(find_other_vertex(muon_2nd, main_vertex)->get_fit_pt(), offset_x) << std::endl;
      }

      //      std::cout << "Xin_C: " << flag_sec << " " << muon_length/units::cm  << " " << flag_inside << " " << flag_weak_dir << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << dir.Angle(dir_vertical)/3.1415926*180. << " " << dir.Angle(dir_drift)/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << " " << n_muon_tracks << " " << total_shower_length/units::cm << std::endl;
      //      std::cout << "Xin_CC: " << flag_sec << " " << muon_length/units::cm  << " " << flag_inside << " " << flag_weak_dir << " " << dir.Theta()/3.1415926*180. << " " << dir.Phi()/3.1415926*180. << " " << dir.Angle(dir_drift)/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;
      
      if (flag_sec && n_muon_tracks <=2 && total_shower_length < 40*units::cm){
	if (( (!flag_inside) && dir.Angle(dir_beam)/3.1415926*180. > 40
	      || (flag_inside  && (flag_weak_dir && (!(dQ_dx_end > 1.4 *43e3/units::cm && dQ_dx_end > 1.2 * dQ_dx_front)) || dir.Angle(dir_beam)/3.1415926*180. > 60)))
	    && (dir.Theta()/3.1415926*180.>=100. || fabs(fabs(dir.Phi()/3.1415926*180.)-90)<=50)
	    ){
	  if (flag_print)       std::cout << "Xin_C: " << flag_sec << " " << muon_length/units::cm  << " " << flag_inside << " " << flag_weak_dir << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << dir.Angle(dir_vertical)/3.1415926*180. << " " << dir.Angle(dir_drift)/3.1415926*180. << " " << dQ_dx_front/(43e3/units::cm) << " " << dQ_dx_end/(43e3/units::cm) << std::endl;
	  flag_cosmic_7 = true;
	}
	/* else if (muon->get_particle_type()==13 && flag_inside && (!muon->is_dir_weak()) && */
	/* 	  (dir.Angle(dir_drift)/3.1415926*180. > 165     // prolonged towards anode ... */
	/* 	   || dir.Angle(dir_vertical)/3.1415926*180. > 150) // vertical down ... */
	/* 	  ){ */
	/*   neutrino_type |= 1UL << 1; */
	/*   return true; */
	/* } */
      }
           
    }else if ((muon !=0 || long_muon !=0) && valid_tracks == 1 ){ // two things ...
      TVector3 dir; 
      double muon_length; 
      if (muon != 0){
      	dir = muon->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      	muon_length = muon->get_length();
      }else if (long_muon != 0){
      	dir = long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm);
      	muon_length = long_muon->get_total_track_length();
      }
      double acc_length = 0;
      bool flag_out = false;
      for (auto it1 = map_vertex_segments[main_vertex].begin(); it1!= map_vertex_segments[main_vertex].end(); it1++){
      	WCPPID::ProtoSegment *sg1 = *it1;
      	if (sg1 == muon) continue;
	if (long_muon !=0)
	  if( sg1 == long_muon->get_start_segment()) continue;
      	TVector3 dir2 = sg1->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
	if (dir.Angle(dir2)/3.1415926*180.>165){
	  WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg1, main_vertex);
	  if (!fid->inside_fiducial_volume(other_vertex->get_fit_pt(), offset_x)) flag_out = true;
	}else{
	  acc_length += sg1->get_length();
	}
	//      	std::cout << dir.Angle(dir2)/3.1415926*180.  << " " << sg1->get_length()/units::cm << " " << muon_length/units::cm<< std::endl;
      }
      if (flag_out && muon_length > 100*units::cm && acc_length < 12*units::cm){
	if (flag_print)       std::cout << "Xin_H: " << std::endl;
	flag_cosmic_8 = true;
      }
      //      std::cout << flag_out << " " << acc_length/units::cm << " " << muon_length/units::cm << std::endl;
    }
  } 
  


  
  {
    bool flag_main_cluster = true;
    
    std::map<int, PointVector> map_cluster_id_points;
    std::map<int, Point> map_cluster_id_high_point;
    std::map<int, int> map_cluster_id_shower_points;
    std::map<int, double> map_cluster_id_length;
    std::set<int> cluster_id_set;
    
    // calculate PCA for large clusters ...
    int num_small_pieces = 0;
    double acc_small_length = 0;
    for (auto it = map_cluster_length.begin(); it != map_cluster_length.end(); it++){
      map_cluster_id_length[it->first->get_cluster_id()] = it->second;
      if (it->second > 3*units::cm){
	cluster_id_set.insert(it->first->get_cluster_id());
	map_cluster_id_shower_points[it->first->get_cluster_id()] = 0;
      }else{
	if (it->first->get_point_cloud()->get_cloud().pts.front().y > 50*units::cm){
	  num_small_pieces ++;
	  acc_small_length += it->second;
	}
      }
    }

    for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      if (cluster_id_set.find(sg->get_cluster_id()) == cluster_id_set.end()) continue;
      PointVector& pts = sg->get_point_vec();
      
      if (pts.size()<=2) continue;
      if (map_cluster_id_high_point.find(sg->get_cluster_id()) == map_cluster_id_high_point.end()) map_cluster_id_high_point[sg->get_cluster_id()] = pts.front();
      for (size_t i=1; i+1<pts.size();i++){
	map_cluster_id_points[sg->get_cluster_id()].push_back(pts.at(i));
	if (pts.at(i).y > map_cluster_id_high_point[sg->get_cluster_id()].y) map_cluster_id_high_point[sg->get_cluster_id()] = pts.at(i);
	  
      }
      if (sg->get_flag_shower())
	map_cluster_id_shower_points[sg->get_cluster_id()] += pts.size();
    }
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
      WCPPID::ProtoVertex *vtx = it->first;
      if (cluster_id_set.find(vtx->get_cluster_id()) == cluster_id_set.end()) continue;
      if (vtx->get_fit_pt().y > map_cluster_id_high_point[vtx->get_cluster_id()].y) map_cluster_id_high_point[vtx->get_cluster_id()] = vtx->get_fit_pt();
      map_cluster_id_points[vtx->get_cluster_id()].push_back(vtx->get_fit_pt());
    }

    int num_cosmic = 0;
    double acc_cosmic_length = 0;
    double acc_total_length = 0;
    double highest_y = -100*units::cm;
    double max_length = 0;
    int num_showers = 0;
    for (auto it = map_cluster_id_points.begin(); it != map_cluster_id_points.end(); it++){
      double angle_cosmic, angle_beam;
       
      if (it->first == main_vertex->get_cluster_id() &&
	  (map_cluster_id_shower_points[it->first] * 1./it->second.size() > 0.7 && map_cluster_id_length[it->first] < 45 *units::cm
	   ||   map_cluster_id_shower_points[it->first] * 1./it->second.size() <0.7 && map_cluster_id_length[it->first] > 40 *units::cm )
	  || map_cluster_id_length[it->first] > 60 *units::cm
	  ){
	Point vector(0,0,0);
	for (size_t i=0; i!= it->second.size();i++){
	  vector.x += it->second.at(i).x - main_vertex->get_fit_pt().x;
	  vector.y += it->second.at(i).y - main_vertex->get_fit_pt().y;
	  vector.z += it->second.at(i).z - main_vertex->get_fit_pt().z;
	}
	TVector3 dir(vector.x, vector.y, vector.z);
	angle_beam = dir.Angle(dir_beam)/3.1415926*180.;
	if (angle_beam > 90) angle_beam = 180 - angle_beam;

	angle_cosmic = 180 - dir.Angle(dir_vertical)/3.1415926*180.;
      }else{
	main_cluster->Calc_PCA(it->second);
	auto vector = main_cluster->get_PCA_axis(0);
	TVector3 dir(vector.x, vector.y, vector.z);
	angle_cosmic = dir.Angle(dir_vertical)/3.1415926*180.;
	angle_beam = dir.Angle(dir_beam)/3.1415926*180.;
	
	if (angle_cosmic > 90) angle_cosmic = 180 - angle_cosmic;
	if (angle_beam > 90) angle_beam = 180 - angle_beam;
      }
      
      // main track is not cosmic like ...
      if ( it->first == main_vertex->get_cluster_id() && map_cluster_id_shower_points[it->first] * 1./it->second.size() <0.3 && map_cluster_id_length[it->first] > 20*units::cm && angle_cosmic > 40){
	flag_main_cluster = false;
      }else if (it->first == main_vertex->get_cluster_id() && map_cluster_id_length[it->first] > 80*units::cm && angle_cosmic > 25){
	flag_main_cluster = false;
      }
      
      //      std::cout << angle_cosmic << " " << angle_beam << " " << map_cluster_id_length[it->first]/units::cm << " " << map_cluster_id_shower_points[it->first] * 1./it->second.size() <<  " " << map_cluster_id_high_point[it->first] << std::endl;
      if (map_cluster_id_high_point[it->first].y > highest_y) highest_y = map_cluster_id_high_point[it->first].y;
      
      if (map_cluster_id_shower_points[it->first] * 1./it->second.size() > 0.7){
	if (angle_cosmic < 30 && angle_beam > 30){
	  acc_cosmic_length += map_cluster_id_length[it->first];
	  num_cosmic ++;
	  if (max_length < map_cluster_id_length[it->first]) max_length = map_cluster_id_length[it->first];
	}else if (angle_cosmic < 35 && angle_beam > 40){
	  acc_cosmic_length += map_cluster_id_length[it->first];
	  num_cosmic ++;
	  if (max_length < map_cluster_id_length[it->first]) max_length = map_cluster_id_length[it->first];
	}
      }else{
	if (angle_cosmic < 20 || angle_cosmic < 30 && highest_y > 100*units::cm && it->first == main_vertex->get_cluster_id() ){
	  acc_cosmic_length += map_cluster_id_length[it->first];
	  num_cosmic ++;
	  if (max_length < map_cluster_id_length[it->first]) max_length = map_cluster_id_length[it->first];
	}
      }
      acc_total_length += map_cluster_id_length[it->first];
      //      std::cout << it->first << " " << map_cluster_id_length[it->first]/units::cm << std::endl;
    }
    
    //    std::cout << num_cosmic << " Xin " << acc_cosmic_length/units::cm << " " << acc_small_length/units::cm << " " << acc_total_length/units::cm << " " << main_vertex->get_fit_pt() << " " << showers.size() << " " << highest_y/units::cm << " " << max_length/units::cm << " " << flag_main_cluster << std::endl;   

    bool flagp_cosmic = false;
    if ( (num_cosmic >2 && acc_cosmic_length + acc_small_length  > 0.55 * (acc_total_length)||
	  num_cosmic >=2 && acc_cosmic_length + acc_small_length  > 0.7 * (acc_total_length) ||
	  num_cosmic ==1 && acc_cosmic_length + acc_small_length  > 0.625 * (acc_total_length) && highest_y > 102*units::cm
	   ) && main_vertex->get_fit_pt().y > 0 && flag_main_cluster && highest_y > 80*units::cm 
	 ) {
      flagp_cosmic = true;
    }
    if (num_cosmic == 1 && acc_cosmic_length > 100*units::cm){
      flagp_cosmic = false;
    }

    if (flagp_cosmic){
      // try to see nueCC events ...
      int n_solid_tracks = 0;
      int n_direct_showers = 0;
      double energy_direct_showers = 0;
      int n_main_showers = 0;
      double energy_main_showers = 0;
      int n_indirect_showers = 0;
      double energy_indirect_showers = 0;
      for (auto it2 = map_vertex_segments[main_vertex].begin(); it2 != map_vertex_segments[main_vertex].end(); it2++){
	WCPPID::ProtoSegment *sg1 = *it2;
	if (sg1->get_flag_shower()) continue;
	double length = sg1->get_length();
	if (sg1->is_dir_weak() && length > 10*units::cm || (!sg1->is_dir_weak())) n_solid_tracks ++;
      }

      for (auto it2 = showers.begin(); it2 != showers.end(); it2++){
	WCPPID::WCShower *shower = *it2;
	if (shower->get_start_segment()->get_particle_type()!=11) continue;
	double Eshower = 0;
	if (shower->get_kine_best() != 0){ 
	  Eshower = shower->get_kine_best();
	}else{
	  Eshower = shower->get_kine_charge();
	}
	if (Eshower > 60*units::MeV){
	  if (shower->get_start_vertex().first == main_vertex && shower->get_start_vertex().second==1) {
	    n_main_showers ++;
	    energy_main_showers += Eshower;
	  
	  }

	  if (shower->get_start_vertex().first->get_cluster_id() == main_vertex->get_cluster_id() && shower->get_start_vertex().second==1){
	    n_direct_showers ++;
	    energy_direct_showers += Eshower;
	  }

	  if (shower->get_start_vertex().second>1){
	    n_indirect_showers ++;
	    energy_indirect_showers += Eshower;
	  }
	}
	
      }

      tagger_info.cosmic_n_solid_tracks = n_solid_tracks;
      tagger_info.cosmic_energy_main_showers = energy_main_showers/units::MeV;
      tagger_info.cosmic_energy_direct_showers = energy_direct_showers/units::MeV;
      tagger_info.cosmic_energy_indirect_showers = energy_indirect_showers/units::MeV;
      tagger_info.cosmic_n_direct_showers = n_direct_showers;
      tagger_info.cosmic_n_indirect_showers = n_indirect_showers;
      tagger_info.cosmic_n_main_showers = n_main_showers;
      tagger_info.cosmic_filled = 1;
      
      if ((n_solid_tracks >0 && energy_main_showers > 80*units::MeV || energy_main_showers > 200*units::MeV) && (energy_indirect_showers < energy_main_showers * 0.6 || n_solid_tracks > 0 && energy_indirect_showers < energy_main_showers)){
	flagp_cosmic = false;
	
	//	std::cout << "kaka: " << n_solid_tracks << " " <<  n_direct_showers << " " << energy_direct_showers << " " << n_main_showers  << " " <<  energy_main_showers << " " << n_indirect_showers << " " <<  energy_indirect_showers  << " " << tagger_info.cosmic_filled << std::endl;

	tagger_info.cosmic_flag = true;
      }else{
	tagger_info.cosmic_flag = false;
      }

      //      std::cout << "cosmic: " << (!flag_cosmic) << " " << n_solid_tracks << " " << energy_main_showers/units::MeV << " " << energy_direct_showers/units::MeV << " " << energy_indirect_showers/units::MeV << " " << n_direct_showers << " " << n_indirect_showers << " "  << n_main_showers << " " << 1 << std::endl;
      
    }
    
    
    if (flagp_cosmic){
      flag_cosmic_9 = true;
      if (flag_print)  std::cout << num_cosmic << " Xin_D: " << acc_cosmic_length/units::cm << " " << acc_small_length/units::cm << " " << acc_total_length/units::cm << " " << main_vertex->get_fit_pt() << " " << showers.size() << " " << highest_y/units::cm << " " << max_length/units::cm << " " << flag_main_cluster << std::endl;     
      
    }
  } // else

  // test front end ...
  for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != main_cluster->get_cluster_id()) continue;

    if (!fid->inside_fiducial_volume(vtx->get_fit_pt(), offset_x) && vtx->get_fit_pt().z < 15*units::cm){
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::ProtoSegment *sg = *it1;
	TVector3 dir = sg->cal_dir_3vector(vtx->get_fit_pt(),15*units::cm);

	double angle_beam = dir.Angle(dir_beam)/3.1415926*180.;
	if (angle_beam > 90) angle_beam = 180 - angle_beam;
	double length = sg->get_length();
	
	if ((!sg->get_flag_shower()) && sg->is_dir_weak() && angle_beam < 25 && length > 10*units::cm){
	  flag_cosmic_10 = true;
	  //      std::cout << true << std::endl;
	  if (flag_print) std::cout << "Xin_E: " << true << std::endl;
	}
	//std::cout << dir.Angle(dir_beam)/3.1415926*180. << " " << sg->get_flag_dir() << " " << sg->is_dir_weak() << std::endl;
      }
    }
    //    std::cout << vtx->get_fit_pt() << " " << fid->inside_fiducial_volume(vtx->get_fit_pt(), offset_x) << std::endl;
  }

  std::cout << flag_cosmic_1 << " " << flag_cosmic_2 << " " << flag_cosmic_3 << " " << flag_cosmic_4 << " " << flag_cosmic_5 << " " << flag_cosmic_6 << " " << flag_cosmic_7 << " " << flag_cosmic_8 << " " << flag_cosmic_9 << " " << flag_cosmic_10 << std::endl;

  flag_cosmic = flag_cosmic_1 || flag_cosmic_2 || flag_cosmic_3 || flag_cosmic_4 || flag_cosmic_5
    || flag_cosmic_6 || flag_cosmic_7 || flag_cosmic_8 || flag_cosmic_9 || flag_cosmic_10;
    
  
  
  if (flag_cosmic){
    neutrino_type |= 1UL << 1;
  }
  
  return flag_cosmic;
}
