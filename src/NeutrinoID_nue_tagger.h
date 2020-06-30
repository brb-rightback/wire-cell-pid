
bool WCPPID::NeutrinoID::nue_tagger(double muon_length){
  bool flag_nue = false;
  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);
  bool flag_print = true;
  bool flag_print_detail = true;
  
  // check main_vertex ...
  {
    WCPPID::WCShower* max_shower = 0;
    double max_energy = 0;
    std::set<WCShower*> good_showers;
        
    auto it = map_vertex_to_shower.find(main_vertex);
    if (it != map_vertex_to_shower.end()){
      bool flag_single_shower = false;
      if (map_vertex_segments[main_vertex].size() == 1) flag_single_shower = true;
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::WCShower *shower = *it1;
	WCPPID::ProtoSegment *sg = shower->get_start_segment();

	//	if (flag_print) std::cout << "Qian_I: " << sg->get_id() << " " << sg->get_particle_type() << " " << shower->get_kine_charge()/units::MeV << " " << shower->get_kine_best()/units::MeV << std::endl;
	
	if (sg->get_particle_type()!=11) continue;

	double tmp_energy =0;
	int n_3seg = 0;
	int n_segs = 0;
	{
	  if (shower->get_kine_best()==0){
	    tmp_energy = shower->get_kine_charge();
	  }else{
	    tmp_energy = shower->get_kine_best();
	  }
	  Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
	  Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
	  n_segs = shower->get_num_main_segments();
	  for (auto it1 = map_vtx_segs.begin(); it1 != map_vtx_segs.end(); it1++){
	    WCPPID::ProtoVertex *vtx1 = it1->first;
	    if (vtx1->get_cluster_id() != sg->get_cluster_id()) continue;
	    if (it1->second.size()>=3) n_3seg ++;
	  }
	}
	
	// higher than 70 MeV, and connect to the vertex ...
	if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end() && (tmp_energy > 80*units::MeV || tmp_energy > 60*units::MeV && n_segs>=3 && n_3seg>0 )){

	  double E_shower = 0;
	  if (shower->get_kine_best() != 0){ 
	    E_shower = shower->get_kine_best();
	  }else{
	    E_shower = shower->get_kine_charge();
	  }
	  if (E_shower > max_energy)  {
	    max_shower = shower;
	    max_energy = E_shower;
	  }


	  good_showers.insert(shower);

	  
	 
	} // 80 MeV, and connected
      } // loop over showers

      // std::cout << max_shower << " " << good_showers.size() << " " << *good_showers.begin() << std::endl;
      
      if (good_showers.find(max_shower)!=good_showers.end() && max_shower !=0){
	Point test_p;
	WCPPID::ProtoSegment *sg = max_shower->get_start_segment();
	if (main_vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
	  test_p = sg->get_point_vec().front();
	}else{
	  test_p = sg->get_point_vec().back();
	}
	TVector3 dir = max_shower->cal_dir_3vector(test_p, 15*units::cm);
	flag_nue = true;

	int num_valid_tracks = 0;
	for (auto it2 = map_vertex_segments[main_vertex].begin(); it2 != map_vertex_segments[main_vertex].end(); it2++){
	  WCPPID::ProtoSegment *sg1 = *it2;
	  if (sg1 == max_shower->get_start_segment()) continue;
	  if ((!sg1->get_flag_shower()) && (sg1->get_length() > 8*units::cm || sg1->is_dir_weak() && sg1->get_length() > 5*units::cm)) num_valid_tracks ++;
	  //	    std::cout << sg1->get_length()/units::cm << " " << sg1->is_dir_weak() << std::endl;
	}
	
	if (flag_nue){
	  WCPPID::ProtoSegment *sg = max_shower->get_start_segment();
	  
	  auto pair_result = gap_identification(main_vertex, sg, flag_single_shower, num_valid_tracks, max_energy);
	  
	  // std::cout << pair_result.first << " " << pair_result.second << std::endl;
	  
	  if (!pair_result.first){ // gap id
	    int mip_id;
	    if (flag_single_shower && max_energy < 400*units::MeV){
	      mip_id = mip_identification(main_vertex, sg, max_shower, flag_single_shower, true, flag_print_detail); // dQ/dx id
	    }else{
	      mip_id = mip_identification(main_vertex, sg, max_shower, flag_single_shower, false, flag_print_detail); // dQ/dx id
	    }
	    
	    if (flag_print) std::cout << "Qian_B: " << max_energy/units::MeV << " " << flag_single_shower << " " << mip_id << std::endl;
	    if (mip_id == -1){
	      flag_nue = false;
	    } else if (mip_id == 0){
	      if (!pi0_identification(main_vertex, sg, max_shower)){ // pi0 identification ...
	      }else{
		flag_nue = false;
		std::cout << "Qian_D: pi0 found " << max_energy/units::MeV << std::endl;
	      }
	    }
	  }else{ // gap id
	    flag_nue = false;
	    if (flag_print) std::cout << "Qian_C: gap founded " << max_energy/units::MeV << " " << pair_result.second << std::endl; 
	  }
	}

	// bad reconstruction
	if (flag_nue && flag_single_shower && stem_direction(max_shower, max_energy, flag_print_detail)){ // single shower ...
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_B: " << max_energy << std::endl;
	} // single shower ...

	if (flag_nue){
	  
	  bool flag1 = bad_reconstruction(max_shower) ;	// bad reconstruction (mis ided track as shower)
	  bool flag2 = bad_reconstruction_1(max_shower, flag_single_shower, num_valid_tracks); // bad reconstruction ... ( stem does not match with shower direction)
	  bool flag3 = low_energy_overlapping(max_shower); // low energy overlapping situation 
	  if ( flag1 || flag2 || flag3 ){  
	    flag_nue = false;
	    if (flag_print) std::cout << "QXin_G: " << max_energy << " " << flag1 << " " << flag2 << " " << flag3 << std::endl;
	  }
	}
	if (flag_nue){
	  bool flag1 = bad_reconstruction_2(main_vertex, max_shower, flag_print_detail);
	  bool flag2 = bad_reconstruction_3(main_vertex, max_shower, flag_print_detail);
	  bool flag3 = high_energy_overlapping(max_shower, flag_print_detail);
	  if ( flag1 || flag2 || flag3
	       // shower cone, and main cluster vs. others ...
	       ) {
	    flag_nue = false;
	    if (flag_print) std::cout << "QXin_H: " << max_energy << " " << flag1 << " " << flag2 << " " << flag3 << std::endl;
	  }
	}
	// vertex inside shower ...
	if (vertex_inside_shower(max_shower, flag_print_detail) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_O: " << max_energy << std::endl;
	}
	// stem length
	if (flag_nue && stem_length(max_shower, max_energy, flag_print_detail)){ // long stem is not preferred ...
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_F: " << max_energy << std::endl;
	}
	if (broken_muon_id(max_shower, flag_print_detail) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_K: " << max_energy << std::endl;
	}
	//overclustering...
	if (track_overclustering(max_shower, flag_print_detail) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_L: " << max_energy << std::endl;
	}
	
	// kinematics ...
	if (flag_nue && compare_muon_energy(max_shower, max_energy, muon_length, flag_print_detail)){ // compare with muon energy ...
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_E: " << max_energy << std::endl;
	}
	if (shower_to_wall(max_shower, max_energy, flag_single_shower) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_J: " << max_energy << std::endl;
	}
	// angular cut ...
	
	if (angular_cut(max_shower, max_energy, dir.Angle(dir_beam)/3.1415926*180.) && flag_nue){
	  //	    max_energy < 650*units::MeV && dir.Angle(dir_beam)/3.1415926*180. > 135 ||
	  //  max_energy>= 650*units::MeV && dir.Angle(dir_beam)/3.1415926*180. > 90){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_M: " << max_energy << " " << dir.Angle(dir_beam)/3.1415926*180. << std::endl;
	}
	

	// single shower
	if (flag_single_shower && single_shower_pio_tagger(max_shower) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_I: " << max_energy << std::endl;
	}
	if (flag_nue && multiple_showers(max_shower, max_energy, flag_print_detail)){ // no multiple EM showers
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_D: " << max_energy << std::endl;
	}
	if (other_showers(max_shower, flag_single_shower, flag_print_detail) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_N: " << max_energy << std::endl;
	}
	if (single_shower(max_shower, flag_single_shower, flag_print_detail) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_P: " << max_energy << std::endl;
	}

	// low-energy michel
	if (low_energy_michel(max_shower) && flag_nue){
	  flag_nue = false;
	  if (flag_print) std::cout << "QXin_Q: " << max_energy << std::endl;
	}
	
	//	std::cout <<  max_shower->get_kine_range() << " " << max_shower->get_kine_dQdx() << " " << max_shower->get_kine_charge() << " " << max_shower->std::endl;
	

	std::cout << "Xin: " << good_showers.size() << " " << flag_nue << " " << " " << flag_single_shower << " " << max_energy/units::MeV << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << fabs(3.1415926/2. - dir.Angle(dir_drift))/3.1415926*180. << " " << dir.Angle(dir_vertical)/3.1415926*180. << " " << max_shower->get_total_length(max_shower->get_start_segment()->get_cluster_id())/units::cm << " " << max_shower->get_total_length()/units::cm << std::endl;
      }
    } // has a shower ...
  } // main vertex
  
  
  if (flag_nue){
    neutrino_type |= 1UL << 5; //nue
  }

  return flag_nue;
}


bool WCPPID::NeutrinoID::low_energy_michel(WCPPID::WCShower* shower, bool flag_print){
  bool flag_bad = false;
  double E_range = shower->get_kine_range();
  double E_dQdx  = shower->get_kine_dQdx();
  double E_charge = shower->get_kine_charge();

  if (E_range == 0){
    E_range = shower->get_start_segment()->cal_kine_range(shower->get_total_length(shower->get_start_segment()->get_cluster_id()));
  }

  int n_3seg = 0;
  int n_segs = 0;
  {
    // 7003_1226_61350  
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
    n_segs = shower->get_num_main_segments();
    for (auto it1 = map_vtx_segs.begin(); it1 != map_vtx_segs.end(); it1++){
      WCPPID::ProtoVertex *vtx1 = it1->first;
      if (vtx1->get_cluster_id() != shower->get_start_segment()->get_cluster_id()) continue;
      if (it1->second.size()>=3) n_3seg ++;
    }
  }
  
  
  if (shower->get_total_length() < 25*units::cm && shower->get_total_length(shower->get_start_segment()->get_cluster_id()) > 0.75 * shower->get_total_length() && n_3seg == 0 ||
      shower->get_total_length() < 18*units::cm && shower->get_total_length(shower->get_start_segment()->get_cluster_id()) > 0.75 * shower->get_total_length() && n_3seg >0 ) flag_bad = true;

  // 7004_1291_64560 + 7026_879_43995
  if (E_charge < 100*units::MeV && E_dQdx < 0.7 * E_charge && shower->get_num_segments() == shower->get_num_main_segments()){
    flag_bad = true;
  }
  
  //  std::cout << "qaqa2: " << E_range << " " << E_dQdx << " " << E_charge << " " << shower->get_start_segment()->get_length()/units::cm << " " << shower->get_total_length(shower->get_start_segment()->get_cluster_id())/units::cm << " " << shower->get_total_length()/units::cm << " " << shower->get_num_segments() << " " << shower->get_num_main_segments() << " " << flag_bad << std::endl;
  
  
  //  std::cout <<  max_shower->get_kine_range() << " " << max_shower->get_kine_dQdx() << " " << max_shower->get_kine_charge() << " " << max_shower->std::endl;
  return flag_bad;
}

bool WCPPID::NeutrinoID::single_shower(WCPPID::WCShower*shower, bool flag_single_shower, bool flag_print){
  bool flag_bad = false;
  
  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);

  WCPPID::ProtoVertex *vertex = shower->get_start_vertex().first;
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  Point vertex_point;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
  }else{
    vertex_point = sg->get_point_vec().back();
  }
  
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }

  std::vector<double> vec_dQ_dx = shower->get_stem_dQ_dx(shower->get_start_vertex().first, shower->get_start_segment(), 20);
  double max_dQ_dx = 0;
  for (size_t i= 0; i!=vec_dQ_dx.size();i++){
    if (vec_dQ_dx.at(i) > max_dQ_dx) max_dQ_dx = vec_dQ_dx.at(i);
    if (i==2) break;
  }
  
  
  
  
  TVector3 dir_shower;
  if (shower->get_start_segment()->get_length() > 12*units::cm){
    dir_shower = shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);
  }else{
    dir_shower = shower->cal_dir_3vector(vertex_point,15*units::cm);
  }
  if (fabs(dir_shower.Angle(dir_drift)/3.1415926*180.-90)<10 || Eshower > 800*units::MeV) dir_shower = shower->cal_dir_3vector(vertex_point,25*units::cm);
  dir_shower = dir_shower.Unit();

  TVector3 dir_shower1 = shower->cal_dir_3vector(vertex_point, 15*units::cm);


  double angle_beam = dir_shower.Angle(dir_beam)/3.1415926*180.;
  double angle_vertical = dir_vertical.Angle(dir_shower)/3.1415926*180.;
  double angle_drift = fabs(3.1415926/2.-dir_shower.Angle(dir_drift))/3.1415926*180.;

  
  if (flag_single_shower){
    // 6572_18_948 + 7020_473_23679
    if (Eshower < 600*units::MeV && (shower->get_total_length(sg->get_cluster_id()) < 0.1 *  shower->get_total_length() && angle_beam > 40 || angle_beam <=40 && shower->get_total_length(sg->get_cluster_id()) < 0.08 *  shower->get_total_length()) ) flag_bad = true;
    if ( (angle_vertical < 20 || angle_vertical > 160) && angle_beam > 80 && max_dQ_dx < 3.0) flag_bad = true;
    if ((angle_beam > 15 || dir_shower1.Angle(dir_beam)/3.1415926*180. > 15) && (angle_drift < 5 || fabs(3.1415926/2.-dir_shower1.Angle(dir_drift))/3.1415926*180. < 5) && Eshower < 1200*units::MeV && max_dQ_dx < 2.4) flag_bad = true;

    if (flag_print && flag_bad) std::cout << "Xin_P_1: " << Eshower << " " << angle_beam << " " << angle_vertical << " " << angle_drift << " " << flag_single_shower << " " << max_dQ_dx << " " << shower->get_total_length(sg->get_cluster_id())/units::cm << " " << shower->get_total_length()/units::cm << std::endl;
  }else{
    int num_valid_tracks = 0;
    double max_length = 0;
    for (auto it2 = map_vertex_segments[vertex].begin(); it2 != map_vertex_segments[vertex].end(); it2++){
      WCPPID::ProtoSegment *sg1 = *it2;
      if (sg1 == shower->get_start_segment()) continue;

      double medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
      double length = sg1->get_length();

      // std::cout << sg1->get_medium_dQ_dx()/(43e3/units::cm) << " " << sg1->get_length()/units::cm << " " << sg1->get_flag_shower() << std::endl;

      // 7022_110_5542
      if ((!sg1->get_flag_shower()) && ((!sg1->is_dir_weak()) || sg1->is_dir_weak() && length > 4.2*units::cm ||
					length > 0.6*units::cm && medium_dQ_dx > 3 ||
					length > 1.6*units::cm && medium_dQ_dx > 2.2 
					)){
	num_valid_tracks ++;
	if (length > max_length) max_length = length;
      }else{
	double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
	if (medium_dQ_dx > dQ_dx_cut){
	  num_valid_tracks ++;
	  if (length > max_length) max_length = length;
	}
      }
    }
    // single shower ...
    if (num_valid_tracks == 0 && angle_beam > 60 && map_vertex_segments[vertex].size() <=3) flag_bad = true;
    // 7017_969_48490
    if (Eshower < 200*units::MeV && angle_vertical < 10 && angle_drift < 5 && max_length < 5*units::cm && num_valid_tracks <=1) flag_bad = true;
    
    if (flag_print && flag_bad) std::cout << "Xin_P_0: " << Eshower << " " << num_valid_tracks << " " << angle_beam << " " << angle_drift << " " << angle_vertical << " " << map_vertex_segments[vertex].size() << std::endl;
  }
  //  std::cout << angle_beam << " " << angle_drift << std::endl;
  
  return flag_bad;
}

bool WCPPID::NeutrinoID::angular_cut(WCPPID::WCShower* shower, double energy, double angle, bool flag_print){
  bool flag_bad = false;

  WCPPID::ProtoVertex *vertex = shower->get_start_vertex().first;
  double acc_forward_length = 0;
  double acc_forward_length1 = 0;
  double acc_backward_length = 0;
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  Point vertex_point;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
  }else{
    vertex_point = sg->get_point_vec().back();
  }
  TVector3 dir_beam(0,0,1);

  TVector3 dir_shower  = shower->cal_dir_3vector(vertex_point, 30*units::cm);
  double max_angle = 0;
  double max_length = 0;
  
  for (auto it = map_vertex_segments[vertex].begin(); it != map_vertex_segments[vertex].end(); it++){
    WCPPID::ProtoSegment *sg1 = *it;
    TVector3 dir1 = sg1->cal_dir_3vector(vertex_point, 15*units::cm);
    double angle = dir1.Angle(dir_beam)/3.1415926*180.;
    
    auto pair_result = calculate_num_daughter_tracks(vertex, sg1, true);
    
    if (angle >90){
      acc_backward_length += pair_result.second;
    }else{
      acc_forward_length += pair_result.second;
      if (angle < 85)
	acc_forward_length1 += pair_result.second;
    }

    
    double angle1 = dir1.Angle(dir_shower)/3.1415926*180.;
    if (angle1 > max_angle) {
      max_angle = angle1;
      max_length = pair_result.second;
    }
    
  }

  //
  /* 7019_194_9724 577.291 170.663 96.7857 130.036 145.594 10.3261 */
/* 5508_69_3463 255.781 144.239 4.52824 83.3248 135.595 4.52824 */
/* 7001_952_47602 432.217 151.1 0 41.7747 6.46575 41.7747 */
/* 7010_1468_73409 133.766 138.688 1.78977 31.7006 167.432 1.78977 */
/* 7014_390_19528 220.937 154.107 49.0435 89.8684 153.851 20.6747 */
/* 7026_446_22327 240.436 141.147 29.2358 86.4908 163.239 8.54359 */
/* 7049_1351_67556 217.105 151.947 19.4925 42.8508 176.869 15.2018 */
/* 7021_461_23078 354.798 139.929 46.9243 16.9505 171.139 46.9243 */
  
  if (energy < 650*units::MeV && angle > 160 || // no matter what ...
      energy < 650*units::MeV && angle > 135 && max_angle > 170 && max_length > 12*units::cm || 
      energy < 650*units::MeV && angle > 135 && (acc_forward_length < 0.8 * acc_backward_length  && acc_forward_length < 15*units::cm || acc_forward_length < 0.6 * acc_backward_length && acc_forward_length >= 15*units::cm && acc_backward_length >= 80*units::cm || acc_forward_length < 0.4 * acc_backward_length && acc_forward_length >= 15*units::cm )|| // 7004_1546_77336
      energy < 650*units::MeV && (acc_forward_length  == 0 || acc_forward_length < 0.03 * acc_backward_length || acc_forward_length1 ==0 || acc_forward_length1 < 0.03*acc_backward_length) && acc_backward_length >0  && (acc_backward_length - shower->get_total_length(shower->get_start_segment()->get_cluster_id()) > acc_forward_length && angle > 90 || angle <=90) || // 7049_1213_60652
      energy>= 650*units::MeV && angle > 90){
    flag_bad = true;
    std::cout << "Xin_M_0: " << energy << " " << angle << " " << acc_forward_length/units::cm << " " << acc_backward_length/units::cm << " " << max_angle << " " << max_length/units::cm << " " << acc_forward_length1/units::cm << " " << shower->get_total_length(shower->get_start_segment()->get_cluster_id())/units::cm << " " << std::endl;
  }


  
  
  // touch boundary the main one
  bool flag_main_outside = false;
  std::vector<double> stm_tol_vec =     {-1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm};
  Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
  for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
    WCPPID::ProtoVertex *vtx1 = it->first;
    if (vtx1->get_cluster_id() != shower->get_start_segment()->get_cluster_id()) continue;
    if (vtx1 == shower->get_start_vertex().first) continue;
    
    if (!fid->inside_fiducial_volume(vtx1->get_fit_pt(), offset_x, &stm_tol_vec)) flag_main_outside = true;
  }
  double main_length = shower->get_total_length(shower->get_start_segment()->get_cluster_id());
  double total_length = shower->get_total_length();

  if ((angle > 90 || energy < 300*units::MeV || angle > 60 && energy < 800*units::MeV) && flag_main_outside) {
    flag_bad = true;

    if (main_length < 0.5 * total_length && total_length > 80*units::cm && angle < 90) flag_bad = false;
    
    std::cout << "Xin_M_1: " << energy << " " << angle << " " << flag_main_outside << " " << main_length/units::cm << " " << total_length/units::cm << std::endl;
  }
  return flag_bad;
}




bool WCPPID::NeutrinoID::track_overclustering(WCPPID::WCShower *shower, bool flag_print){
  bool flag_bad = false;
  double Eshower = 0;
  TVector3 drift_dir(1,0,0);
  
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  WCPPID::ProtoVertex *vertex = shower->get_start_vertex().first;
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  Point vertex_point;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
  }else{
    vertex_point = sg->get_point_vec().back();
  }

  // 6462_58_2942 
  if (sg->get_length() > 80*units::cm) flag_bad = true; // stem too long ...

  // find if there are good tracks inside shower
  if (!flag_bad){
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();

    double total_length = shower->get_total_length();
    double total_length_main = shower->get_total_length(vertex->get_cluster_id());
    for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      if (sg1->get_cluster_id() != vertex->get_cluster_id()) continue;
      auto pair_vertices = find_vertices(sg1);
      if (map_vertex_segments[pair_vertices.first].size()==1 || map_vertex_segments[pair_vertices.second].size()==1){
	if (sg1->get_particle_type()!=11 && (!sg1->is_dir_weak())){
	  double dis1 = sqrt(pow(pair_vertices.first->get_fit_pt().x - vertex_point.x,2) + pow(pair_vertices.first->get_fit_pt().y - vertex_point.y,2) + pow(pair_vertices.first->get_fit_pt().z - vertex_point.z,2));
	  double dis2 = sqrt(pow(pair_vertices.second->get_fit_pt().x - vertex_point.x,2) + pow(pair_vertices.second->get_fit_pt().y - vertex_point.y,2) + pow(pair_vertices.second->get_fit_pt().z - vertex_point.z,2));
	  //
	  // distance cut to avoid close overclustering to save efficiency
	  // 7025_580_29014
	  if (std::min(dis1, dis2) > 10*units::cm && (sg1->get_length() > 0.03 * total_length_main && std::max(map_vertex_segments[pair_vertices.first].size(), map_vertex_segments[pair_vertices.second].size()) < 5 || sg1->get_length() > 0.06 * total_length_main || sg1->get_length() > 3.6*units::cm) ) {
	    flag_bad = true;
	    //std::cout << "qaqa1: " << Eshower/units::MeV << " " << sg1->get_length()/units::cm << " " << sg1->get_particle_type() << " " << sg1->is_dir_weak() << " " << total_length/units::cm << " " << total_length_main/units::cm << " " << std::min(dis1,dis2)/units::cm << " " << std::max(map_vertex_segments[pair_vertices.first].size(), map_vertex_segments[pair_vertices.second].size()) << std::endl;
	  }
	}
	double tmp_length = sg1->get_length();
	// 7055_677_33891
	if (tmp_length > 12*units::cm && (!sg1->get_flag_shower_topology()) && tmp_length > 0.3 * shower->get_total_length(sg->get_cluster_id())){
	  double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/tmp_length, 0.4234);
	  if (sg1->get_medium_dQ_dx()/(43e3/units::cm) > dQ_dx_cut * 1.1) {
	    flag_bad = true;
	    //	    std::cout << "qaqa2:  " << sg1->get_particle_type() << " " << tmp_length/units::cm << " " << sg1->get_medium_dQ_dx()/(43e3/units::cm) << " " << dQ_dx_cut << " " << shower->get_total_length(sg->get_cluster_id())/units::cm << std::endl;
	  }
	}
      }
    }
    if (flag_print && flag_bad) std::cout << "L0: " << flag_bad << std::endl;
  }


  
  if (!flag_bad){
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();

    std::set<WCPPID::ProtoSegment*> muon_segments;
    ProtoSegment *curr_muon_segment = shower->get_start_segment();
    ProtoVertex *curr_muon_vertex = find_other_vertex(curr_muon_segment, vertex);
    bool flag_continue = true;
    muon_segments.insert(curr_muon_segment);
    while (flag_continue){
      flag_continue = false;
      TVector3 dir1 = curr_muon_segment->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);

      for (auto it = map_vtx_segs[curr_muon_vertex].begin(); it != map_vtx_segs[curr_muon_vertex].end(); it++){
	WCPPID::ProtoSegment *sg1 = *it;
	if (muon_segments.find(sg1) != muon_segments.end()) continue;
	TVector3 dir2 = sg1->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);
	//	std::cout << "A: " << dir1.Angle(dir2)/3.1415926*180. << std::endl;
	if (180 - dir1.Angle(dir2)/3.1415926*180. < 15 && sg1->get_length() > 6*units::cm){
	  flag_continue = true;
	  curr_muon_segment = sg1;
	  curr_muon_vertex = find_other_vertex(sg1, curr_muon_vertex);
	  break;
	}
      }
      //std::cout << curr_muon_segment << " " << curr_muon_vertex << std::endl;
      muon_segments.insert(curr_muon_segment);
    }

    double stem_length = 0;
    for (auto it = muon_segments.begin(); it != muon_segments.end(); it++){
      stem_length += (*it)->get_length(); 
    }
    
    TVector3 dir_shower;
    if (shower->get_start_segment()->get_length() > 12*units::cm){
      dir_shower = shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);
    }else{
      dir_shower = shower->cal_dir_3vector(vertex_point,15*units::cm);
    }
    if (fabs(dir_shower.Angle(drift_dir)/3.1415926*180.-90)<10 || Eshower > 800*units::MeV) dir_shower = shower->cal_dir_3vector(vertex_point,25*units::cm);
    dir_shower = dir_shower.Unit();
    
    
    TVector3 dir1 = curr_muon_segment->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);
    dir1 = dir1.Unit();
    //std::vector<std::tuple<double, double, double> > saved_results;
    //double min_angle = 180;
    //std::tuple<double, double, double> min_info;
    for (auto it = map_vtx_segs[curr_muon_vertex].begin(); it != map_vtx_segs[curr_muon_vertex].end(); it++){
      WCPPID::ProtoSegment *sg1 = *it;
      if (muon_segments.find(sg1) != muon_segments.end()) continue;
      TVector3 dir2 = sg1->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);
      dir2 = dir2.Unit();
      double angle = 180 - dir1.Angle(dir2)/3.1415926*180;

      double max_length = 0;
      double max_length1 = 0;
      for (auto it1 = map_vtx_segs.begin(); it1 != map_vtx_segs.end() ; it1++){
	WCPPID::ProtoVertex *vtx1 = it1->first;
	if (vtx1->get_cluster_id() != curr_muon_vertex->get_cluster_id()) continue;
	TVector3 dir3(vtx1->get_fit_pt().x - curr_muon_vertex->get_fit_pt().x, vtx1->get_fit_pt().y - curr_muon_vertex->get_fit_pt().y, vtx1->get_fit_pt().z - curr_muon_vertex->get_fit_pt().z);
	/* double tmp_max_length = 0 */
	/* if (dir3.Mag() > 0){ */
	/*   for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++){ */
	/*     double length = (*it2)->get_length(); */
	/*     if ( */
	/*   } */
	/* } */
	
	if (dir3.Angle(dir2)/3.1415926*180. < 30){
	  if (Eshower > 600*units::MeV){
	    double length = dir3.Cross(dir_shower).Mag();
	    if (length > max_length) {
	      max_length = length;
	    }
	  }else{
	    double length = dir3.Cross(dir1).Mag();
	    //double length = std::min(dir3.Cross(dir1).Mag(), dir3.Cross(dir_shower).Mag());
	    if (length > max_length) {
	      max_length = length;
	    }
	  }
	}
      }

      //      std::cout << "qaqa3: " << Eshower/units::MeV << " " << stem_length/units::cm << " " << max_length/units::cm  << " " << angle << " " << fabs(drift_dir.Angle(dir2)-3.1415926/2.)/3.1415926*180. << " " << fabs(drift_dir.Angle(dir1)-3.1415926/2.)/3.1415926*180. << " " << flag_bad << std::endl;

      if (Eshower < 800*units::MeV && stem_length > 6*units::cm){
	if (stem_length > 40*units::cm || std::max(fabs(drift_dir.Angle(dir2)-3.1415926/2.)/3.1415926*180.,fabs(drift_dir.Angle(dir1)-3.1415926/2.)/3.1415926*180.) < 10 ){
	  // 7014_701_35061
	  if (max_length > 18 * units::cm) flag_bad = true;
	}else{
	  // 7006_462_23130
	  if (max_length > 17 * units::cm) flag_bad = true;
	}
      }
      if (Eshower >=800*units::MeV && stem_length > 40*units::cm && std::max(fabs(drift_dir.Angle(dir2)-3.1415926/2.)/3.1415926*180.,fabs(drift_dir.Angle(dir1)-3.1415926/2.)/3.1415926*180.) > 15 && max_length > 25*units::cm && angle > 30){
	flag_bad = true;
      }
      if (max_length > 30*units::cm && std::max(fabs(drift_dir.Angle(dir2)-3.1415926/2.)/3.1415926*180.,fabs(drift_dir.Angle(dir1)-3.1415926/2.)/3.1415926*180.) > 25 || max_length > 40*units::cm && std::max(fabs(drift_dir.Angle(dir2)-3.1415926/2.)/3.1415926*180.,fabs(drift_dir.Angle(dir1)-3.1415926/2.)/3.1415926*180.) > 20){
	flag_bad = true;
      }
      
    }

    if (flag_print && flag_bad) std::cout << "L1: " << flag_bad << std::endl;


    flag_continue = true;
    while (flag_continue){
      flag_continue = false;

      TVector3 dir1 = curr_muon_segment->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);
      
      // things connected to this vertex
      for (auto it = map_vtx_segs[curr_muon_vertex].begin(); it != map_vtx_segs[curr_muon_vertex].end(); it++){
	WCPPID::ProtoSegment *sg1 = *it;
	if (muon_segments.find(sg1) != muon_segments.end()) continue;
	TVector3 dir2 = sg1->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);
	//	std::cout << "A: " << dir1.Angle(dir2)/3.1415926*180. << std::endl;
	if (180 - dir1.Angle(dir2)/3.1415926*180. < 15 && sg1->get_length() > 6*units::cm){
	  flag_continue = true;
	  curr_muon_segment = sg1;
	  curr_muon_vertex = find_other_vertex(sg1, curr_muon_vertex);
	  break;
	}
      }

      double min_dis=1e9;
      if (!flag_continue){
	// things not connected to this vertex ...
	WCPPID::ProtoSegment *min_seg = 0;
	WCPPID::ProtoVertex *min_vtx = 0;
	for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	  WCPPID::ProtoSegment *sg1 = it->first;
	  bool flag_continue1 = false;
	  for (auto it1 = muon_segments.begin(); it1 != muon_segments.end(); it1++){
	    if (sg1->get_cluster_id() == (*it1)->get_cluster_id()) {
	      flag_continue1 = true;
	      break;
	    }
	  }
	  if (flag_continue1) continue;
	  
	  TVector3 dir2, dir3;
	  double dis1 = sqrt(pow(curr_muon_vertex->get_fit_pt().x - sg1->get_point_vec().front().x,2) + pow(curr_muon_vertex->get_fit_pt().y - sg1->get_point_vec().front().y,2) + pow(curr_muon_vertex->get_fit_pt().z - sg1->get_point_vec().front().z,2));
	  double dis2 = sqrt(pow(curr_muon_vertex->get_fit_pt().x - sg1->get_point_vec().back().x,2) + pow(curr_muon_vertex->get_fit_pt().y - sg1->get_point_vec().back().y,2) + pow(curr_muon_vertex->get_fit_pt().z - sg1->get_point_vec().back().z,2));
	  if (dis1 < dis2){
	    dir2.SetXYZ(sg1->get_point_vec().front().x - curr_muon_vertex->get_fit_pt().x,
			sg1->get_point_vec().front().y - curr_muon_vertex->get_fit_pt().y,
			sg1->get_point_vec().front().z - curr_muon_vertex->get_fit_pt().z);
	    dir3 = sg1->cal_dir_3vector(sg1->get_point_vec().front(), 15*units::cm);
	  }else{
	    dir2.SetXYZ(sg1->get_point_vec().back().x - curr_muon_vertex->get_fit_pt().x,
			sg1->get_point_vec().back().y - curr_muon_vertex->get_fit_pt().y,
			sg1->get_point_vec().back().z - curr_muon_vertex->get_fit_pt().z);
	    dir3 = sg1->cal_dir_3vector(sg1->get_point_vec().back(), 15*units::cm);
	  }
	  double angle1 = 180 - dir1.Angle(dir2)/3.1415926*180.;
	  double angle2 = dir2.Angle(dir3)/3.1415926*180.;
	  double angle3 = 180 - dir1.Angle(dir3)/3.1415926*180. ;

	  //	  std::cout << "angle: " << angle1 << " " << angle2 << " " << angle3 << " " << std::min(dis1, dis2)/units::cm << " " << sg1->get_length()/units::cm << std::endl;
	  
	  if ( (std::min(angle1, angle2) < 10 && angle1 + angle2 < 25 || angle3 < 15 && std::min(dis1, dis2) < 5*units::cm) && std::min(dis1,dis2) <  25*units::cm|| std::min(angle1, angle2) < 15 && angle3 < 30 && std::min(dis1,dis2) > 30*units::cm && sg1->get_length() > 25*units::cm && std::min(dis1,dis2) < 60*units::cm){
	    if (std::min(dis1, dis2) < min_dis){
	      min_dis = std::min(dis1, dis2);
	      min_seg = sg1;
	      auto pair_vertices = find_vertices(min_seg);
	      double dis3 = sqrt(pow(pair_vertices.first->get_fit_pt().x - curr_muon_vertex->get_fit_pt().x, 2) + pow(pair_vertices.first->get_fit_pt().y - curr_muon_vertex->get_fit_pt().y, 2) + pow(pair_vertices.first->get_fit_pt().z - curr_muon_vertex->get_fit_pt().z, 2));
	      double dis4 = sqrt(pow(pair_vertices.second->get_fit_pt().x - curr_muon_vertex->get_fit_pt().x, 2) + pow(pair_vertices.second->get_fit_pt().y - curr_muon_vertex->get_fit_pt().y, 2) + pow(pair_vertices.second->get_fit_pt().z - curr_muon_vertex->get_fit_pt().z, 2));
	      if (dis4 > dis3){
		min_vtx = pair_vertices.second;
	      }else{
		min_vtx = pair_vertices.first;
	      }
	    }
	    
	  } // if satisfy angular cut
	} // loop over segment
	//	std::cout << "min: " << min_dis/units::cm << " " << min_seg << " " << std::endl;
	
	if (min_seg != 0 ){
	  flag_continue = true;
	  curr_muon_segment = min_seg;
	  curr_muon_vertex = min_vtx;
	}
      } // if ...
      if (flag_continue){
	muon_segments.insert(curr_muon_segment);
      }
      //      std::cout << curr_muon_segment << " " << curr_muon_vertex << " " << flag_continue << std::endl;
    } // while loop

    
    stem_length = 0;
    for (auto it = muon_segments.begin(); it != muon_segments.end(); it++){
      stem_length += (*it)->get_length(); 
    }

    if (stem_length > 120*units::cm) flag_bad = true;
    //    std::cout << "qaqa: " << Eshower/units::MeV << " " << stem_length/units::cm << " " << muon_segments.size() << std::endl;
    if (flag_print && flag_bad) std::cout << "L2: " << flag_bad << std::endl;
  }


  
  if (!flag_bad){
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
    TVector3 dir1 = sg->cal_dir_3vector(vertex_point, 15*units::cm);

    for (auto it1 = map_vtx_segs.begin(); it1 != map_vtx_segs.end(); it1++){
      WCPPID::ProtoVertex *vtx1 = it1->first;
      if (vtx1->get_cluster_id()!=vertex->get_cluster_id()) continue;
      if (it1->second.size()!=1) continue;
      if (vtx1 == vertex) continue;
      TVector3 dir2(vtx1->get_fit_pt().x - vertex_point.x, vtx1->get_fit_pt().y - vertex_point.y, vtx1->get_fit_pt().z - vertex_point.z);
      WCPPID::ProtoSegment* sg1 = (*it1->second.begin());
      
      // 7054_155_7797
      if (dir2.Mag() < 10*units::cm && dir1.Angle(dir2)/3.1415926*180. > 15 && (*it1->second.begin())->get_length() > 5*units::cm && (*it1->second.begin())->get_medium_dQ_dx()/(43e3/units::cm) > 1.5 && std::min(fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. , fabs(3.1415926/2.-dir2.Angle(drift_dir))/3.1415926*180.) > 5)
	flag_bad = true;

      double angle = dir1.Angle(dir2)/3.1415926*180.;
      double angle1 = std::max(fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. , fabs(3.1415926/2.-dir2.Angle(drift_dir))/3.1415926*180.);
      double length = sg->get_length();
      double length1 = (*it1->second.begin())->get_length();

      // 7054_767_38376	 + 7020_1327_66376	
      if (dir2.Mag() < 10*units::cm && dir2.Mag() < 0.5 *  length && length > 10*units::cm && ( angle >30 && angle1 > 10 || angle > 60) && length1 > 10*units::cm ||
	  dir2.Mag() < 12*units::cm && dir2.Mag() < 0.75 * length && length > 12.5*units::cm && (angle > 20 && angle1 > 10) && length1 > 20*units::cm ) // 6882_114_5707
	flag_bad = true;

      double end_dQ_dx = 0;
      if (vtx1->get_wcpt().index == sg1->get_wcpt_vec().front().index){
	end_dQ_dx = sg1->get_medium_dQ_dx(0,6)/(43e3/units::cm);
      }else{
	end_dQ_dx = sg1->get_medium_dQ_dx(int(sg1->get_point_vec().size())-1-6, int(sg1->get_point_vec().size())-1)/(43e3/units::cm);
      }
      double medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
	
      // 6649_42_2117
      if (Eshower < 200*units::MeV && end_dQ_dx > 1.6 && sg1->get_length() > shower->get_total_length(sg1->get_cluster_id()) * 0.33 && sg1->get_flag_shower_trajectory()) flag_bad = true;
      //std::cout << "qaqa1: " << Eshower << " " << end_dQ_dx << " " << medium_dQ_dx << " " << sg1->get_length()/units::cm << " " << shower->get_total_length(sg1->get_cluster_id())/units::cm << " " << flag_bad << std::endl;
      
      //if (dir2.Mag() < sg->get_length())
      //      std::cout << "qaqa1: " << Eshower << " " << dir2.Mag()/units::cm << " " << dir1.Angle(dir2)/3.1415926*180. << " " << sg->get_length()/units::cm << " " << sg->is_shower_trajectory() << " " << fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. << " " << fabs(3.1415926/2.-dir2.Angle(drift_dir))/3.1415926*180. << " " << (*it1->second.begin())->get_medium_dQ_dx()/(43e3/units::cm) << " " << (*it1->second.begin())->get_length()/units::cm << " " << flag_bad << std::endl;
    }

    
    WCPPID::ProtoVertex *vtx1 = find_other_vertex(sg, vertex);
    dir1 = sg->cal_dir_3vector(vtx1->get_fit_pt(), 15*units::cm);

    if (map_vtx_segs[vtx1].size()>=2){
      double min_angle = 180;
      double min_length = 0;
      int min_count = 0;
      double min_angle1 = 0;
      
      double max_length = 0;
      int max_count = 0;
      double max_angle = 0;
      double max_angle1 = 0;
    
      for (auto it1 = map_vtx_segs[vtx1].begin(); it1 != map_vtx_segs[vtx1].end(); it1++){
	WCPPID::ProtoSegment *sg1 = *it1;
	if (sg1 == sg) continue;
	TVector3 dir2 = sg1->cal_dir_3vector(vtx1->get_fit_pt(), 6*units::cm);
	auto pair_result = calculate_num_daughter_tracks(vtx1, sg1, true);
	double angle = 180 - dir1.Angle(dir2)/3.1415926*180. ;
	double angle1 = fabs(3.1415926/2.-dir2.Angle(drift_dir))/3.1415926*180.;
	
	if (angle < min_angle){
	  min_angle = angle;
	  min_length = pair_result.second;
	  min_count = pair_result.first;
	  min_angle1 = angle1;
	}
	
	if (pair_result.second > max_length){
	  max_length = pair_result.second;
	  max_angle = angle;
	  max_angle1 = angle1;
	  max_count = pair_result.first;
	}
	
	//      std::cout << "qaqa: " << Eshower << " " << dir1.Angle(dir2)/3.1415926*180. << " " << fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. << " " << fabs(3.1415926/2.-dir2.Angle(drift_dir))/3.1415926*180. << " " << pair_result.first << " " << pair_result.second/units::cm << std::endl;
      }
      double dis1 = sqrt(pow(vtx1->get_fit_pt().x - vertex_point.x,2) + pow(vtx1->get_fit_pt().y - vertex_point.y,2) + pow(vtx1->get_fit_pt().z - vertex_point.z,2));
      
      // 7006_293_14699 
      if (max_angle > 25 && min_angle < max_angle && max_length > 10*units::cm && min_angle < 20 && fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. > 10 && map_vtx_segs[vtx1].size() == 3 && min_count ==1 && (Eshower >= 600*units::MeV && fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. < 40 || Eshower < 600*units::MeV && fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. < 25) ){
	flag_bad = true;
      }
      
      //      std::cout << "qaqa: " << Eshower << " " << map_vtx_segs[vtx1].size() << " " << min_angle << " " << min_length/units::cm << " " << min_count << " " << min_angle1 << " " << max_angle << " " << max_length/units::cm << " " << max_count << " " << max_angle1 << " " << fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180. << " " << flag_bad << " " << dis1/units::cm << std::endl;
    }
    if (flag_print && flag_bad) std::cout << "L3: " << flag_bad << std::endl;
  }


  
   
   
  return flag_bad;
}

bool WCPPID::NeutrinoID::broken_muon_id(WCPPID::WCShower* shower, bool flag_print){
  bool flag_bad = false;
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  
  WCPPID::ProtoVertex *vertex = shower->get_start_vertex().first;
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  Point vertex_point;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
  }else{
    vertex_point = sg->get_point_vec().back();
  }

  //  std::cout << shower->get_total_length()/units::cm << " " << Eshower << std::endl;
  
  if (Eshower < 350*units::MeV  ){
    TVector3 dir = shower->cal_dir_3vector(vertex_point, 15*units::cm);

    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
    
    std::set<WCPPID::ProtoSegment* > muon_segments;
    double add_length = 0;
    ProtoSegment *curr_muon_segment = shower->get_start_segment();
    ProtoVertex *curr_muon_vertex = find_other_vertex(curr_muon_segment, vertex);
    bool flag_continue = true;
    muon_segments.insert(curr_muon_segment);
    while (flag_continue){
      flag_continue = false;

      TVector3 dir1 = curr_muon_segment->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);

      // things connected to this vertex
      for (auto it = map_vtx_segs[curr_muon_vertex].begin(); it != map_vtx_segs[curr_muon_vertex].end(); it++){
	WCPPID::ProtoSegment *sg1 = *it;
	if (muon_segments.find(sg1) != muon_segments.end()) continue;
	TVector3 dir2 = sg1->cal_dir_3vector(curr_muon_vertex->get_fit_pt(), 15*units::cm);


	//std::cout << "A: " << dir1.Angle(dir2)/3.1415926*180. << " "  << " " << sg1->get_length()/units::cm << std::endl;

	if (180 - dir1.Angle(dir2)/3.1415926*180. < 15 && sg1->get_length() > 6*units::cm){
	  flag_continue = true;
	  curr_muon_segment = sg1;
	  curr_muon_vertex = find_other_vertex(sg1, curr_muon_vertex);
	  break;
	}
      }

      double min_dis=1e9;
      if (!flag_continue){
	// things not connected to this vertex ...
	WCPPID::ProtoSegment *min_seg = 0;
	WCPPID::ProtoVertex *min_vtx = 0;
	for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	  WCPPID::ProtoSegment *sg1 = it->first;
	  bool flag_continue1 = false;
	  for (auto it1 = muon_segments.begin(); it1 != muon_segments.end(); it1++){
	    if (sg1->get_cluster_id() == (*it1)->get_cluster_id()) {
	      flag_continue1 = true;
	      break;
	    }
	  }
	  if (flag_continue1) continue;
	  
	  TVector3 dir2, dir3;
	  double dis1 = sqrt(pow(curr_muon_vertex->get_fit_pt().x - sg1->get_point_vec().front().x,2) + pow(curr_muon_vertex->get_fit_pt().y - sg1->get_point_vec().front().y,2) + pow(curr_muon_vertex->get_fit_pt().z - sg1->get_point_vec().front().z,2));
	  double dis2 = sqrt(pow(curr_muon_vertex->get_fit_pt().x - sg1->get_point_vec().back().x,2) + pow(curr_muon_vertex->get_fit_pt().y - sg1->get_point_vec().back().y,2) + pow(curr_muon_vertex->get_fit_pt().z - sg1->get_point_vec().back().z,2));
	  if (dis1 < dis2){
	    dir2.SetXYZ(sg1->get_point_vec().front().x - curr_muon_vertex->get_fit_pt().x,
			sg1->get_point_vec().front().y - curr_muon_vertex->get_fit_pt().y,
			sg1->get_point_vec().front().z - curr_muon_vertex->get_fit_pt().z);
	    dir3 = sg1->cal_dir_3vector(sg1->get_point_vec().front(), 15*units::cm);
	  }else{
	    dir2.SetXYZ(sg1->get_point_vec().back().x - curr_muon_vertex->get_fit_pt().x,
			sg1->get_point_vec().back().y - curr_muon_vertex->get_fit_pt().y,
			sg1->get_point_vec().back().z - curr_muon_vertex->get_fit_pt().z);
	    dir3 = sg1->cal_dir_3vector(sg1->get_point_vec().back(), 15*units::cm);
	  }
	  double angle1 = 180 - dir1.Angle(dir2)/3.1415926*180.;
	  double angle2 = dir2.Angle(dir3)/3.1415926*180.;
	  double angle3 = 180 - dir1.Angle(dir3)/3.1415926*180. ;

	  //std::cout << "angle: " << angle1 << " " << angle2 << " " << angle3 << " " << std::min(dis1, dis2)/units::cm << " " << sg1->get_length()/units::cm << std::endl;
	  
	  if ( (std::min(angle1, angle2) < 10 && angle1 + angle2 < 25 || angle3 < 15 && std::min(dis1, dis2) < 5*units::cm) && std::min(dis1,dis2) <  25*units::cm || std::min(angle1, angle2) < 15 && angle3 < 30 && std::min(dis1,dis2) > 30*units::cm && sg1->get_length() > 25*units::cm && std::min(dis1,dis2) < 60*units::cm ||
	       (std::min(angle1, angle2) < 5 && angle1 + angle2 < 15 || angle3 < 10 && std::min(dis1, dis2) < 5*units::cm) && std::min(dis1,dis2) <  30*units::cm
	       ){
	    if (std::min(dis1, dis2) < min_dis){
	      min_dis = std::min(dis1, dis2);
	      min_seg = sg1;
	      auto pair_vertices = find_vertices(min_seg);
	      double dis3 = sqrt(pow(pair_vertices.first->get_fit_pt().x - curr_muon_vertex->get_fit_pt().x, 2) + pow(pair_vertices.first->get_fit_pt().y - curr_muon_vertex->get_fit_pt().y, 2) + pow(pair_vertices.first->get_fit_pt().z - curr_muon_vertex->get_fit_pt().z, 2));
	      double dis4 = sqrt(pow(pair_vertices.second->get_fit_pt().x - curr_muon_vertex->get_fit_pt().x, 2) + pow(pair_vertices.second->get_fit_pt().y - curr_muon_vertex->get_fit_pt().y, 2) + pow(pair_vertices.second->get_fit_pt().z - curr_muon_vertex->get_fit_pt().z, 2));
	      if (dis4 > dis3){
		min_vtx = pair_vertices.second;
	      }else{
		min_vtx = pair_vertices.first;
	      }
	    }
	    
	  } // if satisfy angular cut
	} // loop over segment
	//	std::cout << "min: " << min_dis/units::cm << " " << min_seg << " " << std::endl;
	
	if (min_seg != 0 ){
	  flag_continue = true;
	  curr_muon_segment = min_seg;
	  curr_muon_vertex = min_vtx;
	}
      } // if ...
      if (flag_continue){
	if (min_dis < 100*units::cm) add_length += min_dis;
	muon_segments.insert(curr_muon_segment);
      }

      
      //      std::cout << curr_muon_segment << " " << curr_muon_vertex << " " << flag_continue << std::endl;
    } // while loop

    
    double acc_length = 0;
    double acc_direct_length = 0;
    std::set<int> tmp_ids;
    for (auto it = muon_segments.begin(); it!= muon_segments.end(); it++){
      acc_length += (*it)->get_length();
      acc_direct_length += (*it)->get_direct_length();
      //      std::cout << (*it)->get_direct_length()/(*it)->get_length() << " " << (*it)->get_length() << std::endl;
      tmp_ids.insert((*it)->get_cluster_id());
    }
    TPCParams& mp = Singleton<TPCParams>::Instance();
    TGraph *g_range = mp.get_muon_r2ke();
    // check muon  ...
    double Ep = g_range->Eval((acc_length)/units::cm) * units::MeV;

    // 7020_348_17421
    double connected_length = 0;
    for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      double tmp_length = sg1->get_length();
      //      std::cout << tmp_length/units::cm << " " << sg1->get_flag_shower_topology() << " " << sg1->get_flag_shower_trajectory() << std::endl;
      if (tmp_ids.find(sg1->get_cluster_id()) != tmp_ids.end() ) connected_length += tmp_length;
    }
    
    {
      // 7022_42_2123
      TVector3 dir = sg->cal_dir_3vector(vertex_point, 15*units::cm);
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (muon_segments.find(sg1)==muon_segments.end() && sg1->get_cluster_id() == sg->get_cluster_id()){
	  auto pair_vertices = find_vertices(sg1);
	  TVector3 dir1 = sg1->cal_dir_3vector(pair_vertices.first->get_fit_pt(),15*units::cm);
	  TVector3 dir2 = sg1->cal_dir_3vector(pair_vertices.second->get_fit_pt(),15*units::cm);
	  double angle1 = std::min(dir1.Angle(dir)/3.1415926*180.,180-dir1.Angle(dir)/3.1415926*180.);
	  double angle2 = std::min(dir2.Angle(dir)/3.1415926*180.,180-dir2.Angle(dir)/3.1415926*180.);
	  if (std::min(angle1, angle2) < 10) muon_segments.insert(sg1);
	}
      }
    }
    int num_muon_main = 0;
    for (auto it = muon_segments.begin(); it != muon_segments.end(); it++){
      if ((*it)->get_cluster_id() == sg->get_cluster_id()) num_muon_main ++;
    }
    
    

    
    
    if (muon_segments.size()>1 && (Ep > Eshower * 0.55 || acc_length > 0.65 * shower->get_total_length() || connected_length > 0.95 * shower->get_total_length()) && tmp_ids.size()>1
	&& (acc_direct_length > 0.94 * acc_length )
	) {
      // 7004_989_49482	 
      if (shower->get_num_main_segments() <=3 && shower->get_num_main_segments() - num_muon_main <2 && (shower->get_num_segments() < muon_segments.size() +6 || acc_length > connected_length * 0.9 ||  acc_length > 0.8 * shower->get_total_length()) ) flag_bad = true; // 6640_173_8673
    }

    if (flag_print && flag_bad)
      std::cout << "Xin_K0: " << muon_segments.size() << " " << acc_length/units::cm << " " << add_length/units::cm << " " << connected_length/units::cm << " " << shower->get_total_length()/units::cm << " " << Ep << " " << Eshower << " " << map_seg_vtxs.size() << " " << acc_direct_length/units::cm << " " << tmp_ids.size() << " " << shower->get_num_main_segments() << " " << num_muon_main << " " << " " << shower->get_num_segments() << std::endl;
    
  } // energy cut



  
  return flag_bad;
}


bool WCPPID::NeutrinoID::shower_to_wall(WCPPID::WCShower* shower, double shower_energy, bool flag_single_shower, bool flag_print){
  //  7054_450_22517,  7049_603_30176,  7017_21_1087
  bool flag_bad = false;
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  WCPPID::ProtoVertex *vertex = shower->get_start_vertex().first;
  Point vertex_point;
  double medium_dQ_dx = 0;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
    medium_dQ_dx = sg->get_medium_dQ_dx(0,6)/(43e3/units::cm);
  }else{
    vertex_point = sg->get_point_vec().back();
    medium_dQ_dx = sg->get_medium_dQ_dx(int(sg->get_point_vec().size())-1-6, int(sg->get_point_vec().size())-1)/(43e3/units::cm);
  }
  TVector3 dir = shower->cal_dir_3vector(vertex_point, 15*units::cm);
  dir = (-1) * dir.Unit();
  double step = 1*units::cm;

  Point test_p;
  test_p.x = vertex_point.x + step * dir.X();
  test_p.y = vertex_point.y + step * dir.Y();
  test_p.z = vertex_point.z + step * dir.Z();

  std::vector<double> stm_tol_vec =     {-1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm, -1.5*units::cm};
  while (fid->inside_fiducial_volume(test_p, offset_x, &stm_tol_vec)){
    test_p.x = test_p.x + step*dir.X();
    test_p.y = test_p.y + step*dir.Y();
    test_p.z = test_p.z + step*dir.Z();
  }

  double dis = sqrt(pow(test_p.x-vertex_point.x,2) + pow(test_p.y-vertex_point.y,2) + pow(test_p.z-vertex_point.z,2));

  if (flag_single_shower){
    std::vector<double> vec_dQ_dx = shower->get_stem_dQ_dx(shower->get_start_vertex().first, shower->get_start_segment(), 20);
    double max_dQ_dx = 0;
    for (size_t i= 0; i!=vec_dQ_dx.size();i++){
      if (vec_dQ_dx.at(i) > max_dQ_dx) max_dQ_dx = vec_dQ_dx.at(i);
      if (i==2) break;
    }
    //    std::cout << max_dQ_dx << std::endl;
    
    if (shower_energy < 300*units::MeV && dis < 15*units::cm && max_dQ_dx < 2.6) flag_bad = true;
    
    // 7018_885_44275
    if ((!flag_bad) && shower_energy < 500*units::MeV){
      for (auto it = map_vertex_to_shower[vertex].begin(); it != map_vertex_to_shower[vertex].end(); it++){
	WCPPID::WCShower *shower1 = *it;
	if (shower1 == shower) continue;
	if (shower1->get_start_vertex().second > 2) continue;
	TVector3 dir1 = shower1->cal_dir_3vector(shower1->get_start_point(), 15*units::cm);
	dir1 = dir1.Unit();
	if (dir1.Angle(dir)/3.1415926*180. < 30){
	  test_p = shower1->get_end_point();
	  while(fid->inside_fiducial_volume(test_p, offset_x, &stm_tol_vec)){
	    test_p.x = test_p.x + step*dir1.X();
	    test_p.y = test_p.y + step*dir1.Y();
	    test_p.z = test_p.z + step*dir1.Z();
	  }
	  double dis1 = sqrt(pow(test_p.x - shower1->get_end_point().x,2) + pow(test_p.y - shower1->get_end_point().y,2) + pow(test_p.z - shower1->get_end_point().z,2));
	  //std::cout << dir1.Angle(dir)/3.1415926*180. << " " << shower1->get_kine_charge() << " " << shower1->get_start_vertex().second << " " << dis1/units::cm << std::endl;
	  if (dis1 < 3*units::cm) {
	    flag_bad = true;
	    std::cout << "Xin_J_4: " << dis1/units::cm << std::endl;
	  }
	}
      }
    }
    
    if ((!flag_bad) && dis < 5*units::cm){
      int n_other_shower = 0;
      int n_pi0 = 0;
      for (auto it = map_vertex_to_shower[vertex].begin(); it != map_vertex_to_shower[vertex].end(); it++){
	WCPPID::WCShower *shower1 = *it;
	if (shower1 == shower) continue;
	if (shower1->get_start_vertex().second > 2) continue;
	if (shower1->get_particle_type()!=11) continue;
	double E_shower = 0;
	if (shower1->get_kine_best() != 0){ 
	  E_shower = shower1->get_kine_best();
	}else{
	  E_shower = shower1->get_kine_charge();
	}
	if (E_shower > 60*units::MeV) {
	  n_other_shower ++;
	}
	if (E_shower > 25*units::MeV && shower1->get_start_vertex().first == vertex && pi0_showers.find(shower1) != pi0_showers.end()) n_pi0 ++;
      }
      
      // 7023_669_33467
      if (n_pi0 <2 || shower_energy > 1000 *units::MeV){ // 7004_1546_77305
	flag_bad = true;
      
	std::cout << "Xin_J_2: " << shower_energy << " " << dis/units::cm << " " << medium_dQ_dx << " " << n_other_shower << " " << max_dQ_dx << " " << n_pi0 << std::endl;
      }
    }
    
    //    if (flag_print) std::cout << "kaka: " << shower_energy << " " << dis/units::cm << " " << flag_bad << std::endl;
    
    if (!flag_bad){
      TVector3 dir2 = shower->cal_dir_3vector(vertex_point, 6*units::cm);
      dir2 = (-1) * dir2.Unit();
      
      for (auto it = showers.begin(); it != showers.end(); it++){
	WCPPID::WCShower *shower1 = *it;
	if (shower1 == shower) continue;
	if (shower1->get_total_length() == 0) continue;
	if (shower1->get_particle_type()!=11){
	  TVector3 dir1(shower1->get_start_point().x - vertex_point.x, shower1->get_start_point().y - vertex_point.y, shower1->get_start_point().z - vertex_point.z);
	  
	  if ((medium_dQ_dx > 1.3  || shower_energy < 300*units::MeV)&& std::min(dir1.Angle(dir)/3.1415926*180.,dir1.Angle(dir2)/3.1415926*180.)<15 && dir1.Mag() < 40*units::cm && max_dQ_dx < 3.0){
	    //	
	    flag_bad = true;
	    std::cout << "Xin_J_1: " << shower_energy << " " << shower1->get_particle_type() << " " << shower1->get_start_segment()->get_particle_type() << " " << dir1.Angle(dir)/3.1415926*180. << " " << dir1.Angle(dir2)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << medium_dQ_dx << " " << shower1->get_total_length()/units::cm << std::endl;
	  }
	}
      }
      
     
      
      for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
	WCPPID::ProtoVertex *vtx1 = it->first;
	if (vtx1->get_cluster_id() == vertex->get_cluster_id() || it->second.size()==1) continue;
	TVector3 dir1(vtx1->get_wcpt().x - vertex_point.x, vtx1->get_wcpt().y - vertex_point.y, vtx1->get_wcpt().z - vertex_point.z);
	if (std::min(dir1.Angle(dir)/3.1415926*180. , dir1.Angle(dir2)/3.1415926*180.) <15)
	  if (dir1.Mag() < 40*units::cm && (shower_energy < 300*units::MeV || medium_dQ_dx > 1.3)){
	    flag_bad = true;
	    std::cout << "Xin_J_0: " << shower_energy << " " << medium_dQ_dx << " " << it->second.size() << " " << dir1.Mag()/units::cm << " " << dir1.Angle(dir)/3.1415926*180. << " " << dir1.Angle(dir2)/3.1415926*180. << std::endl;
	  }
      }
      
    
    }
  }else{
    // not single shower ...
    int num_valid_tracks = 0;
    for (auto it2 = map_vertex_segments[vertex].begin(); it2 != map_vertex_segments[vertex].end(); it2++){
      WCPPID::ProtoSegment *sg1 = *it2;
      if (sg1 == shower->get_start_segment()) continue;
      if ((!sg1->get_flag_shower()) && ((!sg1->is_dir_weak()) || sg1->is_dir_weak() && sg1->get_length() > 5*units::cm)) num_valid_tracks ++;
    }

    if (num_valid_tracks ==0 && dis < 3*units::cm) flag_bad = true;
    //    std::cout << "kaka: " << shower_energy << " " << num_valid_tracks << " " << dis/units::cm << std::endl;
    
  }
  
  
  return flag_bad;
}






std::pair<bool, int> WCPPID::NeutrinoID::gap_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment* sg, bool flag_single_shower, int num_valid_tracks, double E_shower){
  bool flag_gap = false;

  Point vertex_point;
  bool flag_start;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    flag_start = true; // front ...
    vertex_point = sg->get_point_vec().front();
  }else{
    flag_start = false; // back ...
    vertex_point = sg->get_point_vec().back();
  }

  PointVector &pts = sg->get_point_vec();

  int n_points = 0;
  int n_bad = 0;

  bool flag_prolong_u;
  bool flag_prolong_v;
  bool flag_prolong_w;
  bool flag_parallel;
  
  if (flag_start){
    Point closest_p = pts.back();
    double min_dis = 1e9;
    for (int i=0;i<pts.size();i++){
      double dis = fabs(sqrt(pow(pts.at(i).x - vertex_point.x,2) + pow(pts.at(i).y - vertex_point.y,2) + pow(pts.at(i).z - vertex_point.z,2))-3*units::cm);
      if (dis < min_dis){
        min_dis = dis;
	closest_p = pts.at(i);
      }
    }

    TVector3 dir(closest_p.x - vertex_point.x, closest_p.y - vertex_point.y, closest_p.z - vertex_point.z);
    std::vector<bool> flag = main_cluster->check_direction(dir);
    flag_prolong_u = flag.at(0);
    flag_prolong_v = flag.at(1);
    flag_prolong_w = flag.at(2);
    flag_parallel = flag.at(3);

    //    std::cout << flag_prolong_u << " " << flag_prolong_v << " " << flag_prolong_w << " " << flag_parallel << std::endl;
    
    for (int i=0;i+1 < pts.size();i++){
      Point test_p;
      for (int j=0;j!=3;j++){
	test_p.x = pts.at(i).x + j/3.*(pts.at(i+1).x - pts.at(i).x);
	test_p.y = pts.at(i).y + j/3.*(pts.at(i+1).y - pts.at(i).y);
	test_p.z = pts.at(i).z + j/3.*(pts.at(i+1).z - pts.at(i).z);

	//ct_point_cloud->Print(test_p);
	
	int num_bad_ch = 0;
	int num_connect = 0;
	int num_spec = 0;
	// check U
	{
	  WCP::CTPointCloud<double> tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 0);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 0, 0)) num_bad_ch ++;
	    else if (flag_prolong_u) num_spec ++;
	  }
	  //	  std::cout << "U: " << tmp_pts.pts.size() << std::endl;
	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 1);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 1, 0)) num_bad_ch ++;
	    else if (flag_prolong_v) num_spec ++;
	  }
	  //	  std::cout << "V: " << tmp_pts.pts.size() << std::endl;
	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 2);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 2, 0)) num_bad_ch ++;
	  }
	  //	  std::cout << "W: " << tmp_pts.pts.size() << std::endl;
	}
	//	std::cout << num_connect << " " << num_bad_ch << " " << num_spec << std::endl;

	if (num_connect + num_bad_ch + num_spec == 3){
	  if (n_bad == n_points && n_bad <=5) n_bad = 0;
	}else{
	  n_bad ++;
	}
	n_points ++;
	//	if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm,0,0)) n_bad ++;
      }
      double dis = sqrt(pow(pts.at(i+1).x - vertex_point.x,2) + pow(pts.at(i+1).y - vertex_point.y,2) + pow(pts.at(i+1).z - vertex_point.z,2));
      if (dis > 2.4*units::cm) break;
    }
  }else{
    Point closest_p = pts.back();
    double min_dis = 1e9;
    for (int i=0;i<pts.size();i++){
      double dis = fabs(sqrt(pow(pts.at(i).x - vertex_point.x,2) + pow(pts.at(i).y - vertex_point.y,2) + pow(pts.at(i).z - vertex_point.z,2))-3*units::cm);
      if (dis < min_dis){
        min_dis = dis;
	closest_p = pts.at(i);
      }
    }

    TVector3 dir(closest_p.x - vertex_point.x, closest_p.y - vertex_point.y, closest_p.z - vertex_point.z);
    std::vector<bool> flag = main_cluster->check_direction(dir);
    flag_prolong_u = flag.at(0);
    flag_prolong_v = flag.at(1);
    flag_prolong_w = flag.at(2);
    flag_parallel = flag.at(3);
    
    for (int i=int(pts.size())-1;i>0;i--){
      Point test_p;
      for (int j=0;j!=3;j++){
	test_p.x = pts.at(i).x + j/3.*(pts.at(i-1).x - pts.at(i).x);
	test_p.y = pts.at(i).y + j/3.*(pts.at(i-1).y - pts.at(i).y);
	test_p.z = pts.at(i).z + j/3.*(pts.at(i-1).z - pts.at(i).z);

	
	
	int num_bad_ch = 0;
	int num_connect = 0;
	int num_spec = 0;
	
	// check U
	{
	  WCP::CTPointCloud<double> tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 0);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 0, 0)) 	      num_bad_ch ++;
	    else if (flag_prolong_u) num_spec ++;
	  }
	  //	  std::cout << "U: " << tmp_pts.pts.size() << std::endl;
	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 1);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 1, 0)) 	      num_bad_ch ++;
	    else if (flag_prolong_v) num_spec ++;
	  }
	  //	  std::cout << "V: " << tmp_pts.pts.size() << std::endl;
	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 2);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 2, 0)) 	      num_bad_ch ++;
	    else if (flag_prolong_w) num_spec ++;
	  }
	  //	  std::cout << "w: " << tmp_pts.pts.size() << std::endl;
	}
	//	std::cout << num_connect << " " << num_bad_ch << " " << num_spec << std::endl;

	if (num_connect + num_bad_ch + num_spec == 3){
	  if (n_bad == n_points && n_bad <=5) n_bad = 0;
	}else{
	  n_bad ++;
	}
	n_points ++;
	//	if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm,0,0)) n_bad ++;
      }
      double dis = sqrt(pow(pts.at(i-1).x - vertex_point.x,2) + pow(pts.at(i-1).y - vertex_point.y,2) + pow(pts.at(i-1).z - vertex_point.z,2));
      if (dis > 2.4*units::cm) break;
    }
  }


  if (E_shower > 900*units::MeV) { // very high energy ...
    if ((!flag_single_shower) && (!flag_parallel)){
      if (E_shower > 1200*units::MeV){
	if (n_bad > 2./3 * n_points) flag_gap = true;
      }else{
	if (n_bad > 1./3. *n_points) flag_gap = true;
      }
    }
    // 7008_720_36032	 
    if (flag_parallel && (!flag_single_shower)){
      if (n_bad > 1/2. * n_points) flag_gap = true;
    }
  }else if (E_shower > 150*units::MeV){
    if ((!flag_single_shower)){
      if (flag_parallel){
	if (n_bad > 4) flag_gap = true;
      }else{
	if (n_bad > 1) flag_gap = true;
      }
    }else{
      if (n_bad > 2) flag_gap = true;
    }
  }else{
    if (!flag_single_shower){ 
      if (flag_parallel){
	if (n_bad > 3) flag_gap = true;
      }else{ 
	if (n_bad > 1) flag_gap = true;
      } 
    }else{
      if (n_bad > 2) flag_gap = true; 
    }
  }
  // 7021_521_26090 
  if (n_bad >=6 && E_shower < 1000*units::MeV) flag_gap = true;
  if (E_shower <=900*units::MeV && n_bad >1) flag_gap = true;
  
  //  std::cout << "kaka "<< sg->get_id() << " " << E_shower/units::MeV << " " << n_bad << " " << n_points << " " << " " << flag_start << " " << flag_parallel << " " << flag_single_shower << " " << num_valid_tracks << " " << flag_gap << std::endl;
  // hack 
  // if (n_bad >0) flag_gap = true;
    
    
  return std::make_pair(flag_gap, n_bad);
}

int WCPPID::NeutrinoID::mip_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg, WCPPID::WCShower *shower, bool flag_single_shower, bool flag_strong_check, bool flag_print){
  int mip_id = 1; 
  // 1 good, -1 bad, 0 not sure ...
  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }

  Point vertex_point;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
  }else{
    vertex_point = sg->get_point_vec().back();
  }
  
  TVector3 dir_shower;
  if (shower->get_start_segment()->get_length() > 12*units::cm){
    dir_shower = shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);
  }else{
    dir_shower = shower->cal_dir_3vector(vertex_point,15*units::cm);
  }
  if (fabs(dir_shower.Angle(dir_drift)/3.1415926*180.-90)<10 || Eshower > 800*units::MeV) dir_shower = shower->cal_dir_3vector(vertex_point,25*units::cm);
  dir_shower = dir_shower.Unit();
  
  
  double dQ_dx_cut  = 1.45;
  if (Eshower > 1200*units::MeV) dQ_dx_cut = 1.85;
  else if (Eshower > 1000*units::MeV) dQ_dx_cut = 1.6;
  else if (Eshower < 550*units::MeV) dQ_dx_cut = 1.3;
  if (Eshower < 300*units::MeV) dQ_dx_cut = 1.3;
  // hack 
  //  dQ_dx_cut = 1.3;
  //  std::cout << Eshower << " " << dQ_dx_cut << std::endl;
  
  
  std::vector<double> vec_dQ_dx = shower->get_stem_dQ_dx(vertex, sg, 20);
  

  std::vector<int> vec_threshold(vec_dQ_dx.size(), 0);
  for (size_t i=0;i!=vec_dQ_dx.size(); i++){
    if (vec_dQ_dx.at(i)>dQ_dx_cut) vec_threshold.at(i) = 1;
  }

  int n_end_reduction = 0;
  double prev_vec_dQ_dx = vec_dQ_dx.front();
  for (size_t i=1;i<vec_dQ_dx.size();i++){
    if (vec_dQ_dx.at(i) < prev_vec_dQ_dx){
      n_end_reduction = i;
      prev_vec_dQ_dx = vec_dQ_dx.at(i);
      if (vec_dQ_dx.at(i) < dQ_dx_cut) break;
    }
  }
  
  int n_first_mip = 0; // first MIP like ...
  for (size_t i=0; i!= vec_dQ_dx.size();i++){
    n_first_mip = i;
    if (vec_threshold.at(i) ==0 ) break;
  }

  int n_first_non_mip = n_first_mip;
  for (size_t i=n_first_non_mip; i<vec_dQ_dx.size(); i++){
    n_first_non_mip = i;
    if (vec_threshold.at(i)==1) break;
  }

  int n_first_non_mip_1 = n_first_mip;
  for (size_t i=n_first_non_mip; i<vec_dQ_dx.size(); i++){
    n_first_non_mip_1 = i;
    if (vec_threshold.at(i)==1 && i+1 < vec_dQ_dx.size()){
      if (vec_threshold.at(i+1) == 1) break;
    }
  }

  int n_first_non_mip_2 = n_first_mip;
  for (size_t i=n_first_non_mip; i<vec_dQ_dx.size(); i++){
    n_first_non_mip_2 = i;
    if (vec_threshold.at(i)==1 && i+1 < vec_dQ_dx.size()){
      if (vec_threshold.at(i+1) == 1 && i+2 < vec_dQ_dx.size()){
	if (vec_threshold.at(i+2) == 1)    break;
      }
    }
  }
  
  //  for (size_t i=0;i!=vec_dQ_dx.size();i++){
  //  std::cout  << i << " " << vec_dQ_dx.at(i) << std::endl;
  // }

  
  int n_showers = 0;
  int n_protons = 0;
  int n_tracks = 0;
  std::set<WCPPID::WCShower* > connected_showers;
  std::set<WCPPID::WCShower* > tmp_pi0_showers;
  {
    auto it = map_vertex_to_shower.find(vertex);
    if (it != map_vertex_to_shower.end()){
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::WCShower *shower = *it1;
	WCPPID::ProtoSegment *sg1 = shower->get_start_segment();
	if (sg1->get_particle_type()!=11) continue;
	if (map_vertex_segments[main_vertex].find(sg1) != map_vertex_segments[main_vertex].end()){
	  n_showers ++;
	  connected_showers.insert(shower);
	}
	if (pi0_showers.find(shower) != pi0_showers.end())
	  tmp_pi0_showers.insert(shower);
      }
    }
    for (auto it1 = map_vertex_segments[vertex].begin(); it1 != map_vertex_segments[vertex].end(); it1++){
      WCPPID::ProtoSegment *sg1 = *it1;
      if (sg1->get_flag_shower()) continue;
      n_tracks ++;
      if (sg1->get_particle_type()==2212) n_protons ++;
    }
  }

  double lowest_dQ_dx = 100; int n_lowest = 0;
  double highest_dQ_dx = 0; int n_highest = 0;
  int n_below_threshold = 0;
  int n_below_zero = 0;
  for (size_t i=n_first_mip; i < n_first_non_mip_2; i++){
    if (vec_dQ_dx.at(i) < lowest_dQ_dx && i <= 12) {
      lowest_dQ_dx = vec_dQ_dx.at(i);
      n_lowest = i;
    }
    if (vec_dQ_dx.at(i) > highest_dQ_dx) {
      highest_dQ_dx = vec_dQ_dx.at(i);
      n_highest = i;
    }
    if (vec_dQ_dx.at(i) < dQ_dx_cut) n_below_threshold ++;
    if (vec_dQ_dx.at(i) < 0) n_below_zero ++;
  }
  

  

   // quality check ...
  {
    double medium_dQ_dx = 1;
    {
      std::vector<double> tmp_vec_dQ_dx = vec_dQ_dx;
      std::nth_element(tmp_vec_dQ_dx.begin(), tmp_vec_dQ_dx.begin() + tmp_vec_dQ_dx.size()/2, tmp_vec_dQ_dx.end());
      medium_dQ_dx  = *std::next(tmp_vec_dQ_dx.begin(), tmp_vec_dQ_dx.size()/2);
    }
    // low energy only ...
    if (medium_dQ_dx < 0.75 && Eshower < 150*units::MeV) {
      mip_id = -1;
      //      std::cout << "kaka1: " << medium_dQ_dx << " " << Eshower/units::MeV << std::endl;
      return mip_id;
    }
  }
  
  if (Eshower < 800*units::MeV){
    // check overlapping situation (inside shower)
    PointVector test_pts;
    if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
      for (size_t i=0;i!=sg->get_point_vec().size();i++){
	if (i==3) break;
	test_pts.push_back(sg->get_point_vec().at(i));
      }
    }else{
      for (int i=int(sg->get_point_vec().size())-1;i>=0;i--){
	if (i == int(sg->get_point_vec().size())-4) break;
	test_pts.push_back(sg->get_point_vec().at(i));
      }
    }
    WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
    int nconnected = map_vertex_segments[other_vertex].size();
    
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();

    for (size_t i=0; i!= test_pts.size(); i++){
      double min_u = 1e9, min_v = 1e9, min_w = 1e9;
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (sg1 == sg) continue;
	auto tuple_result = sg1->get_closest_2d_dis(test_pts.at(i));
	if (std::get<0>(tuple_result) < min_u) min_u = std::get<0>(tuple_result);
      	if (std::get<1>(tuple_result) < min_v) min_v = std::get<1>(tuple_result);
    	if (std::get<2>(tuple_result) < min_w) min_w = std::get<2>(tuple_result);
	//std::tuple<double, double, double> get_closest_2d_dis(WCP::Point &p);
	//	std::cout << sg1->get_id() << " " << sg->get_id() << " " << min_u << " " << min_v << " " << min_w << std::endl;
      }
      //
      if (min_u < 0.3*units::cm && min_v < 0.3*units::cm && min_w < 0.3*units::cm) {
	// 7017_617_30888 
	if (i==0 && min_u == 0 && min_v ==0 && min_w ==0 || i+1==sg->get_point_vec().size() && min_u ==0 && min_v ==0 && min_w ==0 && nconnected == 2){
	}else{
	  mip_id = -1;
	  // std::cout << "kaka2: " << Eshower << " " << i << " " << min_u/units::cm << " " << min_v/units::cm << " " << min_w/units::cm << " " << sg->get_point_vec().size() << " " << nconnected << std::endl;
	}
	return mip_id;
      }
    }
    
  }

   // shower split ...
  if (n_showers == 2 && n_tracks ==0 && Eshower < 500*units::MeV) {
    mip_id = -1; // bad
    bool flag_inside_pi0 = false;
    double shortest_length = 1e9;
    double shortest_acc_length = 0;
    double shortest_angle = 0;
    TVector3 dir1 = shower->cal_dir_3vector(vertex->get_fit_pt(),6*units::cm);
    for (auto it1 = connected_showers.begin(); it1 != connected_showers.end(); it1++){
      WCPPID::WCShower *shower1 = *it1;
      if (shower1 == shower) continue;
      double medium_dQ_dx = shower1->get_start_segment()->get_medium_dQ_dx()/(43e3/units::cm);
      double length = shower1->get_start_segment()->get_length();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
      if (medium_dQ_dx > dQ_dx_cut) mip_id = 1;
      if (tmp_pi0_showers.find(shower1) != tmp_pi0_showers.end()) flag_inside_pi0 = true;
      
      TVector3 dir2 = shower1->cal_dir_3vector(vertex->get_fit_pt(),6*units::cm);
      if (length < shortest_length) {
	shortest_length = length;
	shortest_angle = dir1.Angle(dir2)/3.1415926*180.;
	shortest_acc_length = shower1->get_total_length(shower1->get_start_segment()->get_cluster_id());
      }
    }
    // 7004_365_18300
    if (mip_id==-1 && (!flag_inside_pi0) && tmp_pi0_showers.size()>0 ) mip_id = 1;
    // 7010_1076_53830
    if (mip_id==-1 && (shortest_angle > 45 && shortest_length > 20*units::cm || shortest_angle > 35 && shortest_acc_length > 40*units::cm)&& shortest_length < 1e9) mip_id = 1;
    
    //    std::cout << "kaka3: " << Eshower << " " << flag_inside_pi0 << " " << tmp_pi0_showers.size() << " " << shortest_length/units::cm << " " << shortest_angle <<" " << short_acc_length/units::cm <<  std::endl;
    return mip_id;
  }

  

  
 
  if (n_first_non_mip_2 - n_first_mip >=2 && // dQ_dx cut ...
      (n_first_mip <=2 || (n_first_mip <= n_end_reduction &&
			   (n_first_mip <=3
			    || n_first_mip <=4 && n_first_non_mip_1 - n_first_mip > 5 && Eshower > 150*units::MeV // 7012_177_8857, 7051_76_3820, fixed 7020_29_1465
			    || n_first_mip <=4 && Eshower > 600*units::MeV
			    || n_first_mip <=5 && Eshower > 800*units::MeV
			    || n_first_mip <=6 && Eshower > 1000*units::MeV
			    || n_first_mip <=10 && Eshower > 1000*units::MeV && n_first_non_mip_1 - n_first_mip > 5
			    || n_first_mip <=10 && Eshower > 1250*units::MeV) )
       //       || (n_first_mip <= n_end_reduction &&
       //   (n_first_non_mip_1 - n_first_mip >= 5 && n_first_mip <= 6))
       )) mip_id = 1;
  else mip_id = -1;

  double max_dQ_dx_sample = 0;
  for (size_t i = n_first_non_mip_2; i < n_first_non_mip_2 + 3;i++){
    if (i >= vec_dQ_dx.size()) break;
    if (vec_dQ_dx.at(i) > max_dQ_dx_sample) max_dQ_dx_sample = vec_dQ_dx.at(i);
  }

  
  
  // 7013_63_3191 + 7004_498_24922 + 7014_786_39346 + 6583_143_7192
  if (mip_id==-1 && n_first_mip <= n_end_reduction && n_first_mip <=5 && ((n_first_non_mip_2 - n_first_mip >= 8 && n_first_non_mip - n_first_mip >=7 || n_first_non_mip_2 - n_first_mip >=5 && max_dQ_dx_sample < 1.6) && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) > 1.75 || n_first_non_mip_2 - n_first_mip >=5 && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) > 3.0)) {
    mip_id = 0;
  }
  // 6640_171_8560 + 7014_954_47722
  if (mip_id == -1 && n_first_mip <= n_end_reduction && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) < 1.45 && (n_first_non_mip_2 - n_first_mip + n_end_reduction >=12 && n_below_threshold + n_end_reduction >=10) && n_first_non_mip_2 - n_first_mip >= 4 && (n_end_reduction <4 && Eshower < 100*units::MeV || n_end_reduction < 7 && Eshower < 200*units::MeV && Eshower >= 100*units::MeV || Eshower >= 200*units::MeV)) {
    if (flag_single_shower)
      mip_id = 0;
    else
      mip_id = 1;
  }
  
  if (flag_strong_check){
    if (!((n_first_mip <=2 || (n_first_mip <= n_end_reduction &&
			       (n_first_mip <=3
				|| n_first_mip <=4 && Eshower > 600*units::MeV
				|| n_first_mip <=5 && Eshower > 800*units::MeV
				|| n_first_mip <=6 && Eshower > 1000*units::MeV
				|| n_first_mip <=10 && Eshower > 1250*units::MeV)))
	  && (n_first_non_mip_2 - n_first_mip > 3 || n_first_non_mip_2 - n_first_mip == 3 && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) > 3.3 ) // 7012_366_18344
	  )) mip_id = -1;
  
    if (mip_id == -1 && n_first_mip <= n_end_reduction && (n_first_mip <=5 && n_first_non_mip_2 - n_first_mip >=7 )&& std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) > 3.3 ) mip_id = 0;
    
    // std::cout << "Xin_A:" << n_first_mip << " " << n_first_non_mip_2 - n_first_mip << std::endl;
  }
  

  if (flag_print && mip_id==-1) std::cout << "Qian_B_4" << std::endl;
  
  if (Eshower < 600*units::MeV){
    int n_good_tracks = 0;
    for (auto it = map_vertex_segments[vertex].begin(); it!= map_vertex_segments[vertex].end(); it++){
      WCPPID::ProtoSegment *sg1 = *it;
      if (sg1 == sg) continue;
      if (!sg1->is_dir_weak() || sg1->get_length() > 10*units::cm) n_good_tracks ++;
    }
    // 6043_4_243 event
    if (n_good_tracks >1 && n_first_non_mip_2 <=2) mip_id = -1;
    //std::cout << n_good_tracks << " " << n_first_non_mip_2 << " " << Eshower << std::endl;
  }

  if (flag_print && mip_id==-1) std::cout << "Qian_B_3" << std::endl;
  
  // single shower situation with energy lower than 500 MeV
  if (mip_id == 1 &&  map_vertex_segments[vertex].size() ==1 && Eshower < 500*units::MeV){
    //    std::cout << vec_dQ_dx.front() << " " << medium_dQ_dx << std::endl;
    if (Eshower < 180*units::MeV || n_first_mip>0 || vec_dQ_dx.front() > 1.15 && n_end_reduction >= n_first_mip && Eshower < 360*units::MeV)
      mip_id = 0;
    // 6058_43_2166, 7003_1636_81828, 7054_364_18210
    if (flag_single_shower && Eshower < 400*units::MeV && n_end_reduction > 0) mip_id = 0;
  }else if (mip_id==1 && map_vertex_segments[vertex].size() > 1 && Eshower < 300*units::MeV){
    if (vec_dQ_dx.size() >=3){ // 7017_482_24127, a dip?
      if (vec_dQ_dx.at(1) < 0.6 || vec_dQ_dx.at(2) < 0.6) mip_id = 0;
    }
    //    std::cout << vec_dQ_dx.at(0) << " " << vec_dQ_dx.at(1) << " " << vec_dQ_dx.at(2) << std::endl;
  }else if (mip_id==1 && map_vertex_segments[vertex].size() > 1 && Eshower < 600*units::MeV){
    TVector3 dir = shower->cal_dir_3vector(vertex_point, 15*units::cm);
    // 7017_579_28979
    if (dir.Angle(dir_beam)/3.1415926*180. > 60 || n_first_non_mip_1 == 1) mip_id = 0;
    
    // 7018_876_43824
    bool flag_all_above = true;
    for (size_t i=0;i!=vec_dQ_dx.size();i++){
      if (vec_dQ_dx.at(i) < 1.2) flag_all_above = false;
      if (i > 5) break;
    }
    if (flag_all_above) mip_id = 0;
  }else if (mip_id == 1 && flag_single_shower && Eshower < 900*units::MeV){
    // 7025_615_30788
    if (flag_single_shower && n_first_mip !=0 ) mip_id = 0;
  }

  if (flag_print && mip_id==-1) std::cout << "Qian_B_2" << std::endl;


  
  TVector3 dir = shower->cal_dir_3vector(vertex_point, 15*units::cm);
  if (Eshower < 300*units::MeV){   
    if (dir.Angle(dir_beam)/3.1415926*180. > 40){
      // 7018_926_46331 // STEM length too short ... 7001_787_39370 ...
      if ((n_first_non_mip_2 - n_first_mip <=3 && n_first_non_mip_2 <=3 || n_first_non_mip_2 - n_first_mip <=2)&& n_first_mip <=1 && max_dQ_dx_sample > 1.9) mip_id = -1;
      // 5337_192_9614 
      if (flag_single_shower && n_first_mip >=3 && n_first_non_mip - n_first_mip <=1  && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) < 2.7) mip_id = -1;
      // 7048_108_5419
      if (flag_single_shower && n_first_mip >=2 && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) < 2.7) mip_id = -1;
    }
    // 7021_586_29303 
    if (dir.Angle(dir_beam)/3.1415926*180. > 30 && Eshower < 200*units::MeV && flag_single_shower){
      if (vec_dQ_dx.at(0)>1.5 && n_first_mip>0 && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) < 2.7) mip_id = -1;
      //std::cout << vec_dQ_dx.at(0) << std::endl;
    }
  }

  if (flag_print && mip_id==-1) std::cout << "Qian_B_1_0" << std::endl;
  
  // 7017_1631_81564
  if (flag_single_shower && Eshower < 500*units::MeV && shower->get_total_length(vertex->get_cluster_id()) > shower->get_total_length() *0.95){
    double min_dQ_dx = 1e9;
    for (size_t i=0;i!=vec_dQ_dx.size();i++){
      if (vec_dQ_dx.at(i) < min_dQ_dx) min_dQ_dx = vec_dQ_dx.at(i);
      if (i>5) break;
    }
    if (n_first_non_mip_2 - n_first_mip<=2 && min_dQ_dx > 1.3){
      mip_id = -1;
      //      std::cout << min_dQ_dx << std::endl;
    }
    //  std::cout << shower->get_total_length(vertex->get_cluster_id()) << " " << shower->get_total_length() << std::endl;
  }


  if (flag_print && mip_id==-1) std::cout << "Qian_B_1" << std::endl;
  
  if (mip_id == 1){
    if (n_below_threshold <=5 && (lowest_dQ_dx < 0.7 || lowest_dQ_dx > 1.1 && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. < 15)) mip_id = 0;
  }
  // 7018_235_11772
  if (lowest_dQ_dx > 1.3 && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. < 15 && Eshower < 1000*units::MeV) mip_id = -1;
  // 7049_1241_62062 + 7017_997_49856
  if (lowest_dQ_dx < 0 && Eshower < 800*units::MeV && n_below_zero > 2) mip_id = -1;
  // 7025_380_19030
  if (lowest_dQ_dx < 0 && Eshower < 800*units::MeV && n_below_zero <=2 && highest_dQ_dx > 1.3) mip_id = -1;
  // 7054_1985_99267
  if (lowest_dQ_dx < 0 && n_lowest <=1 && n_highest < n_lowest && highest_dQ_dx < 0.9) mip_id = -1;
  // 7049_1033_51667
  if (lowest_dQ_dx <0.6 && highest_dQ_dx < 0.8 && n_lowest <=1 && n_highest <=1 && n_first_non_mip_2 - n_first_mip <=2 && max_dQ_dx_sample > 1.8) mip_id = -1;
  
  // 7017_1508_75440, 
  if (lowest_dQ_dx < 0.6 && highest_dQ_dx > 1.3 && n_highest >1 && n_highest < 4 && Eshower < 1000*units::MeV &&fabs(n_lowest - n_highest)>1 ) mip_id = -1;
  //7003_1754_87734 + anti 7026_54_2747
  if (lowest_dQ_dx < 0.9 && n_lowest <=1 && highest_dQ_dx > 1.2 && n_below_threshold <= 4 && Eshower < 1000*units::MeV && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. > 10 && n_first_non_mip_2 <5 && max_dQ_dx_sample > 1.9) mip_id = -1;

  
  // 7012_1370_68520 + anti 7010_451_22560
  if (n_lowest <=2 && n_highest > n_lowest && lowest_dQ_dx > 1.1 && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. <5 && map_vertex_segments[vertex].size() > 1) mip_id = -1;
  // 7055_147_7354
  if (n_lowest <=3 && lowest_dQ_dx < 0.7 && highest_dQ_dx > 1.3 && n_highest< n_lowest && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. <5 ) mip_id = -1;
  // 7012_297_14884
  if (flag_single_shower && n_below_threshold <=3 && highest_dQ_dx > 1.2 && Eshower < 800*units::MeV && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. > 7.5 ) {
    mip_id = -1;
    // 7012_366_18344
    if (n_below_threshold == 3 && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) > 3.5) mip_id = 0;
  }
  
  // 6936_165_8288
  if (Eshower < 800*units::MeV && lowest_dQ_dx < 0.2 && n_lowest <=3 && fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. > 15 && sg->get_length() < 5*units::cm) mip_id = -1;

  if (flag_print && mip_id==-1) std::cout << "Qian_B_0" << std::endl;
  
  double E_direct_max_energy = 0, E_direct_total_energy = 0;
  double E_indirect_max_energy = 0, E_indirect_total_energy = 0;
  int n_direct_showers = 0;
  int n_indirect_showers = 0;
  
  for (auto it1 = showers.begin(); it1 != showers.end() ;it1++){
    WCPPID::WCShower *shower1 = *it1;
    WCPPID::ProtoSegment *sg = shower1->get_start_segment();
    if (sg->get_particle_type()!=11) continue;
    if (shower1 == shower) continue; 
    auto pair_result = shower1->get_start_vertex();
    double E_shower1 = 0;
    if (shower1->get_kine_best() != 0){ 
      E_shower1 = shower1->get_kine_best();
    }else{
      E_shower1 = shower1->get_kine_charge();
    }

    if (pair_result.second == 1){
      E_direct_total_energy += E_shower1;
      if (E_shower1 > E_direct_max_energy) E_direct_max_energy = E_shower1;
      if (E_shower1 > 80*units::MeV) n_direct_showers ++;
    }else if (pair_result.second == 2){
      E_indirect_total_energy += E_shower1;
      if (E_shower1 > E_indirect_max_energy) E_indirect_max_energy = E_shower1;
      if (E_shower1 > 80*units::MeV) n_indirect_showers ++;
    }
  }
  
  //std::cout << "qaqa: " << Eshower << " " << n_direct_showers << " " << E_direct_max_energy << " " << E_direct_total_energy << " " << n_indirect_showers << " " << E_indirect_max_energy << " " << E_indirect_total_energy << std::endl;


  
  // 7049_1070_53534
  if (flag_single_shower && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) > 1.6 && std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) < 3.5 && Eshower < 350*units::MeV && E_indirect_max_energy > 70*units::MeV) mip_id = -1;
  // 7012_1450_72525	
  if (flag_single_shower && E_indirect_max_energy >0.33 * Eshower && mip_id==1) mip_id = 0;

  // 7023_28_1419
  if (mip_id==0 && Eshower < 250*units::MeV && sg->get_flag_shower_trajectory() && sg->get_length() < 5*units::cm) mip_id = -1;

  if (mip_id == 1){
    
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();

    double length1 = shower->get_total_length(sg->get_cluster_id());
    double length2 = shower->get_total_length();
    double min_dis = 1e9;
    for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      if (sg1->get_cluster_id() == sg->get_cluster_id() || sg1->get_length() < 3*units::cm) continue;

      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::ProtoVertex *vtx1 = *it1;
	for (auto it2 = map_vtx_segs.begin(); it2 != map_vtx_segs.end(); it2++){
	  WCPPID::ProtoVertex *vtx2 = it2->first;
	  if (vtx2->get_cluster_id() != sg->get_cluster_id()) continue;
	  double dis = sqrt(pow(vtx1->get_fit_pt().x - vtx2->get_fit_pt().x,2) + pow(vtx1->get_fit_pt().y - vtx2->get_fit_pt().y,2) + pow(vtx1->get_fit_pt().z - vtx2->get_fit_pt().z,2));
	  if (dis < min_dis) min_dis = dis;
	}
      }
    }
    // 7012_1646_82342
    if (length1 < 0.1 * length2 && length1 < 10*units::cm && min_dis > 8*units::cm) mip_id = 0;
    //    std::cout << length1 << " " << length2 << " " << min_dis << std::endl;
  }

  {
    WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
    // 7017_969_48489
    if (map_vertex_segments.find(other_vertex) != map_vertex_segments.end())
      if (map_vertex_segments[other_vertex].size()>2 && sg->get_point_vec().size() <= n_first_mip + 1 && n_first_mip > 2) mip_id = -1;
    //    std::cout << sg->get_point_vec().size() << " " << n_first_mip << " " << map_vertex_segments[other_vertex].size() << std::endl;
  }
  
  if (flag_print) std::cout << "Qian_B0: " << n_lowest << " " << lowest_dQ_dx << " " << n_highest << " " << highest_dQ_dx << " " << n_below_threshold << " " << n_below_zero << " " << fabs(3.1415926/2. - dir_shower.Angle(dir_drift))/3.1415926*180. << " " << mip_id << " " << Eshower << " " << flag_single_shower << " " << sg->get_length()/units::cm << " " << max_dQ_dx_sample << std::endl;
  if (flag_print) std::cout << "Qian_B1: "<< mip_id<< " " << sg->get_id() << " " << n_end_reduction << " " << n_first_mip << " " << n_first_non_mip << " " << n_first_non_mip_1 << " " << n_first_non_mip_2 << " "  << Eshower/units::MeV << " " << n_showers << " " << n_tracks << " " << n_protons << " " << map_vertex_segments[vertex].size() << " " << flag_single_shower << " " << dir.Angle(dir_beam)/3.1415926*180. << " " << std::max(vec_dQ_dx.at(0), vec_dQ_dx.at(1)) << std::endl;
  
  
  
  return mip_id;
}
bool WCPPID::NeutrinoID::high_energy_overlapping(WCPPID::WCShower* shower, bool flag_print){
  bool flag_overlap = false;
  double Eshower = 0 ;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }

  std::vector<double> vec_dQ_dx = shower->get_stem_dQ_dx(shower->get_start_vertex().first, shower->get_start_segment(), 20);
  double max_dQ_dx = 0;
  for (size_t i= 0; i!=vec_dQ_dx.size();i++){
    if (vec_dQ_dx.at(i) > max_dQ_dx) max_dQ_dx = vec_dQ_dx.at(i);
    if (i==2) break;
  }
  
  if (Eshower < 1500*units::MeV){
    WCPPID::ProtoSegment *sg = shower->get_start_segment();
    auto pair_result = shower->get_start_vertex();
    WCPPID::ProtoVertex *vtx = pair_result.first;
    Point vtx_point;
    bool flag_start;
    if (vtx->get_wcpt().index == sg->get_wcpt_vec().front().index){
      vtx_point = sg->get_point_vec().front();
      flag_start = true;
    }else{
      vtx_point = sg->get_point_vec().back();
      flag_start = false;
    }

    bool flag_all_showers = true;
    
    // 7012_1195_59764 + 7017_1158_57929 
    if (pair_result.second == 1){
      TVector3 dir1 = sg->cal_dir_3vector(vtx_point, 15*units::cm);
      int n_valid_tracks = 0;
      double min_angle = 180;
      double min_length = 0;
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	if ((*it)==sg) continue;
	if ((*it)->get_particle_type()==11 || (*it)->get_particle_type()==13 && (*it)->is_dir_weak() && (*it)->get_length() < 6*units::cm){
	  TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
	  if (dir2.Mag() == 0) continue;
	  double angle = dir1.Angle(dir2)/3.1415926*180.;
	  if (angle < min_angle) {
	    min_angle = angle;
	    min_length = (*it)->get_length();
	  }
	}else{
	  flag_all_showers = false;
	}

	double medium_dQ_dx = (*it)->get_medium_dQ_dx()/(43e3/units::cm);

	//	std::cout << (*it)->get_particle_type() << " " << (*it)->is_dir_weak() << " " << (*it)->get_length()/units::cm << " " << medium_dQ_dx << std::endl;
	
	
	if ((!(*it)->is_dir_weak() || (*it)->get_particle_type() == 2212 || ((*it)->get_length() > 20*units::cm))  && (!(*it)->get_flag_shower())) n_valid_tracks ++;
	else if (medium_dQ_dx > 2.0 && (*it)->get_length() > 1.8*units::cm) n_valid_tracks ++; // 7010_20_1012
      }
      
      for (auto it1 = map_vertex_to_shower[vtx].begin(); it1 != map_vertex_to_shower[vtx].end(); it1++){
	WCPPID::WCShower *shower1 = *it1;
	if (shower1->get_start_segment()->get_particle_type()!=11) continue;
	if (shower1 == shower) continue;
	
	auto pair_result1 = shower1->get_start_vertex();
	double Eshower1 = 0;
	if (shower1->get_kine_best() != 0){ 
	  Eshower1 = shower1->get_kine_best();
	}else{
	  Eshower1 = shower1->get_kine_charge();
	}
	//	std::cout << Eshower << std::endl;
	if (pair_result1.second == 1 && Eshower1 > 250*units::MeV) n_valid_tracks ++;
      }
      // 6689_127_6366
      if (max_dQ_dx > 3.6 && n_valid_tracks ==0) n_valid_tracks ++;
      
      if (n_valid_tracks ==0 && min_angle < 30) flag_overlap = true;
      if (n_valid_tracks ==0 && min_angle < 60 && flag_all_showers == 1 && Eshower < 300*units::MeV) flag_overlap = true;
      if (n_valid_tracks ==0 && min_angle < 60 && flag_all_showers == 1 && Eshower < 800*units::MeV && min_length < 5*units::cm) flag_overlap = true;
      
      //      std::cout << "kaka: " << Eshower << " " << n_valid_tracks << " " << min_angle << " " << flag_all_showers << " " << min_length/units::cm << " " << flag_overlap << " " << max_dQ_dx << std::endl;


      
    }

    if (pair_result.second == 1){
      TVector3 dir1 = sg->cal_dir_3vector(vtx_point, 8*units::cm);
      double min_angle = 180;
      WCPPID::ProtoSegment *min_sg = 0;
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	if ((*it)==sg) continue;
	TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
	if (dir2.Mag() == 0) continue;
	double angle = dir1.Angle(dir2)/3.1415926*180.;
	if (angle < min_angle) {
	  min_angle = angle;
	  min_sg = (*it);
	}
      }
      
      int ncount = 0;
      PointVector& pts = sg->get_point_vec();
      double medium_dQ_dx = 0;
      if (flag_start){
	for(auto it1 = pts.begin(); it1 != pts.end(); it1++){
	  double min_dis = 1e9;
	  double min_dQ_dx = 0;
	  for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	    if ((*it)==sg) continue;
	    double dis = (*it)->get_closest_point(*it1).first;
	    if (dis < min_dis) {
	      min_dis = dis;
	    }
	  }
	  if (min_dis < 0.6*units::cm) {
	    ncount ++;
	  }else{
	    break;
	  }
	}
	//	medium_dQ_dx = sg->get_medium_dQ_dx(0, ncount)/(43e3/units::cm);
	if (min_sg != 0)
	  if (sg->get_wcpt_vec().front().index == min_sg->get_wcpt_vec().front().index){
	    medium_dQ_dx += min_sg->get_medium_dQ_dx(0, ncount)/(43e3/units::cm);
	  }else{
	    medium_dQ_dx += min_sg->get_medium_dQ_dx(int(min_sg->get_point_vec().size())-1-ncount, int(min_sg->get_point_vec().size())-1)/(43e3/units::cm);
	  }
      }else{
	for (auto it1 = pts.rbegin(); it1 != pts.rend(); it1++){
	  double min_dis = 1e9;
	  for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	    if ((*it)==sg) continue;
	     double dis = (*it)->get_closest_point(*it1).first;
	     if (dis < min_dis) {
	       min_dis = dis;
	     }
	  }
	  if (min_dis < 0.6*units::cm) {
	    ncount ++;
	  }else{
	    break;
	  }
	}
	//	medium_dQ_dx = sg->get_medium_dQ_dx(int(sg->get_point_vec().size())-1-ncount, int(sg->get_point_vec().size())-1)/(43e3/units::cm);
	if (min_sg !=0)
	  if (sg->get_wcpt_vec().back().index == min_sg->get_wcpt_vec().back().index){
	    medium_dQ_dx += min_sg->get_medium_dQ_dx(int(min_sg->get_point_vec().size())-1-ncount, int(min_sg->get_point_vec().size())-1)/(43e3/units::cm);
	  }else{
	    medium_dQ_dx += min_sg->get_medium_dQ_dx(0, ncount)/(43e3/units::cm);
	  }
      }

     
    
      
      
      if (min_angle < 15 && medium_dQ_dx > 0.95 && ncount > 5) flag_overlap = true;
      if (min_angle < 7.5 && medium_dQ_dx > 0.8 && ncount > 8) flag_overlap = true;
      if (min_angle < 5 && ncount > 12 && medium_dQ_dx > 0.5) flag_overlap = true; 
      
      //      std::cout << "qaqa: " << Eshower << " " << min_angle << " " << ncount << " "  << medium_dQ_dx << " " << flag_overlap << std::endl;
    }
    
  }

  

  return flag_overlap;
}



bool WCPPID::NeutrinoID::low_energy_overlapping(WCPPID::WCShower* shower, bool flag_print){
  bool flag_overlap = false;
  TVector3 dir_beam(0,0,1);
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  auto pair_result = shower->get_start_vertex();
  WCPPID::ProtoVertex *vtx = pair_result.first;
  Point vtx_point;
  if (vtx->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vtx_point = sg->get_point_vec().front();
  }else{
    vtx_point = sg->get_point_vec().back();
  }

  
  Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
  Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
      
  // first case ...  all showers
  if (pair_result.second == 1){
    if (Eshower < 150*units::MeV && map_vertex_segments[vtx].size()==1){
      int nseg = 0;
      for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
	if (it->first->get_cluster_id() == sg->get_cluster_id()) nseg ++;
      }
      //std::cout << nseg << std::endl;
      if (nseg ==2){
	for (auto it = map_vtx_segs.begin(); it!= map_vtx_segs.end(); it++){
	  if (it->first->get_cluster_id() != sg->get_cluster_id()) continue;
	  if (it->second.size()==2){
	    TVector3 dir1 = (*it->second.begin())->cal_dir_3vector(it->first->get_fit_pt(), 5*units::cm);
	    TVector3 dir2 = (*it->second.rbegin())->cal_dir_3vector(it->first->get_fit_pt(), 5*units::cm);
	    //std::cout << "Xin1: " << dir1.Mag() << " " << dir2.Mag() << std::endl;
	    //std::cout << dir1.Angle(dir2)/3.1415926*180. << std::endl;
	    if (dir1.Angle(dir2)/3.1415926*180. < 36){
	      flag_overlap = true;
	      //	      std::cout << "A: " << std::endl;
	    }
	  } // two segment 
	} //loop all vertex
      }
    }else if (map_vertex_segments[vtx].size()>1 && Eshower < 300*units::MeV && shower->get_total_length(sg->get_cluster_id()) < 20*units::cm){
      TVector3 dir1 = sg->cal_dir_3vector(vtx_point, 5*units::cm);
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	if ((*it)==sg) continue;
	TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
	
	if (((*it)->get_length() < 30*units::cm && dir1.Angle(dir2)/3.1415926*180.< 10 || // 7017_1604_80242
	     (*it)->get_length() < 7.5*units::cm && dir1.Angle(dir2)/3.1415926*180.< 17.5 )&& (*it)->get_particle_type() == 13){
	  //	  std::cout << shower->get_total_length(sg->get_cluster_id())/units::cm << std::endl;
	  flag_overlap = true;
	  //	  std::cout << "B: " << std::endl;
	}
	//	std::cout << dir1.Angle(dir2)/3.1415926*180. << " " << (*it)->get_length()/units::cm << " " << (*it)->get_particle_type()<< std::endl;
      }
      if ( dir1.Angle(dir_beam)/3.1415926*180.>60){ // 7010_194_9750 
	int n_valid_tracks = 0;
	double min_angle = 180;
	for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	  if ((*it)==sg) continue;
	  TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
	  double tmp_angle = dir2.Angle(dir1)/3.1415926*180.;
	  if (tmp_angle < min_angle) min_angle = tmp_angle;
	  if ((!(*it)->is_dir_weak() || (*it)->get_particle_type() == 2212 || ((*it)->get_length() > 20*units::cm)) && (!(*it)->get_flag_shower())) n_valid_tracks ++;
	}
	if (n_valid_tracks==0 && min_angle < 80) flag_overlap = true;
	//std::cout << dir1.Angle(dir_beam)/3.1415926*180. << std::endl;
      }
      //std::cout << Eshower/units::MeV << std::endl;
    }else if (map_vertex_segments[vtx].size()==2 && Eshower < 400*units::MeV ){ // 7020_249_12479 raise from 100 to 400 MeV
      TVector3 dir1 = sg->cal_dir_3vector(vtx_point, 5*units::cm);
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	if ((*it)==sg) continue;
	TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
	if ((*it)->is_dir_weak() && (*it)->get_length() < 8*units::cm && dir1.Angle(dir2)/3.1415926*180. < 30){
	  //	  std::cout << "XinQ_A: " << (*it)->get_length() << " " << dir1.Angle(dir2)/3.1415926*180. << std::endl;
	  flag_overlap = true;
	}
	
      }
    }

    //    std::cout << "kaka " << flag_overlap << " " << map_vertex_segments[vtx].size() << " " << shower->get_total_length(sg->get_cluster_id())/units::cm << " " <<  Eshower << std::endl;
    
    if (map_vertex_segments[vtx].size()==1 && shower->get_total_length(sg->get_cluster_id()) < 15*units::cm &&  Eshower > 30*units::MeV && Eshower < 250*units::MeV){
      TVector3 dir1 = sg->cal_dir_3vector(vtx_point, 15*units::cm);
      int n_sum = 0;
      int n_out = 0;
      for (auto it = map_vtx_segs.begin(); it!= map_vtx_segs.end(); it++){
	TVector3 dir2(it->first->get_fit_pt().x - vtx_point.x, it->first->get_fit_pt().y - vtx_point.y, it->first->get_fit_pt().z - vtx_point.z);
	n_sum ++;
	if (dir1.Angle(dir2)/3.1415926*180.>15) n_out ++;
      }
      for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
	PointVector& pts = it->first->get_point_vec();
	for(size_t i=1; i+1 < pts.size();i++){
	  TVector3 dir2(pts.at(i).x - vtx_point.x, pts.at(i).y - vtx_point.y, pts.at(i).z - vtx_point.z);
	  n_sum ++;
	  if (dir1.Angle(dir2)/3.1415926*180.>15) n_out ++;
	}
      }
      //      std::cout << n_out << " " << n_sum << std::endl;
      if (n_out > n_sum/3) {
	flag_overlap = true;
	//s	std::cout << "C: " << std::endl;
      }
      //
    }
    //std::cout << shower->get_total_length(sg->get_cluster_id())/units::cm << " " << Eshower/units::MeV << std::endl;
  } // connection type ...

  return flag_overlap;
}


bool WCPPID::NeutrinoID::pi0_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg, WCPPID::WCShower *shower, double threshold){
  bool flag_pi0 = false;
  
  Point vertex_point = vertex->get_fit_pt();

  std::set<WCPPID::ProtoVertex* > used_vertices;
  for (auto it1 = map_shower_pio_id.begin(); it1 != map_shower_pio_id.end(); it1++){
    WCPPID::WCShower *shower1 = it1->first;
    Map_Proto_Vertex_Segments& map_vtx_segs = shower1->get_map_vtx_segs();
    for (auto it2 = map_vtx_segs.begin(); it2 != map_vtx_segs.end(); it2++){
      WCPPID::ProtoVertex *vtx1 = it2->first;
      used_vertices.insert(vtx1);
    }
  }

  
  auto it = map_shower_pio_id.find(shower);
  if (it != map_shower_pio_id.end()){
    std::vector<WCShower*> tmp_pi0_showers = map_pio_id_showers[it->second];
    auto mass_pair = map_pio_id_mass[it->second];

    if (fabs(mass_pair.first - 135*units::MeV)<35*units::MeV && mass_pair.second == 1 || fabs(mass_pair.first - 135*units::MeV) < 60 * units::MeV && mass_pair.second == 2){
      double Eshower_1 = tmp_pi0_showers.front()->get_kine_charge();
      double Eshower_2 = tmp_pi0_showers.back()->get_kine_charge();
      //      std::cout << tmp_pi0_showers.size() << " " << mass_pair.first << " " << mass_pair.second << " " << Eshower_1/units::MeV << " " << Eshower_2/units::MeV << std::endl;
      double dis1 = sqrt(pow(tmp_pi0_showers.front()->get_start_point().x - vertex_point.x,2) + pow(tmp_pi0_showers.front()->get_start_point().y - vertex_point.y,2) + pow(tmp_pi0_showers.front()->get_start_point().z - vertex_point.z,2));
      double dis2 = sqrt(pow(tmp_pi0_showers.back()->get_start_point().x - vertex_point.x,2) + pow(tmp_pi0_showers.back()->get_start_point().y - vertex_point.y,2) + pow(tmp_pi0_showers.back()->get_start_point().z - vertex_point.z,2));

      //std::cout << Eshower_1 << " " << Eshower_2 << " " << dis1/units::cm << " " << dis2/units::cm << " " << vertex_point << std::endl;
      
      if (std::min(Eshower_1, Eshower_2) > 15*units::MeV && fabs(Eshower_1 - Eshower_2)/(Eshower_1 + Eshower_2) < 0.87)	flag_pi0 = true;
      // 6058_43_2166, 7017_364_18210
      if (std::min(Eshower_1, Eshower_2) > std::max(10*units::MeV, threshold) && std::max(Eshower_1, Eshower_2) < 400*units::MeV) flag_pi0 = true;

      // 7049_875_43775
      if (flag_pi0 && (std::min(Eshower_1, Eshower_2) < 30*units::MeV && std::max(dis1, dis2) > 80*units::cm && fabs(Eshower_1 - Eshower_2)/(Eshower_1 + Eshower_2) > 0.87 || 
		       std::min(Eshower_1, Eshower_2) < 30*units::MeV && std::max(dis1, dis2) > 120*units::cm && fabs(Eshower_1 - Eshower_2)/(Eshower_1 + Eshower_2) > 0.80 ) && std::min(dis1, dis2) == 0) flag_pi0 = false;

    }
  }else{
    TVector3 dir1 = sg->cal_dir_3vector(vertex->get_fit_pt(), 12*units::cm);
    if (dir1.Mag() >0){
      for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
	WCPPID::ProtoVertex *vtx1 = it->first;
	if (vtx1->get_cluster_id() == vertex->get_cluster_id()) continue;
	if (used_vertices.find(vtx1) != used_vertices.end()) continue;

	double acc_length = 0;
	for (auto it1 = map_segment_vertices.begin(); it1 != map_segment_vertices.end(); it1++){
	  WCPPID::ProtoSegment *sg1 = it1->first;
	  if (sg1->get_cluster_id() != vtx1->get_cluster_id()) continue;
	  acc_length += sg1->get_length();
	}
	
	TVector3 dir2(vtx1->get_fit_pt().x - vertex->get_fit_pt().x,vtx1->get_fit_pt().y - vertex->get_fit_pt().y, vtx1->get_fit_pt().z - vertex->get_fit_pt().z);
	if (dir2.Mag()>0){
	  if (dir2.Mag() < 36*units::cm && 180 - dir1.Angle(dir2)/3.1415926*180. < 7.5 && acc_length > 0) {
	    
	    flag_pi0 = true;
	    //	    std::cout << dir1.Angle(dir2)/3.1415926*180. << " " << dir2.Mag()/units::cm << " " << acc_length << std::endl;
	  }
	}
      }
    }
  }
  
  
  return flag_pi0;
}

bool WCPPID::NeutrinoID::single_shower_pio_tagger(WCPPID::WCShower *shower, bool flag_print){
  bool flag_bad = false;
  TVector3 dir_beam(0,0,1);
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  Point vertex_point;
  if (shower->get_start_segment()->get_wcpt_vec().front().index == shower->get_start_vertex().first->get_wcpt().index){
    vertex_point = shower->get_start_segment()->get_point_vec().front();
  }else{
    vertex_point = shower->get_start_segment()->get_point_vec().back();
  }
  TVector3 dir = shower->cal_dir_3vector(vertex_point, 15*units::cm);
  double shower_angle = dir.Angle(dir_beam)/3.1415926*180.;

  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  auto pair_result = shower->get_start_vertex();
  WCPPID::ProtoVertex *vtx = pair_result.first;

  if (Eshower < 250*units::MeV){
    for (auto it1 = map_vertex_to_shower[vtx].begin(); it1 != map_vertex_to_shower[vtx].end(); it1++){
      WCPPID::WCShower *shower1 = *it1;
      if (shower1->get_start_segment()->get_particle_type()!=11) continue;
      if (shower1 == shower) continue;
      if (shower1->get_start_vertex().second >2) continue;
      double Eshower1 = 0;
      if (shower1->get_kine_best() != 0){ 
	Eshower1 = shower1->get_kine_best();
      }else{
	Eshower1 = shower1->get_kine_charge();
      }
      if (Eshower1 > 60*units::MeV && shower1->get_start_vertex().second==2){
	TVector3 dir1(shower1->get_start_point().x - vtx->get_fit_pt().x, shower1->get_start_point().y - vtx->get_fit_pt().y, shower1->get_start_point().z - vtx->get_fit_pt().z);
	TVector3 dir2 = shower1->cal_dir_3vector(shower1->get_start_point(), 15*units::cm);
	//	std::cout << Eshower << " " << Eshower1 << " " << dir1.Angle(dir2)/3.1415926*180. << std::endl;
	if (dir1.Angle(dir2)/3.1415926*180. < 30) flag_bad = true;
      }
    }
  }

  if (flag_bad) std::cout << "Xin_I1: " << flag_bad << std::endl;
  
  if ((!flag_bad)){
    // find the vertex inside main cluster which is furthest away from current vertex
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();

    double max_dis = 0;
    WCPPID::ProtoVertex *max_vtx = 0;
    for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
      WCPPID::ProtoVertex *tmp_vtx = it->first;
      if (tmp_vtx->get_cluster_id() != vtx->get_cluster_id()) continue;
      TVector3 dir1(tmp_vtx->get_fit_pt().x - vertex_point.x, tmp_vtx->get_fit_pt().y - vertex_point.y, tmp_vtx->get_fit_pt().z - vertex_point.z);
      double dis = dir1.Dot(dir);
      if (dis > max_dis){
	max_dis = dis;
	max_vtx = tmp_vtx;
      }
    }
    double max_angle = 0;
    WCPPID::ProtoSegment *max_sg = 0;
    for (auto it = map_vtx_segs[max_vtx].begin(); it != map_vtx_segs[max_vtx].end(); it++){
      WCPPID::ProtoSegment *sg1 = *it;
      TVector3 dir1 = sg1->cal_dir_3vector(max_vtx->get_fit_pt(), 15*units::cm);
      double angle = dir1.Angle(dir)/3.1415926*180.;
      if (angle > max_angle){
	max_angle = angle;
	max_sg = sg1;
      }
    }
    // check dQ/dx there ...
    if (max_vtx!=0 && max_sg !=0){
      double medium_dQ_dx = 0;
      if (max_sg->get_wcpt_vec().front().index == max_vtx->get_wcpt().index){
	medium_dQ_dx = max_sg->get_medium_dQ_dx(0,6)/(43e3/units::cm);
      }else{
	medium_dQ_dx = max_sg->get_medium_dQ_dx(int(max_sg->get_point_vec().size())-7, int(max_sg->get_point_vec().size())-1)/(43e3/units::cm);
      }
      
      double start_dQ_dx = 0;
      std::vector<double>& vec_dQ = sg->get_dQ_vec();
      std::vector<double>& vec_dx = sg->get_dx_vec();
      std::vector<double> vec_dQ_dx;
      for (size_t i=0;i!= vec_dQ.size();i++){
	vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
	//	std::cout << vec_dQ_dx.back() << std::endl;
      }
      int ncount = 0;
      if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index){
	for (auto it2 = vec_dQ_dx.begin(); it2 != vec_dQ_dx.end(); it2++){
	  if (ncount > 2) break;
	  if ((*it2) > start_dQ_dx) start_dQ_dx = *it2;
	  ncount ++;
	}
      }else{
	for (auto it2 = vec_dQ_dx.rbegin(); it2 != vec_dQ_dx.rend(); it2++){
	  if (ncount > 2) break;
	  if ((*it2) > start_dQ_dx) start_dQ_dx = *it2;
	  ncount ++;
	}
      }
      
      
      if ((Eshower < 250*units::MeV || Eshower < 500*units::MeV && shower_angle > 120 || Eshower >=500*units::MeV && shower_angle > 150))
	// 7020_1108_55428
	if (medium_dQ_dx > 1.6 && shower_angle > 60 && start_dQ_dx < 3.6 || medium_dQ_dx > 1.6 && shower_angle < 60 && start_dQ_dx < 2.5) flag_bad = true;

      if (Eshower < 800*units::cm)
	if (medium_dQ_dx > 2.0 && shower_angle > 60 || medium_dQ_dx > 2.0 && start_dQ_dx < 2.5) flag_bad = true;

      if (flag_bad) std::cout << "Xin_I2: " << Eshower << " " << medium_dQ_dx << " " << start_dQ_dx << " " << shower_angle << " " << start_dQ_dx << std::endl;
    }
  }
  
  return flag_bad;
}

bool WCPPID::NeutrinoID::bad_reconstruction_3(WCPPID::ProtoVertex* vertex, WCPPID::WCShower *shower, bool flag_print){
  bool flag_bad = false;
  TVector3 drift_dir(1,0,0);

  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }

  double main_length = shower->get_total_length(vertex->get_cluster_id());
  double total_length = shower->get_total_length();

  Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
  Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();
    
  // find the point in the main vertex that is far away from the main_vertex;
  double max_dis = 0;
  Point max_p = vertex->get_fit_pt();
  for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != vertex->get_cluster_id()) continue;
    double dis = sqrt(pow(vtx->get_fit_pt().x - vertex->get_fit_pt().x,2) + pow(vtx->get_fit_pt().y - vertex->get_fit_pt().y,2) + pow(vtx->get_fit_pt().z - vertex->get_fit_pt().z,2));
    if (dis > max_dis){
      max_dis = dis;
      max_p = vtx->get_fit_pt();
    }
  }

  // find the segment in the shower, not in main cluster, that is close to the max_p;
  double min_dis = 1e9;
  for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() == vertex->get_cluster_id()) continue;
    double dis = sg->get_closest_point(max_p).first;
    if (sg->get_length() < 6 *units::cm) continue;
    //    std::cout << sg->get_length()/units::cm << " " << dis/units::cm << std::endl;
    if (dis < min_dis) min_dis = dis;
  }
  double acc_close_length = 0;
  double min_dis1 = 1e9;
  for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() == vertex->get_cluster_id()) continue;
    double dis = sg->get_closest_point(max_p).first;
    if (dis < min_dis){
      double length1 = sg->get_length();
      acc_close_length += length1;
      if(length1 > 3*units::cm ){
	if (dis < min_dis1) min_dis1 = dis;
      }
    }
  }
  if (min_dis1 > 1e8) min_dis1 = 0;

  //std::cout << acc_close_length/units::cm << " " << min_dis1 << " " << shower->get_start_segment()->get_flag_avoid_muon_check() << std::endl;
  
  //7006_489_24469
  if (acc_close_length > 10*units::cm  || shower->get_start_segment()->get_flag_avoid_muon_check()) min_dis = min_dis1;

  
  if (min_dis < 1e7){
    if (main_length < 0.4*total_length && min_dis > 40*units::cm) flag_bad = true;
    if (main_length < 0.25*total_length && min_dis > 30*units::cm) flag_bad = true;
    if (main_length < 0.16*total_length && min_dis > 20*units::cm) flag_bad = true;
    if (main_length < 0.05*total_length && min_dis > 8*units::cm) flag_bad = true;
    if (main_length < 8*units::cm && main_length < 0.1*total_length && (min_dis > 8*units::cm && Eshower < 300*units::MeV || min_dis > 14*units::cm)) flag_bad = true;

    // 7003_1226_61332
    if (flag_bad && shower->get_start_segment()->get_flag_avoid_muon_check() && main_length > 12*units::cm &&  main_length > 0.1*total_length && min_dis < 40*units::cm) flag_bad = false;
    // 7010_405_20296
    if (map_vertex_segments[shower->get_start_vertex().first].size()==1 && main_length > 20*units::cm && min_dis < 40*units::cm &&  main_length > 0.1*total_length) flag_bad = false;
    // 7049_1874_93741
    if (flag_bad && main_length > 30*units::cm && shower->get_num_main_segments() >=4)     flag_bad = false;
  }
  
  if (flag_print && flag_bad) std::cout << "Xin_H2_1: " << Eshower << " " << main_length/units::cm << " " << total_length/units::cm << " " << min_dis/units::cm << " " << flag_bad << std::endl;
  


  if (!flag_bad){
    Point vertex_point;
    
    if (shower->get_start_segment()->get_wcpt_vec().front().index == shower->get_start_vertex().first->get_wcpt().index){
      vertex_point = shower->get_start_segment()->get_point_vec().front();
    }else{
      vertex_point = shower->get_start_segment()->get_point_vec().back();
    }

    TVector3 dir_sg = shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);

    
    TVector3 dir;
    if (shower->get_start_segment()->get_length() > 12*units::cm){
      dir = shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);
    }else{
      dir = shower->cal_dir_3vector(vertex_point,15*units::cm);
    }

   
    
    if (fabs(dir.Angle(drift_dir)/3.1415926*180.-90)<10) dir = shower->cal_dir_3vector(vertex_point,25*units::cm);

    
    Int_t ncount = 0;     Int_t ncount1 = 0;
    Int_t ncount_15 = 0;  Int_t ncount1_15 = 0;
    Int_t ncount_25 = 0;  Int_t ncount1_25 = 0;
    Int_t ncount_35 = 0;  Int_t ncount1_35 = 0;
    Int_t ncount_45 = 0;  Int_t ncount1_45 = 0;
    
    for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      //      if (sg->get_cluster_id() != vertex->get_cluster_id()) continue;
      PointVector& pts = sg->get_point_vec();
      for (size_t i=1; i+1 < pts.size();i++){
	TVector3 dir1(pts.at(i).x-vertex_point.x, pts.at(i).y-vertex_point.y, pts.at(i).z-vertex_point.z);
	double angle = dir1.Angle(dir)/3.1415926*180.;

	/* TVector3 unit_dir = dir.Unit(); */
	/* TVector3 unit_z(0,0,1); */
	/* unit_z.RotateUz(unit_dir); */
	
	/* TVector3 dir2 = dir1; */
	/* dir2.RotateUz(unit_z); */
	/* std::cout << angle << " " << dir1.Theta()/3.1415926*180. << " " << dir1.Phi()/3.1415926*180. << " " << dir2.Theta()/3.1415926*180. << " " << dir2.Phi()/3.1415926*180. << std::endl; */
	
	if (angle < 15) ncount1_15 ++;
	if (angle < 25) ncount1_25 ++;
	if (angle < 35) ncount1_35 ++;
	if (angle < 45) ncount1_45 ++;
	ncount1 ++;
	
	if (sg->get_cluster_id() == vertex->get_cluster_id()){
	  if (angle < 15) ncount_15 ++;
	  if (angle < 25) ncount_25 ++;
	  if (angle < 35) ncount_35 ++;
	  if (angle < 45) ncount_45 ++;
	  ncount ++;
	}
      }
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	TVector3 dir1((*it1)->get_fit_pt().x-vertex_point.x, (*it1)->get_fit_pt().y-vertex_point.y, (*it1)->get_fit_pt().z-vertex_point.z);
	double angle = dir1.Angle(dir)/3.1415926*180.;

	if (angle < 15) ncount1_15 ++;
	if (angle < 25) ncount1_25 ++;
	if (angle < 35) ncount1_35 ++;
	if (angle < 45) ncount1_45 ++;
	ncount1 ++;
	
	if (sg->get_cluster_id() == vertex->get_cluster_id()){
	  if (angle < 15) ncount_15 ++;
	  if (angle < 25) ncount_25 ++;
	  if (angle < 35) ncount_35 ++;
	  if (angle < 45) ncount_45 ++;
	  ncount ++;
	}
      }
    }

    

    if (ncount_45 < 0.7*ncount || ncount_25 < 0.6*ncount || ncount_25 < 0.8*ncount && ncount_15 < 0.3*ncount
	|| ncount_15 < 0.35*ncount && ncount_25 > 0.9*ncount && Eshower < 1000*units::MeV ) flag_bad = true;
    
    if (ncount1_15 < 0.35 * ncount1 && fabs(dir.Angle(drift_dir)/3.1415926*180.-90) > 15 && (ncount1_25 < 0.95 *ncount1 && Eshower < 1000*units::MeV || Eshower >=1000*units::MeV) || // 7014_1661_83096
	ncount1_15 < 0.2*ncount1 && ncount1_25 < 0.45*ncount1 && Eshower < 600*units::MeV ||
	dir_sg.Angle(dir)/3.1415926*180. > 25 && std::max(fabs(dir.Angle(drift_dir)/3.1415926*180.-90),fabs(dir_sg.Angle(drift_dir)/3.1415926*180.-90)) > 8 && (ncount1_15 < 0.8 * ncount1 && Eshower < 1000*units::MeV || Eshower >= 1000*units::MeV) || // 7001_100_5015
	dir_sg.Angle(dir)/3.1415926*180. > 20 && std::max(fabs(dir.Angle(drift_dir)/3.1415926*180.-90),fabs(dir_sg.Angle(drift_dir)/3.1415926*180.-90)) > 5 && ncount1_15 < 0.5 * ncount1
	) flag_bad = true;

    if (flag_print && flag_bad) std::cout << "Xin_H2_0: " << Eshower << " " << ncount_15/(ncount+1e-9) << " " << ncount_25/(ncount + 1e-9) << " " << ncount_35/(ncount + 1e-9) << " " << ncount_45/(ncount + 1e-9)  << " " << ncount1_15/(ncount1+1e-9) << " " << ncount1_25/(ncount1 + 1e-9) << " " << ncount1_35/(ncount1 + 1e-9) << " " << ncount1_45/(ncount1 + 1e-9) << " " << fabs(dir.Angle(drift_dir)/3.1415926*180.-90) << " " << fabs(dir_sg.Angle(drift_dir)/3.1415926*180.-90) << " " <<   dir_sg.Angle(dir)/3.1415926*180. << " " << flag_bad << std::endl; 

  }
  
  return flag_bad;
}


bool WCPPID::NeutrinoID::bad_reconstruction_2(WCPPID::ProtoVertex* vertex, WCPPID::WCShower* shower, bool flag_print){
  bool flag_bad = false;
  TVector3 dir_drift(1,0,0);
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  double total_length = shower->get_total_length();
  double total_main_length = shower->get_total_length(sg->get_cluster_id());
  double length = sg->get_length();
  double direct_length = sg->get_direct_length();

  TVector3 dir_two_end(sg->get_point_vec().front().x - sg->get_point_vec().back().x,
		       sg->get_point_vec().front().y - sg->get_point_vec().back().y,
		       sg->get_point_vec().front().z - sg->get_point_vec().back().z    );
  
  // low energy, only one segment, and not wiggled ... 7008_184_9232
  if (Eshower < 100*units::MeV && shower->get_num_segments() == 1 && (!sg->get_flag_shower_trajectory()) && direct_length/length > 0.95) flag_bad = true;
  // not wiggle ...  7049_870_43513
  if (Eshower < 100*units::MeV && total_main_length/total_length > 0.95 && length/total_length > 0.85 && (direct_length/length > 0.95 || fabs(dir_two_end.Angle(dir_drift)/3.1415926*180.-90)<5) && sg->get_flag_shower_trajectory()) flag_bad = true;
  // 7008_907_45383
  if (Eshower < 200*units::MeV && total_main_length/total_length > 0.96 && length/total_length > 0.925 && (direct_length/length > 0.95 || fabs(dir_two_end.Angle(dir_drift)/3.1415926*180.-90)<5 && sg->get_flag_shower_trajectory()) && length > 25*units::cm) flag_bad = true;
  
  
  //7014_2_109
  if (Eshower < 100*units::MeV && total_main_length/total_length > 0.95 && length/total_length > 0.95 && direct_length/length > 0.95 && sg->get_flag_shower_topology()) flag_bad = true;

  
  
  if (flag_print && flag_bad) std::cout << "Xin_H4: " << Eshower << " " << total_main_length/total_length << " " << length/total_length << " " << direct_length/length << " " << length/units::cm << " " << sg->get_flag_shower_topology() << sg->get_flag_shower_trajectory() << " " << flag_bad << " " << dir_two_end.Angle(dir_drift)/3.1415926*180. << " " << sg->get_medium_dQ_dx()/(43e3/units::cm) << " " << flag_bad << std::endl;

  if ((!flag_bad)){
    if (Eshower < 150*units::MeV && total_main_length/total_length > 0.95){
      Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
      int n_ele = 0;
      int n_other = 0;
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (sg1->get_cluster_id()!= sg->get_cluster_id()) continue;
	double medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
	
	if (sg1->get_flag_shower_topology() || sg1->get_flag_shower_trajectory() && medium_dQ_dx < 1.3 || sg1->get_direct_length()/sg1->get_length() < 0.92) n_ele ++;
	else if (medium_dQ_dx > 1.3 || sg1->get_direct_length()/sg1->get_length() > 0.95) n_other ++;

	//	std::cout << sg1->get_medium_dQ_dx()/(43e3/units::cm) << " " << sg1->get_direct_length()/sg1->get_length() << std::endl;
      }
       WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
       
      // event 7003_209_10494  & 7020_570_28512
      if (n_ele == 0 && n_other > 0 || n_ele == 1 && n_ele < n_other && n_other <=2) flag_bad = true;
      if (n_ele == 1 && n_other == 0 &&  (!fid->inside_fiducial_volume(other_vertex->get_fit_pt(), offset_x)) ) flag_bad = true;

      if (flag_print && flag_bad) std::cout << "Xin_H1_1: "  << Eshower << " " << total_main_length/total_length << " " << map_vertex_segments[vertex].size() << " " << n_ele << " " << n_other << " " << flag_bad << std::endl;
    }
  }


  
  if ((!flag_bad) && Eshower < 600*units::MeV){
    Point vertex_point;
    Point other_point;
    
    if (sg->get_wcpt_vec().front().index == shower->get_start_vertex().first->get_wcpt().index){
      vertex_point = sg->get_point_vec().front();
      other_point = sg->get_point_vec().back();
    }else{
      vertex_point = sg->get_point_vec().back();
      other_point = sg->get_point_vec().front();
    }
    TVector3 dir = sg->cal_dir_3vector(vertex_point, 15*units::cm);
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    
    double acc_length = 0;
    double total_length = 0;
    for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      if (sg1->get_cluster_id() != sg->get_cluster_id()) continue;
      double dis1 = sqrt(pow(vertex_point.x - sg1->get_point_vec().front().x,2) + pow(vertex_point.y - sg1->get_point_vec().front().y,2) + pow(vertex_point.z - sg1->get_point_vec().front().z,2));
      double dis2 = sqrt(pow(vertex_point.x - sg1->get_point_vec().back().x,2) + pow(vertex_point.y - sg1->get_point_vec().back().y,2) + pow(vertex_point.z - sg1->get_point_vec().back().z,2));
      TVector3 dir1;
      if(dis1 < dis2){
	dir1.SetXYZ(sg1->get_point_vec().back().x - sg1->get_point_vec().front().x,
		    sg1->get_point_vec().back().y - sg1->get_point_vec().front().y,
		    sg1->get_point_vec().back().z - sg1->get_point_vec().front().z);
      }else{
	dir1.SetXYZ(sg1->get_point_vec().front().x - sg1->get_point_vec().back().x,
		    sg1->get_point_vec().front().y - sg1->get_point_vec().back().y,
		    sg1->get_point_vec().front().z - sg1->get_point_vec().back().z);
      }
      double length = sg1->get_length();
      double angle = dir.Angle(dir1)/3.1415926*180.;
      if (angle > 90 && dir1.Mag() > 10*units::cm) acc_length += length;
      total_length += length;
      if (angle > 150 && dir1.Mag() > 10*units::cm) flag_bad = true;
      if (angle > 105 && sg1->get_length() > 15*units::cm) flag_bad = true;

      if (flag_print && flag_bad) std::cout << "Xin_H1_2: " << dir.Angle(dir1)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << sg1->get_length() << " " << flag_bad << std::endl;
    }
    if (acc_length > 0.33 * total_length) flag_bad = true;
    if (flag_print&& flag_bad) std::cout << "Xin_H1_3: " << Eshower << " " << acc_length/units::cm << " " << total_length/units::cm << " " << flag_bad << std::endl;


    
    if ((!flag_bad) && Eshower < 250*units::MeV){  // detection a Michel electron ...
      Point ave_p(0,0,0);
      int num_p = 0;
      double total_length = 0;
      int n_seg = 0;
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (sg1->get_cluster_id() != sg->get_cluster_id()) continue;
	if (sg1 == sg) continue;
	PointVector &pts = sg1->get_point_vec();
	for (size_t i=0;i!=pts.size();i++){
	  ave_p.x += pts.at(i).x;
	  ave_p.y += pts.at(i).y;
	  ave_p.z += pts.at(i).z;
	  num_p++;
	}
	n_seg ++;
	total_length += sg1->get_length();
      }
      if (num_p >0){
	ave_p.x /= num_p;
	ave_p.y /= num_p;
	ave_p.z /= num_p;
      
	TVector3 dir1(ave_p.x - other_point.x, ave_p.y - other_point.y,  ave_p.z - other_point.z);
	double angle1 = std::max(fabs(90-dir1.Angle(dir_drift)/3.1415926*180.), fabs(90 - dir.Angle(dir_drift)/3.1415926*180.));

	// 6090_85_4300
	if ( ((dir1.Mag() > 3*units::cm || total_length > 6*units::cm) && ((!sg->get_flag_avoid_muon_check()) || n_seg > 1)) &&  dir.Angle(dir1)/3.1415926*180. > 60 && sg->get_length() > 10*units::cm) {
	  flag_bad = true;
	}
	// 7018_888_44410
	if (flag_bad && shower->get_num_main_segments() +6 < shower->get_num_segments() && shower->get_total_length(shower->get_start_segment()->get_cluster_id()) < 0.7 * shower->get_total_length()) flag_bad = false;
	//	std::cout << shower->get_num_main_segments() << " " << shower->get_num_segments() << std::endl;
	
	if (flag_print&&flag_bad) std::cout << "Xin_H1_4: " << Eshower << " " << dir.Angle(dir1)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << total_length/units::cm << " " << n_seg << " " << sg->get_flag_avoid_muon_check() << " " << angle1 << std::endl;
      }
    }

    
    
    if ((!flag_bad) &&  Eshower < 600*units::MeV){
      WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
      double min_angle = 180;
      for (auto it = map_vertex_segments[other_vertex].begin(); it != map_vertex_segments[other_vertex].end(); it++){
	WCPPID::ProtoSegment *sg1 = *it;
	if (sg1 == sg)  continue;
	WCPPID::ProtoVertex *vtx_1 = find_other_vertex(sg1, other_vertex);
	TVector3 dir1(vtx_1->get_fit_pt().x - other_vertex->get_fit_pt().x,
		      vtx_1->get_fit_pt().y - other_vertex->get_fit_pt().y,
		      vtx_1->get_fit_pt().z - other_vertex->get_fit_pt().z);
	double angle = dir1.Angle(dir)/3.1415926*180.;
	double angle1 = std::max(fabs(90-dir.Angle(dir_drift)/3.1415926*180.), fabs(90-dir1.Angle(dir_drift)/3.1415926*180.));
	//	std::cout << angle1 << std::endl;
	// 7022_110_5530 + 6787_236_11833
	if ( (angle > 150 && angle1 > 10 && (!sg1->get_flag_shower_trajectory())) && sg1->get_length() > 7.5*units::cm && map_vertex_segments[other_vertex].size()<=4) flag_bad = true; 
	if (angle < min_angle && sg1->get_length() > 6*units::cm) min_angle = angle;

	//	std::cout << Eshower << " " << dir1.Angle(dir)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << sg1->get_length()/units::cm << " " << sg1->get_flag_shower_trajectory() << std::endl;
      }
      //
      // 6952_132_6635 
      if (Eshower < 200*units::MeV && min_angle > 60 && sg->get_length() < 0.2 * shower->get_total_length(vertex->get_cluster_id())) flag_bad = true;
      if (flag_print && flag_bad) std::cout << "Xin_H1_5: " << sg->get_length()/units::cm << " " << shower->get_total_length(vertex->get_cluster_id())/units::cm << std::endl;
    }
  }


  
  if (!flag_bad && Eshower < 150*units::MeV && shower->get_num_main_segments() <=2 && shower->get_total_length(vertex->get_cluster_id()) > shower->get_total_length() * 0.8){
    Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
    double max_dQ_dx = 0;
    for (auto it = map_seg_vtxs.begin(); it!=map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      if (sg1->get_cluster_id() != vertex->get_cluster_id()) continue;
      int ncount = int(sg1->get_point_vec().size())-5;
      for (int i=0;i<ncount;i++){
	double medium_dQ_dx = sg1->get_medium_dQ_dx(i, i+5)/(43e3/units::cm);
	//	std::cout << i << " " << medium_dQ_dx << std::endl;
	if (medium_dQ_dx > max_dQ_dx) max_dQ_dx = medium_dQ_dx;
      }
    }
    if (max_dQ_dx > 1.85) flag_bad = true;
    //    std::cout << max_dQ_dx << std::endl;
  }

  

  if (flag_print && flag_bad) std::cout << "Xin_H1: " << Eshower << " " << direct_length/length << " " << length/total_length << " " << total_main_length/total_length << " " << total_length/units::cm << " " << shower->get_num_main_segments()  << " " << flag_bad << std::endl;

  
  return flag_bad;
}

bool WCPPID::NeutrinoID::bad_reconstruction(WCPPID::WCShower* shower, bool flag_print){
  TVector3 dir_drift(1,0,0);
  bool flag_bad_shower = false;
  
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  auto pair_result = shower->get_start_vertex();
  WCPPID::ProtoVertex *vtx = pair_result.first;
  
  Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
  Map_Proto_Vertex_Segments& map_vtx_segs = shower->get_map_vtx_segs();

  if (pair_result.second == 1 && map_vertex_segments[vtx].size() == 1 && Eshower < 120*units::MeV && map_seg_vtxs.size()<=3){
    if ( (!sg->get_flag_shower_topology()) && (!sg->get_flag_shower_trajectory()) && sg->get_length() > 10*units::cm)
      flag_bad_shower = true;
    //std::cout << sg->get_direct_length()/sg->get_length() << std::endl;
  }

  //  std::cout << flag_bad_shower << std::endl;
  
  if (!flag_bad_shower){
    double max_length = 0;
    double max_dQ_dx = 0;
    double max_length_ratio = 0;
    int n_connected = 0;
    WCPPID::ProtoSegment *max_sg = 0;
    for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      double length = sg1->get_length();
      double direct_length = sg1->get_direct_length();
      double medium_dQ_dx = sg1->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
      
      //  std::cout << length/units::cm << " " << direct_length/units::cm << " " << medium_dQ_dx/(43e3/units::cm) << " " << dQ_dx_cut  << " " << sg1->get_flag_shower_topology() << std::endl;
      
      if (((!sg1->get_flag_shower_topology()) || direct_length > 0.94*length) &&(!sg1->get_flag_avoid_muon_check())){
	auto pair_vertices = find_vertices(sg1);
	double tmp_length = length;
	if (pair_vertices.first != main_vertex){
	  auto results1 = find_cont_muon_segment_nue(sg1, pair_vertices.first, true);
	  if (results1.first != 0){
	    if (((!results1.first->get_flag_shower_topology()) || results1.first->get_direct_length() > 0.94*results1.first->get_length()) &&(!results1.first->get_flag_avoid_muon_check()) )
	      tmp_length += results1.first->get_length();
	  }
	}
	if (pair_vertices.second != main_vertex){
	  auto results2 = find_cont_muon_segment_nue(sg1, pair_vertices.second, true);
	  if (results2.first != 0){
	    if (((!results2.first->get_flag_shower_topology()) || results2.first->get_direct_length() > 0.94*results2.first->get_length()) && (!results2.first->get_flag_avoid_muon_check()))
	      tmp_length += results2.first->get_length();
	  }
	}
	//7001_30_1514 + 7054_1396_69804
	double length_offset = 0;
	if (sg1->get_flag_shower_topology() || sg1->get_cluster_id() != shower->get_start_segment()->get_cluster_id()) length_offset = 6*units::cm;
	
	if (tmp_length -length_offset > max_length) {
	  max_sg = sg1;
	  max_length = tmp_length - length_offset;
	  max_dQ_dx = medium_dQ_dx;
	  max_length_ratio = sg1->get_direct_length()/sg1->get_length();
	  n_connected = 0;
	  if (pair_vertices.first != main_vertex)
	    n_connected += map_vertex_segments[pair_vertices.first].size()-1;
	  if (pair_vertices.second != main_vertex)
	    n_connected += map_vertex_segments[pair_vertices.second].size()-1;
	  // std::cout << sg1->get_id() << " " << map_vertex_segments[pair_vertices.first].size() << " " << map_vertex_segments[pair_vertices.second].size() << " " << tmp_length/units::cm << " " << n_connected << std::endl;
	}	
      }
    }
    
    TVector3 dir_shower = shower->cal_dir_3vector(shower->get_start_point(), 30*units::cm);
    
    
    
    if (Eshower < 200*units::MeV){
      if (n_connected <= 1 && max_length > 38*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 42*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 46*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 50*units::cm){
	flag_bad_shower = true;
      }
    }else if (Eshower < 400*units::MeV){
      if (n_connected <= 1 && max_length > 42*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 49*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 52*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 55*units::cm){
	flag_bad_shower = true;
      }
    }else if (Eshower < 600*units::MeV){
      if (n_connected <= 1 && max_length > 45*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 48*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 54*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 62*units::cm){
	flag_bad_shower = true;
      }
    }else if (Eshower < 800*units::MeV){
      if (n_connected <= 1 && max_length > 51*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 52*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 56*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 62*units::cm){
	flag_bad_shower = true;
	// 7018_530_2652 + 7026_397_19867
	if (map_vertex_segments[main_vertex].size()==1 && max_length < 68*units::cm || n_connected >=6 && max_length < 76*units::cm)
	  flag_bad_shower = false;
      }
    }else if (Eshower < 1500*units::MeV){
      if (n_connected <= 1 && max_length > 55*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 60*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 65*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 75*units::cm){
	flag_bad_shower = true;
      }
    }else{
      if (n_connected <= 1 && max_length > 55*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 65*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 70*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 75*units::cm){
	flag_bad_shower = true;
      }
    }

    if (Eshower > 1000*units::MeV && flag_bad_shower && max_length_ratio < 0.95){
      flag_bad_shower = false;
    }

    
    if (max_length > 0.75*shower->get_total_length() && max_length > 35*units::cm) flag_bad_shower = true;
    //    std::cout << flag_bad_shower << " " << max_length << " " << shower->get_total_length() << std::endl;

    //    std::cout << "kaka: " << Eshower << " " << n_connected << " " << max_length/units::cm << " " << max_dQ_dx/(43e3/units::cm) << " " << max_length_ratio << " " << dir_shower.Angle(dir_drift)/3.1415926*180. << " " << shower->get_start_segment()->get_flag_shower_topology() << " " << map_vertex_segments[main_vertex].size() << " " << flag_bad_shower << std::endl;

    
    if (flag_print){
      std::cout << "XinQ: " << Eshower/units::MeV << " " << max_length/units::cm << " " << n_connected << " " << max_length/shower->get_total_length() << std::endl;
    }

    // 5564_90_4517   
    if ((!flag_bad_shower) ){
      max_length = 0;
      n_connected = 0;
      //std::cout << sg->get_id() << " " << sg->get_flag_shower_topology() << " " << sg->get_flag_shower_trajectory() << " " << shower->get_num_main_segments() << " " << sg->get_length()/units::cm << std::endl;
      if ( ((!sg->get_flag_shower_topology()) || sg->get_flag_shower_topology() && Eshower < 200*units::MeV) && (!sg->get_flag_shower_trajectory()) && shower->get_num_main_segments()  == 1 ){
	double main_length = sg->get_length();
	double tmp_length1 = 0;
	double tmp_n_connected = 0;
	WCPPID::ProtoVertex *other_vtx = find_other_vertex(sg, vtx);
	if (main_length > 10*units::cm && other_vtx !=0){
	  TVector3 dir1 = sg->cal_dir_3vector(other_vtx->get_fit_pt(), 15*units::cm);
	  for (auto it1 = map_seg_vtxs.begin(); it1 != map_seg_vtxs.end(); it1++){
	    WCPPID::ProtoSegment *sg1 = it1->first;
	    if (sg1 == sg) continue;
	    auto pair_vertices = find_vertices(sg1);
	    double tmp_length = sg1->get_length();
	    if (sg1->get_flag_shower_topology() || sg1->get_flag_shower_trajectory() || tmp_length < 10*units::cm ) continue;
	    double dis1 = sqrt(pow(pair_vertices.first->get_fit_pt().x - other_vtx->get_fit_pt().x,2)+ pow(pair_vertices.first->get_fit_pt().y - other_vtx->get_fit_pt().y,2)+ pow(pair_vertices.first->get_fit_pt().z - other_vtx->get_fit_pt().z,2));
	    double dis2 = sqrt(pow(pair_vertices.second->get_fit_pt().x - other_vtx->get_fit_pt().x,2)+ pow(pair_vertices.second->get_fit_pt().y - other_vtx->get_fit_pt().y,2)+ pow(pair_vertices.second->get_fit_pt().z - other_vtx->get_fit_pt().z,2));

	    //	    std::cout << dis1/units::cm << " " << dis2/units::cm << std::endl;
	    
	    if (dis1 < 5*units::cm){
	      TVector3 dir2 = sg1->cal_dir_3vector(pair_vertices.first->get_fit_pt(), 15*units::cm);
	      double angle = dir1.Angle(dir2)/3.1415926*180.;
	      if (angle > 170){
		if (tmp_length + dis1 > tmp_length1){
		  tmp_length1 = tmp_length + dis1;
		  tmp_n_connected = int(map_vertex_segments[pair_vertices.second].size())-1;
		}
	      }
	      //std::cout << angle << std::endl;
	    }else if (dis2 < 5*units::cm){
	      TVector3 dir2 = sg1->cal_dir_3vector(pair_vertices.second->get_fit_pt(), 15*units::cm);
	      double angle = dir1.Angle(dir2)/3.1415926*180.;
	      //std::cout << angle << std::endl;
	      if (angle > 170){
		if (tmp_length + dis2 > tmp_length1){
		  tmp_length1 = tmp_length + dis2;
		  tmp_n_connected = int(map_vertex_segments[pair_vertices.first].size())-1;
		}
	      }
	    }else{
	      TVector3 dir2;
	      Point tmp_point;
	      if (dis1 < dis2){
		dir2 = sg1->cal_dir_3vector(pair_vertices.first->get_fit_pt(), 15*units::cm);
		tmp_point = pair_vertices.first->get_fit_pt();
	      }else{
		dir2 = sg1->cal_dir_3vector(pair_vertices.second->get_fit_pt(), 15*units::cm);
		tmp_point = pair_vertices.second->get_fit_pt();
	      }
	      double angle = dir1.Angle(dir2)/3.1415926*180.;
	      if (angle > 165){
		
		for (auto it2 = map_seg_vtxs.begin(); it2 != map_seg_vtxs.end(); it2++){
		  WCPPID::ProtoSegment *sg2 = it2->first;
		  if (sg2 == sg || sg2 == sg1) continue;
		  auto pair_vertices_1 = find_vertices(sg2);
		  double tmp_length2 = sg2->get_length();
		  if (tmp_length2 < 10*units::cm) continue;
		  
		  double dis3 = sqrt(pow(pair_vertices_1.first->get_fit_pt().x - other_vtx->get_fit_pt().x,2)+ pow(pair_vertices_1.first->get_fit_pt().y - other_vtx->get_fit_pt().y,2)+ pow(pair_vertices_1.first->get_fit_pt().z - other_vtx->get_fit_pt().z,2));
		  double dis4 = sqrt(pow(pair_vertices_1.second->get_fit_pt().x - other_vtx->get_fit_pt().x,2)+ pow(pair_vertices_1.second->get_fit_pt().y - other_vtx->get_fit_pt().y,2)+ pow(pair_vertices_1.second->get_fit_pt().z - other_vtx->get_fit_pt().z,2));

		  double angle1 = 0;
		  if (dis3 < 6*units::cm){
		    TVector3 dir3 = sg2->cal_dir_3vector(pair_vertices_1.first->get_fit_pt(), 15*units::cm);
		    angle1 = dir1.Angle(dir3)/3.1415926*180.;
		  }else if (dis4 < 6*units::cm){
		    TVector3 dir3 = sg2->cal_dir_3vector(pair_vertices_1.second->get_fit_pt(), 15*units::cm);
		    angle1 = dir1.Angle(dir3)/3.1415926*180.;
		  }
		  //5564_90_4517
		  if (angle1 > 170 && tmp_length2 > 0.75*std::min(dis1,dis2)){
		    if (tmp_length + tmp_length2> tmp_length1){
		      tmp_length1 = tmp_length +tmp_length2;
		      if (dis2 < dis1){
			tmp_n_connected = int(map_vertex_segments[pair_vertices.first].size())-1;
		      }else{
			tmp_n_connected = int(map_vertex_segments[pair_vertices.second].size())-1;
		      }
		    }
		  }
		  //		  std::cout << angle1 << " " << dis3/units::cm << " " << dis4/units::cm << " " << sg2->get_length() << " " << std::min(dis1,dis2) << std::endl;
		}
	      }
	    }
	    
	  } // loop over other segment

	  if (tmp_length1 + main_length > max_length){
	    max_length = tmp_length1 + main_length ;
	    n_connected =  tmp_n_connected;
	  }
	} // 
      }

      if (max_length > sg->get_length()) {
	// std::cout << Eshower << " " << max_length/units::cm << " " << n_connected << " " << sg->get_length()/units::cm << std::endl;
	
	if (Eshower < 200*units::MeV){
	  if (n_connected <= 1 && max_length > 36*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 2 && max_length > 42*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 3 && max_length > 48*units::cm){
	    flag_bad_shower = true;
	  }else if (max_length > 54*units::cm){
	    flag_bad_shower = true;
	  }
	}else if (Eshower < 400*units::MeV){
	  if (n_connected <= 1 && max_length > 45*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 2 && max_length > 42*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 3 && max_length > 42*units::cm){
	    flag_bad_shower = true;
	  }else if (max_length > 50*units::cm){
	    flag_bad_shower = true;
	  }
	}else if (Eshower < 800*units::MeV){
	  if (n_connected <= 1 && max_length > 55*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 2 && max_length > 60*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 3 && max_length > 75*units::cm){
	    flag_bad_shower = true;
	  }else if (max_length > 80*units::cm){
	    flag_bad_shower = true;
	  }
	}else if (Eshower < 1500*units::MeV){
	  if (n_connected <= 1 && max_length > 55*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 2 && max_length > 60*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 3 && max_length > 75*units::cm){
	    flag_bad_shower = true;
	  }else if (max_length > 80*units::cm){
	    flag_bad_shower = true;
	  }
	}else{
	  if (n_connected <= 1 && max_length > 50*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 2 && max_length > 60*units::cm){
	    flag_bad_shower = true;
	  }else if (n_connected == 3 && max_length > 75*units::cm){
	    flag_bad_shower = true;
	  }else if (max_length > 80*units::cm){
	    flag_bad_shower = true;
	  }
	}
      }
    }
  }
  
 
 
  
  return flag_bad_shower;
}

bool WCPPID::NeutrinoID::bad_reconstruction_1(WCPPID::WCShower* shower, bool flag_single_shower, int num_valid_tracks){
  TVector3 dir_drift(1,0,0);
  bool flag_bad_shower = false;
  // stem direct does not match with shower direction ...

  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  
  WCPPID::ProtoSegment *sg = shower->get_start_segment();
  PointVector tmp_pts;
  shower->fill_point_vec(tmp_pts, true);
  WCPPID::ProtoVertex *vertex = shower->get_start_vertex().first;
  main_cluster->Calc_PCA(tmp_pts);
  TVector3 dir1(main_cluster->get_PCA_axis(0).x,main_cluster->get_PCA_axis(0).y,main_cluster->get_PCA_axis(0).z);
  
  double angle = 0;
  double angle1 = 0;
  double angle2 = fabs(dir1.Angle(dir_drift)/3.1415926*180.-90);
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    TVector3 dir2 = sg->cal_dir_3vector(sg->get_point_vec().front(), 5*units::cm);
    TVector3 dir3 = shower->cal_dir_3vector(sg->get_point_vec().front(),30*units::cm);
    angle = dir1.Angle(dir2)/3.1415926*180.;
    if (angle > 90) angle = 180-angle;
    angle1 = fabs(dir3.Angle(dir_drift)/3.1415926*180. - 90);
  }else{
    TVector3 dir2 = sg->cal_dir_3vector(sg->get_point_vec().back(), 5*units::cm);
    TVector3 dir3 = shower->cal_dir_3vector(sg->get_point_vec().back(),30*units::cm);
    angle = dir1.Angle(dir2)/3.1415926*180.;
    if (angle > 90) angle = 180-angle;
    angle1 = fabs(dir3.Angle(dir_drift)/3.1415926*180. - 90);
  }

  double max_angle = 0;
  WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
  TVector3 dir_1 = sg->cal_dir_3vector(other_vertex->get_fit_pt(), 10*units::cm);
  for (auto it1 = map_vertex_segments[other_vertex].begin(); it1 != map_vertex_segments[other_vertex].end(); it1++){
    WCPPID::ProtoSegment *sg1 = *it1;
    if (sg1 == sg) continue;
    TVector3 dir_2 = sg1->cal_dir_3vector(other_vertex->get_fit_pt(), 10*units::cm);
    double angle2 = dir_1.Angle(dir_2)/3.1415926*180.;
    // std::cout << angle << std::endl;
    if (angle2 > max_angle) max_angle = angle2;
  }

  
  if (flag_single_shower || num_valid_tracks == 0){
    if (Eshower > 1000*units::MeV){
    }else if (Eshower > 500*units::MeV){ // high energy
      if ((angle1 >10 || angle2 > 10)&& angle > 30){
	flag_bad_shower = true;
      }
    }else{
      if ((angle > 25 && shower->get_num_main_segments() >1 || angle > 30)&& (angle1 > 7.5 || angle2 > 7.5)){
	flag_bad_shower = true;
      }
    }
  }

  // 7012_922_46106
  if (angle > 40 && (angle1 > 7.5 || angle2 > 7.5) && max_angle < 100 ) flag_bad_shower = true; 
  // 7012_785_39252
  if (angle > 20 && (angle1 > 7.5 || angle2 > 7.5) && sg->get_length() > 21 *units::cm && Eshower < 600*units::MeV && sg->get_flag_shower_trajectory()){
    flag_bad_shower = true;
  }
  
  //if (flag_bad_shower)
  //  std::cout << "kaka2: " << Eshower/units::MeV   << " " << angle << " " << angle1 << " " << angle2 << " " << flag_single_shower << " " << num_valid_tracks << " " << sg->get_length()/units::cm << " " << max_angle << " " << shower->get_num_main_segments() << std::endl;
  
  
  return flag_bad_shower;
}



