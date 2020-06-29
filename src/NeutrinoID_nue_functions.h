bool WCPPID::NeutrinoID::stem_direction(WCPPID::WCShower *max_shower, double max_energy, bool flag_print){
  bool flag_bad = false;
  TVector3 dir_drift(1,0,0);
  
  WCPPID::ProtoSegment *sg = max_shower->get_start_segment();
  PointVector tmp_pts;
  max_shower->fill_point_vec(tmp_pts, true);
  main_cluster->Calc_PCA(tmp_pts);
  TVector3 dir1(main_cluster->get_PCA_axis(0).x,main_cluster->get_PCA_axis(0).y,main_cluster->get_PCA_axis(0).z);
  
  double angle = 0;
  double angle1 = 0;
  double ratio = 0;
  double angle2 = fabs(dir1.Angle(dir_drift)/3.1415926*180.-90);
  if (main_vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    TVector3 dir2 = sg->cal_dir_3vector(sg->get_point_vec().front(), 5*units::cm);
    TVector3 dir3 = max_shower->cal_dir_3vector(sg->get_point_vec().front(),30*units::cm);
    angle = dir1.Angle(dir2)/3.1415926*180.;
    if (angle > 90) angle = 180-angle;
    angle1 = fabs(dir3.Angle(dir_drift)/3.1415926*180. - 90);
    ratio = sg->get_direct_length(0,10)/ sg->get_length(0,10);
    //	      std::cout << "Xin_A: " << max_energy/units::MeV << " " << sg->get_direct_length(0,10)/ sg->get_length(0,10) << " " << sg->get_length(0,10)/units::cm << " " << angle << " " << dir3.Angle(dir_drift)/3.1415926*180. << std::endl;
  }else{
    TVector3 dir2 = sg->cal_dir_3vector(sg->get_point_vec().back(), 5*units::cm);
    TVector3 dir3 = max_shower->cal_dir_3vector(sg->get_point_vec().back(),30*units::cm);
    int num = int(sg->get_point_vec().size())-1;
    angle = dir1.Angle(dir2)/3.1415926*180.;
    if (angle > 90) angle = 180-angle;
    angle1 = fabs(dir3.Angle(dir_drift)/3.1415926*180. - 90);
    ratio = sg->get_direct_length(num-10, num)/sg->get_length(num-10,num);
    //	      std::cout << "Xin_A: " << max_energy/units::MeV << " " << sg->get_direct_length(num-10, num)/sg->get_length(num-10,num) << " " << sg->get_length(num-10,num)/units::cm << " " << angle << " " << dir3.Angle(dir_drift)/3.1415926*180. << std::endl;
  }

  // std::cout << "Xin_B: " << max_energy/units::MeV << " " << ratio  << " " << angle << " " << angle1 << " " << angle2 << std::endl;

  if (angle > 18){
    if (max_energy > 1000*units::MeV){
    }else if (max_energy > 500*units::MeV){ // high energy
      if ((angle1 >12.5 || angle2 > 12.5) && angle > 25
	  || (angle1 > 10 || angle2 > 10) && angle > 32){
	flag_bad = true;
	if (flag_print) std::cout << "Xin_B: " << max_energy << " " << angle1 << " " << angle2 << " " << angle << " " << ratio << std::endl;
      }
    }else{
      if (angle > 25 && (angle1 > 7.5 || angle2 > 7.5)){
	flag_bad = true;
	if (flag_print) std::cout << "Xin_B1: " << max_energy << " " << angle1 << " " << angle2 << " " << angle << " " << ratio << std::endl;
      }else if ((angle1 > 7.5 || angle2 > 7.5) && ratio<0.97){
	flag_bad = true;
	if (flag_print) std::cout << "Xin_B2: " << max_energy << " " << angle1 << " " << angle2 << " " << angle << " " << ratio << std::endl;
      }
    }
  }

  return flag_bad;
}

bool WCPPID::NeutrinoID::multiple_showers(WCPPID::WCShower *max_shower, double max_energy, bool flag_print ){
  bool flag_bad = false;

  double E_total = 0;
  int nshowers = 0;
  for (auto it1 = map_vertex_to_shower[main_vertex].begin(); it1 != map_vertex_to_shower[main_vertex].end(); it1++){
    WCPPID::WCShower *shower = *it1;
    WCPPID::ProtoSegment *sg = shower->get_start_segment();
    if (sg->get_particle_type() != 11) continue;
    if (shower == max_shower) continue;
    auto pair_result = shower->get_start_vertex();
    if (pair_result.second>1) continue;
    double E_shower = 0;
    if (shower->get_kine_best() != 0){ 
      E_shower = shower->get_kine_best();
    }else{
      E_shower = shower->get_kine_charge();
    }
    bool flag_pi0 = pi0_identification(shower->get_start_vertex().first, shower->get_start_segment(), shower, 15*units::MeV);
    if (flag_pi0) continue;
    if (shower->get_total_length(shower->get_start_segment()->get_cluster_id()) < shower->get_total_length()*0.1 && shower->get_total_length(shower->get_start_segment()->get_cluster_id()) < 10*units::cm) continue;
    // 7010_532_26643
    if (shower->get_start_segment()->get_length() > 80*units::cm) continue;
    //    if (mip_identification(shower->get_start_vertex().first, shower->get_start_segment(), shower, false, false)==-1 || track_overclustering(shower) || bad_reconstruction_3(shower->get_start_vertex().first, shower) || bad_reconstruction_2(shower->get_start_vertex().first, shower) || gap_identification(shower->get_start_vertex().first, shower->get_start_segment(), false, 1, E_shower).first ) continue;
    
    if ((E_shower > 0.6 * max_energy || E_shower > 0.45 * max_energy && max_energy - E_shower < 150*units::MeV) ){
      flag_bad = true;
      if (flag_print) std::cout << "Xin_D_4: " << max_energy << " " << E_shower << std::endl;
      break;
    }
    if (E_shower > 50*units::MeV){
      E_total += E_shower;
      nshowers ++;
    }
    
    //	    std::cout << "Xin_A: " << max_energy << " " << E_shower << " " << bad_reconstruction(shower) << " " << bad_reconstruction_1(shower) << std::endl;
  }
  if ((E_total > 0.6*max_energy || max_energy < 400*units::MeV && nshowers >=2 && E_total > 0.3 * max_energy ) && (!flag_bad)){
    flag_bad = true;
    if (flag_print) std::cout << "Xin_D_3: " << max_energy << " " << E_total << std::endl;
  }

  double total_other_energy=0;
  double total_other_energy_1=0;
  int total_num_showers = 0;
  for (auto it1 = showers.begin(); it1 != showers.end() ;it1++){
    WCPPID::WCShower *shower = *it1;
    WCPPID::ProtoSegment *sg = shower->get_start_segment();
    if (sg->get_particle_type()!=11) continue;
    if (sg->get_cluster_id() == main_vertex->get_cluster_id()) continue;
    auto pair_result = shower->get_start_vertex();
    double E_shower = 0;
    if (shower->get_kine_best() != 0){ 
      E_shower = shower->get_kine_best();
    }else{
      E_shower = shower->get_kine_charge();
    }
    //  std::cout << sg->get_cluster_id() << " " << pair_result.second << " " << E_shower << " " << total_num_showers << std::endl;
    if (bad_reconstruction(shower)) continue;
    bool flag_pi0 = pi0_identification(shower->get_start_vertex().first, shower->get_start_segment(), shower, 15*units::MeV); 
    if (flag_pi0){ // 7003_1682_84132
      if (map_shower_pio_id.find(shower) != map_shower_pio_id.end()){
	std::vector<WCShower* >& tmp_showers = map_pio_id_showers[map_shower_pio_id[shower]];
	if (find(tmp_showers.begin(), tmp_showers.end(), max_shower)  != tmp_showers.end()){
	  //	  flag_pi0 = false;
	  if (E_shower > max_energy * 0.75) flag_bad = true;
	  //std::cout << "haha: " << std::endl;
	}
      }
    }

    if (flag_pi0) continue;
    //
    if (pair_result.second <=3)  {
      total_other_energy += E_shower;
      if (shower->get_start_vertex().first != main_vertex)
	total_other_energy_1 += E_shower;
      if (E_shower > 50*units::MeV)
	total_num_showers ++;
    }
    
    
    if (pair_result.second > 2) continue;
    if (E_shower > max_energy *1.2  && max_energy < 250*units::MeV){ // not the highest energy ...
      flag_bad = true;
      if (flag_print) std::cout << "Xin_D_2: " << max_energy << " " << E_shower << std::endl;
      break;
    }
  }

  //	  std::cout << "kaka: " << max_energy << " " << total_other_energy << " " << total_other_energy_1 << " "<< total_num_showers << " " << max_shower->get_kine_dQdx() << std::endl;
  
  // 7014_241_12058
  if ((max_energy < 250*units::MeV && (total_other_energy  - max_energy >200*units::MeV 
				       || total_other_energy - max_energy > 60*units::MeV && total_num_showers >=2 )
       || max_energy > 800*units::MeV && total_other_energy_1 > max_energy)
      && (!flag_bad)) {
    flag_bad = true;
    if (flag_print) std::cout << "Xin_D_1: " << max_energy << " " << total_other_energy << std::endl;
  }
	  
  return flag_bad;
}

bool WCPPID::NeutrinoID::other_showers(WCPPID::WCShower *shower, bool flag_single_shower, bool flag_print){
  bool flag_bad = false;
  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    Eshower = shower->get_kine_best();
  }else{
    Eshower = shower->get_kine_charge();
  }
  Point vertex_point;
  if (shower->get_start_vertex().first->get_wcpt().index == shower->get_start_segment()->get_wcpt_vec().front().index){
    vertex_point = shower->get_start_segment()->get_point_vec().front();
  }else{
    vertex_point = shower->get_start_segment()->get_point_vec().back();
  }

  double total_other_energy=0;
  double max_energy = 0;
    
  if (flag_single_shower){
    for (auto it1 = showers.begin(); it1 != showers.end() ;it1++){
      WCPPID::WCShower *shower1 = *it1;
      WCPPID::ProtoSegment *sg1 = shower1->get_start_segment();
      if (sg1->get_particle_type() != 11) continue;
      if (sg1->get_cluster_id() == main_vertex->get_cluster_id()) continue;
      auto pair_result = shower1->get_start_vertex();
      double E_shower1 = 0;
      if (shower1->get_kine_best() != 0){ 
	E_shower1 = shower1->get_kine_best();
      }else{
	E_shower1 = shower1->get_kine_charge();
      }
      
      //std::cout << sg->get_cluster_id() << " " << pair_result.second << " " << E_shower << std::endl;
      if (pair_result.second <=3)
	total_other_energy += E_shower1;
      if (pair_result.second > 2) continue;
	    
      if (E_shower1 > Eshower){ // not the highest energy ...
	flag_bad = true;
	max_energy = E_shower1;
	//      if (flag_print) std::cout << "Xin_A: " << max_energy << " " << E_shower << std::endl;
	break;
      }
    }
    if (Eshower < 150*units::MeV && total_other_energy > 0.27 * Eshower) {
      flag_bad = true;
      //	    if (flag_print) std::cout << "Xin_A1: " << max_energy << " " << total_other_energy << std::endl;
    }
  }
  if (flag_print && flag_bad ) std::cout << "Xin_N1: " << Eshower/units::MeV << " " << total_other_energy/units::MeV << " " << max_energy/units::MeV << " " << flag_single_shower << std::endl;


  double E_direct_max_energy = 0, E_direct_total_energy = 0;
  double E_indirect_max_energy = 0, E_indirect_total_energy = 0;
  int n_direct_showers = 0;
  int n_indirect_showers = 0;
  bool flag_direct_max_pi0 = false;
  bool flag_indirect_max_pi0 = false;
  if (!flag_bad){
    
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

      bool flag_pi0 = pi0_identification(shower1->get_start_vertex().first, shower1->get_start_segment(), shower1, 15*units::MeV);
      if (flag_pi0){ //7003_1682_84132
	if (map_shower_pio_id.find(shower1) != map_shower_pio_id.end()){
	  std::vector<WCShower* >& tmp_showers = map_pio_id_showers[map_shower_pio_id[shower1]];
	  if (find(tmp_showers.begin(), tmp_showers.end(), shower)  != tmp_showers.end()){
	    //	    flag_pi0 = false;
	    
	    //	    std::cout << E_shower1 << " " << max_energy << std::endl;

	    if (E_shower1 > Eshower * 0.75) flag_bad = true;
	    //std::cout << "haha: " << std::endl;
	  }
	}
      }

      //      std::cout << shower1->get_total_length(shower1->get_start_segment()->get_cluster_id())/units::cm << " " << shower1->get_total_length()/units::cm << std::endl;
      
      if (pair_result.second == 1){
	// 7006_387_19382
	if (shower1->get_total_length(shower1->get_start_segment()->get_cluster_id()) < shower1->get_total_length()*0.1 && shower1->get_total_length(shower1->get_start_segment()->get_cluster_id()) < 10*units::cm) continue; 
	// 7021_282_14130 
	//	if (mip_identification(shower1->get_start_vertex().first, shower1->get_start_segment(), shower1, false, false)==-1 || track_overclustering(shower1) || bad_reconstruction_3(shower1->get_start_vertex().first, shower1) || bad_reconstruction_2(shower1->get_start_vertex().first, shower1) ) continue;
	//	if (flag_pi0) continue;
	if (shower1->get_start_segment()->get_length() > 80*units::cm) continue;
	// 6090_89_4498
	if (flag_pi0 && shower1->get_total_length(shower1->get_start_segment()->get_cluster_id()) < 0.12 * shower1->get_total_length()) continue;

	//	std::cout << flag_pi0 << " " << shower1->get_total_length(shower1->get_start_segment()->get_cluster_id()) << " " <<  shower1->get_total_length() << std::endl;
	
	E_direct_total_energy += E_shower1;
	if (E_shower1 > E_direct_max_energy && shower1->get_start_vertex().first == main_vertex) {
	  E_direct_max_energy = E_shower1;
	  if (map_shower_pio_id.find(shower1) != map_shower_pio_id.end()) flag_direct_max_pi0 = true;
	  else flag_direct_max_pi0 = false;
	}
	if (E_shower1 > 80*units::MeV) n_direct_showers ++;
      }else if (pair_result.second == 2){
	
	if (shower1->get_num_segments()<=2){
	  if ( (!shower1->get_start_segment()->get_flag_shower_trajectory()) && (!shower1->get_start_segment()->get_flag_shower_topology()) && shower1->get_start_segment()->get_length() > 45*units::cm && shower1->get_start_segment()->get_length() > 0.95 * shower1->get_total_length()) continue;
	  //	std::cout << shower1->get_start_segment()->get_length()/units::cm << " " << shower1->get_start_segment()->get_flag_shower_trajectory() << " " << shower1->get_start_segment()->get_flag_shower_topology() << std::endl;
	}
	double dis = sqrt(pow(vertex_point.x - shower1->get_start_point().x,2) + pow(vertex_point.y - shower1->get_start_point().y,2) + pow(vertex_point.z - shower1->get_start_point().z,2) );

	//	std::cout << dis/units::cm << std::endl;
	double factor = 1;
	// 6090_89_4498
	if (dis > 72*units::cm) continue;
	else if (dis > 48*units::cm) factor = 0.75;
	
	//std::cout <<  E_shower1 << " " << flag_pi0 << std::endl;
	if (!flag_pi0){
	  E_indirect_total_energy += E_shower1;
	  if (E_shower1 * factor > E_indirect_max_energy) {
	    E_indirect_max_energy = E_shower1 * factor;
	    if (map_shower_pio_id.find(shower1) != map_shower_pio_id.end()) flag_indirect_max_pi0 = true;
	    else flag_indirect_max_pi0 = false;
	  }
	}
	
	if (E_shower1 > 80*units::MeV) n_indirect_showers ++;
      }
    }
    
    if (E_indirect_max_energy > Eshower + 350*units::MeV || E_direct_max_energy > Eshower) flag_bad = true;
    if (Eshower < 1000*units::MeV && n_direct_showers >0 && E_direct_max_energy > 0.33 * Eshower) flag_bad = true;
    if (Eshower >= 1000*units::MeV && n_direct_showers >0 && E_direct_max_energy > 0.33 * Eshower && E_direct_total_energy > 900*units::MeV) flag_bad = true;
    
    // 6748_57_2867 + 7004_1604_80229
    if (flag_indirect_max_pi0){
      if (Eshower < 800*units::MeV && E_indirect_total_energy - E_indirect_max_energy > Eshower && E_indirect_max_energy > 0.5 * Eshower) flag_bad = true;
    }else{
      if (Eshower < 800*units::MeV && E_indirect_total_energy > Eshower *0.6 && E_indirect_max_energy > 0.5 * Eshower) flag_bad = true;
    }    
  }
  //  if (flag_print && flag_bad)
  std::cout << "Xin_N0: " << Eshower << " " << n_direct_showers << " " << E_direct_max_energy << " " << flag_direct_max_pi0 << " " << E_direct_total_energy << " " << n_indirect_showers << " " << E_indirect_max_energy << " " << flag_indirect_max_pi0 << " " << E_indirect_total_energy << std::endl; 
  
  return flag_bad;
}

bool WCPPID::NeutrinoID::stem_length(WCPPID::WCShower *max_shower, double max_energy, bool flag_print){
  bool flag_bad = false;
  WCPPID::ProtoSegment *sg = max_shower->get_start_segment();
  if (max_energy < 500*units::MeV && sg->get_length() > 50*units::cm && (!sg->get_flag_avoid_muon_check())){
    flag_bad = true;
    if (flag_print) std::cout << "Xin_F_0: " << max_energy << " " << sg->get_length()/units::cm << std::endl;
  }
  return flag_bad;
}


bool WCPPID::NeutrinoID::vertex_inside_shower(WCPPID::WCShower *shower, bool flag_print){
  bool flag_bad = false;

  TVector3 drift_dir(1,0,0);
  TVector3 beam_dir(0,0,1);

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


  {
    // shower direction ...
    TVector3 dir1;
    if (vertex->get_wcpt().index == shower->get_start_segment()->get_wcpt_vec().front().index){
      dir1 = shower->cal_dir_3vector(shower->get_start_segment()->get_point_vec().front(),30*units::cm);
    }else{
      dir1= shower->cal_dir_3vector(shower->get_start_segment()->get_point_vec().back(), 30*units::cm);
    }

    double max_angle = 0;
    WCPPID::ProtoSegment *max_sg = 0;
    int num_good_tracks = 0;
    for (auto it = map_vertex_segments[vertex].begin(); it != map_vertex_segments[vertex].end(); it++){
      WCPPID::ProtoSegment *sg1 = *it;
      if (sg1 == shower->get_start_segment()) continue;
      TVector3 dir2 = sg1->cal_dir_3vector(vertex_point,15*units::cm);
      double angle = dir2.Angle(dir1)/3.1415926*180.;
      if ((!sg1->get_flag_shower()) && (!sg1->is_dir_weak())){
	num_good_tracks ++;
      }
      //	    std::cout << sg->get_length()/units::cm << " " << angle << std::endl;
      if (angle > max_angle && sg1->get_length() > 1.0*units::cm){
	max_angle = angle;
	max_sg = sg;
      }
    }

    if (max_sg !=0){
      double tmp_length1 = max_sg->get_length();
      double tmp_length2 = shower->get_start_segment()->get_length();
      if (map_vertex_segments[vertex].size()>=3 && Eshower < 500*units::MeV && num_good_tracks == 0
	  && (max_angle > 150  && (tmp_length1 < 15*units::cm || tmp_length2 < 15*units::cm) && std::max(tmp_length1, tmp_length2) < 25*units::cm
	      || max_angle > 170  && (tmp_length1 < 25*units::cm || tmp_length2 < 25*units::cm) && std::max(tmp_length1, tmp_length2) < 35*units::cm)){
	flag_bad = true;
	if (flag_print) std::cout << "Xin_O_3: " << Eshower << " " << map_vertex_segments[vertex].size() << " " << num_good_tracks << " " << max_angle << " " << max_sg->get_length()/units::cm << " " << shower->get_start_segment()->get_length()/units::cm << std::endl;
      }else if (map_vertex_segments[vertex].size()==2 && Eshower < 500*units::MeV && num_good_tracks == 0 && max_angle > 150 && max_sg->get_particle_type()==13 && (tmp_length1 < 35*units::cm || tmp_length2 < 35*units::cm)){
	flag_bad = true;
	if (flag_print) std::cout << "Xin_O_2: " << Eshower << " " << map_vertex_segments[vertex].size() << " " << num_good_tracks << " " << max_angle << " " << max_sg->get_length()/units::cm << " " << shower->get_start_segment()->get_length()/units::cm << std::endl;
      }
      
      //	    std::cout << "kaka: " << max_energy << " " << max_angle << " " << map_vertex_segments[main_vertex].size() << " " << max_sg->get_particle_type()  << " " << max_sg->is_dir_weak() << " " << max_sg->get_medium_dQ_dx()/(43e3/units::cm) << " " << num_good_tracks << " " << max_sg->get_length()/units::cm << " " << max_shower->get_start_segment()->get_length()/units::cm << " " << std::endl;
    }
  }
  

  if (!flag_bad){
    if (map_vertex_segments[vertex].size()>1){
      TVector3 dir_shower;
      if (shower->get_start_segment()->get_length() > 12*units::cm){
	dir_shower = shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);
      }else{
	dir_shower = shower->cal_dir_3vector(vertex_point,15*units::cm);
      }
      if (fabs(dir_shower.Angle(drift_dir)/3.1415926*180.-90)<10 || Eshower > 800*units::MeV) dir_shower = shower->cal_dir_3vector(vertex_point,25*units::cm);
      dir_shower = dir_shower.Unit();
      
      double max_angle = 0;
      double max_angle1 = 0;
      int max_weak_track = 0;
      double max_length = 0;
      double max_medium_dQ_dx = 0;
      
      double min_angle = 180;
      double min_angle1 = 0;
      int min_weak_track = 0;
      double min_length = 0;
      double min_medium_dQ_dx = 0;
      
      TVector3 dir2 = sg->cal_dir_3vector(vertex_point, 6*units::cm);
      
      
      for (auto it = map_vertex_segments[vertex].begin(); it != map_vertex_segments[vertex].end(); it++){
	WCPPID::ProtoSegment *sg1 = *it;
	if (sg1 == shower->get_start_segment()) continue;
	TVector3 dir1 = sg1->cal_dir_3vector(vertex_point, 15*units::cm);
	
	double angle = std::min(180 - dir1.Angle(dir_shower)/3.1415926*180., 180 - dir1.Angle(dir2)/3.1415926*180.);
	
	if (angle > max_angle){
	  max_angle = angle;
	  max_angle1 = fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180.;
	  max_weak_track = sg1->is_dir_weak();
	  max_length = sg1->get_length();
	  max_medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
	}
	
	if (angle < min_angle){
	  min_angle = angle;
	  min_angle1 = fabs(3.1415926/2.-dir1.Angle(drift_dir))/3.1415926*180.;
	  min_weak_track = sg1->is_dir_weak();
	  min_length = sg1->get_length();
	  min_medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
	}
	
	
      }

      // 6090_85_4300
      if (map_vertex_segments[vertex].size()==2 && (min_angle < 25 && min_weak_track ==1 || min_angle < 20)&& beam_dir.Angle(dir_shower)/3.1415926*180. > 50) flag_bad = true;
      if (map_vertex_segments[vertex].size()==2 && min_angle < 70 && min_angle1 < 10 && fabs(3.1415926/2. - drift_dir.Angle(dir_shower))/3.1415926*180. < 10 && (fabs(3.1415926/2. - drift_dir.Angle(dir_shower))/3.1415926*180. + min_angle1) < 15 && min_medium_dQ_dx > 1.5 && min_medium_dQ_dx < 2.2) {
	flag_bad = true;
	// 7001_100_5003 ... very short ...
	if (min_length < 4*units::cm && min_angle > 45 && sg->get_length() > 30*units::cm) flag_bad = false;
      }
      // 7003_1740_87003
      if (map_vertex_segments[vertex].size()==3 && (min_angle < 15 && min_medium_dQ_dx < 2.1 || min_angle < 17.5 && min_length < 5.0*units::cm && min_medium_dQ_dx < 2.5)&& (min_weak_track == 1 && max_angle > 120 || min_length < 6*units::cm && max_angle > 135 && min_angle < 12.5 && max_weak_track==1)) {
	// 7004_8_428
	if (max_length > 40*units::cm && max_weak_track == 0){
	}else{
	  flag_bad = true;
	}
      }
      
      if (map_vertex_segments[vertex].size()==3 && (min_angle < 5 && min_medium_dQ_dx < 2.1 )&& min_length < 10*units::cm && max_angle > 90 && max_weak_track ==1) flag_bad = true;
      if (map_vertex_segments[vertex].size()==3 && min_angle < 35 && min_angle1 < 10 && fabs(3.1415926/2. - drift_dir.Angle(dir_shower))/3.1415926*180. < 10 && (fabs(3.1415926/2. - drift_dir.Angle(dir_shower))/3.1415926*180. + min_angle1) < 15 && (min_medium_dQ_dx < 2.1 )&& min_weak_track == 1 && max_angle > 120) flag_bad = true;
      
      //if (flag_print && flag_bad)
	std::cout << "Xin_O_1: " << Eshower << " " <<  map_vertex_segments[vertex].size() << " " << max_angle << " " << max_angle1 << " " << max_weak_track << " " << max_length/units::cm << " " << max_medium_dQ_dx << " " << min_angle << " " << min_angle1 << " " << min_weak_track << " " << min_length/units::cm << " " << min_medium_dQ_dx << " " << fabs(3.1415926/2. - drift_dir.Angle(dir_shower))/3.1415926*180. << " " << beam_dir.Angle(dir_shower)/3.1415926*180. << std::endl;
    }
  }
  
  
  return flag_bad;
}


bool WCPPID::NeutrinoID::compare_muon_energy(WCPPID::WCShower *max_shower, double max_energy, double muon_length, bool flag_print ){
  bool flag_bad = false;

  TVector3 dir_drift(1,0,0);
  WCPPID::ProtoVertex *vertex = max_shower->get_start_vertex().first;
  WCPPID::ProtoSegment *sg = max_shower->get_start_segment();
  Point vertex_point;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    vertex_point = sg->get_point_vec().front();
  }else{
    vertex_point = sg->get_point_vec().back();
  }
  TVector3 dir_shower;
  if (max_shower->get_start_segment()->get_length() > 12*units::cm){
    dir_shower = max_shower->get_start_segment()->cal_dir_3vector(vertex_point,15*units::cm);
  }else{
    dir_shower = max_shower->cal_dir_3vector(vertex_point,15*units::cm);
  }
  if (fabs(dir_shower.Angle(dir_drift)/3.1415926*180.-90)<10 || max_energy > 800*units::MeV) dir_shower = max_shower->cal_dir_3vector(vertex_point,25*units::cm);
  dir_shower = dir_shower.Unit();

  TVector3 dir_beam(0,0,1);
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  TGraph *g_range = mp.get_muon_r2ke();
  // check muon  ...
  double E_muon = g_range->Eval(muon_length/units::cm) * units::MeV;
	  
  for (auto it1 = map_vertex_segments[main_vertex].begin(); it1 != map_vertex_segments[main_vertex].end(); it1++){
    WCPPID::ProtoSegment *sg1 = *it1;
    if (sg1->get_particle_type()==13 || sg1->get_particle_type()==2212){ 
      double length = sg1->get_length();
      double medium_dQ_dx = sg1->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
      //	      std::cout << medium_dQ_dx/(43e3/units::cm) << " " << length/units::cm << " " << dQ_dx_cut << std::endl;
      if (medium_dQ_dx < dQ_dx_cut * 43e3/units::cm){
	double tmp_energy =  g_range->Eval(sg1->get_length()/units::cm) * units::MeV;
	if (tmp_energy > E_muon) E_muon = tmp_energy;
      }
    }
  }
  
  //	  std::cout << "kaka: " << E_muon << " " << max_energy << " " << max_shower->get_kine_dQdx() << " " << muon_length/units::cm << " " << max_shower->get_total_length()/units::cm << std::endl;

  double tmp_shower_total_length = max_shower->get_total_length();

  bool flag_numuCC = (neutrino_type >> 2) & 1U;
  
  if (E_muon > max_energy && max_energy < 550*units::MeV ||
      muon_length > tmp_shower_total_length ||
      muon_length > 80*units::cm ||
      muon_length > 0.6 * tmp_shower_total_length && max_energy < 500*units::MeV // round 2 debug
      ) {
    flag_bad = true;
    if (flag_print) std::cout << "Xin_E_1: " << max_energy << " " << E_muon <<  " " << muon_length/units::cm << " " << tmp_shower_total_length/units::cm << " " << dir_beam.Angle(dir_shower)/3.1415926*180. << " " << flag_numuCC << std::endl;
  }
  //	  std::cout << "Xin_A: " << max_energy << " " << E_muon << " " << muon_length/units::cm << std::endl;
  

  return flag_bad;
}
