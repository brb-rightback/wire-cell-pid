bool WCPPID::NeutrinoID::nue_tagger(double muon_length){
  bool flag_nue = false;
  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);
  bool flag_print = true;
  
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
	if (sg->get_particle_type()==13) continue;

	//	if (flag_print) std::cout << "Qian_I: " << sg->get_id() << " " << sg->get_particle_type() << " " << shower->get_kine_charge()/units::MeV << " " << shower->get_kine_best()/units::MeV << std::endl;
	
	// higher than 70 MeV, and connect to the vertex ...
	if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end() && (shower->get_kine_charge() > 80*units::MeV && shower->get_kine_best()==0 || shower->get_kine_best() > 80*units::MeV)){

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

	  if (flag_print) std::cout << "Qian_O: " << shower->get_start_segment()->get_id() << " " << E_shower/units::MeV << " " << bad_reconstruction(shower) << " " << bad_reconstruction_1(shower) << " " << low_energy_overlapping(shower) << std::endl;

	  
	  if (bad_reconstruction(shower)) continue;	// bad reconstruction (mis ided track as shower)
	  if (bad_reconstruction_1(shower)) continue;	// bad reconstruction ... ( stem does not match with shower direction)
	  if (low_energy_overlapping(shower)) continue; // low energy overlapping situation 
	
	  
	 
	  
	  //std::cout << shower->get_kine_charge()/units::MeV  << std::endl;
	  //	  std::cout << ": " << shower->get_total_length(main_vertex->get_cluster_id())/units::cm << " " << muon_length/units::cm << std::endl;
	  if (flag_print) std::cout << "Qian_A: " << E_shower/units::MeV << " " << shower->get_total_length(main_vertex->get_cluster_id())/units::cm << " " << muon_length/units::cm << std::endl;
	  
	  if (shower->get_total_length(main_vertex->get_cluster_id()) < muon_length * 0.5) continue;

	  
	  if (!gap_identification(main_vertex, sg)){ // gap id
	    int mip_id;
	    if (flag_single_shower && E_shower < 400*units::MeV){
	      mip_id = mip_identification(main_vertex,sg, shower, true); // dQ/dx id
	    }else{
	      mip_id = mip_identification(main_vertex,sg, shower); // dQ/dx id
	    }

	    if (flag_print) std::cout << "Qian_B: " << E_shower/units::MeV << " " << mip_id << std::endl;
	    
	    
	    if (mip_id == 1){
	      good_showers.insert(shower);
	      //	      flag_nue = true;
	    }else if (mip_id == 0){
	      if (!pi0_identification(main_vertex, sg, shower)){ // pi0 identification ...
		good_showers.insert(shower);
		//flag_nue = true;
	      }
	    }
	  }else{ // gap id
	    if (flag_print) std::cout << "Qian_C: gap founded " << E_shower/units::MeV << std::endl; 
	  }
	} // 80 MeV, and connected
      } // loop over showers

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
	

	
	if (flag_nue && flag_single_shower){ // single shower ...
	  double total_other_energy=0;
	  for (auto it1 = showers.begin(); it1 != showers.end() ;it1++){
	    WCPPID::WCShower *shower = *it1;
	    WCPPID::ProtoSegment *sg = shower->get_start_segment();
	    if (sg->get_particle_type()==13) continue;
	    if (sg->get_cluster_id() == main_vertex->get_cluster_id()) continue;
	    auto pair_result = shower->get_start_vertex();
	    double E_shower = 0;
	    if (shower->get_kine_best() != 0){ 
	      E_shower = shower->get_kine_best();
	    }else{
	      E_shower = shower->get_kine_charge();
	    }

	    //std::cout << sg->get_cluster_id() << " " << pair_result.second << " " << E_shower << std::endl;
	    if (pair_result.second <=3)
	      total_other_energy += E_shower;
	    if (pair_result.second > 2) continue;
	    

	    if (E_shower > max_energy){ // not the highest energy ...
	      flag_nue = false;
	      if (flag_print) std::cout << "Xin_A: " << max_energy << " " << E_shower << std::endl;
	      break;
	    }
	  }
	  if (max_energy < 150*units::MeV && total_other_energy > 0.27 * max_energy && flag_nue) {
	    flag_nue = false;
	    if (flag_print) std::cout << "Xin_A1: " << max_energy << " " << total_other_energy << std::endl;
	  }
	  //	  std::cout << "Xin_B: " << max_energy << " " << total_other_energy << std::endl;
	  
	  if (flag_nue){ // beginning not consistent with shower itself
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
	    //	    std::cout << "Xin_B: " << max_energy/units::MeV << " " << ratio  << " " << angle << " " << angle1 << " " << angle2 << std::endl;
	    if (angle > 18){
	      if (max_energy > 1000*units::MeV){
	      }else if (max_energy > 500*units::MeV){ // high energy
		if ((angle1 >10 || angle2 > 10)&& angle > 25){
		  flag_nue = false;
		  if (flag_print) std::cout << "Xin_B: " << max_energy << " " << angle1 << " " << angle2 << " " << angle << " " << ratio << std::endl;
		}
	      }else{
		if (angle > 25 && (angle1 > 7.5 || angle2 > 7.5)){
		  flag_nue = false;
		  if (flag_print) std::cout << "Xin_B1: " << max_energy << " " << angle1 << " " << angle2 << " " << angle << " " << ratio << std::endl;
		}else if ((angle1 > 7.5 || angle2 > 7.5) && ratio<0.97){
		  flag_nue = false;
		  if (flag_print) std::cout << "Xin_B2: " << max_energy << " " << angle1 << " " << angle2 << " " << angle << " " << ratio << std::endl;
		}
	      }
	    }
	  }
	}

	if (flag_nue){ // no multiple EM showers
	  double E_total = 0;
	  int nshowers = 0;
	  for (auto it1 = map_vertex_to_shower[main_vertex].begin(); it1 != map_vertex_to_shower[main_vertex].end(); it1++){
	    WCPPID::WCShower *shower = *it1;
	    WCPPID::ProtoSegment *sg = shower->get_start_segment();
	    if (sg->get_particle_type()==13) continue;
	    if (shower == max_shower) continue;
	    auto pair_result = shower->get_start_vertex();
	    if (pair_result.second>1) continue;
	    double E_shower = 0;
	    if (shower->get_kine_best() != 0){ 
	      E_shower = shower->get_kine_best();
	    }else{
	      E_shower = shower->get_kine_charge();
	    }
	    if ((E_shower > 0.6 * max_energy || E_shower > 0.45 * max_energy && max_energy - E_shower < 150*units::MeV) ){
	      flag_nue = false;
	      if (flag_print) std::cout << "Xin_C: " << max_energy << " " << E_shower << std::endl;
	      break;
	    }
	    if (E_shower > 50*units::MeV){
	      E_total += E_shower;
	      nshowers ++;
	    }

	    //	    std::cout << "Xin_A: " << max_energy << " " << E_shower << " " << bad_reconstruction(shower) << " " << bad_reconstruction_1(shower) << std::endl;
	  }
	  if ((E_total > 0.6*max_energy || max_energy < 400*units::MeV && nshowers >=2 && E_total > 0.3 * max_energy ) && flag_nue){
	    flag_nue = false;
	    if (flag_print) std::cout << "Xin_C1: " << max_energy << " " << E_total << std::endl;
	  }

	  double total_other_energy=0;
	  int total_num_showers = 0;
	  for (auto it1 = showers.begin(); it1 != showers.end() ;it1++){
	    WCPPID::WCShower *shower = *it1;
	    WCPPID::ProtoSegment *sg = shower->get_start_segment();
	    if (sg->get_particle_type()==13) continue;
	    if (sg->get_cluster_id() == main_vertex->get_cluster_id()) continue;
	    auto pair_result = shower->get_start_vertex();
	    double E_shower = 0;
	    if (shower->get_kine_best() != 0){ 
	      E_shower = shower->get_kine_best();
	    }else{
	      E_shower = shower->get_kine_charge();
	    }
	    if (bad_reconstruction(shower)) continue;
	  
	    //
	    if (pair_result.second <=3) 	  {
	      total_other_energy += E_shower;
	      if (E_shower > 50*units::MeV)
		total_num_showers ++;
	    }
	    //  std::cout << sg->get_cluster_id() << " " << pair_result.second << " " << E_shower << " " << total_num_showers << std::endl;
	    if (pair_result.second > 2) continue;
	    if (E_shower > max_energy *1.2  && max_energy < 250*units::MeV){ // not the highest energy ...
	      flag_nue = false;
	      if (flag_print) std::cout << "Xin_D: " << max_energy << " " << E_shower << std::endl;
	      break;
	    }
	  }
	  //	  std::cout << max_energy << " " << total_other_energy << " "<< total_num_showers << std::endl;
	  
	  if (max_energy < 250*units::MeV && (total_other_energy  - max_energy >200*units::MeV 
					      || total_other_energy - max_energy > 60*units::MeV && total_num_showers >=2 ) && flag_nue) {
	    flag_nue = false;
	    if (flag_print) std::cout << "Xin_D1: " << max_energy << " " << total_other_energy << std::endl;
	  }
	  
	  

	}

	if (flag_nue){ // compare with muon energy ...
	  TPCParams& mp = Singleton<TPCParams>::Instance();
	  TGraph *g_range = mp.get_muon_r2ke();
	  // check muon  ...
	  double E_muon = g_range->Eval(muon_length/units::cm) * units::MeV;
	  
	  for (auto it1 = map_vertex_segments[main_vertex].begin(); it1 != map_vertex_segments[main_vertex].end(); it1++){
	    WCPPID::ProtoSegment *sg1 = *it1;
	    if (sg1->get_particle_type()==13 || sg1->get_particle_type()==2212){ 
	      double length = sg->get_length();
	      double medium_dQ_dx = sg->get_medium_dQ_dx();
	      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);//0.85+0.95 *sqrt(25./ (length / units::cm));
	      if (medium_dQ_dx < dQ_dx_cut * 43e3/units::cm){
		double tmp_energy =  g_range->Eval(sg1->get_length()/units::cm) * units::MeV;
		if (tmp_energy > E_muon) E_muon = tmp_energy;
	      }
	    }
	  }
	  if (E_muon > max_energy) {
	    flag_nue = false;
	    if (flag_print) std::cout << "Xin_E: " << max_energy << " " << E_muon << std::endl;
	  }
	  //	  std::cout << "Xin_A: " << max_energy << " " << E_muon << " " << muon_length/units::cm << std::endl;
	}

	if (flag_nue){ // long stem is not preferred ...
	  WCPPID::ProtoSegment *sg = max_shower->get_start_segment();
	  if (max_energy < 500*units::MeV && sg->get_length() > 50*units::cm){
	    flag_nue = false;
	    if (flag_print) std::cout << "Xin_F: " << max_energy << " " << sg->get_length()/units::cm << std::endl;
	  }
	  //	  std::cout << "Xin_A: " << sg->get_length()/units::cm << " " << sg->get_direct_length()/units::cm << " " << max_energy << std::endl;
	}

	if (flag_nue){ // vertex inside the shower ...
	  // shower direction ...
	  TVector3 dir1;
	  if (main_vertex->get_wcpt().index == max_shower->get_start_segment()->get_wcpt_vec().front().index){
	    dir1 = max_shower->cal_dir_3vector(max_shower->get_start_segment()->get_point_vec().front(),30*units::cm);
	  }else{
	    dir1= max_shower->cal_dir_3vector(max_shower->get_start_segment()->get_point_vec().back(), 30*units::cm);
	  }

	  double max_angle = 0;
	  WCPPID::ProtoSegment *max_sg = 0;
	  int num_good_tracks = 0;
	  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
	    WCPPID::ProtoSegment *sg = *it;
	    if (sg == max_shower->get_start_segment()) continue;
	    TVector3 dir2 = sg->cal_dir_3vector(main_vertex->get_fit_pt(),15*units::cm);
	    double angle = dir2.Angle(dir1)/3.1415926*180.;
	    if ((!sg->get_flag_shower()) && (!sg->is_dir_weak())){
	      num_good_tracks ++;
	    }
	    //	    std::cout << sg->get_length()/units::cm << " " << angle << std::endl;
	    if (angle > max_angle && sg->get_length() > 1.0*units::cm){
	      max_angle = angle;
	      max_sg = sg;
	    }
	  }

	  if (max_sg !=0){
	    double tmp_length1 = max_sg->get_length();
	    double tmp_length2 = max_shower->get_start_segment()->get_length();
	    if (map_vertex_segments[main_vertex].size()>=3 && num_good_tracks == 0 && max_energy < 500*units::MeV && max_angle > 150  && (tmp_length1 < 15*units::cm || tmp_length2 < 15*units::cm) && std::max(tmp_length1, tmp_length2) < 25*units::cm){
	      flag_nue = false;
	      if (flag_print) std::cout << "Xin_G: " << max_energy << " " << map_vertex_segments[main_vertex].size() << " " << num_good_tracks << " " << max_angle << " " << max_sg->get_length()/units::cm << " " << max_shower->get_start_segment()->get_length()/units::cm << std::endl;
	    }

	    //   std::cout << "Xin_A: " << max_energy << " " << max_angle << " " << map_vertex_segments[main_vertex].size() << " " << max_sg->get_particle_type() << " " << max_sg->get_length()/units::cm << " " << max_sg->is_dir_weak() << " " << max_sg->get_medium_dQ_dx()/(43e3/units::cm) << " " << num_good_tracks << " " << max_shower->get_start_segment()->get_length()/units::cm << " " << std::endl;
	  }
	  
	}

	if (flag_nue){
	  if (bad_reconstruction_2(main_vertex, max_shower)) {
	    flag_nue = false;
	    if (flag_print) std::cout << "Xin_H: " << max_energy << std::endl;
	  }
	}
	
	//	bad_reconstruction(max_shower, true);
	std::cout << "Xin: " << good_showers.size() << " " << flag_nue << " " << max_energy/units::MeV << " " << dir.Angle(dir_beam)/3.1415926*180. << std::endl;
      }
    } // has a shower ...
  } // main vertex
  



  
  if (flag_nue){
    neutrino_type |= 1UL << 5; //nue
  }

  return flag_nue;
}

bool WCPPID::NeutrinoID::gap_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment* sg){
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
    bool flag_prolong_u = flag.at(0);
    bool flag_prolong_v = flag.at(1);
    bool flag_prolong_w = flag.at(2);
    bool flag_parallel = flag.at(3);

    //    std::cout << flag_prolong_u << " " << flag_prolong_v << " " << flag_prolong_w << " " << flag_parallel << std::endl;
    
    for (int i=0;i+1 < pts.size();i++){
      Point test_p;
      for (int j=0;j!=3;j++){
	test_p.x = pts.at(i).x + j/3.*(pts.at(i+1).x - pts.at(i).x);
	test_p.y = pts.at(i).y + j/3.*(pts.at(i+1).y - pts.at(i).y);
	test_p.z = pts.at(i).z + j/3.*(pts.at(i+1).z - pts.at(i).z);

	
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

	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 1);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 1, 0)) num_bad_ch ++;
	    else if (flag_prolong_v) num_spec ++;
	  }

	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 2);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 2, 0)) num_bad_ch ++;
	  }
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
    bool flag_prolong_u = flag.at(0);
    bool flag_prolong_v = flag.at(1);
    bool flag_prolong_w = flag.at(2);
    bool flag_parallel = flag.at(3);
    
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

	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 1);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 1, 0)) 	      num_bad_ch ++;
	    else if (flag_prolong_v) num_spec ++;
	  }

	  tmp_pts = ct_point_cloud->get_closest_points(test_p, 0.2*units::cm, 2);
	  if (tmp_pts.pts.size()>0){
	    num_connect ++;
	  }else{
	    if (ct_point_cloud->get_closest_dead_chs(test_p, 2, 0)) 	      num_bad_ch ++;
	    else if (flag_prolong_w) num_spec ++;
	  }
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
  
  //  std::cout << sg->get_id() << " Xin " << n_bad << " " << n_points << " " << flag_start << std::endl;
  if (n_bad >0) flag_gap = true;
  
  return flag_gap;
}

int WCPPID::NeutrinoID::mip_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg, WCPPID::WCShower *shower, bool flag_strong_check){
  int mip_id = 1; 
  // 1 good, -1 bad, 0 not sure ...

  std::vector<double> vec_dQ_dx;

  std::vector<double>& vec_dQ = sg->get_dQ_vec();
  std::vector<double>& vec_dx = sg->get_dx_vec();

  
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    for (size_t i=0;i!=vec_dQ.size();i++){
      vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
      if (vec_dQ_dx.size() >=20) break;
    }
  }else{
    for (int i=int(vec_dQ.size())-1;i>=0;i--){
      vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
      if (vec_dQ_dx.size() >=20) break;
    }
  }

  std::vector<int> vec_threshold(vec_dQ_dx.size(), 0);
  for (size_t i=0;i!=vec_dQ_dx.size(); i++){
    if (vec_dQ_dx.at(i)>1.3) vec_threshold.at(i) = 1;
  }

  int n_end_reduction = 0;
  double prev_vec_dQ_dx = vec_dQ_dx.front();
  for (size_t i=1;i<vec_dQ_dx.size();i++){
    if (vec_dQ_dx.at(i) < prev_vec_dQ_dx){
      n_end_reduction = i;
      prev_vec_dQ_dx = vec_dQ_dx.at(i);
      if (vec_dQ_dx.at(i) < 1.3) break;
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

  {
    auto it = map_vertex_to_shower.find(vertex);
    if (it != map_vertex_to_shower.end()){
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::WCShower *shower = *it1;
	WCPPID::ProtoSegment *sg = shower->get_start_segment();
	if (sg->get_particle_type()==13) continue;
	if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end()){
	  n_showers ++;
	}
      }
    }
    for (auto it1 = map_vertex_segments[vertex].begin(); it1 != map_vertex_segments[vertex].end(); it1++){
      WCPPID::ProtoSegment *sg = *it1;
      if (sg->get_flag_shower()) continue;
      n_tracks ++;
      if (sg->get_particle_type()==2212) n_protons ++;
    }
  }

  double Eshower = 0;
  if (shower->get_kine_best() != 0){ 
    //std::cout << "Xin: "<< sg->get_id() << " " << n_end_reduction << " " << n_first_mip << " " << n_first_non_mip << " " << n_first_non_mip_1 << " " << n_first_non_mip_2 << " " << shower->get_kine_best()/units::MeV  << " " << n_showers << " " << n_tracks << " " << n_protons << " " << map_vertex_segments[vertex].size() << std::endl;
    Eshower = shower->get_kine_best();
  }else{
    // std::cout << "Xin: "<< sg->get_id() << " " << n_end_reduction << " " << n_first_mip << " " << n_first_non_mip << " " << n_first_non_mip_1 << " " << n_first_non_mip_2 << " "  << shower->get_kine_charge()/units::MeV << " " << n_showers << " " << n_tracks << " " << n_protons << " " << map_vertex_segments[vertex].size() << std::endl;
    Eshower = shower->get_kine_charge();
  }

  //  std::cout << "Xin: "<< sg->get_id() << " " << n_end_reduction << " " << n_first_mip << " " << n_first_non_mip << " " << n_first_non_mip_1 << " " << n_first_non_mip_2 << " "  << shower->get_kine_charge()/units::MeV << " " << n_showers << " " << n_tracks << " " << n_protons << " " << map_vertex_segments[vertex].size() << std::endl;
  
  // quality check ...
  double medium_dQ_dx = 1;
  {
    std::vector<double> tmp_vec_dQ_dx = vec_dQ_dx;
    std::nth_element(tmp_vec_dQ_dx.begin(), tmp_vec_dQ_dx.begin() + tmp_vec_dQ_dx.size()/2, tmp_vec_dQ_dx.end());
    medium_dQ_dx  = *std::next(tmp_vec_dQ_dx.begin(), tmp_vec_dQ_dx.size()/2);
  }
  
  // shower split ...
  if (n_showers == 2 && n_tracks ==0 && Eshower < 500*units::MeV) mip_id = -1; // bad
  else if (n_first_non_mip_2 - n_first_mip >=2 && // dQ_dx cut ...
	   (n_first_mip <=2 || (n_first_mip <= n_end_reduction &&
				(n_first_mip <=3 || n_first_mip <=4 && Eshower > 600*units::MeV || n_first_mip <=6 && Eshower > 1000*units::MeV) ))
	   ) mip_id = 1;
  /* else if (n_first_non_mip_2 - n_first_mip >=1 && // dQ_dx cut ... */
  /* 	   (n_first_mip <=2 || */
  /* 	    n_first_mip <= n_end_reduction && */
  /* 	    (n_first_mip <=3 || n_first_mip <=4 && Eshower > 600*units::MeV || n_first_mip <=6 && Eshower > 1000*units::MeV) )){ */
  /*   if(n_first_mip > 1){  */
  /*     if (vec_dQ_dx.at(n_first_mip-1) < 1.5) mip_id = 1; // nue CC */
  /*     else mip_id = -1; */
  /*   }else{ */
  /*     mip_id = -1; */
  /*   } */
  /* } */
  else mip_id = -1;

  if (Eshower < 600*units::MeV){
    int n_good_tracks = 0;
    for (auto it = map_vertex_segments[vertex].begin(); it!= map_vertex_segments[vertex].end(); it++){
      WCPPID::ProtoSegment *sg1 = *it;
      if (sg1 == sg) continue;
      if (!sg1->is_dir_weak() || sg1->get_length() > 10*units::cm) n_good_tracks ++;
    }
    // 6043_4_243 event
    if (n_good_tracks >1 && n_first_non_mip_2 <=2) mip_id = -1;
    //    std::cout << n_good_tracks << " " << n_first_non_mip_2 << " " << Eshower << std::endl;
  }
  
  if (flag_strong_check){
    if (!((n_first_mip <=2 || (n_first_mip <= n_end_reduction &&
			       (n_first_mip <=3 || n_first_mip <=4 && Eshower > 600*units::MeV || n_first_mip <=6 && Eshower > 1000*units::MeV))) && n_first_non_mip_2 - n_first_mip > 3)) mip_id = -1;
	//  std::cout << "Xin_A:" << n_first_mip << " " << n_first_non_mip_2 - n_first_mip << std::endl;
  }


  
  if (mip_id == 1 &&  map_vertex_segments[vertex].size() ==1 && Eshower < 500*units::MeV){
    //    std::cout << vec_dQ_dx.front() << " " << medium_dQ_dx << std::endl;
    if (Eshower < 180*units::MeV || n_first_mip>0 || vec_dQ_dx.front() > 1.15 && n_end_reduction >= n_first_mip && Eshower < 360*units::MeV)
      mip_id = 0;
  }

  //  std::cout << mip_id << " " << map_vertex_segments[vertex].size() << " " << Eshower/units::MeV << std::endl;
  
  
  
  if (medium_dQ_dx < 0.75 && Eshower < 150*units::MeV) mip_id = -1;
    //std::cout << "Xin: " << medium_dQ_dx << " " << Eshower/units::MeV << std::endl;
  
  {
    // check overlapping situation ...
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
      }
      //std::cout << "Xin_B: " << i << " " << min_u/units::cm << " " << min_v/units::cm << " " << min_w/units::cm << std::endl;
      if (min_u < 0.3*units::cm && min_v < 0.3*units::cm && min_w < 0.3*units::cm) mip_id = -1;
    }
    
    
  }
     
  
  return mip_id;
}


bool WCPPID::NeutrinoID::low_energy_overlapping(WCPPID::WCShower* shower){
  bool flag_overlap = false;

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
    }else if (map_vertex_segments[vtx].size()>1 && Eshower < 250*units::MeV && shower->get_total_length(sg->get_cluster_id()) < 20*units::cm){
      TVector3 dir1 = sg->cal_dir_3vector(vtx_point, 5*units::cm);
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
	if ((*it)==sg) continue;
	TVector3 dir2 = (*it)->cal_dir_3vector(vtx_point, 5*units::cm);
	
	if ((*it)->get_length() < 30*units::cm && dir1.Angle(dir2)/3.1415926*180.< 10 && (*it)->get_particle_type() == 13){
	  //	  std::cout << shower->get_total_length(sg->get_cluster_id())/units::cm << std::endl;
	  flag_overlap = true;
	  //	  std::cout << "B: " << std::endl;
	}
	// std::cout << dir1.Angle(dir2)/3.1415926*180. << " " << (*it)->get_length()/units::cm << " " << (*it)->get_particle_type()<< std::endl;
      }
      //std::cout << Eshower/units::MeV << std::endl;
    }else if (map_vertex_segments[vtx].size()==2 && Eshower < 100*units::MeV ){
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

    if (map_vertex_segments[vtx].size()==1 && shower->get_total_length(sg->get_cluster_id()) < 15*units::cm &&  Eshower > 30*units::MeV ){
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
      //std::cout << n_out << " " << n_sum << std::endl;
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


bool WCPPID::NeutrinoID::pi0_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg, WCPPID::WCShower *shower){
  bool flag_pi0 = false;

  
  auto it = map_shower_pio_id.find(shower);
  if (it != map_shower_pio_id.end()){
    std::vector<WCShower*> tmp_pi0_showers = map_pio_id_showers[it->second];
    auto mass_pair = map_pio_id_mass[it->second];

    if (fabs(mass_pair.first - 135*units::MeV)<35*units::MeV && mass_pair.second == 1 || fabs(mass_pair.first - 135*units::MeV) < 60 * units::MeV && mass_pair.second == 2){
      double Eshower_1 = tmp_pi0_showers.front()->get_kine_charge();
      double Eshower_2 = tmp_pi0_showers.back()->get_kine_charge();
      //      std::cout << tmp_pi0_showers.size() << " " << mass_pair.first << " " << mass_pair.second << " " << Eshower_1/units::MeV << " " << Eshower_2/units::MeV << std::endl;
      if (std::min(Eshower_1, Eshower_2) > 15*units::MeV && fabs(Eshower_1 - Eshower_2)/(Eshower_1 + Eshower_2) < 0.87)
	flag_pi0 = true;
    }
  }else{
    TVector3 dir1 = sg->cal_dir_3vector(vertex->get_fit_pt(), 12*units::cm);
    if (dir1.Mag() >0){
      for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
	WCPPID::ProtoVertex *vtx1 = it->first;
	if (vtx1->get_cluster_id() == vertex->get_cluster_id()) continue;
	TVector3 dir2(vtx1->get_fit_pt().x - vertex->get_fit_pt().x,vtx1->get_fit_pt().y - vertex->get_fit_pt().y, vtx1->get_fit_pt().z - vertex->get_fit_pt().z);
	if (dir2.Mag()>0){
	  if (dir2.Mag() < 36*units::cm && 180 - dir1.Angle(dir2)/3.1415926*180. < 7.5) flag_pi0 = true;
	  //	  std::cout << dir1.Angle(dir2)/3.1415926*180. << " " << dir2.Mag() << std::endl;
	}
      }
    }
  }
  
  
  return flag_pi0;
}


bool WCPPID::NeutrinoID::bad_reconstruction_2(WCPPID::ProtoVertex* vertex, WCPPID::WCShower* shower){
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
  
  //7014_2_109
  if (Eshower < 100*units::MeV && total_main_length/total_length > 0.95 && length/total_length > 0.95 && direct_length/length > 0.95 && sg->get_flag_shower_topology()) flag_bad = true;
  
  // std::cout << Eshower << " " << total_main_length/total_length << " " << length/total_length << " " << direct_length/length << " " << sg->get_flag_shower_trajectory() << " " << flag_bad << " " << dir_two_end.Angle(dir_drift)/3.1415926*180. << " " << sg->get_medium_dQ_dx()/(43e3/units::cm) << std::endl;

  if ((!flag_bad)){
    if (Eshower < 150*units::MeV && total_main_length/total_length > 0.95){
      Map_Proto_Segment_Vertices& map_seg_vtxs = shower->get_map_seg_vtxs();
      int n_ele = 0;
      int n_other = 0;
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (sg1->get_cluster_id()!= sg->get_cluster_id()) continue;
	double medium_dQ_dx = sg1->get_medium_dQ_dx()/(43e3/units::cm);
	
	if (sg1->get_flag_shower_topology() || sg1->get_flag_shower_trajectory() && medium_dQ_dx < 1.3) n_ele ++;
	else if (medium_dQ_dx > 1.3 || sg1->get_direct_length()/sg1->get_length() > 0.95) n_other ++;
	//	std::cout << sg1->get_medium_dQ_dx()/(43e3/units::cm) << " " << sg1->get_direct_length()/sg1->get_length() << std::endl;
      }
       WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
       
      // event 7003 209 10494  & 7020_570_28512
      if (n_ele == 0 && n_other > 0 || n_ele == 1 && n_ele < n_other) flag_bad = true;
      if (n_ele == 1 && n_other == 0 &&  (!fid->inside_fiducial_volume(other_vertex->get_fit_pt(), offset_x)) ) flag_bad = true;
      //std::cout << "XinQ_A: "  << Eshower << " " << total_main_length/total_length << " " << map_vertex_segments[vertex].size() << " " << n_ele << " " << n_other << std::endl;
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
      //      std::cout << dir.Angle(dir1)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << sg1->get_length() << std::endl;
    }
    if (acc_length > 0.33 * total_length) flag_bad = true;
    //std::cout << Eshower << " " << acc_length/units::cm << " " << total_length/units::cm << std::endl;

    
    if ((!flag_bad) && Eshower < 250*units::MeV){
      Point ave_p(0,0,0);
      int num_p = 0;
      double total_length = 0;
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
	total_length += sg1->get_length();
      }
      if (num_p >0){
	ave_p.x /= num_p;
	ave_p.y /= num_p;
	ave_p.z /= num_p;
      
	TVector3 dir1(ave_p.x - other_point.x, ave_p.y - other_point.y,  ave_p.z - other_point.z);
	if ((dir1.Mag() > 3*units::cm || total_length > 6*units::cm)&&  dir.Angle(dir1)/3.1415926*180. > 60 && sg->get_length() > 10*units::cm) {
	  flag_bad = true;
	}
	//	std::cout << Eshower << " " << dir.Angle(dir1)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << total_length/units::cm << std::endl;
      }
    }
    
    if ((!flag_bad) &&  Eshower < 600*units::MeV){
      WCPPID::ProtoVertex *other_vertex = find_other_vertex(sg, vertex);
      for (auto it = map_vertex_segments[other_vertex].begin(); it != map_vertex_segments[other_vertex].end(); it++){
	WCPPID::ProtoSegment *sg1 = *it;
	if (sg1 == sg)  continue;
	WCPPID::ProtoVertex *vtx_1 = find_other_vertex(sg1, other_vertex);
	TVector3 dir1(vtx_1->get_fit_pt().x - other_vertex->get_fit_pt().x,
		      vtx_1->get_fit_pt().y - other_vertex->get_fit_pt().y,
		      vtx_1->get_fit_pt().z - other_vertex->get_fit_pt().z);
	
	if (dir1.Angle(dir)/3.1415926*180. > 150 && sg1->get_length() > 7.5*units::cm && map_vertex_segments[other_vertex].size()<=4) flag_bad = true; 
	//	std::cout << Eshower << " " << dir1.Angle(dir)/3.1415926*180. << " " << dir1.Mag()/units::cm << " " << sg1->get_length()/units::cm << std::endl;
      }      
    }
  }



  //  std::cout << "Xin_A: " << Eshower << " " << direct_length/length << " " << length/total_length << " " << total_main_length/total_length << " " << total_length/units::cm << " " << shower->get_num_main_segments()  << std::endl;

  
  return flag_bad;
}

bool WCPPID::NeutrinoID::bad_reconstruction(WCPPID::WCShower* shower, bool flag_print){
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
  }

  //  std::cout << flag_bad_shower << std::endl;
  
  if (!flag_bad_shower){
    double max_length = 0;
    int n_connected = 0;
    for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      double length = sg1->get_length();
      double direct_length = sg1->get_direct_length();
      double medium_dQ_dx = sg1->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);

      //  std::cout << length/units::cm << " " << direct_length/units::cm << " " << medium_dQ_dx/(43e3/units::cm) << " " << dQ_dx_cut  << " " << sg1->get_flag_shower_topology() << std::endl;
      
      if ((!sg1->get_flag_shower_topology()) || direct_length > 0.94*length){
	auto pair_vertices = find_vertices(sg1);
	double tmp_length = length;
	if (pair_vertices.first != main_vertex){
	  auto results1 = find_cont_muon_segment_nue(sg1, pair_vertices.first, true);
	  if (results1.first != 0){
	    if ((!results1.first->get_flag_shower_topology()) || results1.first->get_direct_length() > 0.94*results1.first->get_length() )
	      tmp_length += results1.first->get_length();
	  }
	}
	if (pair_vertices.second != main_vertex){
	  auto results2 = find_cont_muon_segment_nue(sg1, pair_vertices.second, true);
	  if (results2.first != 0){
	    if ((!results2.first->get_flag_shower_topology()) || results2.first->get_direct_length() > 0.94*results2.first->get_length())
	      tmp_length += results2.first->get_length();
	  }
	}
	
	if (tmp_length > max_length) {
	  max_length = tmp_length;
	  n_connected = 0;
	  if (pair_vertices.first != main_vertex)
	    n_connected += map_vertex_segments[pair_vertices.first].size()-1;
	  if (pair_vertices.second != main_vertex)
	    n_connected += map_vertex_segments[pair_vertices.second].size()-1;
	  //	  std::cout << sg1->get_id() << " " << map_vertex_segments[pair_vertices.first].size() << " " << map_vertex_segments[pair_vertices.second].size() << " " << tmp_length/units::cm << std::endl;
	}	
      }
    }
    
    //    std::cout << Eshower << " " << n_connected << " " << max_length << std::endl;
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
      if (n_connected <= 1 && max_length > 42*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 48*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 54*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 62*units::cm){
	flag_bad_shower = true;
      }
    }else if (Eshower < 800*units::MeV){
      if (n_connected <= 1 && max_length > 46*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 52*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 3 && max_length > 56*units::cm){
	flag_bad_shower = true;
      }else if (max_length > 62*units::cm){
	flag_bad_shower = true;
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

    if (max_length > 0.75*shower->get_total_length()) flag_bad_shower = true;
    //  std::cout << flag_bad_shower << std::endl;
    
    if (flag_print){
      std::cout << "XinQ: " << Eshower/units::MeV << " " << max_length/units::cm << " " << n_connected << " " << max_length/shower->get_total_length() << std::endl;
    }

    
    if (!flag_bad_shower){
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

	    //    std::cout << dis1/units::cm << " " << dis2/units::cm << std::endl;
	    
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
	    }
	  } // loop over other segment

	  if (tmp_length1 + main_length > max_length){
	    max_length = tmp_length1 + main_length ;
	    n_connected =  tmp_n_connected;
	  }
	} // 
      }

      //      std::cout << Eshower << " " << max_length/units::cm << " " << n_connected << std::endl;

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
  
 
 
  
  return flag_bad_shower;
}

bool WCPPID::NeutrinoID::bad_reconstruction_1(WCPPID::WCShower* shower){
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
  //	  std::cout << "Xin_A: " << max_energy/units::MeV   << " " << angle << " " << angle1 << " " << angle2 << std::endl;
  if (Eshower > 1000*units::MeV){
  }else if (Eshower > 500*units::MeV){ // high energy
    if ((angle1 >10 || angle2 > 10)&& angle > 25){
      flag_bad_shower = true;
    }
  }else{
    if (angle > 25 && (angle1 > 7.5 || angle2 > 7.5)){
      flag_bad_shower = true;
    }
  }

  return flag_bad_shower;
}



void WCPPID::NeutrinoID::examine_showers(){

  std::map<WCPPID::ProtoSegment *, WCPPID::WCShower*> map_merge_seg_shower;

  TVector3 drift_dir(1,0,0);
  std::set<WCPPID::WCShower*> del_showers;
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;

    if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) continue;
    
    
    if (sg->get_length() > 45*units::cm) continue;
    //std::cout << sg->get_id() << " " << sg->get_length()/units::cm << " " << sg->get_direct_length()/units::cm << std::endl;
    WCPPID::ProtoVertex *vtx = find_other_vertex(sg, main_vertex);
    
    bool flag_checked = false;
    
    /* for (auto it1 = map_vertex_segments[vtx].begin(); it1 != map_vertex_segments[vtx].end(); it1++){ */
    /*   WCPPID::ProtoSegment *sg1 = *it1; */
    /*   if (sg1 == sg) continue; */
    /*   std::cout << sg1->get_length()/units::cm << " " << sg1->get_particle_type() << " " << sg1->is_dir_weak() << std::endl; */
    /* } */
    
    // case I ...
    if (map_vertex_to_shower.find(vtx) != map_vertex_to_shower.end()){
      TVector3 dir1 = sg->cal_dir_3vector(vtx->get_fit_pt(), 15*units::cm);
      for (auto it1 = map_vertex_to_shower[vtx].begin(); it1 != map_vertex_to_shower[vtx].end(); it1++){
	WCPPID::WCShower *shower = *it1;
	if (shower->get_start_segment()->get_particle_type()==13) continue;
	auto pair_result = shower->get_start_vertex();
	
	TVector3 dir2 = shower->cal_dir_3vector(shower->get_start_point(), 100*units::cm);
	
	double Eshower = 0;
	if (shower->get_kine_best() != 0){ 
	  Eshower = shower->get_kine_best();
	}else{
	  Eshower = shower->get_kine_charge();
	}
	
	std::cout << Eshower/units::MeV << " " << dir1.Angle(dir2)/3.1415926*180. << " " << pair_result.second << " " << dir1.Angle(drift_dir)/3.1415926*180. << " " << dir2.Angle(drift_dir)/3.1415926*180. << " " << sg->get_length()/units::cm << std::endl;
	
	if (Eshower > 100*units::MeV && pair_result.second == 1) flag_checked = true;
	
	if (Eshower > 800*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 30
	    || Eshower > 150*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 10
	    || Eshower > 100*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 10 && sg->get_length() < 25*units::cm
	    || Eshower > 250*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 15
	    || Eshower > 360*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 25 
	    || Eshower > 100*units::MeV && Eshower <= 150*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 15 && sg->get_length() < 15*units::cm && flag_checked
	    ){
	  //std::cout << shower->get_kine_charge()/units::MeV << " " << dir1.Angle(dir2)/3.1415926*180. << " " << pair_result.second << " " << dir1.Angle(drift_dir)/3.1415926*180. << " " << dir2.Angle(drift_dir)/3.1415926*180. << std::endl;
	  map_merge_seg_shower[sg] = shower;
	  continue;
	}
      }
    }
    
    if (map_merge_seg_shower.find(sg) != map_merge_seg_shower.end()) continue;
    
    if (flag_checked) continue;
    if (map_vertex_to_shower.find(main_vertex) != map_vertex_to_shower.end()){
      TVector3 dir1 = sg->cal_dir_3vector(main_vertex->get_fit_pt(), 15*units::cm);
      for (auto it1 = map_vertex_to_shower[main_vertex].begin(); it1 != map_vertex_to_shower[main_vertex].end(); it1++){
	WCPPID::WCShower *shower = *it1;
	if (shower->get_start_segment()->get_particle_type()==13) continue;
	auto pair_result = shower->get_start_vertex();
	if (pair_result.second == 1) continue; //direct connected ...
	TVector3 dir2(shower->get_start_point().x - main_vertex->get_fit_pt().x,
		      shower->get_start_point().y - main_vertex->get_fit_pt().y,
		      shower->get_start_point().z - main_vertex->get_fit_pt().z);

	//	std::cout << shower->get_kine_charge() << " " << dir1.Angle(dir2)/3.1415926*180. << " "<< sg->get_length()/units::cm << std::endl;
	if ((shower->get_kine_charge() > 150*units::MeV &&  dir1.Angle(dir2)/3.1415926*180. < 10 )&& sg->get_length() > 5*units::cm ){
	  map_merge_seg_shower[sg] = shower;
	  continue;
	}
      }
    } // find shower ...
  
  } // loop over segment


  // merge ...
  //  std::cout << "Xin: " << map_merge_seg_shower.size() << std::endl;
  for (auto it = map_merge_seg_shower.begin(); it != map_merge_seg_shower.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    WCPPID::WCShower *shower = it->second;

    std::cout << "EM shower modification: " << shower->get_start_segment()->get_id() << " -> " << sg->get_id() << std::endl;
    
    auto pair_result = shower->get_start_vertex();
    
    for (auto it1 = map_vertex_to_shower[pair_result.first].begin(); it1!= map_vertex_to_shower[pair_result.first].end(); it1++){
      WCPPID::WCShower *shower1 = *it1;
      if (shower == shower1) continue;
      if (shower1->get_start_vertex().second==1){
	del_showers.insert(shower1);
      }
    }
    
    if (pair_result.second != 1){
      shower->add_segment(sg, map_segment_vertices);
      shower->set_start_vertex(main_vertex, 1);
      shower->set_start_segment(sg);
      shower->set_start_point(main_vertex->get_fit_pt());
      std::set<WCPPID::ProtoSegment* > tmp_used_segments;
      shower->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, tmp_used_segments);
      sg->set_particle_type(11);
      shower->update_particle_type();
      shower->calculate_kinematics();
      //std::cout << shower->get_kine_charge()/units::MeV << std::endl;
      double kine_charge = cal_kine_charge(shower);
      shower->set_kine_charge(kine_charge);
      shower->set_flag_kinematics(true);
      //      std::cout << kine_charge/units::MeV << std::endl;
    }else{
      shower->add_segment(sg, map_segment_vertices);
      shower->set_start_vertex(main_vertex, 1);
      shower->set_start_segment(sg);
      shower->set_start_point(main_vertex->get_fit_pt());
      std::set<WCPPID::ProtoSegment* > tmp_used_segments;
      shower->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, tmp_used_segments);
      sg->set_particle_type(11);
      shower->update_particle_type();
      shower->calculate_kinematics();
      double kine_charge = cal_kine_charge(shower);
      shower->set_kine_charge(kine_charge);
      shower->set_flag_kinematics(true);
      //std::cout << kine_charge/units::MeV << std::endl;
    }
  }
  
  //std::cout << del_showers.size() << std::endl;
  for (auto it1 = del_showers.begin(); it1!=del_showers.end(); it1 ++){
    WCPPID::WCShower *shower1 = *it1;
    if (shower1->get_start_vertex().first != main_vertex){
      showers.erase(find(showers.begin(), showers.end(), shower1));
      delete shower1;
    }
  }
  
  if (map_merge_seg_shower.size()>0)
    update_shower_maps();
}
