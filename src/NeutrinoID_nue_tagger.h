bool WCPPID::NeutrinoID::nue_tagger(double muon_length){
  bool flag_nue = false;
  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);

 
  
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
	
	// higher than 70 MeV, and connect to the vertex ...
	if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end() && (shower->get_kine_charge() > 70*units::MeV && shower->get_kine_best()==0 || shower->get_kine_best() > 70*units::MeV)){

	  if (bad_reconstruction(shower)) continue;	// bad reconstruction ...
	  if (bad_reconstruction_1(shower)) continue;	// bad reconstruction ...
	  if (low_energy_overlapping(shower)) continue; // low energy overlapping 
	
	  
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
	  
	  //  std::cout << shower->get_kine_charge()/units::MeV << " " << muon_kine_energy/units::MeV << std::endl;
	  
	  // std::cout << "Xin: " << shower->get_total_length(main_vertex->get_cluster_id())/units::cm << " " << muon_length/units::cm << std::endl;
	  if (shower->get_total_length(main_vertex->get_cluster_id()) < muon_length * 0.5) continue;
	  
	  if (!gap_identification(main_vertex, sg)){ // gap id
	    int mip_id;
	    if (flag_single_shower && E_shower < 400*units::MeV){
	      mip_id = mip_identification(main_vertex,sg, shower, true); // dQ/dx id
	    }else{
	      mip_id = mip_identification(main_vertex,sg, shower); // dQ/dx id
	    }
	    
	    if (mip_id == 1){
	      good_showers.insert(shower);
	      //	      flag_nue = true;
	    }else if (mip_id == 0){
	      if (!pi0_identification(main_vertex, sg, shower)){ // pi0 identification ...
		good_showers.insert(shower);
		//flag_nue = true;
	      }
	    }
	  } // gap id
 	  
	} // 70 MeV, and connected
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
	  for (auto it1 = showers.begin(); it1 != showers.end() ;it1++){
	    WCPPID::WCShower *shower = *it1;
	    WCPPID::ProtoSegment *sg = shower->get_start_segment();
	    if (sg->get_particle_type()==13) continue;
	    if (sg->get_cluster_id() == main_vertex->get_cluster_id()) continue;
	    auto pair_result = shower->get_start_vertex();
	    if (pair_result.second > 2) continue;
	    double E_shower = 0;
	    if (shower->get_kine_best() != 0){ 
	      E_shower = shower->get_kine_best();
	    }else{
	      E_shower = shower->get_kine_charge();
	    }
	    if (E_shower > max_energy){ // not the highest energy ...
	      flag_nue = false;
	      break;
	    }
	  }
	  if (flag_nue){
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
		}
	      }else{
		if (angle > 25 && (angle1 > 7.5 || angle2 > 7.5)){
		  flag_nue = false;
		}else if ((angle1 > 7.5 || angle2 > 7.5) && ratio<0.97){
		  flag_nue = false;
		}
	      }
	    }
	  }
	}

	if (flag_nue){
	  double E_total = 0;
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
	      break;
	    }
	    if (E_shower > 70*units::MeV)
	      E_total += E_shower;
	    //	    std::cout << "Xin_A: " << max_energy << " " << E_shower << " " << bad_reconstruction(shower) << " " << bad_reconstruction_1(shower) << std::endl;
	  }
	  if (E_total > 0.6*max_energy) flag_nue = false;
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
	  if (E_muon > max_energy) flag_nue = false;
	  //	  std::cout << "Xin_A: " << max_energy << " " << E_muon << " " << muon_length/units::cm << std::endl;
	}
	
	std::cout << "Xin: " << good_showers.size() << " " << flag_nue << " " << max_energy/units::MeV << " " << dir.Angle(dir_beam)/3.1415926*180. << std::endl;
      }
    } // has a shower ...
  } // main vertex
  

  if (flag_nue){
    neutrino_type |= 1UL << 5; //numu
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
	    //	    std::cout << dir1.Angle(dir2)/3.1415926*180. << std::endl;
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

bool WCPPID::NeutrinoID::bad_reconstruction(WCPPID::WCShower* shower){
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

  if (!flag_bad_shower){
    double max_length = 0;
    int n_connected = 0;
    for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      double length = sg1->get_length();
      double direct_length = sg1->get_direct_length();
      double medium_dQ_dx = sg1->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);

      //   std::cout << length << " " << direct_length << " " << medium_dQ_dx/(43e3/units::cm) << " " << dQ_dx_cut  << " " << sg1->get_flag_shower_topology() << std::endl;
      
      if ((!sg1->get_flag_shower_topology()) || direct_length > 0.94*length){
	auto pair_vertices = find_vertices(sg1);
	double tmp_length = length;
	auto results1 = find_cont_muon_segment(sg1, pair_vertices.first, true);
	if (results1.first != 0){
	  tmp_length += results1.first->get_length();
	}
	auto results2 = find_cont_muon_segment(sg1, pair_vertices.second, true);
	if (results2.first != 0){
	  tmp_length += results2.first->get_length();
	}
	if (tmp_length > max_length) {
	  max_length = tmp_length;
	  n_connected = map_vertex_segments[pair_vertices.first].size() + map_vertex_segments[pair_vertices.second].size()-2;
	  //	  std::cout << sg1->get_id() << " " << map_vertex_segments[pair_vertices.first].size() << " " << map_vertex_segments[pair_vertices.second].size() << " " << tmp_length/units::cm << std::endl;
	}
	
      }
    }

    if (Eshower < 400*units::MeV){
      if (n_connected == 1 && max_length > 45*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 42*units::cm){
	flag_bad_shower = true;
      }
    }else if (Eshower < 800*units::MeV){
      if (n_connected <= 1 && max_length > 55*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 60*units::cm){
	flag_bad_shower = true;
      }
    }else if (Eshower < 1500*units::MeV){
      if (n_connected <= 1 && max_length > 55*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 60*units::cm){
	flag_bad_shower = true;
      }
    }else{
      if (n_connected <= 1 && max_length > 55*units::cm){
	flag_bad_shower = true;
      }else if (n_connected == 2 && max_length > 60*units::cm){
	flag_bad_shower = true;
      }
    }
    //    std::cout << "Xin: " << Eshower/units::MeV << " " << max_length/units::cm << " " << n_connected << std::endl;
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
  
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) continue;
    WCPPID::ProtoVertex *vtx = find_other_vertex(sg, main_vertex);

    bool flag_checked = false;
    
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
	
	if (Eshower > 100*units::MeV && pair_result.second == 1) flag_checked = true;
	
	if (Eshower > 800*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 30
	    || Eshower > 150*units::MeV && 180 - dir1.Angle(dir2)/3.1415926*180. < 10){
	  //	  std::cout << shower->get_kine_charge()/units::MeV << " " << dir1.Angle(dir2)/3.1415926*180. << " " << pair_result.second << " " << dir1.Angle(drift_dir)/3.1415926*180. << " " << dir2.Angle(drift_dir)/3.1415926*180. << std::endl;
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
	if (shower->get_kine_charge() > 150*units::MeV &&  dir1.Angle(dir2)/3.1415926*180. < 10 && sg->get_length() > 5*units::cm ){
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
      //      std::cout << shower->get_kine_charge()/units::MeV << std::endl;
      double kine_charge = cal_kine_charge(shower);
      shower->set_kine_charge(kine_charge);
      shower->set_flag_kinematics(true);
      //std::cout << kine_charge/units::MeV << std::endl;
    }else{
      shower->add_segment(sg, map_segment_vertices);
      shower->set_start_vertex(main_vertex, 1);
      shower->set_start_segment(sg);
      shower->set_start_point(main_vertex->get_fit_pt());
      sg->set_particle_type(11);
      shower->update_particle_type();
      shower->calculate_kinematics();
      double kine_charge = cal_kine_charge(shower);
      shower->set_kine_charge(kine_charge);
      shower->set_flag_kinematics(true);
    }
  }
  if (map_merge_seg_shower.size()>0)
    update_shower_maps();
}
