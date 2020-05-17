bool WCPPID::NeutrinoID::nue_tagger(double muon_length){
  bool flag_nue = false;

  // check main_vertex ...
  {
    auto it = map_vertex_to_shower.find(main_vertex);
    if (it != map_vertex_to_shower.end()){
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::WCShower *shower = *it1;
	WCPPID::ProtoSegment *sg = shower->get_start_segment();
	if (sg->get_particle_type()==13) continue;
	
	if (bad_reconstruction(shower)) continue;
	
	// higher than 70 MeV, and connect to the vertex ...
	if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end() && (shower->get_kine_charge() > 70*units::MeV && shower->get_kine_best()==0 || shower->get_kine_best() > 70*units::MeV)){
	  double E_shower = shower->get_kine_charge();
	  //  std::cout << shower->get_kine_charge()/units::MeV << " " << muon_kine_energy/units::MeV << std::endl;
	  
	  // std::cout << "Xin: " << shower->get_total_length(main_vertex->get_cluster_id())/units::cm << " " << muon_length/units::cm << std::endl;
	  if (shower->get_total_length(main_vertex->get_cluster_id()) < muon_length * 0.5) continue;
	  
	  if (!gap_identification(main_vertex, sg)){ // gap id
	    int mip_id = mip_identification(main_vertex,sg, shower);
	    if (mip_id == 1){
	      flag_nue = true;
	    }else if (mip_id == 0){
	      if (!pi0_identification(main_vertex, sg)){
		flag_nue = true;
	      }
	    }
	  } // gap id
	} // 70 MeV, and connected
      } // loop over showers
    } // has a shower ...
  } // main vertex
  

  if (flag_nue){
    neutrino_type |= 1UL << 5; //numu
  }

  return flag_nue;
}

bool WCPPID::NeutrinoID::gap_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment* sg){
  bool flag_gap = false;

  bool flag_start;
  if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
    flag_start = true; // front ... 
  }else{
    flag_start = false; // back ...
  }

  PointVector &pts = sg->get_point_vec();

  int n_points = 0;
  int n_bad = 0;
  
  if (flag_start){
    Point closest_p = pts.back();
    double min_dis = 1e9;
    for (int i=0;i<pts.size();i++){
      double dis = fabs(sqrt(pow(pts.at(i).x - vertex->get_fit_pt().x,2) + pow(pts.at(i).y - vertex->get_fit_pt().y,2) + pow(pts.at(i).z - vertex->get_fit_pt().z,2))-3*units::cm);
      if (dis < min_dis){
        min_dis = dis;
	closest_p = pts.at(i);
      }
    }

    TVector3 dir(closest_p.x - vertex->get_fit_pt().x, closest_p.y - vertex->get_fit_pt().y, closest_p.z - vertex->get_fit_pt().z);
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
      double dis = sqrt(pow(pts.at(i+1).x - vertex->get_fit_pt().x,2) + pow(pts.at(i+1).y - vertex->get_fit_pt().y,2) + pow(pts.at(i+1).z - vertex->get_fit_pt().z,2));
      if (dis > 2.4*units::cm) break;
    }
  }else{
    Point closest_p = pts.back();
    double min_dis = 1e9;
    for (int i=0;i<pts.size();i++){
      double dis = fabs(sqrt(pow(pts.at(i).x - vertex->get_fit_pt().x,2) + pow(pts.at(i).y - vertex->get_fit_pt().y,2) + pow(pts.at(i).z - vertex->get_fit_pt().z,2))-3*units::cm);
      if (dis < min_dis){
        min_dis = dis;
	closest_p = pts.at(i);
      }
    }

    TVector3 dir(closest_p.x - vertex->get_fit_pt().x, closest_p.y - vertex->get_fit_pt().y, closest_p.z - vertex->get_fit_pt().z);
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
      double dis = sqrt(pow(pts.at(i-1).x - vertex->get_fit_pt().x,2) + pow(pts.at(i-1).y - vertex->get_fit_pt().y,2) + pow(pts.at(i-1).z - vertex->get_fit_pt().z,2));
      if (dis > 2.4*units::cm) break;
    }
  }
  
  //  std::cout << sg->get_id() << " Xin " << n_bad << " " << n_points << " " << flag_start << std::endl;
  if (n_bad >0) flag_gap = true;
  
  return flag_gap;
}

int WCPPID::NeutrinoID::mip_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg, WCPPID::WCShower *shower){
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
  //}

  
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
    //std::cout << "Xin: "<< sg->get_id() << " " << n_end_reduction << " " << n_first_mip << " " << n_first_non_mip << " " << n_first_non_mip_1 << " " << n_first_non_mip_2 << " " << shower->get_kine_best()/units::MeV  << " " << n_showers << " " << n_tracks << " " << n_protons << std::endl;
    Eshower = shower->get_kine_best();
  }else{
    //    std::cout << "Xin: "<< sg->get_id() << " " << n_end_reduction << " " << n_first_mip << " " << n_first_non_mip << " " << n_first_non_mip_1 << " " << n_first_non_mip_2 << " "  << shower->get_kine_charge()/units::MeV << " " << n_showers << " " << n_tracks << " " << n_protons << std::endl;
    Eshower = shower->get_kine_charge();
  }

  if (n_showers == 2 && n_tracks ==0 && Eshower < 500*units::MeV) mip_id = -1; // bad
  else if (n_first_non_mip_2 - n_first_mip >=2 &&
	   (n_first_mip <=2 || (n_first_mip == n_end_reduction &&
				(n_first_mip <=3 || n_first_mip <=4 && Eshower > 600*units::MeV || n_first_mip <=6 && Eshower > 1500*units::MeV) ))
	   ) mip_id = 1;
  else mip_id = -1;

  double medium_dQ_dx = 1;
  {
    std::nth_element(vec_dQ_dx.begin(), vec_dQ_dx.begin() + vec_dQ_dx.size()/2, vec_dQ_dx.end());
    medium_dQ_dx  = *std::next(vec_dQ_dx.begin(), vec_dQ_dx.size()/2);
    if (medium_dQ_dx < 0.75 && Eshower < 150*units::MeV) mip_id = -1;
    //std::cout << "Xin: " << medium_dQ_dx << " " << Eshower/units::MeV << std::endl;
  }
 
  
  
  return mip_id;
}

bool WCPPID::NeutrinoID::pi0_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg){
  bool flag_pi0 = false;

  return flag_pi0;
}

bool WCPPID::NeutrinoID::bad_reconstruction(WCPPID::WCShower* shower){
  bool flag_bad_shower = false;
  
  return flag_bad_shower;
}
