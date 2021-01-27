 void WCPPID::NeutrinoID::collect_2D_charges(){
  double min_ts=1e9, max_ts=-1e9;
  double min_uch=1e9, max_uch=-1e9;
  double min_vch=1e9, max_vch=-1e9;
  double min_wch=1e9, max_wch=-1e9;

  std::vector<float> results = main_cluster->get_time_ch_range();
  if (results.at(0) < min_ts) min_ts = results.at(0);
  if (results.at(1) > max_ts) max_ts = results.at(1);
  if (results.at(2) < min_uch) min_uch = results.at(2);
  if (results.at(3) > max_uch) max_uch = results.at(3);
  if (results.at(4) < min_vch) min_vch = results.at(4);
  if (results.at(5) > max_vch) max_vch = results.at(5);
  if (results.at(6) < min_wch) min_wch = results.at(6);
  if (results.at(7) > max_wch) max_wch = results.at(7);

  for (auto it = other_clusters.begin(); it!= other_clusters.end(); it++){
    results = (*it)->get_time_ch_range();
    if (results.at(0) < min_ts) min_ts = results.at(0)-1;
    if (results.at(1) > max_ts) max_ts = results.at(1)+1;
    if (results.at(2) < min_uch) min_uch = results.at(2)-1;
    if (results.at(3) > max_uch) max_uch = results.at(3)+1;
    if (results.at(4) < min_vch) min_vch = results.at(4)-1;
    if (results.at(5) > max_vch) max_vch = results.at(5)+1;
    if (results.at(6) < min_wch) min_wch = results.at(6)-1;
    if (results.at(7) > max_wch) max_wch = results.at(7)+1;
  }
  charge_2d_u = ct_point_cloud->get_overlap_good_ch_charge(min_ts, max_ts, min_uch, max_uch, 0);
  charge_2d_v = ct_point_cloud->get_overlap_good_ch_charge(min_ts, max_ts, min_vch, max_vch, 1);
  charge_2d_w = ct_point_cloud->get_overlap_good_ch_charge(min_ts, max_ts, min_wch, max_wch, 2);

  //   std::cout << min_ts << " " << max_ts << " " << min_uch << " " << max_uch << " " << min_vch << " " << max_vch << " " << min_wch << " " << max_wch << " " << charge_2d_u.size() << " " << charge_2d_v.size() << " " << charge_2d_w.size() << std::endl;
  main_cluster->fill_2d_charge_dead_chs(charge_2d_u, charge_2d_v, charge_2d_w);
  for (auto it = other_clusters.begin(); it!= other_clusters.end(); it++){
    (*it)->fill_2d_charge_dead_chs(charge_2d_u, charge_2d_v, charge_2d_w);
  }
  // std::cout << min_ts << " " << max_ts << " " << min_uch << " " << max_uch << " " << min_vch << " " << max_vch << " " << min_wch << " " << max_wch << " " << charge_2d_u.size() << " " << charge_2d_v.size() << " " << charge_2d_w.size() << std::endl;
  //  for (auto it = charge_2d_w.begin(); it!=charge_2d_w.end(); it++){
  // if (it->first.first > 2400  || it->first.second < 4800 && it->first.second > 8256) std::cout << "Something wrong! " << it->first.first << " " << it->first.second << std::endl;
  //}
} 


double WCPPID::NeutrinoID::cal_kine_charge(WCPPID::WCShower *shower){

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_offset = mp.get_time_offset();
  int nrebin = mp.get_nrebin();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w + 0.5;
  double offset_u = -first_u_dis/pitch_u + 0.5;
  double offset_v = -first_v_dis/pitch_v + 0.5;
  
  double first_t_dis = main_cluster->get_point_cloud()->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - main_cluster->get_point_cloud()->get_cloud().pts[0].x;
  double slope_xt = 1./time_slice_width;
  double offset_t =  first_t_dis/time_slice_width + 0.5;



  
  double kine_energy = 0;
  
  // to be improved ...
  double fudge_factor = 0.95;
  double recom_factor  = 0.7;

  if (shower->get_flag_shower()) {
    recom_factor = 0.5; // assume shower
    fudge_factor = 0.8; // shower ...
  }else if (fabs(shower->get_particle_type())==2212){
    recom_factor = 0.35;
  }

  //std::cout << fudge_factor << " " << recom_factor << std::endl;

  
  double sum_u_charge = 0;
  double sum_v_charge = 0;
  double sum_w_charge = 0;

  shower->rebuild_point_clouds();
  WCP::ToyPointCloud* pcloud1 = shower->get_associated_pcloud();
  WCP::ToyPointCloud* pcloud2 = shower->get_fit_pcloud();
  
  if (pcloud1 == 0 && pcloud2 == 0) return 0;
  else if (pcloud1 ==0 && pcloud2 !=0) pcloud1 = pcloud2;
  else if (pcloud2 ==0 && pcloud1 !=0) pcloud2 = pcloud1;
  
  for (auto it = charge_2d_u.begin(); it!= charge_2d_u.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 0);
    double dis = 1e9;
    int point_index = -1;
    if (pcloud1!=0) {
      auto results = pcloud1->get_closest_2d_dis_index(p2d.first, p2d.second, 0);
      dis = results.first;
      point_index = results.second;
    }
    if (dis < 0.6*units::cm) {
      Point test_p(pcloud1->get_cloud().pts[point_index].x, pcloud1->get_cloud().pts[point_index].y, pcloud1->get_cloud().pts[point_index].z);
      double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
      sum_u_charge += it->second.first*factor;
    }else{
      dis = 1e9;
      point_index = -1;
      if (pcloud2!=0) {
	auto results = pcloud2->get_closest_2d_dis_index(p2d.first, p2d.second, 0);
	dis = results.first;
	point_index = results.second;
      }
      if (dis < 0.6*units::cm) {
	Point test_p(pcloud2->get_cloud().pts[point_index].x, pcloud2->get_cloud().pts[point_index].y, pcloud2->get_cloud().pts[point_index].z);
	double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
	sum_u_charge += it->second.first*factor;
      }
    }
  }

  for (auto it = charge_2d_v.begin(); it!= charge_2d_v.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 1);
    double dis = 1e9;
    int point_index = -1;
    if (pcloud1!=0) {
      auto results = pcloud1->get_closest_2d_dis_index(p2d.first, p2d.second, 1);
      dis = results.first;
      point_index = results.second;
    }
    if (dis < 0.6*units::cm) {
      Point test_p(pcloud1->get_cloud().pts[point_index].x, pcloud1->get_cloud().pts[point_index].y, pcloud1->get_cloud().pts[point_index].z);
      double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
      sum_v_charge += it->second.first*factor;
    }else{
      dis = 1e9;
      point_index = -1;
      if (pcloud2!=0) {
	auto results = pcloud2->get_closest_2d_dis_index(p2d.first, p2d.second, 1);
	dis = results.first;
	point_index = results.second;
      }
      if (dis < 0.6*units::cm) {
	Point test_p(pcloud2->get_cloud().pts[point_index].x, pcloud2->get_cloud().pts[point_index].y, pcloud2->get_cloud().pts[point_index].z);
	double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
	sum_v_charge += it->second.first * factor;
      }
    }
  }

  
  for (auto it = charge_2d_w.begin(); it!= charge_2d_w.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 2);
    double dis = 1e9;
    int point_index = -1;
    if (pcloud1!=0) {
      auto results = pcloud1->get_closest_2d_dis_index(p2d.first, p2d.second, 2);
      dis = results.first;
      point_index = results.second;
    }
    if (dis < 0.6*units::cm) {
      Point test_p(pcloud1->get_cloud().pts[point_index].x, pcloud1->get_cloud().pts[point_index].y, pcloud1->get_cloud().pts[point_index].z);
      double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
      sum_w_charge += it->second.first*factor;
    }else{
      dis = 1e9;
      point_index = -1;
      if (pcloud2!=0) {
       auto results = pcloud2->get_closest_2d_dis_index(p2d.first, p2d.second, 2);
       dis = results.first;
       point_index = results.second;
      }
      if (dis < 0.6*units::cm) {
	Point test_p(pcloud2->get_cloud().pts[point_index].x, pcloud2->get_cloud().pts[point_index].y, pcloud2->get_cloud().pts[point_index].z);
	double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
	sum_w_charge += it->second.first*factor;
      }
    }
  }

  
  //  std::cout << sum_u_charge << " " << sum_v_charge << " " << sum_w_charge << std::endl;
  double charge[3];
  charge[0] = sum_u_charge;
  charge[1] = sum_v_charge;
  charge[2] = sum_w_charge;
  double weight[3]={0.25,0.25,1};

  int min_index = 0, max_index = 0, med_index = 0;
  double min_charge = 1e9, max_charge = -1e9;
  for (int i=0;i!=3;i++){
    if (min_charge > charge[i]){
      min_charge = charge[i];
      min_index = i;
    }
    if (max_charge < charge[i]){
      max_charge = charge[i];
      max_index = i;
    }
  }
  if (min_index != max_index){
    for (int i=0;i!=3;i++){
      if (i==min_index) continue;
      if (i==max_index) continue;
      med_index = i;
    }
  }else{
    min_index = 0;
    med_index = 1;
    max_index = 2;
  }
  //  std::cout << min_index << " " << med_index << " " << max_index << std::endl;

  double min_asy = fabs(charge[med_index] - charge[min_index])/(charge[med_index] + charge[min_index]);
  double max_asy = fabs(charge[med_index] - charge[max_index])/(charge[med_index] + charge[max_index]);
  
  //  std::cout << med_index << " " << min_index << " " << max_index << " " << charge[0] << " " << charge[1] << " " << charge[2] << " " << min_asy << " " << max_asy << " " << (0.25 * sum_u_charge + 0.25*sum_v_charge + sum_w_charge)/1.5 << std::endl;

  // default case ...
  double overall_charge = (weight[0]*charge[0] + weight[1]*charge[1]+weight[2]*charge[2])/(weight[0] + weight[1] + weight[2]);
  
  if (max_asy > 0.04){ // exclude the maximal charge ... 
    // if (min_asy < 0.04){
    overall_charge =(weight[med_index] * charge[med_index] + weight[min_index]*charge[min_index])/(weight[med_index] + weight[min_index]);
    // }else{
    // overall_charge = charge[med_index];
    // }
  }
    

  
  
  // default ...
  kine_energy = overall_charge/recom_factor/fudge_factor*23.6/1e6 * units::MeV;
  
  
  return kine_energy;
}


double WCPPID::NeutrinoID::cal_corr_factor(WCP::Point& p, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw){
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double factor = 1;
  double central_U = offset_u + (slope_yu * p.y + slope_zu * p.z);
  if (central_U >=296 && central_U <=327 ||
      central_U >=336 && central_U <=337 ||
      central_U >=343 && central_U <=351 ||
      central_U >=376 && central_U <=400 ||
      central_U >=410 && central_U <=484 ||
      central_U >=501 && central_U <=524 ||
      central_U >=536 && central_U <=671)
    factor = factor/0.7;
  
  if (mp.get_flag_corr()){
    factor *= mp.get_corr_factor(p, offset_u,  slope_yu,  slope_zu,  offset_v,  slope_yv,  slope_zv,  offset_w,  slope_yw,  slope_zw);
  }
  return factor;
}


double WCPPID::NeutrinoID::cal_kine_charge(WCPPID::ProtoSegment *sg){
   TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_offset = mp.get_time_offset();
  int nrebin = mp.get_nrebin();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w + 0.5;
  double offset_u = -first_u_dis/pitch_u + 0.5;
  double offset_v = -first_v_dis/pitch_v + 0.5;
  
  double first_t_dis = main_cluster->get_point_cloud()->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - main_cluster->get_point_cloud()->get_cloud().pts[0].x;
  double slope_xt = 1./time_slice_width;
  double offset_t =  first_t_dis/time_slice_width + 0.5;

  
  double kine_energy = 0;
  
  // to be improved ...
  double fudge_factor = 0.95;
  double recom_factor  = 0.7;
  
  if (sg->get_flag_shower_topology()) {
    recom_factor = 0.5; // assume shower
    fudge_factor = 0.8; // shower ...
  }else if (fabs(sg->get_particle_type())==2212){
    recom_factor = 0.35;
  }

  double sum_u_charge = 0;
  double sum_v_charge = 0;
  double sum_w_charge = 0;

  WCP::ToyPointCloud* pcloud1 = sg->get_associated_pcloud();
  WCP::ToyPointCloud* pcloud2 = sg->get_fit_pcloud();

  if (pcloud1 == 0 && pcloud2 == 0) return 0;
  else if (pcloud1 ==0 && pcloud2 !=0) pcloud1 = pcloud2;
  else if (pcloud2 ==0 && pcloud1 !=0) pcloud2 = pcloud1;
  
  //  std::cout << sg->get_cluster_id() << " " << sg->get_id() << " " << pcloud1 << " " << pcloud2 << " " << pcloud1->get_num_points() << " " << pcloud2->get_num_points() << std::endl;
  
  for (auto it = charge_2d_u.begin(); it!= charge_2d_u.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 0);

    auto results = pcloud1->get_closest_2d_dis_index(p2d.first, p2d.second, 0);
    double dis = results.first;
    int point_index = results.second;
    if (dis < 0.6*units::cm) {
      Point test_p(pcloud1->get_cloud().pts[point_index].x, pcloud1->get_cloud().pts[point_index].y, pcloud1->get_cloud().pts[point_index].z);
      double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
      sum_u_charge += it->second.first*factor;
    }else{
      results = pcloud2->get_closest_2d_dis_index(p2d.first, p2d.second, 0);
      dis = results.first;
      point_index = results.second;
      if (dis < 0.6*units::cm) {
	Point test_p(pcloud2->get_cloud().pts[point_index].x, pcloud2->get_cloud().pts[point_index].y, pcloud2->get_cloud().pts[point_index].z);
	double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
	sum_u_charge += it->second.first*factor;
      }
    }
  }

  for (auto it = charge_2d_v.begin(); it!= charge_2d_v.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 1);
    auto results = pcloud1->get_closest_2d_dis_index(p2d.first, p2d.second, 1);
    double dis = results.first;
    int point_index = results.second;
    if (dis < 0.6*units::cm) {
      Point test_p(pcloud1->get_cloud().pts[point_index].x, pcloud1->get_cloud().pts[point_index].y, pcloud1->get_cloud().pts[point_index].z);
      double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
      sum_v_charge += it->second.first * factor;
    }else{
      results = pcloud2->get_closest_2d_dis_index(p2d.first, p2d.second, 1);
      dis = results.first;
      point_index = results.second;
      if (dis < 0.6*units::cm) {
	Point test_p(pcloud2->get_cloud().pts[point_index].x, pcloud2->get_cloud().pts[point_index].y, pcloud2->get_cloud().pts[point_index].z);
	double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
	sum_v_charge += it->second.first*factor;
      }
    }
  }

  
  for (auto it = charge_2d_w.begin(); it!= charge_2d_w.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 2);
    auto results = pcloud1->get_closest_2d_dis_index(p2d.first, p2d.second, 2);
    double dis = results.first;
    int point_index = results.second;
    if (dis < 0.6*units::cm) {
      Point test_p(pcloud1->get_cloud().pts[point_index].x, pcloud1->get_cloud().pts[point_index].y, pcloud1->get_cloud().pts[point_index].z);
      double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
      sum_w_charge += it->second.first * factor;
    }else{
      results = pcloud2->get_closest_2d_dis_index(p2d.first, p2d.second, 2);
      dis = results.first;
      point_index = results.second;
      if (dis < 0.6*units::cm) {
	Point test_p(pcloud2->get_cloud().pts[point_index].x, pcloud2->get_cloud().pts[point_index].y, pcloud2->get_cloud().pts[point_index].z);
	double factor = cal_corr_factor(test_p, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
	sum_w_charge += it->second.first*factor;
      }
    }
  }


  //std::cout << sum_u_charge << " " << sum_v_charge << " " << sum_w_charge << std::endl;

   double charge[3];
  charge[0] = sum_u_charge;
  charge[1] = sum_v_charge;
  charge[2] = sum_w_charge;
  double weight[3]={0.25,0.25,1};

  int min_index = 0, max_index = 0, med_index = 0;
  double min_charge = 1e9, max_charge = -1e9;
  for (int i=0;i!=3;i++){
    if (min_charge > charge[i]){
      min_charge = charge[i];
      min_index = i;
    }
    if (max_charge < charge[i]){
      max_charge = charge[i];
      max_index = i;
    }
  }
  if (min_index != max_index){
    for (int i=0;i!=3;i++){
      if (i==min_index) continue;
      if (i==max_index) continue;
      med_index = i;
    }
  }else{
    min_index = 0;
    med_index = 1;
    max_index = 2;
  }
  //  std::cout << min_index << " " << med_index << " " << max_index << std::endl;

  double min_asy = fabs(charge[med_index] - charge[min_index])/(charge[med_index] + charge[min_index]);
  double max_asy = fabs(charge[med_index] - charge[max_index])/(charge[med_index] + charge[max_index]);
  
  //  std::cout << med_index << " " << min_index << " " << max_index << " " << charge[0] << " " << charge[1] << " " << charge[2] << " " << min_asy << " " << max_asy << " " << (0.25 * sum_u_charge + 0.25*sum_v_charge + sum_w_charge)/1.5 << std::endl;

  // default case ...
  double overall_charge = (weight[0]*charge[0] + weight[1]*charge[1]+weight[2]*charge[2])/(weight[0] + weight[1] + weight[2]);
  
  if (max_asy > 0.04){ // exclude the maximal charge ... 
    // if (min_asy < 0.04){
    overall_charge =(weight[med_index] * charge[med_index] + weight[min_index]*charge[min_index])/(weight[med_index] + weight[min_index]);
    // }else{
    // overall_charge = charge[med_index];
    // }
  }


  
  kine_energy = overall_charge/recom_factor/fudge_factor*23.6/1e6 * units::MeV;
  

  
  return kine_energy;
}
