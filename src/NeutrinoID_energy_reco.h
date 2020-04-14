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

  double sum_u_charge = 0;
  double sum_v_charge = 0;
  double sum_w_charge = 0;

  shower->build_point_clouds();
  WCP::ToyPointCloud* pcloud1 = shower->get_associated_pcloud();
  WCP::ToyPointCloud* pcloud2 = shower->get_fit_pcloud();
  
  
  for (auto it = charge_2d_u.begin(); it!= charge_2d_u.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 0);
    double dis = 1e9;
    if (pcloud1!=0) dis = pcloud1->get_closest_2d_dis(p2d.first, p2d.second, 0);
    if (dis < 0.6*units::cm) sum_u_charge += it->second.first;
    else{
      dis = 1e9;
      if (pcloud2!=0) dis = pcloud2->get_closest_2d_dis(p2d.first, p2d.second, 0);
      if (dis < 0.6*units::cm) sum_u_charge += it->second.first;
    }
  }

  for (auto it = charge_2d_v.begin(); it!= charge_2d_v.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 1);
    double dis = 1e9;
    if (pcloud1!=0) dis = pcloud1->get_closest_2d_dis(p2d.first, p2d.second, 1);
    if (dis < 0.6*units::cm) sum_v_charge += it->second.first;
    else{
      dis = 1e9;
      if (pcloud2!=0) dis = pcloud2->get_closest_2d_dis(p2d.first, p2d.second, 1);
      if (dis < 0.6*units::cm) sum_v_charge += it->second.first;
    }
  }

  
  for (auto it = charge_2d_w.begin(); it!= charge_2d_w.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 2);
    double dis = 1e9;
    if (pcloud1!=0) dis = pcloud1->get_closest_2d_dis(p2d.first, p2d.second, 2);
    if (dis < 0.6*units::cm) sum_w_charge += it->second.first;
    else{
      dis = 1e9;
      if (pcloud2!=0) dis = pcloud2->get_closest_2d_dis(p2d.first, p2d.second, 2);
      if (dis < 0.6*units::cm) sum_w_charge += it->second.first;
    }
  }


  //  std::cout << sum_u_charge << " " << sum_v_charge << " " << sum_w_charge << std::endl;

  kine_energy = (0.25 * sum_u_charge + 0.25*sum_v_charge + sum_w_charge)/1.5/recom_factor/fudge_factor*23.6/1e6 * units::MeV;

  
  return kine_energy;
}

double WCPPID::NeutrinoID::cal_kine_charge(WCPPID::ProtoSegment *sg){
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
  
  for (auto it = charge_2d_u.begin(); it!= charge_2d_u.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 0);
    double dis = pcloud1->get_closest_2d_dis(p2d.first, p2d.second, 0);
    if (dis < 0.6*units::cm) sum_u_charge += it->second.first;
    else{
      dis = pcloud2->get_closest_2d_dis(p2d.first, p2d.second, 0);
      if (dis < 0.6*units::cm) sum_u_charge += it->second.first;
    }
  }

  for (auto it = charge_2d_v.begin(); it!= charge_2d_v.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 1);
    double dis = pcloud1->get_closest_2d_dis(p2d.first, p2d.second, 1);
    if (dis < 0.6*units::cm) sum_v_charge += it->second.first;
    else{
      dis = pcloud2->get_closest_2d_dis(p2d.first, p2d.second, 1);
      if (dis < 0.6*units::cm) sum_v_charge += it->second.first;
    }
  }

  
  for (auto it = charge_2d_w.begin(); it!= charge_2d_w.end(); it++){
    std::pair<double, double> p2d = ct_point_cloud->convert_time_ch_2Dpoint(it->first.first, it->first.second, 2);
    double dis = pcloud1->get_closest_2d_dis(p2d.first, p2d.second, 2);
    if (dis < 0.6*units::cm) sum_w_charge += it->second.first;
    else{
      dis = pcloud2->get_closest_2d_dis(p2d.first, p2d.second, 2);
      if (dis < 0.6*units::cm) sum_w_charge += it->second.first;
    }
  }


  //std::cout << sum_u_charge << " " << sum_v_charge << " " << sum_w_charge << std::endl;

  kine_energy = (0.25 * sum_u_charge + 0.25*sum_v_charge + sum_w_charge)/1.5/recom_factor/fudge_factor*23.6/1e6 * units::MeV;
  

  
  return kine_energy;
}
