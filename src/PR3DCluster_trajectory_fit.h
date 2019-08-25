void WireCellPID::PR3DCluster::prepare_data(WireCell::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge){
  
  std::vector<int> proj_channel;
  std::vector<int> proj_timeslice;
  std::vector<int> proj_charge;
  std::vector<int> proj_charge_err;
  std::vector<int> proj_flag;
  get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);

  int min_time = 1e9;
  int max_time = -1e9;
  int min_uch = 1e9;
  int max_uch = -1e9;
  int min_vch = 1e9;
  int max_vch = -1e9;
  int min_wch = 1e9;
  int max_wch = -1e9;
  
  // std::cout << proj_charge.size() << " " << proj_flag.size() << std::endl;
  for (size_t i=0;i!=proj_channel.size();i++){
    if (min_time > proj_timeslice.at(i)) min_time = proj_timeslice.at(i);
    if (max_time < proj_timeslice.at(i)) max_time = proj_timeslice.at(i);
	
    if (proj_channel.at(i)<2400){
      map_2D_ut_charge[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));

      if (min_uch > proj_channel.at(i)) min_uch = proj_channel.at(i);
      if (max_uch < proj_channel.at(i)) max_uch = proj_channel.at(i);
    }else if (proj_channel.at(i)<4800){
      map_2D_vt_charge[std::make_pair(proj_channel.at(i)-2400,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));

      if (min_vch > proj_channel.at(i)) min_vch = proj_channel.at(i);
      if (max_vch < proj_channel.at(i)) max_vch = proj_channel.at(i);
    }else{
      map_2D_wt_charge[std::make_pair(proj_channel.at(i)-4800,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));
      if (min_wch > proj_channel.at(i)) min_wch = proj_channel.at(i);
      if (max_wch < proj_channel.at(i)) max_wch = proj_channel.at(i);
    }
  }

  // flag 0 for dead or overlapping channels
  // flag 1 for good channels
  // flag 2 for isolated channels
  // flag 3 for additional channels ...
  
  // add the rest of live channels within range???
  std::map<std::pair<int,int>, std::pair<double,double> > map_u_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_uch, max_uch, 0);
  std::map<std::pair<int,int>, std::pair<double,double> > map_v_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_vch, max_vch, 1);
  std::map<std::pair<int,int>, std::pair<double,double> > map_w_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_wch, max_wch, 2);

  for (auto it = map_u_tcc.begin(); it!=map_u_tcc.end(); it++){
    if (map_2D_ut_charge.find(std::make_pair(it->first.second, it->first.first))==map_2D_ut_charge.end()){
      map_2D_ut_charge[std::make_pair(it->first.second, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }

  for (auto it = map_v_tcc.begin(); it!=map_v_tcc.end(); it++){
    if (map_2D_vt_charge.find(std::make_pair(it->first.second-2400, it->first.first))==map_2D_vt_charge.end()){
      map_2D_vt_charge[std::make_pair(it->first.second-2400, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }

  for (auto it = map_w_tcc.begin(); it!=map_w_tcc.end(); it++){
    if (map_2D_wt_charge.find(std::make_pair(it->first.second-4800, it->first.first))==map_2D_wt_charge.end()){
      map_2D_wt_charge[std::make_pair(it->first.second, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //      std::cout << it->first.first << " " << it->first.second << std::endl;
    }else{
      //      std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }
  
  
  //  for (auto it = map_w_tcc.begin(); it!=map_w_tcc.end(); it++){
  // std::cout << it->first.first << " " << it->first.second << std::endl;
  // }
  
}

WireCell::PointVector WireCellPID::PR3DCluster::organize_wcps_path(std::list<WCPointCloud<double>::WCPoint>& path_wcps_list,  double low_dis_limit){

  PointVector pts;
  
  std::vector<WCPointCloud<double>::WCPoint> temp_wcps_vec(path_wcps_list.begin(), path_wcps_list.end());

  // fill in the beginning point ...
  {
    Point p1(temp_wcps_vec.front().x, temp_wcps_vec.front().y, temp_wcps_vec.front().z);
    Point p2(temp_wcps_vec.front().x, temp_wcps_vec.front().y, temp_wcps_vec.front().z);
    double dis1 = 0;
    for (auto it = temp_wcps_vec.begin(); it!=temp_wcps_vec.end(); it++){
      p2.x = (*it).x;
      p2.y = (*it).y;
      p2.z = (*it).z;
      dis1 = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis1 > low_dis_limit) break;
    }
    if (dis1!=0){
      p1.x += (p1.x - p2.x)/dis1 * low_dis_limit/2.;
      p1.y += (p1.y - p2.y)/dis1 * low_dis_limit/2.;
      p1.z += (p1.z - p2.z)/dis1 * low_dis_limit/2.;
      pts.push_back(p1);
    }
  }

  // fill in the middle part
  for (size_t i=0;i!=temp_wcps_vec.size(); i++){
    Point p1(temp_wcps_vec.at(i).x, temp_wcps_vec.at(i).y, temp_wcps_vec.at(i).z);
    if (i==0) {
      pts.push_back(p1);
    }else{
      double dis = sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2));
    
      if (dis < low_dis_limit * 0.6 ){
	continue;
      }else if (dis < low_dis_limit * 1.6){
	pts.push_back(p1);
      }else{
	int npoints = std::round(dis/low_dis_limit);

	for (int j=0;j!=npoints;j++){
	  Point p(pts.back().x + (p1.x-pts.back().x) / npoints * (j+1),
		  pts.back().y + (p1.y-pts.back().y) / npoints * (j+1),
		  pts.back().z + (p1.z-pts.back().z) / npoints * (j+1));
	  pts.push_back(p);
	}
      }

    }
  }
  

  // fill in the end part
  {
    Point p1(temp_wcps_vec.back().x, temp_wcps_vec.back().y, temp_wcps_vec.back().z);
    Point p2(temp_wcps_vec.back().x, temp_wcps_vec.back().y, temp_wcps_vec.back().z);
    double dis1 = 0;
    for (auto it = temp_wcps_vec.rbegin(); it!=temp_wcps_vec.rend(); it++){
      p2.x = (*it).x;
      p2.y = (*it).y;
      p2.z = (*it).z;
      dis1 = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis1 > low_dis_limit) break;
    }
    if (dis1!=0){
      p1.x += (p1.x - p2.x)/dis1 * low_dis_limit/2.;
      p1.y += (p1.y - p2.y)/dis1 * low_dis_limit/2.;
      p1.z += (p1.z - p2.z)/dis1 * low_dis_limit/2.;
      pts.push_back(p1);
    }
  }

    
  
  
  return pts;
}
