std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> WireCellPID::PR3DCluster::get_highest_lowest_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  for (size_t i=1;i<cloud.pts.size();i++){
    if (cloud.pts[i].y > highest_wcp.y)
      highest_wcp = cloud.pts[i];
    if (cloud.pts[i].y < lowest_wcp.y)
      lowest_wcp = cloud.pts[i];
  }
  return std::make_pair(highest_wcp,lowest_wcp);
}

std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> WireCellPID::PR3DCluster::get_front_back_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  for (size_t i=1;i<cloud.pts.size();i++){
    if (cloud.pts[i].z > highest_wcp.z)
      highest_wcp = cloud.pts[i];
    if (cloud.pts[i].z < lowest_wcp.z)
      lowest_wcp = cloud.pts[i];
  }
  return std::make_pair(highest_wcp,lowest_wcp);
}


std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> WireCellPID::PR3DCluster::get_earliest_latest_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  for (size_t i=1;i<cloud.pts.size();i++){
    if (cloud.pts[i].x > highest_wcp.x)
      highest_wcp = cloud.pts[i];
    if (cloud.pts[i].x < lowest_wcp.x)
      lowest_wcp = cloud.pts[i];
  }
  return std::make_pair(lowest_wcp,highest_wcp);
}

std::vector<std::vector<WCPointCloud<double>::WCPoint>> WireCellPID::PR3DCluster::get_extreme_wcps(){
  
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  Calc_PCA();
  WCPointCloud<double>::WCPoint wcps[8];
  for (int i=0;i!=8;i++){
    wcps[i] = cloud.pts[0];
  }
  Vector main_axis = get_PCA_axis(0);
  if (main_axis.y <0){
    main_axis.x = -main_axis.x;
    main_axis.y = -main_axis.y;
    main_axis.z = -main_axis.z;
  }
  double high_value = wcps[0].x*main_axis.x + wcps[0].y*main_axis.y + wcps[0].z*main_axis.z;
  double low_value = wcps[1].x * main_axis.x + wcps[1].y * main_axis.y + wcps[1].z * main_axis.z ;

  for (size_t i=1;i<cloud.pts.size();i++){
    double value = cloud.pts[i].x*main_axis.x + cloud.pts[i].y*main_axis.y + cloud.pts[i].z*main_axis.z;
    if (value > high_value){
      wcps[0] = cloud.pts[i];
      high_value = value;
    }
    if (value < low_value){
      wcps[1] = cloud.pts[i];
      low_value = value;
    }
    // top down
    if (cloud.pts[i].y > wcps[2].y)
      wcps[2] = cloud.pts[i];
    if (cloud.pts[i].y < wcps[3].y)
      wcps[3] = cloud.pts[i];
    
    // front back
    if (cloud.pts[i].z > wcps[4].z)
      wcps[4] = cloud.pts[i];
    if (cloud.pts[i].z < wcps[5].z)
      wcps[5] = cloud.pts[i];
    
    // early late
    if (cloud.pts[i].x > wcps[6].x)
      wcps[6] = cloud.pts[i];
    if (cloud.pts[i].x < wcps[7].x)
      wcps[7] = cloud.pts[i];
  }

  
  

  
  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps;

  {
    // first extreme along the main axis
    std::vector<WCPointCloud<double>::WCPoint> saved_wcps;
    saved_wcps.push_back(wcps[0]);
    out_vec_wcps.push_back(saved_wcps);
  }

  {
    // second extreme along the main axis
    std::vector<WCPointCloud<double>::WCPoint> saved_wcps;
    saved_wcps.push_back(wcps[1]);
    out_vec_wcps.push_back(saved_wcps);
  }
  
  // std::cout << std::endl;
  for (int i=2;i!=8;i++){
    //if (cluster_id==16)
    //std::cout << i << " C " << wcps[i].x/units::cm << " " << wcps[i].y/units::cm << " " << wcps[i].z/units::cm << std::endl;
    
    bool flag_save = true;
    for (size_t j=0;j!=out_vec_wcps.size(); j++){
      double dis = sqrt(pow(out_vec_wcps[j].at(0).x-wcps[i].x,2) + pow(out_vec_wcps[j].at(0).y - wcps[i].y,2) + pow(out_vec_wcps[j].at(0).z - wcps[i].z,2));
      if (dis < 5*units::cm){
	out_vec_wcps.at(j).push_back(wcps[i]);
	flag_save = false;
	break;
      }
    }
    
    if (flag_save){
      std::vector<WCPointCloud<double>::WCPoint> saved_wcps;
      saved_wcps.push_back(wcps[i]);
      out_vec_wcps.push_back(saved_wcps);
    }
  }

  return out_vec_wcps;  
}

std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> WireCellPID::PR3DCluster::get_main_axis_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  Calc_PCA();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  Vector main_axis = get_PCA_axis(0);
  if (main_axis.y <0){
    main_axis.x = -main_axis.x;
    main_axis.y = -main_axis.y;
    main_axis.z = -main_axis.z;
  }
  
  double high_value = highest_wcp.x*main_axis.x + highest_wcp.y*main_axis.y + highest_wcp.z*main_axis.z;
  double low_value = lowest_wcp.x * main_axis.x + lowest_wcp.y * main_axis.y + lowest_wcp.z * main_axis.z ;
  
  for (size_t i=1;i<cloud.pts.size();i++){
    double value = cloud.pts[i].x*main_axis.x + cloud.pts[i].y*main_axis.y + cloud.pts[i].z*main_axis.z;
    if (value > high_value){
      highest_wcp = cloud.pts[i];
      high_value = value;
    }
    
    if (value < low_value){
      lowest_wcp = cloud.pts[i];
      low_value = value;
    }
  }
  return std::make_pair(highest_wcp,lowest_wcp);
}

std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> WireCellPID::PR3DCluster::get_two_boundary_wcps(){
  Create_point_cloud();
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint extreme_wcp[14];
  
  
  Calc_PCA();
  Vector main_axis = get_PCA_axis(0);
  Vector second_axis = get_PCA_axis(1);
  double values[4];
  
  
  bool flag_init = false;
  for (size_t i=0;i<cloud.pts.size();i++){
    if (cloud.pts[i].mcell->Estimate_total_charge() < 1500) continue;
    if (!flag_init){
      flag_init = true;
      for (int j=0;j!=14;j++){
	extreme_wcp[j] = cloud.pts[i];
      }
      values[0] = extreme_wcp[0].x*main_axis.x + extreme_wcp[0].y*main_axis.y + extreme_wcp[0].z*main_axis.z;
      values[1] = extreme_wcp[1].x*main_axis.x + extreme_wcp[1].y*main_axis.y + extreme_wcp[1].z*main_axis.z;
      values[2] = extreme_wcp[2].x*second_axis.x + extreme_wcp[2].y*second_axis.y + extreme_wcp[2].z*second_axis.z;
      values[3] = extreme_wcp[3].x*second_axis.x + extreme_wcp[3].y*second_axis.y + extreme_wcp[3].z*second_axis.z;
    }else{
      // main axis points
      double value = cloud.pts[i].x*main_axis.x + cloud.pts[i].y*main_axis.y + cloud.pts[i].z*main_axis.z;
      if (value > values[0]){
	values[0] = value;
	extreme_wcp[0] = cloud.pts[i];
      }
      if (value < values[1]){
	values[1] = value;
	extreme_wcp[1] = cloud.pts[i];
      }
      // second axis points
      value = cloud.pts[i].x*second_axis.x + cloud.pts[i].y*second_axis.y + cloud.pts[i].z*second_axis.z;
      if (value > values[2]){
	values[2] = value;
	extreme_wcp[2] = cloud.pts[i];
      }
      if (value < values[3]){
	values[3] = value;
	extreme_wcp[3] = cloud.pts[i];
      }
      // early-late  x
      if (cloud.pts[i].x > extreme_wcp[4].x){
	extreme_wcp[4] = cloud.pts[i];
      }
      if (cloud.pts[i].x < extreme_wcp[5].x){
	extreme_wcp[5] = cloud.pts[i];
      }
      // top-bottom  y 
      if (cloud.pts[i].y > extreme_wcp[6].y){
	extreme_wcp[6] = cloud.pts[i];
      }
      if (cloud.pts[i].y < extreme_wcp[7].y){
	extreme_wcp[7] = cloud.pts[i];
      }
      // left-right  z or w
      if (cloud.pts[i].z > extreme_wcp[8].z){
	extreme_wcp[8] = cloud.pts[i];
      }
      if (cloud.pts[i].z < extreme_wcp[9].z){
	extreme_wcp[9] = cloud.pts[i];
      }
      // U range
      if (cloud.pts[i].index_u > extreme_wcp[10].index_u){
	extreme_wcp[10] = cloud.pts[i];
      }
      if (cloud.pts[i].index_u < extreme_wcp[11].index_u){
	extreme_wcp[11] = cloud.pts[i];
      }
      // V range
      if (cloud.pts[i].index_v > extreme_wcp[12].index_v){
	extreme_wcp[12] = cloud.pts[i];
      }
      if (cloud.pts[i].index_v < extreme_wcp[13].index_v){
	extreme_wcp[13] = cloud.pts[i];
      }
    }
  }
  
  if (!flag_init){
    flag_init = true;
    for (int i=0;i!=14;i++){
      extreme_wcp[i] = cloud.pts[0];
    }
    extreme_wcp[1] = cloud.pts[1];
  }

  
  std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> boundary_points;
  boundary_points.first = extreme_wcp[0];
  boundary_points.second = extreme_wcp[1];
  double boundary_value = fabs(boundary_points.first.x-boundary_points.second.x)/(2.22*units::mm)
    + fabs(boundary_points.first.index_u - boundary_points.second.index_u)
    + fabs(boundary_points.first.index_v - boundary_points.second.index_v)
    + fabs(boundary_points.first.index_w - boundary_points.second.index_w);
  
  for (int i=0;i<14;i++){
    for (int j=i+1;j<14;j++){
      double value = fabs(extreme_wcp[i].x - extreme_wcp[j].x)/(2.22*units::mm)
	+ fabs(extreme_wcp[i].index_u - extreme_wcp[j].index_u)
	+ fabs(extreme_wcp[i].index_v - extreme_wcp[j].index_v)
	+ fabs(extreme_wcp[i].index_w - extreme_wcp[j].index_w);
      if (value > boundary_value){
	boundary_value = value;
	if (extreme_wcp[i].y > extreme_wcp[j].y){
	  boundary_points.first = extreme_wcp[i];
	  boundary_points.second = extreme_wcp[j];
	}else{
	  boundary_points.first = extreme_wcp[j];
	  boundary_points.second = extreme_wcp[i];
	}
      }
    }
  }
  return boundary_points;
}


std::pair<Point,Point> WireCellPID::PR3DCluster::get_two_extreme_points(){
  Create_point_cloud();
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint extreme_wcp[6];
  for (int i=0;i!=6;i++){
    extreme_wcp[i] = cloud.pts[0];
  }
  for (size_t i=1;i< cloud.pts.size();i++){
    if (cloud.pts[i].y > extreme_wcp[0].y)
      extreme_wcp[0] = cloud.pts[i];
    if (cloud.pts[i].y < extreme_wcp[1].y)
      extreme_wcp[1] = cloud.pts[i];
    
    if (cloud.pts[i].x > extreme_wcp[2].x)
      extreme_wcp[2] = cloud.pts[i];
    if (cloud.pts[i].x < extreme_wcp[3].x)
      extreme_wcp[3] = cloud.pts[i];
    
    if (cloud.pts[i].z > extreme_wcp[4].z)
      extreme_wcp[4] = cloud.pts[i];
    if (cloud.pts[i].z < extreme_wcp[5].z)
      extreme_wcp[5] = cloud.pts[i];
  }

  double max_dis = -1;
  WCPointCloud<double>::WCPoint wcp1, wcp2;
  for (int i=0;i!=6;i++){
    for (int j=i+1;j!=6;j++){
      double dis = sqrt(pow(extreme_wcp[i].x - extreme_wcp[j].x,2)+pow(extreme_wcp[i].y - extreme_wcp[j].y,2)+pow(extreme_wcp[i].z - extreme_wcp[j].z,2));
      if (dis > max_dis){
	max_dis = dis;
	wcp1 = extreme_wcp[i];
	wcp2 = extreme_wcp[j];
      }
    }
  }
  Point p1(wcp1.x,wcp1.y,wcp1.z);
  Point p2(wcp2.x,wcp2.y,wcp2.z);
  p1 = calc_ave_pos(p1,5*units::cm);
  p2 = calc_ave_pos(p2,5*units::cm);
  
  return std::make_pair(p1,p2);
}


void WireCellPID::PR3DCluster::dijkstra_shortest_paths(WCPointCloud<double>::WCPoint& wcp){
  if (graph==(MCUGraph*)0)
    Create_graph();
  if (wcp.index==source_wcp_index)
    return;
  source_wcp_index = wcp.index;
  parents.resize(num_vertices(*graph));
  distances.resize(num_vertices(*graph));
  
  auto v0 = vertex(wcp.index,*graph);
  boost::dijkstra_shortest_paths(*graph, v0,
				 weight_map(get(edge_weight, *graph))
				 .predecessor_map(&parents[0])
				 .distance_map(&distances[0])
				 );
}


void WireCellPID::PR3DCluster::cal_shortest_path(WCPointCloud<double>::WCPoint& wcp_target){
  dest_wcp_index = wcp_target.index;
  path_wcps.clear();
  path_mcells.clear();

  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  int prev_i = -1;
  for(int i = dest_wcp_index; i!=source_wcp_index; i = parents[i]) {
    if (path_wcps.size()==0){
      path_wcps.push_front(cloud.pts[i]);
      path_mcells.push_front(cloud.pts[i].mcell);
    }else{
      path_wcps.push_front(cloud.pts[i]);
      if (cloud.pts[i].mcell!=path_mcells.front())
	path_mcells.push_front(cloud.pts[i].mcell);
    }
    if (i==prev_i) break;
    prev_i = i;
  }
  path_wcps.push_front(cloud.pts[source_wcp_index]);
  if (cloud.pts[source_wcp_index].mcell!=path_mcells.front())
    path_mcells.push_front(cloud.pts[source_wcp_index].mcell);
  
}


