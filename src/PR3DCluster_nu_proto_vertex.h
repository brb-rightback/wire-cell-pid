bool WCPPID::PR3DCluster::proto_break_tracks(WCP::WCPointCloud<double>::WCPoint& first_wcp, WCP::WCPointCloud<double>::WCPoint& curr_wcp, WCP::WCPointCloud<double>::WCPoint& last_wcp, std::list<WCP::WCPointCloud<double>::WCPoint>& wcps_list1, std::list<WCP::WCPointCloud<double>::WCPoint>& wcps_list2){

  double dis1 = sqrt(pow(curr_wcp.x-first_wcp.x,2) + pow(curr_wcp.y-first_wcp.y,2) + pow(curr_wcp.z-first_wcp.z,2));
  double dis2 = sqrt(pow(curr_wcp.x-last_wcp.x,2) + pow(curr_wcp.y-last_wcp.y,2) + pow(curr_wcp.z-last_wcp.z,2));

  if (dis1 > 1*units::cm && dis2 > 1*units::cm){
    dijkstra_shortest_paths(first_wcp,2);
    cal_shortest_path(curr_wcp,2);

    wcps_list1 = path_wcps;

    dijkstra_shortest_paths(curr_wcp,2);
    cal_shortest_path(last_wcp,2);
    wcps_list2 = path_wcps;

    
    auto it1 = wcps_list1.rbegin();
    int count = 0;
    for (auto it = wcps_list2.begin(); it!=wcps_list2.end(); it++){
      if ( (*it).index==(*it1).index){
	count ++;
	it1++;
      }
    }
    for (int i=0;i!=count;i++){
      if (i!=count-1){
	wcps_list1.pop_back();
	wcps_list2.pop_front();
      }
    }
    
    
    return true;
  }else {
    return false;
  }
}



WCP::WCPointCloud<double>::WCPoint WCPPID::PR3DCluster::proto_extend_point(WCP::Point& p, TVector3& dir){

  float step_dis = 1*units::cm;
  
  WCP::WCPointCloud<double>::WCPoint curr_wcp = point_cloud_steiner->get_closest_wcpoint(p);
  WCP::WCPointCloud<double>::WCPoint next_wcp = curr_wcp;

  WCP::Point test_p;
  
  bool flag_continue = true;
  while(flag_continue){
    flag_continue = false;
    
    for (int i=0; i!=3; i++){
      test_p.x = curr_wcp.x + dir.X() * step_dis * (i+1);
      test_p.y = curr_wcp.y + dir.Y() * step_dis * (i+1);
      test_p.z = curr_wcp.z + dir.Z() * step_dis * (i+1);
      
      next_wcp = point_cloud->get_closest_wcpoint(test_p);
      TVector3 dir1(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z);
      
      /* std::cout << dir.X() << " " << dir.Y() << " " << dir.Z() << " " << dir1.X() << " " << dir1.Y() << " " << dir1.Z() << " " << drift_dir.Angle(dir)/3.1415926*180. << " " << drift_dir.Angle(dir1)/3.1415926*180. << std::endl; */
      /* std::cout <<curr_wcp.x << " " << curr_wcp.y << " " << curr_wcp.z << " " << i << " " << dir1.Mag() << " " << dir1.Angle(dir)/3.14151926*180. << " " << test_p << " " << next_wcp.x << " " << next_wcp.y << " " << next_wcp.z << std::endl; */
      
      if (dir1.Mag()!=0 && (dir1.Angle(dir)/3.1415926*180. < 30)){
	//|| (fabs(drift_dir.Angle(dir)/3.1415926*180.-90.)<10. && fabs(drift_dir.Angle(dir1)/3.1415926*180.-90.)<10.) && dir1.Angle(dir)/3.1415926*180. < 60)){
	flag_continue = true;
	curr_wcp = next_wcp;
	dir = dir1 + dir * 5 * units::cm; // momentum trick ...
	dir = dir.Unit();
	break;
      }
      
      next_wcp = point_cloud_steiner->get_closest_wcpoint(test_p);
      TVector3 dir2(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z);
      
      //	std::cout << dir2.Angle(dir)/3.1415926*180. << " " << i << " " << drift_dir.Angle(dir2)/3.1415926*180. << " "<< drift_dir.Angle(dir)/3.1415926*180. << std::endl;
      
      if (dir2.Mag()!=0 && (dir2.Angle(dir)/3.1415926*180. < 30)){
	//|| (fabs(drift_dir.Angle(dir)/3.1415926*180.-90.)<10. && fabs(drift_dir.Angle(dir2)/3.1415926*180.-90.)<10.) && dir2.Angle(dir)/3.1415926*180. < 60)){
	flag_continue = true;
	curr_wcp = next_wcp;
	dir = dir1 + dir * 5 * units::cm; // momentum trick ...
	dir = dir.Unit();
	break;
      }
      
    }
  }

  
  test_p.x = curr_wcp.x;
  test_p.y = curr_wcp.y;
  test_p.z = curr_wcp.z;
  curr_wcp = point_cloud_steiner->get_closest_wcpoint(test_p);

  return curr_wcp;
  
}
    
