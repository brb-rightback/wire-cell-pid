void WireCellPID::PR3DCluster::do_stm_crawl(WireCell::WCPointCloud<double>::WCPoint& first_wcp, WireCell::WCPointCloud<double>::WCPoint& last_wcp){
  path_wcps.clear();
  
  //Crawl in point cloud 
  float step_dis = 1*units::cm;
  
  // initialization
  WireCell::WCPointCloud<double>::WCPoint curr_wcp = first_wcp;
  WireCell::WCPointCloud<double>::WCPoint next_wcp = first_wcp;
  Point p(first_wcp.x, first_wcp.y, first_wcp.z);
  TVector3 dir = VHoughTrans(p, 15*units::cm); // start direction
  Point test_p;
  
  bool flag_continue = true;
  while(flag_continue){
    flag_continue = false;
    for (int i=0; i!=3; i++){
      test_p.x = curr_wcp.x + dir.X() * step_dis * (i+1);
      test_p.y = curr_wcp.y + dir.Y() * step_dis * (i+1);
      test_p.z = curr_wcp.z + dir.Z() * step_dis * (i+1);

      next_wcp = point_cloud->get_closest_wcpoint(test_p);
      TVector3 dir1(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z);
      //std::cout <<curr_wcp.x << " " << curr_wcp.y << " " << curr_wcp.z << " " << i << " " << dir1.Mag() << " " << dir1.Angle(dir)/3.14151926*180. << " " << test_p << std::endl;
      if (dir1.Mag()!=0 && dir1.Angle(dir)/3.1415926*180. < 30){
	flag_continue = true;
	curr_wcp = next_wcp;
      	dir = dir1 + dir * 5 * units::cm; // momentum trick ...
	dir = dir.Unit();
	break;
      }
      
      next_wcp = point_cloud_steiner->get_closest_wcpoint(test_p);
      TVector3 dir2(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z);
      if (dir2.Mag()!=0 && dir2.Angle(dir)/3.1415926*180. < 30){
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

  
  //std::cout << first_wcp.x << " " << first_wcp.y << " " << first_wcp.z << " " << curr_wcp.x << " " << curr_wcp.y << " " << curr_wcp.z << std::endl;
    
  //Find the shortest path between first and middle point
  dijkstra_shortest_paths(first_wcp,2);
  cal_shortest_path(curr_wcp,2);

  //If middle point is far away from last_wcp, find the shortest path
  double dis = sqrt(pow(curr_wcp.x-last_wcp.x,2) + pow(curr_wcp.y-last_wcp.y,2) + pow(curr_wcp.z-last_wcp.z,2));
  if (dis > 1*units::cm){
    std::list<WireCell::WCPointCloud<double>::WCPoint> temp_path_wcps = path_wcps;
    dijkstra_shortest_paths(curr_wcp,2);
    cal_shortest_path(last_wcp,2);
    for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
      if (it == path_wcps.begin()) continue;
      temp_path_wcps.push_back(*it);
    }
    path_wcps = temp_path_wcps;
  }
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_mcells.size()==0){
      path_mcells.push_back( (*it).mcell);
    }else{
      if ( (*it).mcell!=path_mcells.back())
	path_mcells.push_back( (*it).mcell);
    }
  }
       
  // for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
  //  std::cout << (*it).x << " " << (*it).y << " " << (*it).z << std::endl;
  // }
  
  
  
  
  
 
}
