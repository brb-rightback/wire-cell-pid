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
