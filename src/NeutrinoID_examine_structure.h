void WCPPID::NeutrinoID::examine_structure(WCPPID::PR3DCluster *temp_cluster){
  examine_structure_1(temp_cluster);
}

bool WCPPID::NeutrinoID::examine_structure_1(WCPPID::PR3DCluster *temp_cluster){
  // look at the each segment, if the more straight one is better
  // change ...
  bool flag_update = false;
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
     WCPPID::ProtoSegment *sg = it->first;
     if (sg->get_cluster_id()!=temp_cluster->get_cluster_id()) continue;
     //     std::cout << sg->get_length()/units::cm << std::endl;
     double length = sg->get_length();
     double medium_dQ_dx = sg->get_medium_dQ_dx() / (43e3/units::cm);
     //     std::cout << length/units::cm << " " << medium_dQ_dx << std::endl;
     if (length < 5*units::cm || length < 8*units::cm && medium_dQ_dx > 1.5){
       std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> pair_vertices = find_vertices(sg);
       PointVector& pts = sg->get_point_vec();
       std::vector<WCP::WCPointCloud<double>::WCPoint >& wcps = sg->get_wcpt_vec();
       // check the track
       double step_size = 0.6*units::cm;
       Point start_p = pts.front();
       Point end_p = pts.back();
       int ncount = std::round(sqrt(pow(start_p.x-end_p.x,2)+pow(start_p.y-end_p.y,2)+pow(start_p.z-end_p.z,2))/step_size);
       PointVector new_pts;
       bool flag_replace = true;
       int n_bad = 0;
       for (int i=1;i<ncount;i++){
	 Point test_p;
	 test_p.x = start_p.x + (end_p.x - start_p.x)/ncount*i;
	 test_p.y = start_p.y + (end_p.y - start_p.y)/ncount*i;
	 test_p.z = start_p.z + (end_p.z - start_p.z)/ncount*i;
	 new_pts.push_back(test_p);
	 if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
	 if (n_bad>1) flag_replace = false;
	   //	 std::cout << i << " " << test_p << " " << ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0) << " " << sg->get_closest_point(test_p).first/units::cm << std::endl;
       }

       if (flag_replace){
	 WCP::ToyPointCloud* pcloud_steiner = temp_cluster->get_point_cloud_steiner();
	 WCP::WCPointCloud<double>::WCPoint start_wcp = wcps.front();
	 WCP::WCPointCloud<double>::WCPoint end_wcp = wcps.back();
	 wcps.clear();
	 wcps.push_back(start_wcp);
	 for (size_t i=0; i!=new_pts.size(); i++){
	   WCP::WCPointCloud<double>::WCPoint& wcp =  pcloud_steiner->get_closest_wcpoint(new_pts.at(i));
	   if (wcp.index != wcps.back().index) wcps.push_back(wcp);
	 }
	 if (end_wcp.index != wcps.back().index) wcps.push_back(end_wcp);
	 flag_update = true;
	 // std::cout << "Updated Straight Track" << std::endl;
       } // replace
     } // if length cut and dQ/dx cut
  } // loop over all the segments ...
  return flag_update;
}
