void WCPPID::NeutrinoID::examine_structure(WCPPID::PR3DCluster *temp_cluster){
  // change 2 to 1
  
  

  
  if (examine_structure_2(temp_cluster) )
     temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
  // straighten 1
  if ( examine_structure_1(temp_cluster) )
     temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
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
	 //	 std::cout << start_wcp.index << std::endl;
	 wcps.push_back(start_wcp);
	 for (size_t i=0; i!=new_pts.size(); i++){
	   WCP::WCPointCloud<double>::WCPoint& wcp =  pcloud_steiner->get_closest_wcpoint(new_pts.at(i));
	   if (wcp.index != wcps.back().index ) {
	     //std::cout << new_pts.at(i) << " " << wcp.index << std::endl;
	     wcps.push_back(wcp);
	   }
	   if (wcps.back().index == end_wcp.index) break;
	 }
	 if (end_wcp.index != wcps.back().index) {
	   //std::cout << end_wcp.index << std::endl;
	   wcps.push_back(end_wcp);
	 }
	 flag_update = true;
	 std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " id " << sg->get_id() << " replace Track Content with Straight Line" << std::endl;
       } // replace
     } // if length cut and dQ/dx cut
  } // loop over all the segments ...
  return flag_update;
}


bool WCPPID::NeutrinoID::examine_structure_2(WCPPID::PR3DCluster *temp_cluster){
  bool flag_update = false;

  bool flag_continue = true;
  while(flag_continue){
    flag_continue = false;

    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
      WCPPID::ProtoVertex *vtx = it->first;
      //      std::cout <<"A: " << vtx->get_id() << " " << it->second.size() << std::endl;
      if (vtx->get_cluster_id()!=temp_cluster->get_cluster_id()) continue;
      if (it->second.size()!=2) continue;
      WCPPID::ProtoSegment *sg1 = *(it->second.begin());
      WCPPID::ProtoSegment *sg2 = *(it->second.rbegin());
      double length1 = sg1->get_length();
      double length2 = sg2->get_length();


      
      //if (length1 > 5*units::cm || length2 > 5*units::cm) continue;

      //      std::cout << "asd " << std::endl;
      WCPPID::ProtoVertex *vtx1 = find_other_vertex(sg1, vtx);
      WCPPID::ProtoVertex *vtx2 = find_other_vertex(sg2, vtx);


      
      double step_size = 0.6*units::cm;
      Point start_p = vtx1->get_fit_pt();
      Point end_p = vtx2->get_fit_pt();
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
	//	std::cout << i << " " << test_p << " " << ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0) << std::endl;
      }

      // std::cout << flag_replace << " " << sg1->get_id() << " " << sg2->get_id() << " " << length1/units::cm << " " << length2/units::cm << " " << n_bad << " " << ncount * step_size << std::endl;
      
      if (flag_replace){
	// form a new segment
	std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Merge two short segments: " << sg1->get_id() << " " << sg2->get_id() << " with a straight one, vtx id "  << vtx->get_id() << std::endl;
	
	WCP::ToyPointCloud* pcloud_steiner = temp_cluster->get_point_cloud_steiner();
	std::list<WCP::WCPointCloud<double>::WCPoint > wcps;
	WCP::WCPointCloud<double>::WCPoint start_wcp = vtx1->get_wcpt();
	WCP::WCPointCloud<double>::WCPoint end_wcp = vtx2->get_wcpt();
	wcps.push_back(start_wcp);
	for (size_t i=0; i!=new_pts.size(); i++){
	  WCP::WCPointCloud<double>::WCPoint& wcp =  pcloud_steiner->get_closest_wcpoint(new_pts.at(i));
	  if (wcp.index != wcps.back().index) wcps.push_back(wcp);
	}
	if (end_wcp.index != wcps.back().index) wcps.push_back(end_wcp);

	if (vtx1->get_wcpt().index != vtx2->get_wcpt().index){
	  WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps, temp_cluster->get_cluster_id()); acc_segment_id++;
	  add_proto_connection(vtx1,sg3,temp_cluster);
	  add_proto_connection(vtx2,sg3,temp_cluster);
	}else{
	  for (auto it5 = map_vertex_segments[vtx2].begin(); it5 != map_vertex_segments[vtx2].end(); it5++){
	    add_proto_connection(vtx1, *it5, temp_cluster);
	  }
	  del_proto_vertex(vtx2);
	}

	del_proto_segment(sg1);
	del_proto_segment(sg2);
	del_proto_vertex(vtx);

	
	flag_update = true;
	flag_continue = true;
	break;
      }
    }
  } // continue

  
  
  return flag_update;
}


bool WCPPID::NeutrinoID::examine_structure_3(WCPPID::PR3DCluster *temp_cluster){
  bool flag_update = false;
  bool flag_continue = true;
  while(flag_continue){
    flag_continue = false;
    
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
      WCPPID::ProtoVertex *vtx = it->first;
      //      std::cout <<"A: " << vtx->get_id() << " " << it->second.size() << std::endl;
      if (vtx->get_cluster_id()!=temp_cluster->get_cluster_id()) continue;
      if (it->second.size()!=2) continue;
      WCPPID::ProtoSegment *sg1 = *(it->second.begin());
      WCPPID::ProtoSegment *sg2 = *(it->second.rbegin());
      WCPPID::ProtoVertex *vtx1 = find_other_vertex(sg1, vtx);
      WCPPID::ProtoVertex *vtx2 = find_other_vertex(sg2, vtx);

      TVector3 dir1 = sg1->cal_dir_3vector(vtx->get_fit_pt(),10*units::cm);
      TVector3 dir2 = sg2->cal_dir_3vector(vtx->get_fit_pt(),10*units::cm);

      TVector3 dir3 = sg1->cal_dir_3vector(vtx->get_fit_pt(),3*units::cm);
      TVector3 dir4 = sg2->cal_dir_3vector(vtx->get_fit_pt(),3*units::cm);

      
      if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180. < 18 && (3.1415926 - dir3.Angle(dir4))/3.1415926*180. <27 ){
	std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Merge two segments: " << sg1->get_id() << " " << sg2->get_id() << " into one according to angle "  << std::endl;
	
	WCP::ToyPointCloud* pcloud_steiner = temp_cluster->get_point_cloud_steiner();

	std::list<WCP::WCPointCloud<double>::WCPoint > wcps(sg1->get_wcpt_vec().begin(), sg1->get_wcpt_vec().end());
	std::vector<WCP::WCPointCloud<double>::WCPoint > wcps_2 = sg2->get_wcpt_vec();

	if (wcps.front().index == wcps_2.front().index){
	  wcps.pop_front();
	  for (auto it1 = wcps_2.begin(); it1 != wcps_2.end(); it1++){
	    wcps.push_front(*it1);
	  }
	}else if (wcps.front().index == wcps_2.back().index){
	  wcps.pop_front();
	  for (auto it1 = wcps_2.rbegin(); it1 != wcps_2.rend(); it1++){
	    wcps.push_front(*it1);
	  }
	}else if (wcps.back().index == wcps_2.front().index){
	  wcps.pop_back();
	  for (auto it1 = wcps_2.begin(); it1 != wcps_2.end(); it1++){
	    wcps.push_back(*it1);
	  }
	}else if (wcps.back().index == wcps_2.back().index){
	  wcps.pop_back();
	  for (auto it1 = wcps_2.rbegin(); it1 != wcps_2.rend(); it1++){
	    wcps.push_back(*it1);
	  }
	}
	
	
	WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps, temp_cluster->get_cluster_id()); acc_segment_id++; 

	//delete the old segments and old vertices
	add_proto_connection(vtx1,sg3,temp_cluster);
	add_proto_connection(vtx2,sg3,temp_cluster);

	del_proto_segment(sg1);
	del_proto_segment(sg2);
	del_proto_vertex(vtx);

	
	flag_update = true;
	flag_continue = true;
	break;
      }
    }
  }


  
  
  
  return flag_update;
}


bool WCPPID::NeutrinoID::examine_structure_4(WCPPID::ProtoVertex* vertex, WCPPID::PR3DCluster *temp_cluster){
  bool flag_update = false;
  for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (it->second.size()<2) continue;
    if (vtx != vertex) continue;

    ToyPointCloud* pcloud = temp_cluster->get_point_cloud_steiner();
    std::vector<bool>& flag_terminals = temp_cluster->get_flag_steiner_terminal();
    
    std::vector<TVector3> saved_dirs;
    int nshowers = 0;
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      TVector3 dir = get_dir(vtx, *it1);
      if (dir.Mag()!=0) saved_dirs.push_back(dir);
    }

    std::vector<WCP::WCPointCloud<double>::WCPoint > candidate_wcps = pcloud->get_closest_wcpoints(vtx->get_fit_pt(), 6*units::cm);

    double max_dis = 0;
    WCP::WCPointCloud<double>::WCPoint max_wcp;

    for (size_t i=0; i!=candidate_wcps.size(); i++){
      if (flag_terminals.at(candidate_wcps.at(i).index)){
	double dis = sqrt(pow(candidate_wcps.at(i).x - vtx->get_fit_pt().x ,2) + pow(candidate_wcps.at(i).y - vtx->get_fit_pt().y,2) + pow(candidate_wcps.at(i).z - vtx->get_fit_pt().z,2));
	double min_dis = 1e9;
	double min_dis_u = 1e9;
	double min_dis_v = 1e9;
	double min_dis_w = 1e9;
	Point test_p(candidate_wcps.at(i).x, candidate_wcps.at(i).y, candidate_wcps.at(i).z);
      

	for (auto it1 = map_segment_vertices.begin(); it1!=map_segment_vertices.end(); it1++){
	  WCPPID::ProtoSegment *sg = it1->first;
	  //if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
	  WCP::WCPointCloud<double>::WCPoint closest_wcp = sg->get_closest_wcpt(test_p);
	  std::pair<double, WCP::Point> results = sg->get_closest_point(test_p);
	  if (results.first < min_dis) min_dis = results.first;
	  double tmp_dis = sqrt(pow(closest_wcp.x - test_p.x,2)+pow(closest_wcp.y-test_p.y,2)+pow(closest_wcp.z - test_p.z,2));
	  if (tmp_dis < min_dis) min_dis = tmp_dis;
	
	  std::tuple<double, double, double> results_2d = sg->get_closest_2d_dis(test_p);
	  if (std::get<0>(results_2d) < min_dis_u) min_dis_u = std::get<0>(results_2d);
	  if (std::get<1>(results_2d) < min_dis_v) min_dis_v = std::get<1>(results_2d);
	  if (std::get<2>(results_2d) < min_dis_w) min_dis_w = std::get<2>(results_2d);
	}

	//	std::cout << " " << candidate_wcps.at(i).index_u << " " << candidate_wcps.at(i).index_v << " " << candidate_wcps.at(i).index_w << " " << dis/units::cm << " " << min_dis/units::cm << " " << min_dis_u/units::cm << " " <<  min_dis_v/units::cm << " " <<  min_dis_w/units::cm << std::endl;
	
	if (min_dis > 0.9*units::cm   && min_dis_u + min_dis_v + min_dis_w > 1.8*units::cm &&
	    (min_dis_u > 0.8*units::cm && min_dis_v > 0.8*units::cm
	     || min_dis_u > 0.8*units::cm && min_dis_w > 0.8*units::cm
	     || min_dis_v > 0.8*units::cm && min_dis_w > 0.8*units::cm
	     )){
	  double step_size = 0.6*units::cm;
	  Point start_p = vtx->get_fit_pt();
	  Point end_p = test_p;
	  int ncount = std::round(sqrt(pow(start_p.x-end_p.x,2)+pow(start_p.y-end_p.y,2)+pow(start_p.z-end_p.z,2))/step_size);
	  bool flag_pass = true;
	  int n_bad = 0;
	  for (int j=1;j<ncount+1;j++){
	    Point test_p1;
	    test_p1.x = start_p.x + (end_p.x - start_p.x)/ncount*j;
	    test_p1.y = start_p.y + (end_p.y - start_p.y)/ncount*j;
	    test_p1.z = start_p.z + (end_p.z - start_p.z)/ncount*j;
	    if (!ct_point_cloud->is_good_point(test_p1, 0.2*units::cm, 0, 0)) n_bad ++;
	    if (n_bad>0) {
	      flag_pass = false;
	      break;
	    }
	  }

	  //	  std::cout << dis/units::cm << " " << min_dis/units::cm << " " << min_dis_u/units::cm << " " <<  min_dis_v/units::cm << " " <<  min_dis_w/units::cm << std::endl;
	  
	  if (flag_pass){
	    
	    if (max_dis < dis){
	      max_dis = dis;
	      max_wcp = candidate_wcps.at(i);
	    }
	  }
	}
      } // must be steiner terminal
    } // loop over points ...

    //    std::cout << max_dis/units::cm << std::endl;
    
    if (max_dis > 1.6*units::cm){
      WCPPID::ProtoVertex *v1 = new WCPPID::ProtoVertex(acc_vertex_id, max_wcp, temp_cluster->get_cluster_id()); acc_vertex_id++;
      std::list<WCP::WCPointCloud<double>::WCPoint> wcp_list;
      wcp_list.push_back(vtx->get_wcpt());
      
      double step_size = 1*units::cm;
      Point start_p = vtx->get_fit_pt();
      Point end_p(max_wcp.x, max_wcp.y, max_wcp.z);
      int ncount = std::round(sqrt(pow(start_p.x-end_p.x,2)+pow(start_p.y-end_p.y,2)+pow(start_p.z-end_p.z,2))/step_size);
      
      for (size_t j=1; j<ncount;j++){
	Point tmp_p;
	tmp_p.x = start_p.x + (end_p.x - start_p.x)/ncount*j;
	tmp_p.y = start_p.y + (end_p.y - start_p.y)/ncount*j;
	tmp_p.z = start_p.z + (end_p.z - start_p.z)/ncount*j;
	WCP::WCPointCloud<double>::WCPoint& wcp =  pcloud->get_closest_wcpoint(tmp_p);
	if (wcp.index != wcp_list.back().index) wcp_list.push_back(wcp);
      }
      if (max_wcp.index != wcp_list.back().index)
	wcp_list.push_back(max_wcp);
      
      std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Add a track to the main vertex" << " " << acc_segment_id << " " << wcp_list.size() << " " << v1->get_wcpt().index << " " << vtx->get_wcpt().index << std::endl;
      WCPPID::ProtoSegment* sg1 = new WCPPID::ProtoSegment(acc_segment_id, wcp_list, temp_cluster->get_cluster_id()); acc_segment_id++;
      add_proto_connection(v1,sg1,temp_cluster);
      add_proto_connection(vtx,sg1,temp_cluster);
      flag_update = true;
    }
  } // loop over vertices ...
  return flag_update;
}
