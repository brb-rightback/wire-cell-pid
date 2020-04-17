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
	 wcps.push_back(start_wcp);
	 for (size_t i=0; i!=new_pts.size(); i++){
	   WCP::WCPointCloud<double>::WCPoint& wcp =  pcloud_steiner->get_closest_wcpoint(new_pts.at(i));
	   if (wcp.index != wcps.back().index) wcps.push_back(wcp);
	 }
	 if (end_wcp.index != wcps.back().index) wcps.push_back(end_wcp);
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


      
      if (length1 > 5*units::cm || length2 > 5*units::cm) continue;

      //      std::cout << "asd " << std::endl;
      WCPPID::ProtoVertex *vtx1 = find_other_vertex(sg1, vtx);
      WCPPID::ProtoVertex *vtx2 = find_other_vertex(sg2, vtx);

      //      std::cout << vtx1 << " " << vtx2 << std::endl;
      
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
	//	 std::cout << i << " " << test_p << " " << ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0) << " " << sg->get_closest_point(test_p).first/units::cm << std::endl;
      }
      //      std::cout << flag_replace << std::endl;
      
      if (flag_replace){
	// form a new segment
	std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Merge two short segments: " << sg1->get_id() << " " << sg2->get_id() << " with a straight one "  << std::endl;
	
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
	WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps, temp_cluster->get_cluster_id()); acc_segment_id++;
	//	std::cout << wcps.size() << std::endl;
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
