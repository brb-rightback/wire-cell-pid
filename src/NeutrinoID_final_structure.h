
bool WCPPID::NeutrinoID::examine_structure_final(WCPPID::PR3DCluster* temp_cluster){
  // this is after the main vertex determination, but before the other clustering
  
  // results need to repeat the fit and the PID related work ...
  examine_structure_final_1(temp_cluster);

  
}

bool WCPPID::NeutrinoID::examine_structure_final_1(WCPPID::PR3DCluster* temp_cluster){
  // merge two segments if a direct connection is better ...
  bool flag_update = false;
  bool flag_continue = true;

  std::set<WCPPID::ProtoSegment*> new_segments; // need to fit and then PID ...
  
  while(flag_continue){
    flag_continue = false;
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
      WCPPID::ProtoVertex *vtx = it->first;
      if (vtx->get_cluster_id()!=temp_cluster->get_cluster_id()) continue;
      if (it->second.size()!=2) continue;
      WCPPID::ProtoSegment *sg1 = *(it->second.begin());
      WCPPID::ProtoSegment *sg2 = *(it->second.rbegin());
      double length1 = sg1->get_length();
      double length2 = sg2->get_length();
      
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
	std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Final stage merge two short segments: " << sg1->get_id() << " " << sg2->get_id() << " with a straight one, vtx id "  << vtx->get_id() << std::endl;
	
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
  if (flag_update)
    temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);

  return flag_update;
}




/* bool WCPPID::NeutrinoID::examine_vertices_5(WCPPID::PR3DCluster* temp_cluster){ */

/*   std::cout << "haha1 " << std::endl; */


  
/*   for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){ */
/*     WCPPID::ProtoVertex *vtx = it->first; */
/*     if (vtx->get_cluster_id() != temp_cluster->get_cluster_id()) continue; */
/*     if (it->second.size()!=2) continue; // require two segments ... */
    
/*     for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){ */
/*       WCPPID::ProtoSegment *sg1 = (*it1); */
/*       WCPPID::ProtoVertex *vtx1 = find_other_vertex(sg1, vtx); */
      
/*       auto it2 = map_vertex_segments.find(vtx1); */
/*       if ((it2->second).size()!=2) continue; // require two segments ... */
/*       WCPPID::ProtoSegment *sg2  = 0; */
/*       for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++){ */
/* 	if (*it3 != sg1){ */
/* 	  sg2 = *it3; */
/* 	  break; */
/* 	} */
/*       } */

/*       WCPPID::ProtoVertex *vtx2 = find_other_vertex(sg2, vtx1); */
/*       double dis = sqrt(pow(vtx->get_fit_pt().x - vtx2->get_fit_pt().x,2) + pow(vtx->get_fit_pt().y - vtx2->get_fit_pt().y,2) + pow(vtx->get_fit_pt().z - vtx2->get_fit_pt().z,2)); */
/*       if (dis < 5*units::cm){ */
/* 	WCPPID::ProtoSegment *sg = 0; */
/* 	for (auto it3 = it->second.begin(); it3 != it->second.end(); it3++){ */
/* 	  if ( *it3 != sg1) { */
/* 	    sg = *it3; */
/* 	    break; */
/* 	  } */
/* 	} */
/* 	std::cout << sg1->get_id() << " " << sg2->get_id() << std::endl; */
/* 	examine_vertices_5(sg, vtx, vtx2, sg1, vtx1, sg2); */
/*       } */
/*     } */
/*   } */

/*     std::cout << "haha2 " << std::endl; */
/* } */


/* bool WCPPID::NeutrinoID::examine_vertices_5(WCPPID::ProtoSegment *sg, WCPPID::ProtoVertex* vtx, WCPPID::ProtoVertex *vtx2, WCPPID::ProtoSegment *sg1, WCPPID::ProtoVertex *vtx1, WCPPID::ProtoSegment *sg2){ */
/*   TVector3 dir1(vtx2->get_fit_pt().x - vtx->get_fit_pt().x, vtx2->get_fit_pt().y - vtx->get_fit_pt().y, vtx2->get_fit_pt().z - vtx->get_fit_pt().z); */
/*   TVector3 dir2 = sg->cal_dir_3vector(vtx->get_fit_pt(), 6*units::cm); */

/*   std::cout << dir1.Angle(dir2)/3.1415926*180. << std::endl; */
/* } */
