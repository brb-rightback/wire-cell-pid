
bool WCPPID::NeutrinoID::find_proto_vertex(WCPPID::PR3DCluster *temp_cluster, bool flag_break_track, int nrounds_find_other_tracks, bool flag_back_search){

  if (temp_cluster->get_point_cloud_steiner()==0) return false;
  if (temp_cluster->get_point_cloud_steiner()->get_num_points()<2) return false;
  
  WCPPID::ProtoSegment* sg1 = init_first_segment(temp_cluster, flag_back_search);

 
  
  if (sg1 == 0) return false;
  //  std::cout << "haha1 " << std::endl;
  if (temp_cluster == main_cluster) main_cluster_initial_pair_vertices = find_vertices(sg1);
  
  if (sg1->get_wcpt_vec().size()>1){
    // break tracks ... 
    if (flag_break_track){
      std::vector<WCPPID::ProtoSegment*> remaining_segments;
      remaining_segments.push_back(sg1);
      break_segments(remaining_segments, temp_cluster);

      
      
      //      if (flag_check_end_segments)	check_end_segments(temp_cluster);
      // if a straight length is better for a segment ...
      examine_structure(temp_cluster);
      
    }else{
      
      temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);

      
    }
    // std::cout << "haha2 " << std::endl;
    
    /* for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){ */
    /*   WCPPID::ProtoVertex *vtx = it->first; */
    /*   if (vtx->get_cluster_id() == 50){ */
    /* 	std::cout << vtx->get_id() << " " << vtx->get_fit_pt() << " " << vtx->get_wcpt().x << " " << vtx->get_wcpt().y << " " << vtx->get_wcpt().z << std::endl; */
    /*   } */
    /* } */
    
    // find other segments ...
    for (size_t i=0;i!=nrounds_find_other_tracks;i++){
      find_other_segments(temp_cluster, flag_break_track);
    }

    if (temp_cluster == main_cluster){
      // merge two tracks if their angles are consistent
      if ( examine_structure_3(temp_cluster) )
	temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
      
    }
    
    
    
    
    // examine the vertices ...
    examine_vertices(temp_cluster);
    
    // examine the two initial points ...
    if (temp_cluster == main_cluster && main_cluster_initial_pair_vertices.first!=0)
      examine_vertices_3();

    
    
    return true;
  }
  
  return false;

  

  // for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
  //   WCPPID::ProtoSegment *sg = it->first;
  //   std::cout << sg->get_id() << " " << sg->get_point_vec().size() << " " << sg->get_point_vec().front() << " " << sg->get_point_vec().back() << std::endl;
  //   for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  //     std::cout << (*it1)->get_id() << " " << (*it1) << std::endl;
  //   }
  // }
  // for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
  //   std::cout << (it->first)->get_id() << std::endl;
  //   for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  //     std::cout << (*it1)->get_id() << " " << (*it1) << std::endl;
  //   }
  // }
  
  //  temp_cluster->set_fit_parameters(map_vertex_segments, map_segment_vertices);
  // practice other components ...
}


void WCPPID::NeutrinoID::check_end_segments(WCPPID::PR3DCluster* temp_cluster){
  bool flag_check = true;
  while(flag_check){
    flag_check = false;
    for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> tmp_vtxs = find_vertices(sg);
      if (map_vertex_segments[tmp_vtxs.first].size() == 1 && map_vertex_segments[tmp_vtxs.second].size() > 1){
	
	if (sg->get_length() < 5*units::cm){
	  del_proto_vertex(tmp_vtxs.first);
	  del_proto_segment(sg);
	  flag_check = true;
	  break;
	}
	//	std::cout << sg->get_length()/units::cm << " " << sg->get_medium_dQ_dx() << std::endl;
      }else if (map_vertex_segments[tmp_vtxs.first].size() > 1 && map_vertex_segments[tmp_vtxs.second].size() == 1){
	//	std::cout << sg->get_length()/units::cm << " " << sg->get_medium_dQ_dx()<< std::endl;
	if (sg->get_length() < 5*units::cm){
	  del_proto_vertex(tmp_vtxs.second);
	  del_proto_segment(sg);
	  flag_check = true;
	  break;
	}
      }
    }
  }
}



void WCPPID::NeutrinoID::init_point_segment(WCPPID::PR3DCluster *temp_cluster){
  // do the first search of the trajectory ...
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(1);
  // good for the first track
  temp_cluster->dijkstra_shortest_paths(wcps.first,1); 
  temp_cluster->cal_shortest_path(wcps.second,1);
  WCPPID::ProtoVertex *v1=0;
  WCPPID::ProtoVertex *v2=0;
  WCPPID::ProtoSegment *sg1=0;

  if (temp_cluster->get_path_wcps().size()>1){
    v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.first, temp_cluster->get_cluster_id()); acc_vertex_id++;
    v2 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.second, temp_cluster->get_cluster_id()); acc_vertex_id++;
    sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
    add_proto_connection(v1,sg1,temp_cluster);
    add_proto_connection(v2,sg1,temp_cluster);

    temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
  }
}

WCPPID::ProtoSegment* WCPPID::NeutrinoID::init_first_segment(WCPPID::PR3DCluster *temp_cluster, bool flag_back_search){
  // do the first search of the trajectory ...
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2);

  if (temp_cluster == main_cluster){
    // main cluster, start from the downstream point ...
    if (flag_back_search){
      if (wcps.first.z < wcps.second.z){
	auto wcp1 = wcps.first;
	auto wcp2 = wcps.second;
	wcps.first = wcp2;
	wcps.second = wcp1;
      }
    }else{
      if (wcps.first.z > wcps.second.z){
	auto wcp1 = wcps.first;
	auto wcp2 = wcps.second;
	wcps.first = wcp2;
	wcps.second = wcp1;
      }
    }
  }else{
    Point tp1(wcps.first.x, wcps.first.y, wcps.first.z);
    Point tp2(wcps.second.x, wcps.second.y, wcps.second.z);
    double dis1 = main_cluster->get_point_cloud_steiner()->get_closest_dis(tp1);
    double dis2 = main_cluster->get_point_cloud_steiner()->get_closest_dis(tp2);
    if (dis2 < dis1){
      auto wcp1 = wcps.first;
      auto wcp2 = wcps.second;
      wcps.first = wcp2;
      wcps.second = wcp1;
    }
  }
  
  
  // good for the first track
  temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
  temp_cluster->cal_shortest_path(wcps.second,2);

  
  /* std::cout << temp_cluster->get_cluster_id() << " " << wcps.first.index << " " << wcps.first.x << " " << wcps.first.y << " " << wcps.first.z << " " << wcps.first.index_u << " " << wcps.first.index_v << " " << wcps.first.index_w << " " << wcps.second.index << " " << wcps.second.x << " " << wcps.second.y << " " << wcps.second.z << " " << wcps.second.index_u << " " << wcps.second.index_v << " " << wcps.second.index_w << " " << std::endl; */
  /* { */
  /*   Point test_p(wcps.first.x, wcps.first.y, wcps.first.z); */
  /*   std::vector<int> results = ct_point_cloud->convert_3Dpoint_time_ch(test_p); */
  /*   std::cout << results.at(0) << " " << results.at(1) << " " << results.at(2) << " " << results.at(3) << std::endl; */
  /* } */
  //std::cout << temp_cluster->get_path_wcps().front().index << " " << temp_cluster->get_path_wcps().back().index << std::endl;
  
  WCPPID::ProtoVertex *v1=0;
  WCPPID::ProtoVertex *v2=0;
  WCPPID::ProtoSegment *sg1=0;
  
  if (temp_cluster->get_path_wcps().size()>1){
    v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.first, temp_cluster->get_cluster_id()); acc_vertex_id++;
    v2 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.second, temp_cluster->get_cluster_id()); acc_vertex_id++;
    sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
    add_proto_connection(v1,sg1,temp_cluster);
    add_proto_connection(v2,sg1,temp_cluster);

    temp_cluster->collect_charge_trajectory(*ct_point_cloud);

    //    std::cout << "abc1" << std::endl;
    // fit dQ/dx and everything ...
    temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true);

    //    std::cout << "abc" << std::endl;
    
    v1->set_fit(temp_cluster->get_fine_tracking_path().front(), temp_cluster->get_dQ().front(), temp_cluster->get_dx().front(), temp_cluster->get_pu().front(), temp_cluster->get_pv().front(), temp_cluster->get_pw().front(), temp_cluster->get_pt().front(), temp_cluster->get_reduced_chi2().front());
    v2->set_fit(temp_cluster->get_fine_tracking_path().back(), temp_cluster->get_dQ().back(), temp_cluster->get_dx().back(), temp_cluster->get_pu().back(), temp_cluster->get_pv().back(), temp_cluster->get_pw().back(), temp_cluster->get_pt().back(), temp_cluster->get_reduced_chi2().back());
    sg1->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2());

    //    std::cout << temp_cluster->get_fine_tracking_path().size() << " " << temp_cluster->get_dx().size() << " " << temp_cluster->get_fine_tracking_path().front() << " " << temp_cluster->get_fine_tracking_path().back() << " " << v1->get_pw() << " " << wcps.first.index_w << std::endl;
  }
  //else{
  //  v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.first, temp_cluster->get_cluster_id()); acc_vertex_id++;
  //  sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
  //   add_proto_connection(v1,sg1,temp_cluster);
  //}  
  
  // std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg1].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
  //  del_proto_vertex(v1); // del_proto_vertex(v1);
  // std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg1].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
   // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v1].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;
   // del_proto_segment(sg1); // del_proto_segment(sg1);
   // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v1].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;

  // temporary fit ...
  
  // finish the temporary fit ...

  //  std::cout << temp_cluster->get_fine_tracking_path().size() << " " << temp_cluster->get_fine_tracking_path().at(40) << std::endl;

  
  // sg1->print_dis();
  // std::cout << v1->get_fit_init_dis()/units::cm << " " << v2->get_fit_init_dis()/units::cm << std::endl;
  return sg1;
}




void WCPPID::NeutrinoID::break_segments(std::vector<WCPPID::ProtoSegment*>& remaining_segments, WCPPID::PR3DCluster* temp_cluster){
  bool flag_print = false;
  
  int count = 0;
  
  while(remaining_segments.size()!=0 && count < 2){
    WCPPID::ProtoSegment* curr_sg = remaining_segments.back();
    remaining_segments.pop_back();
    WCPPID::ProtoVertex *start_v=0, *end_v=0;
    for (auto it = map_segment_vertices[curr_sg].begin(); it!=map_segment_vertices[curr_sg].end(); it++){
      if ((*it)->get_wcpt().index == curr_sg->get_wcpt_vec().front().index) start_v = *it;
      if ((*it)->get_wcpt().index == curr_sg->get_wcpt_vec().back().index) end_v = *it;
    }
    if (start_v==0 || end_v==0){
      std::cout << "Error in finding vertices for a segment" << std::endl; 
      //std::cout << "Vertex: " << start_v->get_wcpt().index << " " << end_v->get_wcpt().index << std::endl;
    }
    
    // initialize the start points
    WCP::WCPointCloud<double>::WCPoint break_wcp = start_v->get_wcpt();
    Point test_start_p = curr_sg->get_point_vec().front();
    
    while(sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2) +
	       pow(start_v->get_wcpt().y - break_wcp.y,2) +
	       pow(start_v->get_wcpt().z - break_wcp.z,2)) <= 1*units::cm &&
	  sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2) +
	       pow(end_v->get_wcpt().y - break_wcp.y,2) +
	       pow(end_v->get_wcpt().z - break_wcp.z,2)) > 1*units::cm){
      std::tuple<Point, TVector3, TVector3, bool> kink_tuple = curr_sg->search_kink(test_start_p);
      if (std::get<1>(kink_tuple).Mag()!=0 ){
	// find the extreme point ... PR3DCluster function
	break_wcp = temp_cluster->proto_extend_point(std::get<0>(kink_tuple), std::get<1>(kink_tuple), std::get<2>(kink_tuple), std::get<3>(kink_tuple));
	//std::cout << "Break point: " << remaining_segments.size() << " " << std::get<0>(kink_tuple) << " " << std::get<1>(kink_tuple).Mag() << " " << std::get<3>(kink_tuple) << " " << break_wcp.x/units::cm << " " << break_wcp.y/units::cm << " " << break_wcp.z/units::cm << " " << break_wcp.index << std::endl;
        if (sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2) +
	       pow(start_v->get_wcpt().y - break_wcp.y,2) +
	       pow(start_v->get_wcpt().z - break_wcp.z,2)) <= 1*units::cm &&
	    sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2) +
		 pow(end_v->get_wcpt().y - break_wcp.y,2) +
		 pow(end_v->get_wcpt().z - break_wcp.z,2)) > 1*units::cm)
	  test_start_p = std::get<0>(kink_tuple);
      }else{
	break;
      }
    }

    if (sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2) +
	     pow(start_v->get_wcpt().y - break_wcp.y,2) +
	     pow(start_v->get_wcpt().z - break_wcp.z,2)) > 1*units::cm){
      // find the new path ... PR3DCluster function
      // break into two tracks and adding a ProtoVertex in this function ...
      std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list1;
      std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list2;
      bool flag_break = temp_cluster->proto_break_tracks(start_v->get_wcpt(), break_wcp, end_v->get_wcpt(), wcps_list1, wcps_list2);

      //      std::cout << break_wcp.index_u << " " << break_wcp.index_v << " " << break_wcp.index_w << " " << break_wcp.mcell->GetTimeSlice() << std::endl;
      std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " breaking point: " << flag_break << " " << wcps_list1.front().index << " " << wcps_list1.back().index << " " << wcps_list2.front().index << " " << wcps_list2.back().index << std::endl;
	  
      if (flag_break){
	TVector3 tv1(end_v->get_wcpt().x - start_v->get_wcpt().x, end_v->get_wcpt().y - start_v->get_wcpt().y, end_v->get_wcpt().z - start_v->get_wcpt().z);
	TVector3 tv2(end_v->get_wcpt().x - break_wcp.x, end_v->get_wcpt().y - break_wcp.y, end_v->get_wcpt().z - break_wcp.z);

	double min_dis = 1e9;
	for (auto it5 = wcps_list1.begin(); it5!=wcps_list1.end(); it5++){
	  double dis = sqrt(pow((*it5).x - end_v->get_wcpt().x,2) + pow((*it5).y - end_v->get_wcpt().y,2) + pow((*it5).z - end_v->get_wcpt().z,2));
	  if (dis < min_dis) min_dis = dis;
	}
	if (min_dis/units::cm <1.5 && tv1.Angle(tv2)/3.1415926*180.> 120){
	  WCPPID::ProtoVertex *v3 = new WCPPID::ProtoVertex(acc_vertex_id, break_wcp, temp_cluster->get_cluster_id()); acc_vertex_id ++;
	  WCPPID::ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list1, temp_cluster->get_cluster_id()); acc_segment_id++;
	  add_proto_connection(start_v, sg2, temp_cluster);
	  add_proto_connection(v3, sg2, temp_cluster);
	  del_proto_vertex(end_v);
	  del_proto_segment(curr_sg);
	  flag_print = true;
	  // fit dQ/dx here, do not exclude others
	  temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, false);	
	}else{
	  //	std::cout << tv1.Angle(tv2)/3.1415926*180. << " " << min_dis/units::cm << std::endl;

	 
	  
	  WCPPID::ProtoVertex *v3 = new WCPPID::ProtoVertex(acc_vertex_id, break_wcp, temp_cluster->get_cluster_id()); acc_vertex_id ++;
	  WCPPID::ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list1, temp_cluster->get_cluster_id()); acc_segment_id++;
	  WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list2, temp_cluster->get_cluster_id()); acc_segment_id++;

	  //	  std::cout <<"break " << curr_sg->get_id() << " " << sg2->get_id() << " " << sg3->get_id() << std::endl;
	  
	  add_proto_connection(start_v, sg2, temp_cluster);
	  add_proto_connection(v3, sg2, temp_cluster);
	  add_proto_connection(v3, sg3, temp_cluster);
	  add_proto_connection(end_v, sg3, temp_cluster);
	  
	  del_proto_segment(curr_sg);
	  //delete curr_sg;
	  
	  flag_print = true;
	  // fit dQ/dx here, do not exclude others
	  temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, false);	
	  remaining_segments.push_back(sg3);
	}
      }
    }
    
    //    std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg2].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
    // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v3].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;
  }

  if (flag_print)
    std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Break tracks  -- # of Vertices: " << map_vertex_segments.size() << "; # of Segments: " << map_segment_vertices.size() << std::endl;
}




void WCPPID::NeutrinoID::find_other_segments(WCPPID::PR3DCluster* temp_cluster, bool flag_break_track, double search_range, double scaling_2d){
  

  const int N = temp_cluster->get_point_cloud_steiner()->get_num_points();
  std::vector<bool>& flag_tagged = temp_cluster->get_flag_tagged_steiner_graph();
  flag_tagged.resize(N, false);
  int num_tagged = 0;
  WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud_steiner()->get_cloud();
  for (size_t i=0;i!=N;i++){
    Point p(cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z);
    double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
    double min_3d_dis = 1e9;
    for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      if (sg1->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
      std::pair<double, WCP::Point> closest_dis_point = sg1->get_closest_point(p);
      std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p);
      if (closest_dis_point.first < min_3d_dis) min_3d_dis = closest_dis_point.first;
      if (closest_dis_point.first < search_range){
	flag_tagged[i] = true;
	num_tagged ++;
	break;
      }
      if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis);
      if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis);
      if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis);
    }

    if (!flag_tagged[i]){
      if ((min_dis_u < scaling_2d * search_range || ct_point_cloud->get_closest_dead_chs(p, 0) ) &&
	  (min_dis_v < scaling_2d * search_range || ct_point_cloud->get_closest_dead_chs(p, 1) ) &&
	  (min_dis_w < scaling_2d * search_range || ct_point_cloud->get_closest_dead_chs(p, 2) ) )
	flag_tagged[i] = true;
	  // std::cout << min_dis_u/units::cm << " " << min_dis_v/units::cm << " " << min_dis_w/units::cm << std::endl;
    }
    
    //    if (cloud.pts[i].index_u< 820 && cloud.pts[i].index_u >800 && cloud.pts[i].index_w < 5525-4800 && cloud.pts[i].index_w > 5505-4800 && p.x < 200*units::cm && p.x > 193*units::cm && min_dis_v < 0.6*units::cm){
    //std::cout << p << " " << cloud.pts[i].index << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v+2400 << " " << cloud.pts[i].index_w+4800 << " " << min_dis_u/units::cm << " " << min_dis_v/units::cm << " " << min_dis_w/units::cm << " " <<  min_3d_dis/units::cm << " " << flag_tagged[i] << std::endl;
    //}
    
    
  }



  // Now figure out the Terminal Graph ...
  std::vector<int> terminals;
  std::map<int, int> map_oindex_tindex;
  std::vector<bool>& flag_steiner_terminal = temp_cluster->get_flag_steiner_terminal();
  for (size_t i = 0;i!=flag_steiner_terminal.size();i++){
    if (flag_steiner_terminal[i]){
      map_oindex_tindex[i] = terminals.size();
      terminals.push_back(i);
    }
  }
  //  std::cout << N << " abc " << terminals.size() << std::endl;

  using Vertex = typename boost::graph_traits<WCPPID::MCUGraph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<WCPPID::MCUGraph>::edge_descriptor;
  using Base = typename boost::property<edge_base_t, Edge>;
  using EdgeWeightMap = typename boost::property_map<WCPPID::MCUGraph, boost::edge_weight_t>::type;
  using Weight = typename boost::property_traits<EdgeWeightMap>::value_type;
  using WeightProperty =
    typename boost::property<boost::edge_weight_t, Weight, Base>;
  using TerminalGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, WeightProperty>;
  using EdgeTerminal =
    typename boost::graph_traits<TerminalGraph>::edge_descriptor;
  
  WCPPID::MCUGraph* graph_steiner = temp_cluster->get_graph_steiner();
  EdgeWeightMap edge_weight = boost::choose_pmap(boost::get_param(boost::no_named_parameters(), boost::edge_weight), *graph_steiner, boost::edge_weight);
  
  // distance array used in the dijkstra runs
  std::vector<Weight> distance(N);

  std::vector<Vertex> nearest_terminal(num_vertices(*graph_steiner));
  auto index = get(boost::vertex_index, *graph_steiner);
  auto nearest_terminal_map = boost::make_iterator_property_map(nearest_terminal.begin(), get(boost::vertex_index, *graph_steiner));
  for (auto terminal : terminals) {
    nearest_terminal_map[terminal] = terminal;
  }
  // compute voronoi diagram each vertex get nearest terminal and last edge on
  // path to nearest terminal
  auto distance_map = make_iterator_property_map(distance.begin(), index);
  std::vector<Edge> vpred(N);
  auto last_edge = boost::make_iterator_property_map(vpred.begin(), get(boost::vertex_index, *graph_steiner));
  boost::dijkstra_shortest_paths(*graph_steiner, terminals.begin(), terminals.end(), boost::dummy_property_map(),
				 distance_map, edge_weight, index, paal::utils::less(),
				 boost::closed_plus<Weight>(), std::numeric_limits<Weight>::max(), 0,
				 boost::make_dijkstra_visitor(paal::detail::make_nearest_recorder(
												  nearest_terminal_map, last_edge, boost::on_edge_relaxed{})));


  //  std::cout << "abc1 " << std::endl; 
  // computing distances between terminals
  // creating terminal_graph
  TerminalGraph terminal_graph(N);
  std::map<std::pair<int, int>, std::pair<Weight, Edge> > map_saved_edge;

  for (auto w : boost::as_array(edges(*graph_steiner))) {
    auto const &nearest_to_source = nearest_terminal_map[source(w, *graph_steiner)];
    auto const &nearest_to_target = nearest_terminal_map[target(w, *graph_steiner)];
    if (nearest_to_source != nearest_to_target) {
      Weight temp_weight = distance[source(w, *graph_steiner)] + distance[target(w, *graph_steiner)] + edge_weight[w];
      if (map_saved_edge.find(std::make_pair(nearest_to_source, nearest_to_target))!=map_saved_edge.end()){
	if (temp_weight < map_saved_edge[std::make_pair(nearest_to_source, nearest_to_target)].first)
	  map_saved_edge[std::make_pair(nearest_to_source, nearest_to_target)] = std::make_pair(temp_weight,w);
      }else if (map_saved_edge.find(std::make_pair(nearest_to_target, nearest_to_source))!=map_saved_edge.end()){
	if (temp_weight < map_saved_edge[std::make_pair(nearest_to_target, nearest_to_source)].first)
	  map_saved_edge[std::make_pair(nearest_to_target, nearest_to_source)] = std::make_pair(temp_weight,w);
      }else{
	map_saved_edge[std::make_pair(nearest_to_source, nearest_to_target)] = std::make_pair(temp_weight,w);
      }
    }
  }

  std::vector<Edge> terminal_edge;
  for (auto it = map_saved_edge.begin(); it!=map_saved_edge.end(); it++){
    std::pair<Edge, bool> p = add_edge(it->first.first, it->first.second,
				       WeightProperty(it->second.first, Base(it->second.second)),terminal_graph);
    //      terminal_edge.push_back(p.first);
  }
  // minimal spanning tree ...
  boost::kruskal_minimum_spanning_tree(terminal_graph,
				       std::back_inserter(terminal_edge));
  
  // new Graph ...
  std::map<int, std::set<int> > map_connection;
  TerminalGraph terminal_graph_cluster(terminals.size());
  for (auto it = terminal_edge.begin(); it!=terminal_edge.end(); it++){
    if (flag_tagged[source(*it, terminal_graph)] == flag_tagged[target(*it, terminal_graph)]){
      add_edge(map_oindex_tindex[source(*it, terminal_graph)],map_oindex_tindex[target(*it, terminal_graph)], terminal_graph_cluster);
    }else{
      if (map_connection.find(source(*it, terminal_graph))==map_connection.end()){
	std::set<int> temp_results;
	temp_results.insert(target(*it,terminal_graph));
	map_connection[source(*it, terminal_graph)] = temp_results;
      }else{
	map_connection[source(*it, terminal_graph)].insert(target(*it,terminal_graph));
      }
      if (map_connection.find(target(*it, terminal_graph))==map_connection.end()){
	std::set<int> temp_results;
	temp_results.insert(source(*it,terminal_graph));
	map_connection[target(*it, terminal_graph)] = temp_results;
      }else{
	map_connection[target(*it, terminal_graph)].insert(source(*it,terminal_graph));
      }
    }
    //      std::cout << source(*it, terminal_graph) << " " << target(*it, terminal_graph) << std::endl;
  }
  //std::cout << map_connection.size() << std::endl;
  //    std::cout << map_saved_edge.size() << " " << terminal_edge.size() << std::endl;

  //  std::cout << "abc2 " << std::endl; 
  
  std::vector<int> component(num_vertices(terminal_graph_cluster));
  const int num = connected_components(terminal_graph_cluster,&component[0]);
  //  std::cout << "jaja: " << num << std::endl;
  int ncounts[num];
  std::vector<std::vector<int> > sep_clusters(num);
  for (size_t i=0;i!=num;i++){
    ncounts[i] = 0;
    std::vector<int> temp;
    sep_clusters[i] = temp;
  }  
  {
    std::vector<int>::size_type i;
    for (i=0;i!=component.size(); ++i){
      ncounts[component[i]]++;
      sep_clusters[component[i]].push_back(terminals.at(i));
    }
  }

  //test
  //  for (size_t j=0;j!=sep_clusters.size();j++){
  // for(size_t i=0;i!=sep_clusters[j].size();i++){
  //  if (sep_clusters[j].at(i) > 530 && sep_clusters[j].at(i) < 570)
  //	std::cout << j << " " << i << " " << sep_clusters[j].at(i) << " " << flag_tagged[sep_clusters[j].at(i)] << std::endl;
  //}
  //}

  


  std::vector<WCPPID::Res_proto_segment> temp_segments(num);
  std::set<int> remaining_segments;
  for (int i=0;i!=num;i++){
    remaining_segments.insert(i);
  }

  for (size_t i=0;i!=num;i++){
    // inside the original track 
    if (flag_tagged[sep_clusters[i].front()] ) {
      remaining_segments.erase(i);
      continue;
    }
    temp_segments.at(i).group_num = i;
    temp_segments.at(i).number_points = ncounts[i];
    
    int special_A = -1;
    for (size_t j=0;j!=ncounts[i];j++){
      if (map_connection.find(sep_clusters[i].at(j))!=map_connection.end()){
	special_A = sep_clusters[i].at(j);
	break;
      }
    }

    //    if (special_A !=-1){
    //  std::cout << i << " " << special_A << " " << sep_clusters[i].size() << std::endl;
    //  for (auto it1 = map_connection[special_A].begin(); it1 != map_connection[special_A].end(); it1++){
    // 	std::cout << *it1 << std::endl;
    //  }
    // }
    
    int special_B = special_A;
    double min_dis = 0;
    
    int number_not_faked = 0;
    double max_dis_u = 0, max_dis_v = 0, max_dis_w = 0;
    for (size_t j=0;j!=ncounts[i];j++){
      double dis = sqrt(pow(cloud.pts[sep_clusters[i].at(j)].x - cloud.pts[special_A].x,2)
			+ pow(cloud.pts[sep_clusters[i].at(j)].y - cloud.pts[special_A].y,2)
			+ pow(cloud.pts[sep_clusters[i].at(j)].z - cloud.pts[special_A].z,2));
      if (dis > min_dis){
	min_dis = dis;
	special_B = sep_clusters[i].at(j); // furthest away from special_A
      }

      // also judge whether this track is fake ...
      WCP::Point p(cloud.pts[sep_clusters[i].at(j)].x, cloud.pts[sep_clusters[i].at(j)].y, cloud.pts[sep_clusters[i].at(j)].z);
      double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
      for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (sg1->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
	std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p);
	if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis);
	if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis);
	if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis);
      }
      int flag_num = 0;
      if (min_dis_u > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 0))) flag_num ++;
      if (min_dis_v > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 1))) flag_num ++;
      if (min_dis_w > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 2))) flag_num ++;

      //      std::cout << cloud.pts[sep_clusters[i].at(j)].index << " " << min_dis_u/units::cm << " " << min_dis_v/units::cm << " " << min_dis_w/units::cm << " " << ct_point_cloud->get_closest_dead_chs(p, 0) << " " << ct_point_cloud->get_closest_dead_chs(p, 1) << " " << ct_point_cloud->get_closest_dead_chs(p, 2) << std::endl;

      if (min_dis_u > max_dis_u && (!ct_point_cloud->get_closest_dead_chs(p, 0))) max_dis_u = min_dis_u;
      if (min_dis_v > max_dis_v && (!ct_point_cloud->get_closest_dead_chs(p, 1))) max_dis_v = min_dis_v;
      if (min_dis_w > max_dis_w && (!ct_point_cloud->get_closest_dead_chs(p, 2))) max_dis_w = min_dis_w;
      if (flag_num>=2) number_not_faked ++;
      
    }

    
    double length = sqrt(pow(cloud.pts[special_A].x - cloud.pts[special_B].x,2) + pow(cloud.pts[special_A].y - cloud.pts[special_B].y,2) + pow(cloud.pts[special_A].z - cloud.pts[special_B].z,2)); 
    // check ...
    
    if ( length < 3*units::cm){
      int save_index = special_A;
      double save_dis = 1e9;
      for (auto it1 = map_connection[special_A].begin(); it1!=map_connection[special_A].end();it1++){
    	double temp_dis = sqrt(pow(cloud.pts[special_A].x - cloud.pts[*it1].x,2) + pow(cloud.pts[*it1].y - cloud.pts[*it1].y,2) + pow(cloud.pts[*it1].z - cloud.pts[*it1].z,2));
    	if (temp_dis < save_dis){
    	  save_index = *it1;
    	  save_dis = temp_dis;
    	}
      }
      special_A = save_index;
      length = sqrt(pow(cloud.pts[special_A].x - cloud.pts[special_B].x,2) + pow(cloud.pts[special_A].y - cloud.pts[special_B].y,2) + pow(cloud.pts[special_A].z - cloud.pts[special_B].z,2));
    }
    //    std::cout << special_A << " " << special_B << " " << length/units::cm << std::endl;
    
    temp_segments.at(i).special_A = special_A;
    temp_segments.at(i).special_B = special_B;
    temp_segments.at(i).length = length;
    temp_segments.at(i).number_not_faked = number_not_faked;
    temp_segments.at(i).max_dis_u = max_dis_u;
    temp_segments.at(i).max_dis_v = max_dis_v;
    temp_segments.at(i).max_dis_w = max_dis_w;

    // std::cout << special_A << " " << special_B << " " << length/units::cm << " " << temp_segments.at(i).number_points << " "  << number_not_faked << std::endl;
    
    if (temp_segments.at(i).number_points ==1  //  only one point 
    	|| number_not_faked == 0 &&
	(length < 3.5*units::cm  // very short & fake
	 || (number_not_faked < 0.25 * temp_segments.at(i).number_points || number_not_faked < 0.4 * temp_segments.at(i).number_points && length < 7 * units::cm) && max_dis_u/units::cm < 3 && max_dis_v/units::cm < 3 && max_dis_w/units::cm < 3 && max_dis_u + max_dis_v + max_dis_w < 6*units::cm)  // many fake things and very close to each other ...
     	){
      remaining_segments.erase(i);
    }
  }

  /* for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){  */
  /*   std::cout << "jaja1: " << *it << " " << temp_segments.at(*it).number_not_faked << " " << temp_segments.at(*it).number_points << " " <<  temp_segments.at(*it).max_dis_u/units::cm << " " << temp_segments.at(*it).max_dis_v/units::cm << " " << temp_segments.at(*it).max_dis_w/units::cm << " " << temp_segments.at(*it).length/units::cm << std::endl;  */
  /* }  */
  // std::cout << remaining_segments.size() << std::endl;

  // plan to examine things ... 
  std::vector<int> saved_clusters;
  std::vector<std::pair<int, int> > saved_cluster_points(num);

  /* saved_clusters.push_back(8); */
  /* saved_cluster_points.at(8) = std::make_pair(temp_segments.at(8).special_A, temp_segments.at(8).special_B); */
  
  ProtoSegmentSelection temp_segments_1, temp_segments_2;
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    if (it->first->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    temp_segments_1.push_back(it->first);
  }

  while (remaining_segments.size()>0){
    //    std::cout << num << " " << remaining_segments.size() << std::endl;
    // find the maximal length cluster
    double max_length = 0; int max_length_cluster = -1;
    double max_number_not_faked = 0;
    for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){
      if (temp_segments.at(*it).number_not_faked > max_number_not_faked){
  	max_length_cluster = *it;
  	max_number_not_faked = temp_segments.at(*it).number_not_faked;
  	max_length = temp_segments.at(*it).length;
      }else if (temp_segments.at(*it).number_not_faked == max_number_not_faked){
  	if (temp_segments.at(*it).length > max_length){
  	  max_length_cluster = *it;
  	  max_number_not_faked = temp_segments.at(*it).number_not_faked;
  	  max_length = temp_segments.at(*it).length;
  	}
      }
    }
    // save things ...
    saved_clusters.push_back(max_length_cluster);
    saved_cluster_points.at(max_length_cluster) = std::make_pair(temp_segments.at(max_length_cluster).special_A, temp_segments.at(max_length_cluster).special_B);
    remaining_segments.erase(max_length_cluster);

    // form a new segment with this cluster ...
    // good for a new segement ...
    temp_cluster->dijkstra_shortest_paths(cloud.pts[temp_segments.at(max_length_cluster).special_A],2);
    temp_cluster->cal_shortest_path(cloud.pts[temp_segments.at(max_length_cluster).special_B],2);
    WCPPID::ProtoSegment *tmp_sg = new WCPPID::ProtoSegment(-1, temp_cluster->get_path_wcps() ,temp_cluster->get_cluster_id());
    temp_segments_1.push_back(tmp_sg);
    temp_segments_2.push_back(tmp_sg);

    // examine the rest of points ...
    for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){
      //reset
      temp_segments.at(*it).number_not_faked = 0;
      temp_segments.at(*it).max_dis_u = 0;
      temp_segments.at(*it).max_dis_v = 0;
      temp_segments.at(*it).max_dis_w = 0;

      for (size_t j=0;j!=ncounts[*it];j++){
  	// also judge whether this track is fake ...
  	WCP::Point p(cloud.pts[sep_clusters[*it].at(j)].x, cloud.pts[sep_clusters[*it].at(j)].y, cloud.pts[sep_clusters[*it].at(j)].z);
  	double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
  	for (auto it1 = temp_segments_1.begin(); it1 != temp_segments_1.end(); it1++){
  	  WCPPID::ProtoSegment *sg1 = *it1;
  	  std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p);
  	  if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis);
  	  if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis);
  	  if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis);
  	}
  	int flag_num = 0;
  	if (min_dis_u > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 0))) flag_num ++;
  	if (min_dis_v > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 1))) flag_num ++;
  	if (min_dis_w > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 2))) flag_num ++;
  	if (flag_num>=2) temp_segments.at(*it).number_not_faked ++;
  	if (min_dis_u > temp_segments.at(*it).max_dis_u && (!ct_point_cloud->get_closest_dead_chs(p, 0))) temp_segments.at(*it).max_dis_u = min_dis_u;
  	if (min_dis_v > temp_segments.at(*it).max_dis_v && (!ct_point_cloud->get_closest_dead_chs(p, 1))) temp_segments.at(*it).max_dis_v = min_dis_v;
  	if (min_dis_w > temp_segments.at(*it).max_dis_w && (!ct_point_cloud->get_closest_dead_chs(p, 2))) temp_segments.at(*it).max_dis_w = min_dis_w;
      }

       if (temp_segments.at(*it).number_points ==1  //  only one point
       	|| temp_segments.at(*it).number_not_faked ==0 &&
       	   (temp_segments.at(*it).length < 3.5*units::cm  // very short & fake
       	    || (temp_segments.at(*it).number_not_faked < 0.25 * temp_segments.at(*it).number_points || temp_segments.at(*it).number_not_faked < 0.4 * temp_segments.at(*it).number_points && temp_segments.at(*it).length < 7 * units::cm) && temp_segments.at(*it).max_dis_u/units::cm < 3 && temp_segments.at(*it).max_dis_v/units::cm < 3 && temp_segments.at(*it).max_dis_w/units::cm < 3 && temp_segments.at(*it).max_dis_u + temp_segments.at(*it).max_dis_v + temp_segments.at(*it).max_dis_w < 6*units::cm)  // many fake things and very close to each other ...
       	)
       	 remaining_segments.erase(*it);
      
    }

    /* for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){ */
    /*   std::cout << "jaja1: " << *it << " " << temp_segments.at(*it).number_not_faked << " " << temp_segments.at(*it).number_points << " " <<  temp_segments.at(*it).max_dis_u/units::cm << " " << temp_segments.at(*it).max_dis_v/units::cm << " " << temp_segments.at(*it).max_dis_w/units::cm << " " << temp_segments.at(*it).length/units::cm << std::endl; */
    /* } */
    /* std::cout << std::endl; */
    //    break;
  }
  for (auto it = temp_segments_2.begin(); it != temp_segments_2.end(); it++){
    delete *it;
  }

  //  std::cout << "abc3 " << std::endl;
  
  bool flag_final_fit = true;
  
  std::vector<WCPPID::ProtoSegment*> new_segments;
  // save segments ...
  for (auto it = saved_clusters.begin(); it!= saved_clusters.end(); it++){
    //    std::cout << (saved_cluster_points.at(*it)).first << " " << (saved_cluster_points.at(*it)).second << " " << cloud.pts.size() << std::endl;
    // within the new segment, OK ...
    temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).first],2); 
    temp_cluster->cal_shortest_path(cloud.pts[(saved_cluster_points.at(*it)).second],2);


    
    std::list<WCP::WCPointCloud<double>::WCPoint> temp_path_list = temp_cluster->get_path_wcps();
    temp_cluster->collect_charge_trajectory(*ct_point_cloud);
    // do not fit dQ/dx at beginning ...
    temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false, false);

    //    std::cout << "abc3.2 " <<  " " << temp_cluster->get_fine_tracking_path().size() << std::endl; 
    
    if (temp_cluster->get_fine_tracking_path().size() >1){      
      WCPPID::ProtoVertex *v1 = find_vertex_other_segment(temp_cluster, true, cloud.pts[(saved_cluster_points.at(*it)).first]);
      //      std::cout << v1 << std::endl;
      WCPPID::ProtoVertex *v2 = find_vertex_other_segment(temp_cluster, false, cloud.pts[(saved_cluster_points.at(*it)).second]);
      //      std::cout << "abc3.3 " <<  v1 << " " << v2 << std::endl;
      // protection against a corner case crashing the code
      if (v1==v2){
      	double tmp_dis1 = sqrt(pow(v1->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).first].x, 2) + pow(v1->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).first].y,2) + pow(v1->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).first].z,2));
      	double tmp_dis2 = sqrt(pow(v1->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).second].x, 2) + pow(v1->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).second].y,2) + pow(v1->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).second].z,2));
      	if (tmp_dis1 < tmp_dis2){
      	  // change v2
      	  v2 = new WCPPID::ProtoVertex(acc_vertex_id, cloud.pts[(saved_cluster_points.at(*it)).second], temp_cluster->get_cluster_id()); acc_vertex_id++;
      	}else{
      	  // change v1
      	  v1 = new WCPPID::ProtoVertex(acc_vertex_id, cloud.pts[(saved_cluster_points.at(*it)).first], temp_cluster->get_cluster_id()); acc_vertex_id++;
      	}
      }
      
      //      std::cout << v1 << " " << v2 << " " << sqrt(pow(cloud.pts[(saved_cluster_points.at(*it)).first].x - cloud.pts[(saved_cluster_points.at(*it)).second].x,2) + pow(cloud.pts[(saved_cluster_points.at(*it)).first].y - cloud.pts[(saved_cluster_points.at(*it)).second].y,2) + pow(cloud.pts[(saved_cluster_points.at(*it)).first].z - cloud.pts[(saved_cluster_points.at(*it)).second].z,2)) << std::endl;
      
      if (v1->get_wcpt().index == cloud.pts[(saved_cluster_points.at(*it)).first].index && v2->get_wcpt().index == cloud.pts[(saved_cluster_points.at(*it)).second].index){
	// this is between v1 and v2 ...
      }else{	
	//std::cout << temp_path_list.size() << std::endl;
	if (v1->get_wcpt().index != cloud.pts[(saved_cluster_points.at(*it)).first].index){
	  double tmp_length = sqrt(pow(v1->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).first].x,2) + pow(v1->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).first].y,2) + pow(v1->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).first].z,2));
	  int tmp_count = std::round(tmp_length/(0.6*units::cm));
	  for (int k=1;k<tmp_count;k++){
	    Point tmp_p;
	    tmp_p.x = cloud.pts[(saved_cluster_points.at(*it)).first].x + (v1->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).first].x)*1./tmp_count*k;
	    tmp_p.y = cloud.pts[(saved_cluster_points.at(*it)).first].y + (v1->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).first].y)*1./tmp_count*k;
	    tmp_p.z = cloud.pts[(saved_cluster_points.at(*it)).first].z + (v1->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).first].z)*1./tmp_count*k;
	    
	    WCP::WCPointCloud<double>::WCPoint& tmp_wcp = temp_cluster->get_point_cloud_steiner()->get_closest_wcpoint(tmp_p);
	    if (tmp_wcp.index != temp_path_list.front().index) temp_path_list.push_front(tmp_wcp);
	  }
	  if (v1->get_wcpt().index != temp_path_list.front().index) temp_path_list.push_front(v1->get_wcpt());

	  
	}

	//for (auto it1 = temp_path_list.begin(); it1!=temp_path_list.end(); it1++){
	//   std::cout <<  " B: " << (*it1).index << std::endl;
	// }
	
	//std::cout << temp_path_list.size() << std::endl;
	
	if (v2->get_wcpt().index != cloud.pts[(saved_cluster_points.at(*it)).second].index){
	  double tmp_length = sqrt(pow(v2->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).second].x,2) + pow(v2->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).second].y,2) + pow(v2->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).second].z,2));
	  int tmp_count = std::round(tmp_length/(0.6*units::cm));
	  for (int k=1;k<tmp_count;k++){
	    Point tmp_p;
	    tmp_p.x = cloud.pts[(saved_cluster_points.at(*it)).second].x + (v2->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).second].x)*1./tmp_count*k;
	    tmp_p.y = cloud.pts[(saved_cluster_points.at(*it)).second].y + (v2->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).second].y)*1./tmp_count*k;
	    tmp_p.z = cloud.pts[(saved_cluster_points.at(*it)).second].z + (v2->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).second].z)*1./tmp_count*k;
	    
	    WCP::WCPointCloud<double>::WCPoint& tmp_wcp = temp_cluster->get_point_cloud_steiner()->get_closest_wcpoint(tmp_p);
	    if (tmp_wcp.index != temp_path_list.back().index) temp_path_list.push_back(tmp_wcp);
	  }
	  if (v2->get_wcpt().index != temp_path_list.back().index) temp_path_list.push_back(v2->get_wcpt());
	  
	  
	}

	/* std::cout << v1->get_wcpt().index << " " << temp_path_list.front().index << " " << temp_path_list.size() << " " << temp_path_list.back().index << " " << v2->get_wcpt().index << " " << temp_cluster->get_path_wcps().size() << " "<< temp_cluster->get_path_wcps().back().index << " " << temp_cluster->get_cluster_id() << " " << v2->get_cluster_id() << std::endl; */
	/* std::cout << cloud.pts[(saved_cluster_points.at(*it)).second].index << " " << cloud.pts[(saved_cluster_points.at(*it)).first].index << std::endl; */
	//temp_cluster->dijkstra_shortest_paths(v1->get_wcpt(),2);
	/* temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).first],2); */
	/* temp_cluster->cal_shortest_path(v2->get_wcpt(),2); */
	/* std::cout << " A: " << temp_cluster->get_path_wcps().size() << std::endl; */
	/* for (auto it1 = temp_cluster->get_path_wcps().begin(); it1!= temp_cluster->get_path_wcps().end();it1++){ */
	/*   std::cout << it1->index << std::endl; */
	/* } */
	/* temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).first],2); */
	/* temp_cluster->cal_shortest_path(cloud.pts[(saved_cluster_points.at(*it)).second],2); */
	/* std::cout << " A: " << temp_cluster->get_path_wcps().size() << std::endl; */
	/* for (auto it1 = temp_cluster->get_path_wcps().begin(); it1!= temp_cluster->get_path_wcps().end();it1++){ */
	/*   std::cout << it1->index << std::endl; */
	/* } */
	
	/* temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).second],2); */
	/* temp_cluster->cal_shortest_path(v2->get_wcpt(),2); */
	/* std::cout << " A: " << temp_cluster->get_path_wcps().size() << std::endl; */
	/* for (auto it1 = temp_cluster->get_path_wcps().begin(); it1!= temp_cluster->get_path_wcps().end();it1++){ */
	/*   std::cout << it1->index << std::endl; */
	/* } */
	
	//	std::cout << temp_path_list.size() << std::endl;
	temp_cluster->set_path_wcps(temp_path_list);
      }

      if (map_vertex_segments.find(v1)!=map_vertex_segments.end() || map_vertex_segments.find(v2) != map_vertex_segments.end()){
	WCPPID::ProtoSegment *sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
	//std::cout << v1 << " " << v2 << " " << sg1 << std::endl;
	add_proto_connection(v1,sg1,temp_cluster);
	add_proto_connection(v2,sg1,temp_cluster);

	//	std::cout << sg1->get_id() << std::endl;

	if (sg1->get_length() > 30*units::cm)	new_segments.push_back(sg1);

	// do dQ/dx fitting and exclude others ...
	temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
	flag_final_fit = false;
	// update output ...
	std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Other tracks -- # of Vertices: " << map_vertex_segments.size() << "; # of Segments: " << map_segment_vertices.size() << std::endl;
      }else{
	// hack

	TVector3 drift_dir(1,0,0);
	WCPPID::ProtoSegment* new_sg = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
	TVector3 dir(v1->get_wcpt().x - v2->get_wcpt().x, v1->get_wcpt().y - v2->get_wcpt().y, v1->get_wcpt().z - v2->get_wcpt().z);
	// Now need a parallel search ...
	bool flag_parallel = false;
	if (dir.Mag() > 10*units::cm){
	  for (auto it1 = map_vertex_segments.begin(); it1!=map_vertex_segments.end();it1++){
	    WCPPID::ProtoVertex *vtx = it1->first;
	    if (vtx == v1 || vtx == v2) continue;
	    TVector3 dir1(vtx->get_fit_pt().x - v1->get_wcpt().x, vtx->get_fit_pt().y - v1->get_wcpt().y, vtx->get_fit_pt().z - v1->get_wcpt().z);
	    TVector3 dir2(vtx->get_fit_pt().x - v2->get_wcpt().x, vtx->get_fit_pt().y - v2->get_wcpt().y, vtx->get_fit_pt().z - v2->get_wcpt().z);
	    // check one vertex
	    if (dir1.Mag() < 6*units::cm && fabs(drift_dir.Angle(dir1)/3.1415926*180.-90.)< 15. ) flag_parallel = modify_vertex_isochronous(vtx, v1, new_sg, v2, temp_cluster);
	    //check the other vertex
	    if ((!flag_parallel) && dir2.Mag() < 6*units::cm && fabs(drift_dir.Angle(dir2)/3.1415926*180.-90.) < 15. ) flag_parallel = modify_vertex_isochronous(vtx, v2, new_sg, v1, temp_cluster);
	    if (flag_parallel) break;
	  }
	  if (!flag_parallel){
	    for (auto it1 = map_segment_vertices.begin(); it1 != map_segment_vertices.end(); it1++){
	      WCPPID::ProtoSegment *sg1 = it1->first;
	      if (sg1->get_closest_point(v1->get_fit_pt()).first < 6*units::cm){
		flag_parallel = modify_segment_isochronous(sg1, v1, new_sg, v2, temp_cluster);
	      }
	      if ((!flag_parallel) && sg1->get_closest_point(v2->get_fit_pt()).first < 6*units::cm){
		flag_parallel = modify_segment_isochronous(sg1, v2, new_sg, v1, temp_cluster);
	      }
	      if (flag_parallel) break;
	    }
	  }
	}
	
	if (!flag_parallel){
	  residual_segment_candidates.push_back(std::make_tuple(temp_cluster, v1->get_wcpt().index, v2->get_wcpt().index));
	  std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Isolated residual segment found: " << sqrt(pow(v1->get_fit_pt().x - v2->get_fit_pt().x,2) + pow(v1->get_fit_pt().y - v2->get_fit_pt().y,2) + pow(v1->get_fit_pt().z - v2->get_fit_pt().z,2))/units::cm << " cm " << v1->get_fit_pt() << " " << v2->get_fit_pt() << std::endl;
	  delete v1;
	  delete v2;
	  delete new_sg;
	}else{
	  temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
	}
      }

      //hack ...
      //break;
    } // path > 1
  } // loop over saved stuff

  //  std::cout << "abc4 " << std::endl;
  
  if (flag_break_track)
    break_segments(new_segments, temp_cluster);

  //  if (flag_final_fit)
  // temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
}

WCPPID::ProtoVertex* WCPPID::NeutrinoID::find_vertex_other_segment(WCPPID::PR3DCluster* temp_cluster, bool flag_forward, WCP::WCPointCloud<double>::WCPoint& wcp){
  WCPPID::ProtoVertex *v1 = 0;
  
  std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, Point> check_results = check_end_point(temp_cluster,temp_cluster->get_fine_tracking_path(), flag_forward);
  //  std::cout << "A: " << flag_forward << " " << std::get<0>(check_results) << " " << std::get<1>(check_results) << " " << std::get<2>(check_results) << std::endl;
  if (std::get<0>(check_results)==0 && std::get<1>(check_results)==0){
    check_results = check_end_point(temp_cluster,temp_cluster->get_fine_tracking_path(), flag_forward, 1.2*units::cm, 2.5*units::cm);
    //std::cout << "B: " << flag_forward << " " <<  std::get<0>(check_results) << " " << std::get<1>(check_results) << " " << std::get<2>(check_results) << std::endl;
  }
  if (std::get<0>(check_results)==0 && std::get<1>(check_results)==0){
    check_results = check_end_point(temp_cluster, temp_cluster->get_fine_tracking_path(), flag_forward, 1.5*units::cm, 3.0*units::cm);
  }

  // std::cout << std::get<0>(check_results) << " " << std::get<1>(check_results) << " " << std::get<2>(check_results) << std::endl;
 // hack ...
 // std::get<0>(check_results) = 0;
 //std::get<1>(check_results) = 0;

 
 if (std::get<0>(check_results) != 0){
   v1 = std::get<0>(check_results);
 }else if (std::get<1>(check_results) != 0){
   
   WCPPID::ProtoSegment *sg1 = std::get<1>(check_results);
   Point test_p = std::get<2>(check_results);
   WCP::WCPointCloud<double>::WCPoint break_wcp = sg1->get_closest_wcpt(test_p);
   WCPPID::ProtoVertex *start_v = 0, *end_v = 0;
   for (auto it2 = map_segment_vertices[sg1].begin(); it2!=map_segment_vertices[sg1].end(); it2++){
     if ((*it2)->get_wcpt().index == sg1->get_wcpt_vec().front().index) start_v = *it2;
     if ((*it2)->get_wcpt().index == sg1->get_wcpt_vec().back().index) end_v = *it2;
   }
   if (start_v ==0 || end_v ==0) std::cout << "Error in finding vertices for a sgement" << std::endl;

 
   
   // std::cout << sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2)+pow(start_v->get_wcpt().y - break_wcp.y,2) + pow(start_v->get_wcpt().z - break_wcp.z,2))/units::cm << " " << sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2)+pow(end_v->get_wcpt().y - break_wcp.y,2) + pow(end_v->get_wcpt().z - break_wcp.z,2))/units::cm << std::endl;
   if (sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2)+pow(start_v->get_wcpt().y - break_wcp.y,2) + pow(start_v->get_wcpt().z - break_wcp.z,2)) < 0.9*units::cm){
     v1 = start_v;
   }else if (sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2)+pow(end_v->get_wcpt().y - break_wcp.y,2) + pow(end_v->get_wcpt().z - break_wcp.z,2)) < 0.9*units::cm){
     v1 = end_v;
   }else{
     std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list1;
     std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list2;
     //std::cout << start_v->get_wcpt().index << " " << break_wcp.index << " " << end_v->get_wcpt().index << std::endl;
     temp_cluster->proto_break_tracks(start_v->get_wcpt(), break_wcp, end_v->get_wcpt(), wcps_list1, wcps_list2, true);
     //     std::cout << "xin " << wcps_list1.size() << " " << wcps_list2.size() << std::endl;
   
     
     v1 = new WCPPID::ProtoVertex(acc_vertex_id, break_wcp, temp_cluster->get_cluster_id()); acc_vertex_id++;
     WCPPID::ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list1, temp_cluster->get_cluster_id()); acc_segment_id++;
     WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list2, temp_cluster->get_cluster_id()); acc_segment_id++;

     //     std::cout << "Other: " << sg1->get_id() << " " << sg2->get_id() << " " << sg3->get_id() << std::endl;
     //	  std::cout << start_v->get_wcpt().index << " " << break_wcp.index << " " << end_v->get_wcpt().index << " " << wcps_list1.size() << " " << wcps_list2.size() << std::endl;
     //std::cout << map_vertex_segments.size() << " " << map_segment_vertices.size() << std::endl;
     add_proto_connection(start_v, sg2, temp_cluster);
     add_proto_connection(v1, sg2, temp_cluster);
     add_proto_connection(v1, sg3, temp_cluster);
     add_proto_connection(end_v, sg3, temp_cluster);
     del_proto_segment(sg1);
   } 
 }else{
   v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcp, temp_cluster->get_cluster_id()); acc_vertex_id++;
 }
 
 
 return v1;
}

bool WCPPID::NeutrinoID::modify_segment_isochronous(WCPPID::ProtoSegment* sg1, WCPPID::ProtoVertex *v1, WCPPID::ProtoSegment* sg, WCPPID::ProtoVertex *v2, WCPPID::PR3DCluster* temp_cluster){
  bool flag = false;

  TVector3 dir1 = sg->cal_dir_3vector(v1->get_fit_pt(), 15*units::cm) * (-1);
  if (dir1.X()==0) return flag;

  PointVector& pts = sg1->get_point_vec();
  WCP::ToyPointCloud *pcloud = temp_cluster->get_point_cloud_steiner();
  WCP::WCPointCloud<double>::WCPoint vtx_new_wcp;
  
  
  TVector3 drift_dir(1,0,0);
  for (size_t i=0;i!=pts.size();i++){
    TVector3 dir(pts.at(i).x - v1->get_fit_pt().x, pts.at(i).y - v1->get_fit_pt().y, pts.at(i).z - v1->get_fit_pt().z);

    if (dir.Mag() > 6*units::cm || fabs(drift_dir.Angle(dir)/3.1415926*180.-90.)>= 15.) continue;
    
    Point test_p = v1->get_fit_pt();
    test_p.x += pts.at(i).x - v1->get_fit_pt().x;
    test_p.y += dir1.Y()/dir1.X() *(pts.at(i).x - v1->get_fit_pt().x);
    test_p.z += dir1.Z()/dir1.X() *(pts.at(i).x - v1->get_fit_pt().x);
    vtx_new_wcp = pcloud->get_closest_wcpoint(test_p);
    

    double tmp_dis = sqrt(pow(vtx_new_wcp.x - v1->get_fit_pt().x,2)+pow(vtx_new_wcp.y - v1->get_fit_pt().y,2)+pow(vtx_new_wcp.z - v1->get_fit_pt().z,2));

    if (tmp_dis > 5*units::cm) continue;
  
    // check connectivity between vtx_new_wcp and the vtx one ...
    double step_size = 0.6*units::cm;
    Point start_p = pts.at(i);
    Point end_p(vtx_new_wcp.x, vtx_new_wcp.y, vtx_new_wcp.z);
    int ncount = std::round(sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/step_size);
    int n_bad = 0;
    for (int i=1;i<ncount;i++){
      Point test_p(start_p.x + (end_p.x - start_p.x)/ncount*i,
		   start_p.y + (end_p.y - start_p.y)/ncount*i,
		   start_p.z + (end_p.z - start_p.z)/ncount*i);
      if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
    }
    
    if (n_bad ==0){
      flag = true;
      break;
    }
  }
  
  
  if (!flag) return flag;

  // modify the new segment ...
  v1->set_wcpt(vtx_new_wcp);
  temp_cluster->dijkstra_shortest_paths(vtx_new_wcp, 2);
  temp_cluster->cal_shortest_path(v2->get_wcpt(), 2);
  {
    std::vector<WCP::WCPointCloud<double>::WCPoint >& wcpt_vec = sg->get_wcpt_vec();
    wcpt_vec.clear();
    wcpt_vec.reserve(temp_cluster->get_path_wcps().size());
    std::copy(std::begin(temp_cluster->get_path_wcps()), std::end(temp_cluster->get_path_wcps()), std::back_inserter(wcpt_vec));
  }
  add_proto_connection(v1, sg, temp_cluster);
  add_proto_connection(v2, sg,  temp_cluster);

  // create two new segments
  auto pair_vertices = find_vertices(sg1);
  {
    temp_cluster->cal_shortest_path(pair_vertices.first->get_wcpt(),2); 
    ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;  
  
    add_proto_connection(pair_vertices.first, sg2, temp_cluster); 
    add_proto_connection(v1, sg2, temp_cluster);
  }
  {
    temp_cluster->cal_shortest_path(pair_vertices.second->get_wcpt(),2); 
    ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;  
  
    add_proto_connection(pair_vertices.second, sg2, temp_cluster); 
    add_proto_connection(v1, sg2, temp_cluster);
  }
  del_proto_segment(sg1);
  
  return flag;
}


bool WCPPID::NeutrinoID::modify_vertex_isochronous(WCPPID::ProtoVertex* vtx, WCPPID::ProtoVertex *v1, WCPPID::ProtoSegment* sg, WCPPID::ProtoVertex *v2, WCPPID::PR3DCluster* temp_cluster){
  bool flag = false;

  TVector3 dir = sg->cal_dir_3vector(v1->get_fit_pt(), 15*units::cm) * (-1);
  if (dir.X()==0) return flag;
  
  Point test_p = v1->get_fit_pt();
  test_p.x += vtx->get_fit_pt().x - v1->get_fit_pt().x;
  test_p.y += dir.Y()/dir.X() *(vtx->get_fit_pt().x - v1->get_fit_pt().x);
  test_p.z += dir.Z()/dir.X() *(vtx->get_fit_pt().x - v1->get_fit_pt().x);

  WCP::ToyPointCloud *pcloud = temp_cluster->get_point_cloud_steiner();
  WCP::WCPointCloud<double>::WCPoint& vtx_new_wcp = pcloud->get_closest_wcpoint(test_p);

  double tmp_dis = sqrt(pow(vtx_new_wcp.x - v1->get_fit_pt().x,2)+pow(vtx_new_wcp.y - v1->get_fit_pt().y,2)+pow(vtx_new_wcp.z - v1->get_fit_pt().z,2));

  if (tmp_dis > 5*units::cm) return flag;
  
  // check connectivity between vtx_new_wcp and the vtx one ...
  double step_size = 0.6*units::cm;
  Point start_p = vtx->get_fit_pt(); 
  Point end_p(vtx_new_wcp.x, vtx_new_wcp.y, vtx_new_wcp.z);
  int ncount = std::round(sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/step_size);
  int n_bad = 0;
  for (int i=1;i<ncount;i++){
    Point test_p(start_p.x + (end_p.x - start_p.x)/ncount*i,
		 start_p.y + (end_p.y - start_p.y)/ncount*i,
		 start_p.z + (end_p.z - start_p.z)/ncount*i);
    if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
  }
  if (n_bad >0) return flag;
  //  std::cout << tmp_dis/units::cm << " " << n_bad << std::endl;

  
  // Now Update everything ...
  for (auto it = map_vertex_segments[vtx].begin(); it!= map_vertex_segments[vtx].end(); it++){
    WCPPID::ProtoSegment *sg1 = *it;
    std::vector<WCP::WCPointCloud<double>::WCPoint >& wcpt_vec = sg1->get_wcpt_vec();
      
    WCP::WCPointCloud<double>::WCPoint vtx_other_wcp;
    if (wcpt_vec.front().index == vtx->get_wcpt().index){
      vtx_other_wcp = wcpt_vec.back();
    }else{
      vtx_other_wcp = wcpt_vec.front();
    }
    temp_cluster->dijkstra_shortest_paths(vtx_new_wcp, 2);
    temp_cluster->cal_shortest_path(vtx_other_wcp, 2);
    //  std::cout << temp_cluster->get_path_wcps().front().index << " " << temp_cluster->get_path_wcps().back().index << " " << vtx_new_wcp.index << " " << vtx_other_wcp.index << std::endl;
    wcpt_vec.clear();
    wcpt_vec.reserve(temp_cluster->get_path_wcps().size());
    std::copy(std::begin(temp_cluster->get_path_wcps()), std::end(temp_cluster->get_path_wcps()), std::back_inserter(wcpt_vec));
    //    std::cout << sg1->get_id() << " " << sg->get_id() << std::endl;
  }
  vtx->set_wcpt(vtx_new_wcp);
  
  temp_cluster->dijkstra_shortest_paths(vtx_new_wcp, 2);
  temp_cluster->cal_shortest_path(v2->get_wcpt(), 2);
  {
    std::vector<WCP::WCPointCloud<double>::WCPoint >& wcpt_vec = sg->get_wcpt_vec();
    wcpt_vec.clear();
    wcpt_vec.reserve(temp_cluster->get_path_wcps().size());
    std::copy(std::begin(temp_cluster->get_path_wcps()), std::end(temp_cluster->get_path_wcps()), std::back_inserter(wcpt_vec));
  }
  delete v1;
  add_proto_connection(vtx, sg, temp_cluster);
  add_proto_connection(v2, sg,  temp_cluster);
  flag = true;
  std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " shift a isochronous vertex with adding a segment " << sg->get_id() << std::endl;
  
  return flag;
}
    



std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, WCP::Point> WCPPID::NeutrinoID::check_end_point(WCPPID::PR3DCluster* temp_cluster, WCP::PointVector& tracking_path, bool flag_front, double vtx_cut1, double vtx_cut2, double sg_cut1, double sg_cut2 ){

  if (tracking_path.size()<2) std::cout << "vector size wrong!" << std::endl;
  int ncount = 5;
  Point test_p;
  if (flag_front){
    test_p = tracking_path.front();
  }else{
    test_p = tracking_path.back();
  }
  WCPPID::ProtoVertex *vtx=0;
  WCPPID::ProtoSegment *seg=0;
  
  Point dir_p(0,0,0);
  for (int i=0;i!=ncount;i++){
    if (flag_front){
      if (i+1<tracking_path.size()){
	dir_p.x = test_p.x - tracking_path.at(i+1).x;
	dir_p.y = test_p.y - tracking_path.at(i+1).y;
	dir_p.z = test_p.z - tracking_path.at(i+1).z;
      }else{
	dir_p.x = test_p.x - tracking_path.back().x;
	dir_p.y = test_p.y - tracking_path.back().y;
	dir_p.z = test_p.z - tracking_path.back().z;
      }
    }else{
      if (tracking_path.size() >= i+2){
	dir_p.x = test_p.x - tracking_path.at(tracking_path.size()-2-i).x;
	dir_p.y = test_p.y - tracking_path.at(tracking_path.size()-2-i).y;
	dir_p.z = test_p.z - tracking_path.at(tracking_path.size()-2-i).z;
      }else{
	dir_p.x = test_p.x - tracking_path.front().x;
	dir_p.x = test_p.y - tracking_path.front().y;
	dir_p.x = test_p.z - tracking_path.front().z;
      }
    }
    TVector3 temp_dir(dir_p.x, dir_p.y, dir_p.z);
    if (temp_dir.Mag() == 0 ) continue;
    
    temp_dir = temp_dir.Unit();
    Line temp_l(test_p, temp_dir);

    // check vertex
    for (auto it1 = map_vertex_segments.begin(); it1!= map_vertex_segments.end(); it1++){
      WCPPID::ProtoVertex *test_v = it1->first;
      if (test_v->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
      Point p1(test_v->get_wcpt().x, test_v->get_wcpt().y, test_v->get_wcpt().z);
      Point p2 = test_v->get_fit_pt();
      
      double dis1 = sqrt(pow(test_p.x - p1.x,2)
			 + pow(test_p.y - p1.y,2)
			 + pow(test_p.z - p1.z,2));
      double dis2 = temp_l.closest_dis(p1);
      double dis3 = sqrt(pow(test_p.x - p2.x,2)
			 + pow(test_p.y - p2.y,2)
			 + pow(test_p.z - p2.z,2));
      double dis4 = temp_l.closest_dis(p2);

      //      if (dis1 < 5*units::cm && dis2 < 5*units::cm)
      //	std::cout << i << " " << test_p << " " << temp_dir.X() << " " << temp_dir.Y() << " " << temp_dir.Z() << " " << dis1/units::cm << " " << dis2/units::cm << " " << dis3/units::cm << " " << dis4/units::cm << " " << it1->second.size() << " "  << std::endl;

      
      
      if ( std::max(dis1,dis2) < 5*units::cm && (std::min(dis1,dis2) < vtx_cut1 && std::max(dis1,dis2) < vtx_cut2 ||
						 std::min(dis3,dis4) < vtx_cut1 * 1.3 && std::max(dis1,dis2) < vtx_cut2 * 2 && it1->second.size()==1) ){
	TVector3 test_dir(p1.x-test_p.x,p1.y-test_p.y,p1.z-test_p.z);
	//	std::cout << test_dir.Angle(temp_dir)/3.1415926*180. << std::endl;
	if (test_dir.Angle(temp_dir)/3.1415926*180. < 90){
	  vtx = test_v;
	  break;
	}
      }
      if ( std::max(dis3,dis4) < 5*units::cm && (std::min(dis3,dis4) < vtx_cut1 && std::max(dis3,dis4) < vtx_cut2 ||
						 std::min(dis3,dis4) < vtx_cut1 * 1.3 && std::max(dis3,dis4) < vtx_cut2 * 2 && it1->second.size()==1 )){
	TVector3 test_dir(p2.x-test_p.x,p2.y-test_p.y,p2.z-test_p.z);
	//	std::cout << test_dir.Angle(temp_dir)/3.1415926*180. << std::endl;
	if (test_dir.Angle(temp_dir)/3.1415926*180. < 90){
	  vtx = test_v;
	  break;
	}
      }
    }

    if (vtx!=0) break;
  }
   
  // check segment 
  if (vtx == 0){
    WCPPID::ProtoSegment* min_sg=0;
    Point min_point(0,0,0);
    double min_dis = 1e9;
      
    for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
      if (it->first->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
      std::pair<double, WCP::Point> results = it->first->get_closest_point(test_p);
      if (results.first < sg_cut1){
	PointVector& points = it->first->get_point_vec();

	Point dir_p(0,0,0);
	for (int i=0;i!=ncount;i++){
	  if (flag_front){
	    if (i+1<tracking_path.size()){
	      dir_p.x = test_p.x - tracking_path.at(i+1).x;
	      dir_p.y = test_p.y - tracking_path.at(i+1).y;
	      dir_p.z = test_p.z - tracking_path.at(i+1).z;
	    }else{
	      dir_p.x = test_p.x - tracking_path.back().x;
	      dir_p.y = test_p.y - tracking_path.back().y;
	      dir_p.z = test_p.z - tracking_path.back().z;
	    }
	  }else{
	    if (tracking_path.size() >= i+2){
	      dir_p.x = test_p.x - tracking_path.at(tracking_path.size()-2-i).x;
	      dir_p.y = test_p.y - tracking_path.at(tracking_path.size()-2-i).y;
	      dir_p.z = test_p.z - tracking_path.at(tracking_path.size()-2-i).z;
	    }else{
	      dir_p.x = test_p.x - tracking_path.front().x;
	      dir_p.x = test_p.y - tracking_path.front().y;
	      dir_p.x = test_p.z - tracking_path.front().z;
	    }
	  }
	  TVector3 temp_dir(dir_p.x, dir_p.y, dir_p.z);
	  if (temp_dir.Mag()==0) continue;
	  temp_dir = temp_dir.Unit();
	  Line temp_l(test_p, temp_dir);
	
	  for (size_t i= 0; i!=points.size();i++){
	    double dis1 = temp_l.closest_dis(points.at(i) );
	    if (dis1 < sg_cut2 ){
	      TVector3 test_dir(points.at(i).x-test_p.x,points.at(i).y-test_p.y,points.at(i).z-test_p.z);
	      //	      std::cout << test_dir.Angle(temp_dir)/3.1415926*180. << " " << test_dir.Mag()/units::cm << " " << dis1/units::cm << std::endl;
	      dis1 = sqrt(pow(dis1,2)+pow(test_dir.Mag()/4.,2));
	      if (test_dir.Angle(temp_dir)/3.1415926*180. < 90 && dis1 < min_dis){
		min_dis = dis1;
		min_point = points.at(i);
		min_sg = it->first;
	      }
	    }
	  }
	}
      }
    }
    
    if (min_sg!=0){
      seg = min_sg;
      test_p = min_point;
    }
  }
    
  
  // return ...
  return std::make_tuple(vtx, seg, test_p);
}



bool WCPPID::NeutrinoID::add_proto_connection(WCPPID::ProtoVertex *pv, WCPPID::ProtoSegment *ps, WCPPID::PR3DCluster* cluster){

  if (pv->get_wcpt().index != ps->get_wcpt_vec().front().index && pv->get_wcpt().index != ps->get_wcpt_vec().back().index){
    std::cout << "Error! Vertex and Segment does not match " << pv->get_wcpt().index << " " << ps->get_wcpt_vec().front().index << " " << ps->get_wcpt_vec().back().index << std::endl;
    return false;
  }

  if (map_vertex_cluster.find(pv)==map_vertex_cluster.end())
    proto_vertices.push_back(pv);

  map_vertex_cluster[pv] = cluster;
  if (map_cluster_vertices.find(cluster)!=map_cluster_vertices.end()){
    map_cluster_vertices[cluster].insert(pv);
  }else{
    ProtoVertexSet temp_vertices;
    temp_vertices.insert(pv);
    map_cluster_vertices[cluster] = temp_vertices;
  }

  if (map_segment_cluster.find(ps)==map_segment_cluster.end())
    proto_segments.push_back(ps);
  
  map_segment_cluster[ps] = cluster;
  if (map_cluster_segments.find(cluster)!=map_cluster_segments.end()){
    map_cluster_segments[cluster].insert(ps);
  }else{
    ProtoSegmentSet temp_segments;
    temp_segments.insert(ps);
    map_cluster_segments[cluster] = temp_segments;
  }

  if (map_vertex_segments.find(pv)!=map_vertex_segments.end()){
    map_vertex_segments[pv].insert(ps);
  }else{
    ProtoSegmentSet temp_segments;
    temp_segments.insert(ps);
    map_vertex_segments[pv] = temp_segments;
  }

  if (map_segment_vertices.find(ps)!=map_segment_vertices.end()){
    map_segment_vertices[ps].insert(pv);
  }else{
    ProtoVertexSet temp_vertices;
    map_segment_vertices[ps] = temp_vertices;
    map_segment_vertices[ps].insert(pv);
  }

  return true;
}


bool WCPPID::NeutrinoID::del_proto_vertex(WCPPID::ProtoVertex *pv){
  if (map_vertex_segments.find(pv)!=map_vertex_segments.end()){
    // overall
    //    proto_vertices.erase(proto_vertices.find(pv));
    // map between vertex and segment
    for (auto it = map_vertex_segments[pv].begin(); it!= map_vertex_segments[pv].end(); it++){ // loop very semeng *it
      map_segment_vertices[*it].erase(map_segment_vertices[*it].find(pv));
    }
    map_vertex_segments.erase(pv);
    // map between cluster and vertices
    map_cluster_vertices[map_vertex_cluster[pv]].erase(map_cluster_vertices[map_vertex_cluster[pv]].find(pv));
    map_vertex_cluster.erase(pv);
    
    return true;
  }else{
    return false;
  }
}

bool WCPPID::NeutrinoID::del_proto_segment(WCPPID::ProtoSegment *ps){
  if (map_segment_vertices.find(ps)!=map_segment_vertices.end()){
    // overall
    //    proto_segments.erase(proto_segments.find(ps));
    //map between segment and vertex
    for (auto it = map_segment_vertices[ps].begin(); it!=map_segment_vertices[ps].end(); it++){
      map_vertex_segments[*it].erase(map_vertex_segments[*it].find(ps));
    }
    map_segment_vertices.erase(ps);
    //map between cluster and segment
    map_cluster_segments[map_segment_cluster[ps]].erase(map_cluster_segments[map_segment_cluster[ps]].find(ps));
    map_segment_cluster.erase(ps);
    return true;
  }else{
    return false;
  }
  
}

void WCPPID::NeutrinoID::examine_segment(WCPPID::PR3DCluster* temp_cluster){
  std::set<WCPPID::ProtoSegment*> segments_to_be_removed;
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    std::vector<WCPPID::ProtoSegment*> tmp_segments(it->second.begin(), it->second.end());
    for (size_t i=0;i!=tmp_segments.size();i++){
      auto it1 = find_vertices(tmp_segments.at(i));
      for (size_t j=i+1;j<tmp_segments.size();j++){
	auto it2 = find_vertices(tmp_segments.at(j));
	if (it1.first->get_wcpt().index == it2.first->get_wcpt().index &&
	    it1.second->get_wcpt().index == it2.second->get_wcpt().index ||
	    it1.first->get_wcpt().index == it2.second->get_wcpt().index &&
	    it1.second->get_wcpt().index == it2.first->get_wcpt().index){
	  segments_to_be_removed.insert(tmp_segments.at(j));
	}
      }
    }
  }
  for (auto it = segments_to_be_removed.begin(); it != segments_to_be_removed.end(); it++){
    del_proto_segment(*it);
  }
  
 
}

void WCPPID::NeutrinoID::examine_vertices(WCPPID::PR3DCluster* temp_cluster){

  bool flag_continue = true;
  while(flag_continue){
    flag_continue = false;
    
    examine_segment(temp_cluster);
  
  
    // merge vertex if the kink is not at right location
    flag_continue = flag_continue || examine_vertices_1(temp_cluster);
  
    //  std::cout << "haha2 " << std::endl;
    if (find_vertices(temp_cluster).size() > 2)
      // merge vertices if they are too close ...
      flag_continue = flag_continue || examine_vertices_2(temp_cluster);

    flag_continue = flag_continue || examine_vertices_4(temp_cluster);
  }

}





// examine if the initial vertices are not at the extreme location ...
void WCPPID::NeutrinoID::examine_vertices_3(){
  // examine main_cluster_initial_pair_vertices ...
  std::vector<WCPPID::ProtoVertex* > temp_vertices;
  if (main_cluster_initial_pair_vertices.first!=0) temp_vertices.push_back(main_cluster_initial_pair_vertices.first);
  if (main_cluster_initial_pair_vertices.second!=0) temp_vertices.push_back(main_cluster_initial_pair_vertices.second);
  bool flag_refit = false;
  for (size_t i=0;i!=temp_vertices.size();i++){
    WCPPID::ProtoVertex *vtx = temp_vertices.at(i);
    auto it = map_vertex_segments.find(vtx);
    if (it == map_vertex_segments.end()) continue;
    if (it->second.size()>1) continue; // more than one track ...
    WCPPID::ProtoSegment *sg = *(it->second.begin());

    bool flag_start;
    if (sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
      flag_start = true;
    else if (sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
      flag_start = false;

    std::vector<WCP::WCPointCloud<double>::WCPoint >& wcps = sg->get_wcpt_vec();
    
    WCP::WCPointCloud<double>::WCPoint wcp2;
    if (flag_start) wcp2 = wcps.back();
    else wcp2 = wcps.front();
    
    auto wcp1 = main_cluster->get_local_extension(vtx->get_wcpt(),2);

    if (wcp1.index == vtx->get_wcpt().index) continue;
    
    if (flag_start){
      main_cluster->dijkstra_shortest_paths(wcp1,2);  
      main_cluster->cal_shortest_path(wcp2,2);
    }else{
      main_cluster->dijkstra_shortest_paths(wcp2,2);  
      main_cluster->cal_shortest_path(wcp1,2);
    }
    if (main_cluster->get_path_wcps().size() < wcps.size()*2){
      vtx->set_wcpt(wcp1);
      wcps.clear();
      for (auto it = main_cluster->get_path_wcps().begin(); it!=main_cluster->get_path_wcps().end(); it++){
	wcps.push_back(*it);
      }
      flag_refit = true;
    }
    
    //    std::cout << wcps.size() << " " << main_cluster->get_path_wcps().size() << std::endl;
  }
  if (flag_refit)
    main_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);

  //  std::cout << "haha " << std::endl;
  
  
  std::set<WCPPID::ProtoSegment*> segments_to_be_removed;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != main_cluster->get_cluster_id()) continue;
    auto pair_vertices = find_vertices(sg);
    //    std::cout << sg->get_id() << " " << sg->get_direct_length()/units::cm << std::endl;
    if ((map_vertex_segments[pair_vertices.first].size()==1 ||
  	 map_vertex_segments[pair_vertices.second].size()==1)&&sg->get_direct_length() < 5*units::cm){
      PointVector& pts = sg->get_point_vec();
      int num_unique = 0;
      for (size_t i=0;i!=pts.size();i++){
  	double min_u = 1e9;
  	double min_v = 1e9;
  	double min_w = 1e9;
  	for (auto it1 = map_segment_vertices.begin(); it1!= map_segment_vertices.end(); it1++){
  	  WCPPID::ProtoSegment *sg1 = it1->first;
  	  if (sg == sg1) continue;
  	  std::tuple<double, double, double> results = sg1->get_closest_2d_dis(pts.at(i));
  	  if (std::get<0>(results) < min_u) min_u = std::get<0>(results);
  	  if (std::get<1>(results) < min_v) min_v = std::get<1>(results);
  	  if (std::get<2>(results) < min_w) min_w = std::get<2>(results);
  	}
  	if (min_u > 0.6*units::cm || min_v > 0.6*units::cm || min_w > 0.6*units::cm)
  	  num_unique ++;
  	//	if (min_u > 0.6*units::cm
	//	std::cout << sg->get_id() << " " << i << " " << min_u/units::cm << " " << min_v/units::cm << " " << min_w/units::cm << std::endl;
      }
      if (num_unique == 0){
  	segments_to_be_removed.insert(sg);
      }
    }
  }

  std::set<WCPPID::ProtoVertex*> can_vertices;
  
  for (auto it = segments_to_be_removed.begin(); it != segments_to_be_removed.end(); it++){
    auto pair_vertices = find_vertices(*it);
    if (map_vertex_segments[pair_vertices.first].size()>1) can_vertices.insert(pair_vertices.first);
    if (map_vertex_segments[pair_vertices.second].size()>1) can_vertices.insert(pair_vertices.second);
    del_proto_segment(*it);
  }
  for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    if (it->second.size()==0)
      del_proto_vertex(it->first);
  }
  //std::cout << can_vertices.size() << std::endl;
  for (auto it = can_vertices.begin(); it!=can_vertices.end();it++){
    examine_structure_4(*it, main_cluster);
  }
  
  if (segments_to_be_removed.size()>0)
    main_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
  
  
  
}

bool WCPPID::NeutrinoID::examine_vertices_2(WCPPID::PR3DCluster* temp_cluster){
  // search for merge ...
  bool flag_continue = false;
  
  //std::map<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> map_vertex_replace;
  WCPPID::ProtoVertex *v1 = 0;
  WCPPID::ProtoVertex *v2 = 0;
  WCPPID::ProtoSegment *sg = 0;
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    auto tmp_vertices = find_vertices(sg);
    v1 = tmp_vertices.first;
    v2 = tmp_vertices.second;
    if (v1 !=0 &&  v2 !=0){
      double dis = sqrt(pow(v1->get_fit_pt().x-v2->get_fit_pt().x,2) + pow(v1->get_fit_pt().y-v2->get_fit_pt().y,2) + pow(v1->get_fit_pt().z-v2->get_fit_pt().z,2));
      
      //std::cout << sg->get_cluster_id() << " " << sg->get_id() << " " << dis/units::cm << std::endl;
      if ( dis < 0.45*units::cm){
	flag_continue = true;
	break;
      }else{
	v1 = 0;
	v2 = 0;
	sg = 0;
      }
    }
  }
  
  // merge and clean up map ...
  if (v1!=0 && v2!=0){
    std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Merge Vertices Type II" << std::endl;
    // delete the segment in between
    del_proto_segment(sg);
    WCPPID::ProtoSegmentSelection tmp_segments;
    for (auto it1 = map_vertex_segments[v2].begin(); it1 != map_vertex_segments[v2].end(); it1++){
      tmp_segments.push_back(*it1);
    }
    for (auto it1 = tmp_segments.begin(); it1!=tmp_segments.end(); it1++){
      WCPPID::ProtoVertex *v3 = find_other_vertex(*it1, v2);
      temp_cluster->dijkstra_shortest_paths(v3->get_wcpt(),2);  
      temp_cluster->cal_shortest_path(v1->get_wcpt(),2);
      ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++; 
      del_proto_segment(*it1);
      add_proto_connection(v3, sg2, temp_cluster);
      add_proto_connection(v1, sg2, temp_cluster);
    }
    del_proto_vertex(v2); 
    
    WCPPID::ProtoVertexSelection tmp_vertices;
    for (auto it1 = map_vertex_segments.begin(); it1!=map_vertex_segments.end(); it1++){
      if (it1->second.size()==0) tmp_vertices.push_back(it1->first);
    }
    for (auto it1 = tmp_vertices.begin(); it1!=tmp_vertices.end(); it1++){
      del_proto_vertex(*it1);
    }
    
    temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
  }
  
  return flag_continue;
  
}

bool WCPPID::NeutrinoID::examine_vertices_4(WCPPID::PR3DCluster *temp_cluster){
  bool flag_continue = false;
  WCP::ToyPointCloud *pcloud = temp_cluster->get_point_cloud_steiner();
 
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    if (sg->get_direct_length() < 2.0*units::cm){
      auto pair_vertices = find_vertices(sg);
      WCPPID::ProtoVertex *v1 = pair_vertices.first;
      WCPPID::ProtoVertex *v2 = pair_vertices.second;
      //if (map_vertex_segments[v1].size() <= 2 || map_vertex_segments[v2].size() <=2) continue;
      
      if (map_vertex_segments[v1].size()>=2 && examine_vertices_4(v1, v2) ){

	for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	  ProtoSegment *sg1 = *it1;
	  if (sg1 == sg) continue;
	  
	  WCP::WCPointCloud<double>::WCPoint& vtx_new_wcp = v2->get_wcpt();
	  
	  std::vector<WCP::WCPointCloud<double>::WCPoint >& vec_wcps = sg1->get_wcpt_vec();
	  bool flag_front = false;
	  if (vec_wcps.front().index == v1->get_wcpt().index) flag_front = true;

	  WCP::WCPointCloud<double>::WCPoint min_wcp;
	  double min_dis = 1e9;
	  double max_dis = std::max(sqrt(pow(vec_wcps.front().x -v1->get_fit_pt().x,2)+pow(vec_wcps.front().y -v1->get_fit_pt().y,2)+pow(vec_wcps.front().z -v1->get_fit_pt().z,2)), sqrt(pow(vec_wcps.back().x -v1->get_fit_pt().x,2)+pow(vec_wcps.back().y -v1->get_fit_pt().y,2)+pow(vec_wcps.back().z -v1->get_fit_pt().z,2)));
	  double dis_cut = 0;
	  double default_dis_cut = 2.5*units::cm; 
	  if (max_dis > 2 * default_dis_cut) dis_cut = default_dis_cut;
	  for (size_t j=0;j!=vec_wcps.size();j++){
	    double dis1 = sqrt(pow(vec_wcps.at(j).x -v1->get_fit_pt().x,2)+pow(vec_wcps.at(j).y -v1->get_fit_pt().y,2)+pow(vec_wcps.at(j).z -v1->get_fit_pt().z,2));
	    double dis = fabs(dis1-3*units::cm);
	    if (dis < min_dis && dis1 > dis_cut){
	      min_wcp = vec_wcps.at(j);
	      min_dis = dis;
	    }
	  }
	    // establish the shortest path ...
	  std::list<WCP::WCPointCloud<double>::WCPoint> new_list;
	  
	  new_list.push_back(vtx_new_wcp);
	  {
	    double dis_step = 2.0*units::cm;
	    int ncount = std::round(sqrt(pow(vtx_new_wcp.x - min_wcp.x,2) + pow(vtx_new_wcp.y - min_wcp.y,2) + pow(vtx_new_wcp.z - min_wcp.z,2))/dis_step);
	    if (ncount <2) ncount = 2;
	    
	    for (int qx=1;qx<ncount;qx++){
	      Point tmp_p(vtx_new_wcp.x + (min_wcp.x-vtx_new_wcp.x)/ncount*qx, vtx_new_wcp.y + (min_wcp.y-vtx_new_wcp.y)/ncount*qx, vtx_new_wcp.z + (min_wcp.z-vtx_new_wcp.z)/ncount*qx);
	      WCP::WCPointCloud<double>::WCPoint& tmp_wcp = pcloud->get_closest_wcpoint(tmp_p);
	      //std::cout << qx << " " << sqrt(pow(tmp_wcp.x - tmp_p.x,2) + pow(tmp_wcp.y - tmp_p.y,2) + pow(tmp_wcp.z - tmp_p.z,2))/units::cm << std::endl;
	      if (sqrt(pow(tmp_wcp.x - tmp_p.x,2) + pow(tmp_wcp.y - tmp_p.y,2) + pow(tmp_wcp.z - tmp_p.z,2)) > 0.3*units::cm) continue; // too large ...
	      if (tmp_wcp.index != new_list.back().index && tmp_wcp.index != min_wcp.index)
		new_list.push_back(tmp_wcp);
	    }
	  }
	  new_list.push_back(min_wcp);
	  
	  std::list<WCP::WCPointCloud<double>::WCPoint> old_list;
	  std::copy(vec_wcps.begin(), vec_wcps.end(), std::back_inserter(old_list));
	  // append the list back ...
	  //std::cout << min_wcp.index << " " << temp_cluster->get_path_wcps().size() << std::endl;
	  if (flag_front){
	    //std::cout << "f1 " << old_list.size() << std::endl;
	    
	    while(old_list.front().index != min_wcp.index && old_list.size()>0 ){
	      old_list.pop_front();
	    }
	    old_list.pop_front();
	    //      std::cout << "f2 " << old_list.size() << std::endl;
	    for (auto it = new_list.rbegin(); it!=new_list.rend(); it++){
	      old_list.push_front(*it);
	    }
	    //      std::cout << "f3 " << old_list.size() << std::endl;
	  }else{
	    //      std::cout << "b1 " << old_list.size() << std::endl;
	    while(old_list.back().index != min_wcp.index && old_list.size()>0 ){
	      old_list.pop_back();
	    }
	    old_list.pop_back();
	    
	    //      std::cout << "b2 " << old_list.size() << std::endl;
	    for (auto it = new_list.rbegin(); it!=new_list.rend(); it++){
	      old_list.push_back(*it);
	    }
	    //      std::cout << "b3 " << old_list.size() << std::endl;
	  }
	  vec_wcps.clear();
	  vec_wcps.reserve(old_list.size());
	  std::copy(std::begin(old_list), std::end(old_list), std::back_inserter(vec_wcps));

	  add_proto_connection(v2, sg1, temp_cluster);	  
	}
	
	del_proto_vertex(v1);
	del_proto_segment(sg);
	
	flag_continue = true;
	temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);

      }else if (map_vertex_segments[v2].size()>=2 && examine_vertices_4(v2, v1)){
	// merge v2's segments to v1 ...

	for (auto it1 = map_vertex_segments[v2].begin(); it1 != map_vertex_segments[v2].end(); it1++){
	  WCPPID::ProtoSegment *sg1 = *it1;
	  if (sg1 == sg) continue;
	
	  WCP::WCPointCloud<double>::WCPoint& vtx_new_wcp = v1->get_wcpt();
	  
	  std::vector<WCP::WCPointCloud<double>::WCPoint >& vec_wcps = sg1->get_wcpt_vec();
	  bool flag_front = false;
	  if (vec_wcps.front().index == v2->get_wcpt().index) flag_front = true;
	  
	  WCP::WCPointCloud<double>::WCPoint min_wcp;
	  double min_dis = 1e9;
	  double max_dis = std::max(sqrt(pow(vec_wcps.front().x -v2->get_fit_pt().x,2)+pow(vec_wcps.front().y -v2->get_fit_pt().y,2)+pow(vec_wcps.front().z -v2->get_fit_pt().z,2)), sqrt(pow(vec_wcps.back().x -v2->get_fit_pt().x,2)+pow(vec_wcps.back().y -v2->get_fit_pt().y,2)+pow(vec_wcps.back().z -v2->get_fit_pt().z,2)));
	  double dis_cut = 0;
	  double default_dis_cut = 2.5*units::cm; 
	  if (max_dis > 2 * default_dis_cut) dis_cut = default_dis_cut;
	  for (size_t j=0;j!=vec_wcps.size();j++){
	    double dis1 = sqrt(pow(vec_wcps.at(j).x -v2->get_fit_pt().x,2)+pow(vec_wcps.at(j).y -v2->get_fit_pt().y,2)+pow(vec_wcps.at(j).z -v2->get_fit_pt().z,2));
	    double dis = fabs(dis1 - 3.0*units::cm);
	    if (dis < min_dis && dis1 > dis_cut){
	      min_wcp = vec_wcps.at(j);
	      min_dis = dis;
	    }
	  }
	  // establish the shortest path ...
	  std::list<WCP::WCPointCloud<double>::WCPoint> new_list;
	  
	  new_list.push_back(vtx_new_wcp);
	  {
	    double dis_step = 2.0*units::cm;
	    int ncount = std::round(sqrt(pow(vtx_new_wcp.x - min_wcp.x,2) + pow(vtx_new_wcp.y - min_wcp.y,2) + pow(vtx_new_wcp.z - min_wcp.z,2))/dis_step);
	    if (ncount <2) ncount = 2;
	    
	    for (int qx=1;qx<ncount;qx++){
	      Point tmp_p(vtx_new_wcp.x + (min_wcp.x-vtx_new_wcp.x)/ncount*qx, vtx_new_wcp.y + (min_wcp.y-vtx_new_wcp.y)/ncount*qx, vtx_new_wcp.z + (min_wcp.z-vtx_new_wcp.z)/ncount*qx);
	      WCP::WCPointCloud<double>::WCPoint& tmp_wcp = pcloud->get_closest_wcpoint(tmp_p);
	      //std::cout << qx << " " << sqrt(pow(tmp_wcp.x - tmp_p.x,2) + pow(tmp_wcp.y - tmp_p.y,2) + pow(tmp_wcp.z - tmp_p.z,2))/units::cm << std::endl;
	      if (sqrt(pow(tmp_wcp.x - tmp_p.x,2) + pow(tmp_wcp.y - tmp_p.y,2) + pow(tmp_wcp.z - tmp_p.z,2)) > 0.3*units::cm) continue; // too large ...
	      if (tmp_wcp.index != new_list.back().index && tmp_wcp.index != min_wcp.index)
		new_list.push_back(tmp_wcp);
	    }
	  }
	  new_list.push_back(min_wcp);
	  
	  std::list<WCP::WCPointCloud<double>::WCPoint> old_list;
	  std::copy(vec_wcps.begin(), vec_wcps.end(), std::back_inserter(old_list));
	  // append the list back ...
	  //std::cout << min_wcp.index << " " << temp_cluster->get_path_wcps().size() << std::endl;
	  if (flag_front){
	    //std::cout << "f1 " << old_list.size() << std::endl;
	    
	    while(old_list.front().index != min_wcp.index && old_list.size()>0 ){
	      old_list.pop_front();
	    }
	    old_list.pop_front();
	    //      std::cout << "f2 " << old_list.size() << std::endl;
	    for (auto it = new_list.rbegin(); it!=new_list.rend(); it++){
	      old_list.push_front(*it);
	    }
	    //      std::cout << "f3 " << old_list.size() << std::endl;
	  }else{
	    //      std::cout << "b1 " << old_list.size() << std::endl;
	    while(old_list.back().index != min_wcp.index && old_list.size()>0 ){
	      old_list.pop_back();
	    }
	    old_list.pop_back();
	    
	    //      std::cout << "b2 " << old_list.size() << std::endl;
	    for (auto it = new_list.rbegin(); it!=new_list.rend(); it++){
	      old_list.push_back(*it);
	    }
	    //      std::cout << "b3 " << old_list.size() << std::endl;
	  }
	  vec_wcps.clear();
	  vec_wcps.reserve(old_list.size());
	  std::copy(std::begin(old_list), std::end(old_list), std::back_inserter(vec_wcps));

	  add_proto_connection(v1, sg1, temp_cluster);	  
	}
	
	del_proto_vertex(v2);
	del_proto_segment(sg);
	
	flag_continue = true;
	temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
	
      }
      if (flag_continue){
	std::cout << "Merge Vertices Type III: segment id: " << sg->get_id() << std::endl;
	break;
      }
    }
  }
  
  return flag_continue;
}

bool WCPPID::NeutrinoID::examine_vertices_4(WCPPID::ProtoVertex *v1, WCPPID::ProtoVertex *v2){

  WCPPID::ProtoSegment *sg1 = find_segment(v1,v2);
  
  bool flag = true;
  // check segments of v1 with respect v2 ...
  for (auto it = map_vertex_segments[v1].begin(); it != map_vertex_segments[v1].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    if (sg == sg1) continue;
    
    PointVector& pts = sg->get_point_vec();
    Point min_point = pts.front();
    double min_dis = 1e9;
    for (size_t i=0;i!=pts.size();i++){
      double dis = fabs(sqrt(pow(pts.at(i).x - v1->get_fit_pt().x,2) + pow(pts.at(i).y - v1->get_fit_pt().y,2) + pow(pts.at(i).z - v1->get_fit_pt().z,2))-3*units::cm);
      if (dis < min_dis){
	min_dis = dis;
	min_point = pts.at(i);
      }
    }

    //
    double step_size = 0.6*units::cm;
    Point start_p = min_point;
    Point end_p = v2->get_fit_pt();
    int ncount = std::round(sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/step_size);
    int n_bad = 0;
    for (int i=1;i<ncount;i++){
      Point test_p(start_p.x + (end_p.x - start_p.x)/ncount*i,
		   start_p.y + (end_p.y - start_p.y)/ncount*i,
		   start_p.z + (end_p.z - start_p.z)/ncount*i);
      if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
    }
    //    std::cout << n_bad << " " << ncount << " " << sg->get_id() << " " << sg1->get_id() << " " << sg1->get_length()/units::cm << " " << sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/units::cm << std::endl;
    if (n_bad !=0) {
      flag = false;
      break;
    }
  }
  //    std::cout << std::endl;
  
  return flag;
}



bool WCPPID::NeutrinoID::examine_vertices_1(WCPPID::PR3DCluster* temp_cluster){
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_offset = mp.get_time_offset();
  int nrebin = mp.get_nrebin();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w + 0.5;
  double offset_u = -first_u_dis/pitch_u + 0.5;
  double offset_v = -first_v_dis/pitch_v + 0.5;
  
  double first_t_dis = temp_cluster->get_point_cloud()->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - temp_cluster->get_point_cloud()->get_cloud().pts[0].x;
  double slope_xt = 1./time_slice_width;
  double offset_t =  first_t_dis/time_slice_width + 0.5;

  /* for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){ */
  /*   WCPPID::ProtoSegment *sg = it->first; */
  /*   if (sg->get_cluster_id()!=temp_cluster->get_cluster_id()) continue; */
  /*   auto it1 = find_vertices(sg); */
  /*   std::cout << it1.first->get_wcpt().index << " " << it1.second->get_wcpt().index << std::endl; */
  /* } */


  
  bool flag_continue = false;


  //std::map<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> map_vertex_replace;
  WCPPID::ProtoVertex *v1 = 0;
  WCPPID::ProtoVertex *v2 = 0;
  WCPPID::ProtoVertex *v3 = 0;
  
  for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    if (it->second.size()==2){ // potential check
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::ProtoSegment* sg = (*it1);
	
	if (sg->get_length() > 4*units::cm) continue;
	for (auto it2 = map_segment_vertices[sg].begin(); it2 != map_segment_vertices[sg].end(); it2++){
	  WCPPID::ProtoVertex *vtx1 = (*it2);
	  if (vtx1 == vtx) continue;
	  if (map_vertex_segments[vtx1].size() >= 2){ // the other track needs to be larger than 2
	    //	    std::cout << sg->get_id() << " " << vtx->get_id() << " " << vtx1->get_id() << " " <<  sg->get_length()/units::cm <<std::endl;
	    
	    if (examine_vertices_1(vtx, vtx1, offset_t, slope_xt, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw)){
	      v1 = vtx;
	      v2 = vtx1;
	      
	      for (auto it3 = it->second.begin(); it3 != it->second.end(); it3++){
		if (*it3 == sg) continue;
		for (auto it4 = map_segment_vertices[*it3].begin(); it4!=map_segment_vertices[*it3].end();it4++){
		  if (*it4 == v1 ) continue;
		  v3 = *it4;
		}
	      }
	      
	      flag_continue = true;
	      break;
	    }
	  }
	} // loop over vertices ..
	if (flag_continue) break;
      } // search for segments
      if (flag_continue) break;
    } // require two segment
  } // loop over vertices ...
  
  //    std::cout << map_vertex_replace.size() << std::endl;
  // replace ...
  if (v1!=0 && v2!=0){
    
    //      std::cout << v1 << " " << v2 << " " << v3 << std::endl;
    ProtoSegment *sg = find_segment(v1,v2);
    ProtoSegment *sg1 = find_segment(v1,v3);
    
    //      std::cout << v2->get_cluster_id() << " " << v3->get_cluster_id() << std::endl;
    
    //  std::cout << sg<< " " << sg1 << std::endl;
    temp_cluster->dijkstra_shortest_paths(v3->get_wcpt(),2); 
    temp_cluster->cal_shortest_path(v2->get_wcpt(),2);
    ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
    
    std::cout << "Cluster: " << temp_cluster->get_cluster_id() << " Merge Vertices Type I " << sg->get_id() << " + " << sg1->get_id() << " -> " << sg2->get_id() << std::endl;
    add_proto_connection(v2, sg2, temp_cluster);
    add_proto_connection(v3, sg2, temp_cluster);
    
    del_proto_vertex(v1);
    del_proto_segment(sg);
    del_proto_segment(sg1);
    
    temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
    
  }
  
  // hack for now
  // flag_continue = false;

  return flag_continue;
}

bool WCPPID::NeutrinoID::examine_vertices_1(WCPPID::ProtoVertex* v1, WCPPID::ProtoVertex *v2, double offset_t, double slope_xt, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw){
  // U, V, W ...
  Point& v1_p = v1->get_fit_pt();
  double v1_t = offset_t + slope_xt * v1_p.x;
  double v1_u = offset_u + slope_yu * v1_p.y + slope_zu * v1_p.z;
  double v1_v = offset_v + slope_yv * v1_p.y + slope_zv * v1_p.z;
  double v1_w = offset_w + slope_yw * v1_p.y + slope_zw * v1_p.z;
  
  Point& v2_p = v2->get_fit_pt();
  double v2_t = offset_t + slope_xt * v2_p.x;
  double v2_u = offset_u + slope_yu * v2_p.y + slope_zu * v2_p.z;
  double v2_v = offset_v + slope_yv * v2_p.y + slope_zv * v2_p.z;
  double v2_w = offset_w + slope_yw * v2_p.y + slope_zw * v2_p.z;

  // find the track segment in between
  WCPPID::ProtoSegment *sg = find_segment(v1,v2);
  if (sg == 0 || map_vertex_segments[v1].size()!=2) return false;
  
  
  int ncount_close = 0;
  int ncount_dead = 0;
  int ncount_line = 0;

  // one view must be close
  if (sqrt(pow(v1_u - v2_u, 2) + pow(v1_t - v2_t,2)) < 2.5){
    ncount_close ++;
  }else{
    // one view must be dead? or shorter distance
    bool flag_dead = true;
    for (size_t i=0;i!=sg->get_point_vec().size();i++){
      if (!ct_point_cloud->get_closest_dead_chs(sg->get_point_vec().at(i),0))
	flag_dead = false;
    }
    if (flag_dead) ncount_dead++;
    else{
      // third view must be like a line
      PointVector& pts_1 = sg->get_point_vec();
      WCPPID::ProtoSegment *sg1 = 0;
      for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	// sg is one track
	if (*it1 == sg) continue;
	sg1 = *it1;
      }

      PointVector& pts_2 = sg1->get_point_vec();
      TVector3 v1(v2_u-v1_u, v2_t-v1_t, 0);
      TVector3 v2(0,0,0);
      double min_dis = 1e9;
      Point start_p = v2_p;
      Point end_p ;
      
      for (size_t i=0;i!=pts_2.size();i++){
	double p_t = offset_t + slope_xt * pts_2.at(i).x;
	double p_u = offset_u + slope_yu * pts_2.at(i).y + slope_zu * pts_2.at(i).z;
	TVector3 v3(p_u-v1_u, p_t-v1_t,0);
	double dis = fabs(v3.Mag()-9);
	if (dis < min_dis){
	  min_dis = dis;
	  v2 = v3;
	  end_p = pts_2.at(i);
	}
      }

      if (180-v1.Angle(v2)/3.1415926*180. < 30 || v1.Mag()< 8 && 180-v1.Angle(v2)/3.1415926*180. < 35) ncount_line++;
      else{
	double step_size = 0.6*units::cm;
	int ncount = std::round(sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/step_size);
	int n_bad = 0;
	for (int i=1;i!=ncount;i++){
	  Point test_p(start_p.x + (end_p.x - start_p.x)/ncount*i,
		       start_p.y + (end_p.y - start_p.y)/ncount*i,
		       start_p.z + (end_p.z - start_p.z)/ncount*i);
	  if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
	}
	if (n_bad<=1) ncount_line ++;
      }
      //      std::cout << v1.Mag() << " " << v2.Mag() << " " << v1.Angle(v2)/3.1415926*180. << std::endl;
    }
  }

  // std::cout << ncount_close << " " << ncount_dead << " " << ncount_line << std::endl;
  
   // one view must be close
  if (sqrt(pow(v1_v - v2_v, 2) + pow(v1_t - v2_t,2)) < 2.5){
    ncount_close ++;
  }else{
    // one view must be dead? or shorter distance
    bool flag_dead = true;
    for (size_t i=0;i!=sg->get_point_vec().size();i++){
      if (!ct_point_cloud->get_closest_dead_chs(sg->get_point_vec().at(i),1))
	flag_dead = false;
    }
    if (flag_dead) ncount_dead++;
    else{
      // third view must be like a line
      PointVector& pts_1 = sg->get_point_vec();
      WCPPID::ProtoSegment *sg1 = 0;
      for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	// sg is one track
	if (*it1 == sg) continue;
	sg1 = *it1;
      }

      PointVector& pts_2 = sg1->get_point_vec();
      TVector3 v1(v2_v-v1_v, v2_t-v1_t, 0);
      TVector3 v2(0,0,0);
      double min_dis = 1e9;

      Point start_p = v2_p;
      Point end_p ;
      
      
      for (size_t i=0;i!=pts_2.size();i++){
	double p_t = offset_t + slope_xt * pts_2.at(i).x;
	double p_v = offset_v + slope_yv * pts_2.at(i).y + slope_zv * pts_2.at(i).z;
	TVector3 v3(p_v-v1_v, p_t-v1_t,0);
	double dis = fabs(v3.Mag()-9);
	if (dis < min_dis){
	  min_dis = dis;
	  v2 = v3;
	  end_p = pts_2.at(i);
	}
      }

      if (180-v1.Angle(v2)/3.1415926*180. < 30 || v1.Mag()< 8 && 180-v1.Angle(v2)/3.1415926*180. < 35) ncount_line++;
      else{
	double step_size = 0.6*units::cm;
	int ncount = std::round(sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/step_size);
	int n_bad = 0;
	for (int i=1;i!=ncount;i++){
	  Point test_p(start_p.x + (end_p.x - start_p.x)/ncount*i,
		       start_p.y + (end_p.y - start_p.y)/ncount*i,
		       start_p.z + (end_p.z - start_p.z)/ncount*i);
	  if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
	}
	if (n_bad<=1) ncount_line ++;
      }
      //      std::cout << v1.Mag() << " " << v2.Mag() << " " << v1.Angle(v2)/3.1415926*180. << std::endl;
    }
  }

  // std::cout << ncount_close << " " << ncount_dead << " " << ncount_line << std::endl;

   // one view must be close
  if (sqrt(pow(v1_w - v2_w, 2) + pow(v1_t - v2_t,2)) < 2.5){
    //    std::cout << "haha " << sqrt(pow(v1_w - v2_w, 2) + pow(v1_t - v2_t,2)) << std::endl;
    ncount_close ++;
  }else{
    // one view must be dead? or shorter distance
    bool flag_dead = true;
    for (size_t i=0;i!=sg->get_point_vec().size();i++){
      if (!ct_point_cloud->get_closest_dead_chs(sg->get_point_vec().at(i),2))
	flag_dead = false;
    }

    
    
    if (flag_dead) ncount_dead++;
    else{
      // third view must be like a line
      PointVector& pts_1 = sg->get_point_vec();
      WCPPID::ProtoSegment *sg1 = 0;
      for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	// sg is one track
	if (*it1 == sg) continue;
	sg1 = *it1;
      }

      PointVector& pts_2 = sg1->get_point_vec();
      TVector3 v1(v2_w-v1_w, v2_t-v1_t, 0);
      TVector3 v2(0,0,0);
      double min_dis = 1e9;
      Point start_p = v2_p;
      Point end_p ;
      for (size_t i=0;i!=pts_2.size();i++){
	double p_t = offset_t + slope_xt * pts_2.at(i).x;
	double p_w = offset_w + slope_yw * pts_2.at(i).y + slope_zw * pts_2.at(i).z;
	TVector3 v3(p_w-v1_w, p_t-v1_t,0);
	double dis = fabs(v3.Mag()-9);
	if (dis < min_dis){
	  min_dis = dis;
	  v2 = v3;
	  end_p = pts_2.at(i);
	}
      }

      if (180-v1.Angle(v2)/3.1415926*180. < 30 || v1.Mag()< 8 && 180-v1.Angle(v2)/3.1415926*180. < 35) ncount_line++;
      else{
	double step_size = 0.6*units::cm;
	int ncount = std::round(sqrt(pow(start_p.x - end_p.x,2) + pow(start_p.y - end_p.y,2) + pow(start_p.z - end_p.z,2))/step_size);
	int n_bad = 0;
	for (int i=1;i!=ncount;i++){
	  Point test_p(start_p.x + (end_p.x - start_p.x)/ncount*i,
		       start_p.y + (end_p.y - start_p.y)/ncount*i,
		       start_p.z + (end_p.z - start_p.z)/ncount*i);
	  if (!ct_point_cloud->is_good_point(test_p, 0.2*units::cm, 0, 0)) n_bad ++;
	}
	if (n_bad<=1) ncount_line ++;
      }
      //      std::cout << v1.Mag() << " " << v2.Mag() << " " << v1.Angle(v2)/3.1415926*180. << std::endl;
    }
  }


  //  std::cout << ncount_close << " " << ncount_dead << " " << ncount_line << std::endl;
  
  if (ncount_close >=2 ||
      ncount_close ==1 && ncount_dead ==1 & ncount_line>=1 ||
      ncount_close ==1 && ncount_dead == 2 ||
      ncount_close ==1 && ncount_line>=2 || 
      ncount_line >= 3)
    return true;

  
  //  std::cout <<  << " " << sqrt(pow(v1_u - v2_u, 2) + pow(v1_t - v2_t,2)) << " " << sqrt(pow(v1_v - v2_v, 2) + pow(v1_t - v2_t,2)) << std::endl;

  
  return false;
}

std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> WCPPID::NeutrinoID::find_vertices(WCPPID::ProtoSegment* sg){
  if (map_segment_vertices.find(sg) == map_segment_vertices.end())
    return std::make_pair( (WCPPID::ProtoVertex*)0, (WCPPID::ProtoVertex*)0);

  WCPPID::ProtoVertex *v1 = 0;
  WCPPID::ProtoVertex *v2 = 0;
  
  for (auto it = map_segment_vertices[sg].begin(); it != map_segment_vertices[sg].end(); it++){
    if (v1==0) {
      v1 = *it;
    }else{
      v2 = *it;
    }
  }

  return std::make_pair(v1,v2);
}

WCPPID::ProtoVertex* WCPPID::NeutrinoID::find_other_vertex(WCPPID::ProtoSegment *sg, WCPPID::ProtoVertex* v1){
  auto results = find_vertices(sg);
  if (v1 == results.first){
    return results.second;
  }else if (v1 == results.second){
    return results.first;
  }else{
    return (WCPPID::ProtoVertex*)0;
  }
}

WCPPID::ProtoSegment* WCPPID::NeutrinoID::find_segment(WCPPID::ProtoVertex *v1, WCPPID::ProtoVertex *v2){

  WCPPID::ProtoSegment* sg = 0;
  for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
    sg = *it1;
    if (map_segment_vertices[sg].find(v2) == map_segment_vertices[sg].end()){
      sg = 0;
    }else{
      break;
    }
  }
  return sg;
}


WCPPID::ProtoVertexSelection WCPPID::NeutrinoID::find_vertices(WCPPID::PR3DCluster* temp_cluster){
  WCPPID::ProtoVertexSelection tmp_vertices;
  for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
    if (temp_cluster->get_cluster_id() == (it->first)->get_cluster_id()) tmp_vertices.push_back(it->first);
  }
  return tmp_vertices;
}

WCPPID::ProtoSegmentSelection WCPPID::NeutrinoID::find_segments(WCPPID::PR3DCluster* temp_cluster){
  WCPPID::ProtoSegmentSelection tmp_segments;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    if (temp_cluster->get_cluster_id() == (it->first)->get_cluster_id()) tmp_segments.push_back(it->first);
  }
  return tmp_segments;
}
