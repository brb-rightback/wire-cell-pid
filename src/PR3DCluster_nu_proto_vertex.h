bool WCPPID::PR3DCluster::proto_break_tracks(WCP::WCPointCloud<double>::WCPoint& first_wcp, WCP::WCPointCloud<double>::WCPoint& curr_wcp, WCP::WCPointCloud<double>::WCPoint& last_wcp, std::list<WCP::WCPointCloud<double>::WCPoint>& wcps_list1, std::list<WCP::WCPointCloud<double>::WCPoint>& wcps_list2, bool flag_pass_check){

  double dis1 = sqrt(pow(curr_wcp.x-first_wcp.x,2) + pow(curr_wcp.y-first_wcp.y,2) + pow(curr_wcp.z-first_wcp.z,2));
  double dis2 = sqrt(pow(curr_wcp.x-last_wcp.x,2) + pow(curr_wcp.y-last_wcp.y,2) + pow(curr_wcp.z-last_wcp.z,2));

  //  std::cout << first_wcp.index << " " << curr_wcp.index << " " << last_wcp.index << " " << dis1/units::cm << " " << dis2/units::cm << std::endl;
  
  if (dis1 > 1*units::cm && dis2 > 1*units::cm || flag_pass_check){
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
    curr_wcp = wcps_list1.back();
    
    return true;
  }else {
    return false;
  }
}



WCP::WCPointCloud<double>::WCPoint WCPPID::PR3DCluster::proto_extend_point(WCP::Point& p, TVector3& dir, TVector3& dir_other, bool flag_continue){

  float step_dis = 1*units::cm;
  
  WCP::WCPointCloud<double>::WCPoint curr_wcp = point_cloud_steiner->get_closest_wcpoint(p);
  WCP::WCPointCloud<double>::WCPoint next_wcp = curr_wcp;

  WCP::WCPointCloud<double>::WCPoint saved_start_wcp = curr_wcp; // saved start position
  TVector3 saved_dir = dir; // saved start direction
  
  WCP::Point test_p;
  while(flag_continue){
    flag_continue = false;
    
    for (int i=0; i!=3; i++){
      test_p.x = curr_wcp.x + dir.X() * step_dis * (i+1);
      test_p.y = curr_wcp.y + dir.Y() * step_dis * (i+1);
      test_p.z = curr_wcp.z + dir.Z() * step_dis * (i+1);
      
      next_wcp = point_cloud_steiner->get_closest_wcpoint(test_p);
      TVector3 dir2(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z);
      
      if (dir2.Mag()!=0 && (dir2.Angle(dir)/3.1415926*180. < 25)){
	//|| (fabs(drift_dir.Angle(dir)/3.1415926*180.-90.)<10. && fabs(drift_dir.Angle(dir2)/3.1415926*180.-90.)<10.) && dir2.Angle(dir)/3.1415926*180. < 60)){
	flag_continue = true;
	curr_wcp = next_wcp;
	dir = dir2 + dir * 5 * units::cm; // momentum trick ...
	dir = dir.Unit();
	break;
      }
            
      next_wcp = point_cloud->get_closest_wcpoint(test_p);
      TVector3 dir1(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z);
      
      /* std::cout << dir.X() << " " << dir.Y() << " " << dir.Z() << " " << dir1.X() << " " << dir1.Y() << " " << dir1.Z() << " " << drift_dir.Angle(dir)/3.1415926*180. << " " << drift_dir.Angle(dir1)/3.1415926*180. << std::endl; */
      // std::cout <<curr_wcp.x << " " << curr_wcp.y << " " << curr_wcp.z << " " << i << " " << dir1.Mag() << " " << dir1.Angle(dir)/3.14151926*180. << " " << test_p << " " << next_wcp.x << " " << next_wcp.y << " " << next_wcp.z << std::endl; 
      
      if (dir1.Mag()!=0 && (dir1.Angle(dir)/3.1415926*180. < 17.5)){
      	//|| (fabs(drift_dir.Angle(dir)/3.1415926*180.-90.)<10. && fabs(drift_dir.Angle(dir1)/3.1415926*180.-90.)<10.) && dir1.Angle(dir)/3.1415926*180. < 60)){
      	flag_continue = true;
      	curr_wcp = next_wcp;
      	dir = dir1 + dir * 5 * units::cm; // momentum trick ...
      	dir = dir.Unit();
      	break;
      }
    }
  }

  /* if (curr_wcp.index != saved_start_wcp.index || sqrt(pow(curr_wcp.x - saved_start_wcp.x,2) + pow(curr_wcp.y - saved_start_wcp.y,2) + pow(curr_wcp.z - saved_start_wcp.z,2)) > 0.01*units::cm){ */
    // forward point ...
    //    std::cout << "Forward search" << std::endl;
  test_p.x = curr_wcp.x;
  test_p.y = curr_wcp.y;
  test_p.z = curr_wcp.z;
  curr_wcp = point_cloud_steiner->get_closest_wcpoint(test_p);
  /* }else{ */
  /*   // backward search ... */
  /*   //std::cout << "Backward search" << std::endl; */
    
  /*   while(flag_continue){ */
  /*     flag_continue = false; */
      
  /*     for (int i=0; i!=3; i++){ */
  /* 	test_p.x = curr_wcp.x + dir_other.X() * step_dis * (i+1); */
  /* 	test_p.y = curr_wcp.y + dir_other.Y() * step_dis * (i+1); */
  /* 	test_p.z = curr_wcp.z + dir_other.Z() * step_dis * (i+1); */
      
  /* 	next_wcp = point_cloud_steiner->get_closest_wcpoint(test_p); */
  /* 	TVector3 dir2(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z); */
  /* 	TVector3 dir3(saved_start_wcp.x - next_wcp.x, saved_start_wcp.y - next_wcp.y, saved_start_wcp.z - next_wcp.z);   */
  /* 	if (dir2.Mag()!=0 && (dir2.Angle(dir_other)/3.1415926*180. < 25) && (dir3.Angle(saved_dir)/3.1415926*180.) > 10.) { */
  /* 	  flag_continue = true; */
  /* 	  curr_wcp = next_wcp; */
  /* 	  dir_other = dir2 + dir_other * 5 * units::cm; // momentum trick ... */
  /* 	  dir_other = dir_other.Unit(); */
  /* 	  break; */
  /* 	} */
	
  /* 	next_wcp = point_cloud->get_closest_wcpoint(test_p); */
  /* 	TVector3 dir1(next_wcp.x - curr_wcp.x, next_wcp.y - curr_wcp.y, next_wcp.z - curr_wcp.z); */
            
  /* 	if (dir1.Mag()!=0 && (dir1.Angle(dir_other)/3.1415926*180. < 17.5) && (dir3.Angle(saved_dir)/3.1415926*180.) > 10.){ */
  /* 	  flag_continue = true; */
  /* 	  curr_wcp = next_wcp; */
  /* 	  dir_other = dir1 + dir_other * 5 * units::cm; // momentum trick ... */
  /* 	  dir_other = dir_other.Unit(); */
  /* 	  break; */
  /* 	} */
  /*     } */
  /*   } */
  /*   test_p.x = curr_wcp.x; */
  /*   test_p.y = curr_wcp.y; */
  /*   test_p.z = curr_wcp.z; */
  /*   curr_wcp = point_cloud_steiner->get_closest_wcpoint(test_p); */
  /* } */
  
  
  
  return curr_wcp;
}

void WCPPID::PR3DCluster::set_fit_parameters(WCPPID::ProtoVertex* vtx){
  if (vtx->get_cluster_id()==cluster_id) {
    fine_tracking_path.push_back((vtx)->get_fit_pt());
    dQ.push_back((vtx)->get_dQ());
    dx.push_back((vtx)->get_dx());
    pu.push_back((vtx)->get_pu());
    pv.push_back((vtx)->get_pv());
    pw.push_back((vtx)->get_pw());
    pt.push_back((vtx)->get_pt());
    reduced_chi2.push_back((vtx)->get_reduced_chi2());
    
    sub_cluster_rr.push_back(-1);
    
    flag_vertex.push_back(true);
    sub_cluster_id.push_back(-1);
    flag_shower.push_back(false);
  }
}

void WCPPID::PR3DCluster::set_fit_parameters(WCPPID::Map_Proto_Vertex_Segments& map_vertex_segments, WCPPID::Map_Proto_Segment_Vertices& map_segment_vertices){
  fine_tracking_path.clear();
  dQ.clear();
  dx.clear();
  pu.clear();
  pv.clear();
  pw.clear();
  pt.clear();
  reduced_chi2.clear();
  flag_vertex.clear();
  sub_cluster_id.clear();
  flag_shower.clear();
  sub_cluster_rr.clear();

  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
    set_fit_parameters(it->first);
  }

  //  int tmp_id = cluster_id*1000 + 1; // hack ...
  for (auto it=map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoVertex *start_v=0, *end_v=0;
    for (auto it1 = map_segment_vertices[it->first].begin(); it1!=map_segment_vertices[it->first].end(); it1++){
      if ((*it1)->get_wcpt().index == it->first->get_wcpt_vec().front().index) start_v = *it1;
      if ((*it1)->get_wcpt().index == it->first->get_wcpt_vec().back().index) end_v = *it1;
    }
    int start_n = map_vertex_segments[start_v].size();
    int end_n = map_vertex_segments[end_v].size();
    
    
    set_fit_parameters(it->first, start_n, end_n);
  }
  //
}

void WCPPID::PR3DCluster::set_fit_parameters(WCPPID::ProtoSegment* seg, int start_n, int end_n){
  if ((seg)->get_cluster_id()==cluster_id){
    fine_tracking_path.insert(fine_tracking_path.end(),(seg)->get_point_vec().begin(), (seg)->get_point_vec().end());
    dQ.insert(dQ.end(),(seg)->get_dQ_vec().begin(), (seg)->get_dQ_vec().end());
    dx.insert(dx.end(),(seg)->get_dx_vec().begin(), (seg)->get_dx_vec().end());
    pu.insert(pu.end(),(seg)->get_pu_vec().begin(), (seg)->get_pu_vec().end());
    pv.insert(pv.end(),(seg)->get_pv_vec().begin(), (seg)->get_pv_vec().end());
    pw.insert(pw.end(),(seg)->get_pw_vec().begin(), (seg)->get_pw_vec().end());
    pt.insert(pt.end(),(seg)->get_pt_vec().begin(), (seg)->get_pt_vec().end());
    reduced_chi2.insert(reduced_chi2.end(),(seg)->get_reduced_chi2_vec().begin(), (seg)->get_reduced_chi2_vec().end());

    // calculate rr vector ...
    WCP::PointVector& pts = (seg)->get_point_vec();
    std::vector<double> L(pts.size(),0);
    std::vector<double> rr(pts.size(),0);
    double acc_length = 0;
    for (size_t i=0;i+1<pts.size();i++){
      acc_length += sqrt(pow(pts.at(i+1).x-pts.at(i).x,2) + pow(pts.at(i+1).y - pts.at(i).y ,2) + pow(pts.at(i+1).z - pts.at(i).z,2)); 
      L.at(i+1) = acc_length;
    }
    
   
    if ((seg)->get_flag_dir()==1){// forward direction
      for (size_t i=0;i!=pts.size();i++){
      	rr.at(pts.size()-1-i) = L.back() - L.at(pts.size()-1-i);
      }
    }else if ((seg)->get_flag_dir()==-1){ // reverse direction
      rr = L;
    }else{
      rr = L; // quick
    }
    if (start_n>1) rr.front() = -1;
    if (end_n >1) rr.back() = -1;
    sub_cluster_rr.insert(sub_cluster_rr.end(), rr.begin(), rr.end());

    //    std::cout << seg->get_flag_dir() << " " << reduced_chi2.size() << " " << sub_cluster_rr.size() << std::endl;

    
    bool is_shower = (seg)->get_flag_shower();
    for (size_t i=0;i!=(seg)->get_point_vec().size();i++){
      flag_vertex.push_back(false);
      sub_cluster_id.push_back(cluster_id*1000 + (seg)->get_id());
      flag_shower.push_back(is_shower);
    }
  }
}


void WCPPID::PR3DCluster::set_fit_parameters(WCPPID::ProtoVertexSelection& temp_vertices, WCPPID::ProtoSegmentSelection& temp_segments, Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices){
  fine_tracking_path.clear();
  dQ.clear();
  dx.clear();
  pu.clear();
  pv.clear();
  pw.clear();
  pt.clear();
  reduced_chi2.clear();
  flag_vertex.clear();
  sub_cluster_id.clear();
  flag_shower.clear();

  for (auto it = temp_vertices.begin(); it!=temp_vertices.end(); it++){
    set_fit_parameters(*it);
  }
  
  for (auto it=temp_segments.begin(); it!=temp_segments.end(); it++){
     WCPPID::ProtoVertex *start_v=0, *end_v=0;
    for (auto it1 = map_segment_vertices[*it].begin(); it1!=map_segment_vertices[*it].end(); it1++){
      if ((*it1)->get_wcpt().index == (*it)->get_wcpt_vec().front().index) start_v = *it1;
      if ((*it1)->get_wcpt().index == (*it)->get_wcpt_vec().back().index) end_v = *it1;
    }
    int start_n = map_vertex_segments[start_v].size();
    int end_n = map_vertex_segments[end_v].size();
    set_fit_parameters(*it, start_n, end_n);
  }
  //
}
