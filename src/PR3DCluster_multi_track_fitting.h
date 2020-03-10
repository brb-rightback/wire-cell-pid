
void WCPPID::PR3DCluster::do_multi_tracking(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, double time, bool flag_dQ_dx_fit_reg){
  // collect charge
  collect_charge_multi_trajectory(map_segment_vertices, ct_point_cloud);

  // clear things
  fine_tracking_path.clear();
  dQ.clear();
  dx.clear();
  pu.clear();
  pv.clear();
  pw.clear();
  pt.clear();
  reduced_chi2.clear();

  bool flag_1st_tracking = true;
  bool flag_2nd_tracking = true;
  bool flag_dQ_dx = true;

  // prepare the data for the fit, do not contain everything ...
  // form from the cluster ...
  std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_ut_charge;
  std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_vt_charge;
  std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_wt_charge;
  prepare_data(ct_point_cloud, global_wc_map, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

  // first round of organizing the path from the path_wcps (shortest path)
  double low_dis_limit = 1.2*units::cm;
  double end_point_limit = 0.6*units::cm;
  organize_segments_path(map_vertex_segments, map_segment_vertices, low_dis_limit, end_point_limit);
  
  // form association ...
  // map 3D index to set of 2D points
  std::map<int,std::pair<std::set<std::pair<int,int>>, float> > map_3D_2DU_set;
  std::map<int,std::pair<std::set<std::pair<int,int>>, float> > map_3D_2DV_set;
  std::map<int,std::pair<std::set<std::pair<int,int>>, float> > map_3D_2DW_set;
  std::map<int,std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, int> > map_3D_tuple;
  // map 2D points to 3D indices
  std::map<std::pair<int,int>,std::set<int>> map_2DU_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DV_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DW_3D_set;
  
  if (flag_1st_tracking){
    form_map_multi_segments(map_vertex_segments, map_segment_vertices, ct_point_cloud,
			    map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,
			    map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set, map_3D_tuple,
			    map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set);
    
  }
  
}

void WCPPID::PR3DCluster::form_map_multi_segments(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud,
						  std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge,
						  std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<int, std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, int> >& map_3D_tuple,
						  std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,
						  double end_point_factor, double mid_point_factor, int nlevel, double time_cut, double charge_cut){
}

void WCPPID::PR3DCluster::organize_segments_path(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, double low_dis_limit, double end_point_limit){
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    WCPPID::ProtoVertex *start_v = 0, *end_v = 0;
    for (auto it1=it->second.begin(); it1!=it->second.end(); it1++){
      WCPPID::ProtoVertex *vt = *it1;
      if ( vt->get_wcpt().index == sg->get_wcpt_vec().front().index){
	start_v = vt;
      }else if ( vt->get_wcpt().index == sg->get_wcpt_vec().back().index){
	end_v = vt;
      }
    }
    bool flag_startv_end = true;
    bool flag_endv_end = true;
    if (map_vertex_segments[start_v].size()>1) flag_startv_end = false;
    if (map_vertex_segments[end_v].size()>1) flag_endv_end = false;

    PointVector pts;
    std::vector<WCPointCloud<double>::WCPoint> temp_wcps_vec = sg->get_wcpt_vec();
    Point start_p, end_p;
    
    // start_ point
    start_p.x = temp_wcps_vec.front().x;
    start_p.y = temp_wcps_vec.front().y;
    start_p.z = temp_wcps_vec.front().z;

    if (flag_startv_end){
      Point p2(temp_wcps_vec.front().x, temp_wcps_vec.front().y, temp_wcps_vec.front().z);
      double dis1 = 0;
      for (auto it = temp_wcps_vec.begin(); it!=temp_wcps_vec.end(); it++){
	p2.x = (*it).x;
	p2.y = (*it).y;
	p2.z = (*it).z;
	dis1 = sqrt(pow(start_p.x-p2.x,2)+pow(start_p.y-p2.y,2)+pow(start_p.z-p2.z,2));
	if (dis1 > low_dis_limit) break;
      }
      if (dis1!=0){
	start_p.x += (start_p.x-p2.x)/dis1*end_point_limit;
	start_p.y += (start_p.y-p2.y)/dis1*end_point_limit;
	start_p.z += (start_p.z-p2.z)/dis1*end_point_limit;
      }
    }
    
    // end_point
    end_p.x = temp_wcps_vec.back().x;
    end_p.y = temp_wcps_vec.back().y;
    end_p.z = temp_wcps_vec.back().z;
    if (flag_endv_end){
      Point p2(temp_wcps_vec.back().x, temp_wcps_vec.back().y, temp_wcps_vec.back().z);
      double dis1 = 0;
      for (auto it = temp_wcps_vec.rbegin(); it!=temp_wcps_vec.rend(); it++){
	p2.x = (*it).x;
	p2.y = (*it).y;
	p2.z = (*it).z;
	dis1 = sqrt(pow(end_p.x-p2.x,2)+pow(end_p.y-p2.y,2)+pow(end_p.z-p2.z,2));
	if (dis1 > low_dis_limit) break;
      }
      if (dis1!=0){
	end_p.x += (end_p.x - p2.x)/dis1 * end_point_limit;
	end_p.y += (end_p.y - p2.y)/dis1 * end_point_limit;
	end_p.z += (end_p.z - p2.z)/dis1 * end_point_limit;
      }
    }
    start_v->set_fit_pt(start_p);
    end_v->set_fit_pt(end_p);

    // middle points
    pts.push_back(start_p);
    for (size_t i=0;i!=temp_wcps_vec.size(); i++){
      Point p1(temp_wcps_vec.at(i).x, temp_wcps_vec.at(i).y, temp_wcps_vec.at(i).z);
      double dis = low_dis_limit;
      double dis1 = sqrt(pow(p1.x-end_p.x,2) + pow(p1.y-end_p.y,2) + pow(p1.z-end_p.z,2));
      if (pts.size()>0)
	dis = sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2));

      if (dis1 < low_dis_limit * 0.8){
	continue;
      }else if (dis < low_dis_limit * 0.8 ){
	continue;
      }else if (dis < low_dis_limit * 1.6){
	pts.push_back(p1);
	//std::cout << p1 << " " << sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2))/units::cm << std::endl;
      }else{
	int npoints = std::round(dis/low_dis_limit);
	
	Point p_save = pts.back();
	for (int j=0;j!=npoints;j++){
	  Point p(p_save.x + (p1.x-p_save.x) / npoints * (j+1),
		  p_save.y + (p1.y-p_save.y) / npoints * (j+1),
		  p_save.z + (p1.z-p_save.z) / npoints * (j+1));
	  pts.push_back(p);
	  //	std::cout << p << " " << sqrt(pow(p.x-pts.back().x,2)+pow(p.y-pts.back().y,2)+pow(p.z-pts.back().z,2))/units::cm << std::endl;
	}
	
      }
    }
    pts.push_back(end_p);

    sg->set_point_vec(pts);
    //std::cout << sg << " " << start_p << " " << end_p << " " << temp_wcps_vec.size() << " " << pts.size() << std::endl;
    //    std::cout << sg << flag_startv_end << " " << flag_endv_end << std::endl;
    //    std::cout << sg << " " << start_v << " " << end_v << std::endl;
  }
}
    


void WCPPID::PR3DCluster::collect_charge_multi_trajectory(std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet>& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud, double dis_cut, double range_cut){
  //clear up ...
  collected_charge_map.clear();
  
  std::set<std::pair<int,int>> existing_tcs;

  // form a set cotaining everything inside the cluster
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    GeomWireSelection& uwires = mcell->get_uwires();
    GeomWireSelection& vwires = mcell->get_vwires();
    GeomWireSelection& wwires = mcell->get_wwires();
    for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
      const GeomWire *wire = (*it1);
      existing_tcs.insert(std::make_pair(time_slice,wire->channel()));
    }
    for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
      const GeomWire *wire = (*it1);
      existing_tcs.insert(std::make_pair(time_slice,wire->channel()));
    }
    for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
      const GeomWire *wire = (*it1);
      existing_tcs.insert(std::make_pair(time_slice,wire->channel()));
    }
  }


  for (auto it2 = map_segment_vertices.begin(); it2 != map_segment_vertices.end(); it2++){
    WCPPID::ProtoSegment* sg1 = it2->first;
    
    // form a trajectory according to dis and fine tracking?
    PointVector traj_pts;
    //  PointVector& pts = get_fine_tracking_path();
    
    // not using fine tracking results ... 
    PointVector pts;
    std::vector<WCPointCloud<double>::WCPoint>& path_wcps = sg1->get_wcpt_vec();
    for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
      Point p((*it).x, (*it).y, (*it).z);
      pts.push_back(p);
    }

     //std::cout << "trajectory points " << pts.size() << std::endl;

    // if (cluster_id==13){
    //   for (int i=0; i!=pts.size(); i++){
    //     std::cout << i << " " << pts.at(i).x/units::cm << " " << pts.at(i).y/units::cm << " "<< pts.at(i).z/units::cm << std::endl;
    //   }
    // }

    for (int i=0; i!=pts.size(); i++){
      if (pts.at(i).y <-120*units::cm || pts.at(i).y > 120*units::cm ||
	  pts.at(i).z < -5*units::cm || pts.at(i).z > 1070*units::cm) continue;
      
      if (traj_pts.size()==0){
	traj_pts.push_back(pts.at(i));
      }else{
	double dis = sqrt(pow(pts.at(i).x-pts.at(i-1).x,2) +
			  pow(pts.at(i).y-pts.at(i-1).y,2) +
			  pow(pts.at(i).z-pts.at(i-1).z,2));
	if (dis <= dis_cut){
	  traj_pts.push_back(pts.at(i));
	}else{
	  int nseg = dis / dis_cut + 1;
	  //	std::cout << i << " " << pts.at(i).x/units::cm << " " << pts.at(i).y/units::cm << " "<< pts.at(i).z/units::cm << " " << dis << " " << dis_cut << " " << nseg << std::endl;
	  for (int j=0; j!=nseg;j++){
	    Point temp_pt;
	    temp_pt.x = pts.at(i-1).x + (pts.at(i).x-pts.at(i-1).x) *(j+1.)/nseg;
	    temp_pt.y = pts.at(i-1).y + (pts.at(i).y-pts.at(i-1).y) *(j+1.)/nseg;
	    temp_pt.z = pts.at(i-1).z + (pts.at(i).z-pts.at(i-1).z) *(j+1.)/nseg;
	    traj_pts.push_back(temp_pt);
	  }
	}
      }
    }

    // collect the nearby points, and compare with existing maps
    for (size_t i=0;i!=traj_pts.size();i++){
      //std::cout << i << " " << traj_pts.at(i).x/units::cm << " " << traj_pts.at(i).y/units::cm << " " << traj_pts.at(i).z/units::cm << " " << range_cut << std::endl;
      WCP::CTPointCloud<double> nearby_points = ct_point_cloud.get_closest_points(traj_pts.at(i),range_cut,0);
      
      //    std::cout << "0 " << nearby_points.pts.size() << std::endl;
      
      for (size_t j=0;j!=nearby_points.pts.size();j++){
	if (existing_tcs.find(std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)) == existing_tcs.end()){
	  collected_charge_map[std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)] = std::make_pair(nearby_points.pts.at(j).charge,nearby_points.pts.at(j).charge_err);
	}
      }
      nearby_points = ct_point_cloud.get_closest_points(traj_pts.at(i),range_cut,1);
      //std::cout << "1 " << nearby_points.pts.size() << std::endl;
      
      for (size_t j=0;j!=nearby_points.pts.size();j++){
	if (existing_tcs.find(std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)) == existing_tcs.end()){
	  collected_charge_map[std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)] = std::make_pair(nearby_points.pts.at(j).charge,nearby_points.pts.at(j).charge_err);
	}
      }
      nearby_points = ct_point_cloud.get_closest_points(traj_pts.at(i),range_cut,2);
      
      //    std::cout << "2 " << nearby_points.pts.size() << std::endl;
      
      for (size_t j=0;j!=nearby_points.pts.size();j++){
	if (existing_tcs.find(std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)) == existing_tcs.end()){
	  collected_charge_map[std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)] = std::make_pair(nearby_points.pts.at(j).charge,nearby_points.pts.at(j).charge_err);
	}
      }
    }
    
    //    std::cout << collected_charge_map.size() << std::endl;
    
  }

  
}
