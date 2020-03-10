
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
  //std::cout << path_wcps.size() << std::endl;

  
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
