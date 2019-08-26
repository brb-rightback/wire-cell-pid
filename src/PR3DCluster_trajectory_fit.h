std::vector<int> WireCellPID::PR3DCluster::examine_point_association(std::set<std::pair<int,int> >& temp_2dut, std::set<std::pair<int,int> >& temp_2dvt, std::set<std::pair<int,int> >& temp_2dwt,
								     std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge, double charge_cut){

  std::set<int> temp_types_u;
  std::set<int> temp_types_v;
  std::set<int> temp_types_w;
    
  std::set<std::pair<int,int> > saved_2dut;
  std::set<std::pair<int,int> > saved_2dvt;
  std::set<std::pair<int,int> > saved_2dwt;

  std::vector<int> results;
  results.resize(3,0);
  
  for (auto it = temp_2dut.begin(); it!=temp_2dut.end(); it++){
    auto it1 = map_2D_ut_charge.find(*it);
    if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > charge_cut ){
      temp_types_u.insert(std::get<2>(it1->second));
      if (std::get<2>(it1->second)==0) results.at(0)++;
      saved_2dut.insert(*it);
    }
  }

  for (auto it = temp_2dvt.begin(); it!=temp_2dvt.end(); it++){
    auto it1 = map_2D_vt_charge.find(*it);
    if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > charge_cut ){
      temp_types_v.insert(std::get<2>(it1->second));
      if (std::get<2>(it1->second)==0) results.at(1)++;
      saved_2dvt.insert(*it);
    }
  }

  for (auto it = temp_2dwt.begin(); it!=temp_2dwt.end(); it++){
    auto it1 = map_2D_wt_charge.find(*it);
    if (it1!=map_2D_wt_charge.end() && std::get<0>(it1->second) > charge_cut ){
      temp_types_w.insert(std::get<2>(it1->second));
      if (std::get<2>(it1->second)==0) results.at(2)++;
      saved_2dwt.insert(*it);
    }
  }

  if (temp_types_u.find(0)!=temp_types_u.end() && temp_types_u.size()==1){
    saved_2dut.clear();
    results.at(0) = 0;
  }
  if (temp_types_v.find(0)!=temp_types_v.end() && temp_types_v.size()==1){
    saved_2dvt.clear();
    results.at(1) = 0;
  }
  if (temp_types_w.find(0)!=temp_types_w.end() && temp_types_w.size()==1){
    saved_2dwt.clear();
    results.at(2) = 0;
  }
  temp_2dut = saved_2dut;
  temp_2dvt = saved_2dvt;
  temp_2dwt = saved_2dwt;

  return results;
}

void WireCellPID::PR3DCluster::form_point_association(WireCell::Point &p, std::set<std::pair<int,int> >& temp_2dut, std::set<std::pair<int,int> >& temp_2dvt, std::set<std::pair<int,int> >& temp_2dwt, WireCell::ToyCTPointCloud& ct_point_cloud, double dis_cut, int nlevel, double time_cut ){
  // global information
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  
  float coef1 = 2 * pow(sin(angle_u),2);
  float coef2 = 2 * (pow(sin(angle_u),2) - pow(cos(angle_u),2));

  typedef boost::property_map<MCUGraph, boost::vertex_index_t>::type IndexMap;
  typedef boost::graph_traits<MCUGraph>::adjacency_iterator adjacency_iterator;
  
  // original point cloud ...
  if (point_cloud!=0 && graph!=0){
    WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
    IndexMap index = get(boost::vertex_index,*graph);
    if (cloud.pts.size()>0){
      WireCell::WCPointCloud<double>::WCPoint wcp = point_cloud->get_closest_wcpoint(p);
      double temp_dis = sqrt(pow(wcp.x-p.x,2)+pow(wcp.y-p.y,2)+pow(wcp.z-p.z,2));
      
      //std::cout << temp_dis/units::cm << " " << dis_cut/units::cm << std::endl;
      if (temp_dis < dis_cut){
	std::set<int> total_vertices_found;
	std::set<int> vertices_to_be_examined;
	std::set<int> vertices_saved_for_next;
	total_vertices_found.insert(wcp.index);
	vertices_to_be_examined.insert(wcp.index);

	//	std::cout << cloud.pts.size() << std::endl;
	
	for (int j=0;j!=nlevel;j++){
	  for (auto it = vertices_to_be_examined.begin(); it!=vertices_to_be_examined.end(); it++){
	    int temp_current_index = (*it);

	    //	    std::cout << temp_current_index << " " << cloud.pts.size() << std::endl;

	    std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(temp_current_index,*graph),*graph);
	    for (; neighbors.first!=neighbors.second; ++neighbors.first){
	      //std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
	      if (total_vertices_found.find(index(*neighbors.first))==total_vertices_found.end()){
		total_vertices_found.insert(index(*neighbors.first));
		vertices_saved_for_next.insert(index(*neighbors.first));
	      }
	    }
	  }
	  vertices_to_be_examined = vertices_saved_for_next;
	}
	
	SMGCSet nearby_mcells_set;
	for (auto it = total_vertices_found.begin(); it!=total_vertices_found.end(); it++){
	  SlimMergeGeomCell *mcell = cloud.pts[*it].mcell;
	  nearby_mcells_set.insert(mcell);
	}
	
	int cur_time_slice = cloud.pts[wcp.index].mcell->GetTimeSlice();
	int cur_wire_u = cloud.pts[wcp.index].index_u;
	int cur_wire_v = cloud.pts[wcp.index].index_v;
	int cur_wire_w = cloud.pts[wcp.index].index_w;
	
	double dis_cut_u = dis_cut;
	double dis_cut_v = dis_cut;
	double dis_cut_w = dis_cut;
	
	
	double max_time_slice_u = 0;
	double max_time_slice_v = 0;
	double max_time_slice_w = 0;
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_u)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
		max_time_slice_u = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_v)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
		max_time_slice_v = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_w)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
		max_time_slice_w = fabs(this_time_slice-cur_time_slice);
	  }
	  
	}
	
	if (max_time_slice_u * time_slice_width*1.2 < dis_cut_u)
	  dis_cut_u = max_time_slice_u * time_slice_width*1.2;
	if (max_time_slice_v * time_slice_width*1.2 < dis_cut_v)
	  dis_cut_v = max_time_slice_v * time_slice_width*1.2;
	if (max_time_slice_w * time_slice_width*1.2 < dis_cut_w)
	  dis_cut_w = max_time_slice_w * time_slice_width*1.2;
	
	
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  
	  double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    float min_u_dis;
	    if (cur_wire_u < uwires.front()->index()){
	      min_u_dis = uwires.front()->index()-cur_wire_u;
	    }else if (cur_wire_u >= uwires.front()->index() &&
		      cur_wire_u <= uwires.back()->index()){
	      min_u_dis = 0;
	    }else{
	      min_u_dis = cur_wire_u-uwires.back()->index();
	    }
	    float min_v_dis;
	    if (cur_wire_v < vwires.front()->index()){
	      min_v_dis = vwires.front()->index()-cur_wire_v;
	    }else if (cur_wire_v >= vwires.front()->index() &&
		      cur_wire_v <= vwires.back()->index()){
	      min_v_dis = 0;
	    }else{
	      min_v_dis = cur_wire_v-vwires.back()->index();
	    }
	    float min_w_dis;
	    if (cur_wire_w < wwires.front()->index()){
	      min_w_dis = wwires.front()->index()-cur_wire_w;
	    }else if (cur_wire_w >= wwires.front()->index() &&
		      cur_wire_w <= wwires.back()->index()){
	      min_w_dis = 0;
	    }else{
	      min_w_dis = cur_wire_w-wwires.back()->index();
	    }
	    float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	    
	    if ( range_u > 0 && range_v >0 && range_w > 0){
	      float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
	      float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
	      float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
	      float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
	      float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
	      float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;
	      for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
		//auto it1 = map_2D_ut_charge.find(std::make_pair(j,this_time_slice));
		//if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > 0 )
		temp_2dut.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
		//auto it1 = map_2D_vt_charge.find(std::make_pair(j,this_time_slice));
		//if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > 0 )
		temp_2dvt.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
		//auto it1 = map_2D_wt_charge.find(std::make_pair(j,this_time_slice));
		//if (it1!=map_2D_wt_charge.end() && std::get<0>(it1->second) > 0 )
		temp_2dwt.insert(std::make_pair(j,this_time_slice));
	      }
	    } // cuts
	  } // cuts
	} // loop over all mcells
      } // find a good closest point?
    } // point cloud exist
  }
    
  // steiner tree point cloud ...
  if (point_cloud_steiner!=0 && graph_steiner!=0){
    WireCell::WCPointCloud<double>& cloud = point_cloud_steiner->get_cloud();
    IndexMap index = get(boost::vertex_index,*graph_steiner);

    if (cloud.pts.size()>0){
      WireCell::WCPointCloud<double>::WCPoint wcp = point_cloud_steiner->get_closest_wcpoint(p);
      double temp_dis = sqrt(pow(wcp.x-p.x,2)+pow(wcp.y-p.y,2)+pow(wcp.z-p.z,2));
      
      //std::cout << temp_dis/units::cm << " " << dis_cut/units::cm << std::endl;
      if (temp_dis < dis_cut){
	std::set<int> total_vertices_found;
	std::set<int> vertices_to_be_examined;
	std::set<int> vertices_saved_for_next;
	total_vertices_found.insert(wcp.index);
	vertices_to_be_examined.insert(wcp.index);
	
	for (int j=0;j!=nlevel;j++){
	  for (auto it = vertices_to_be_examined.begin(); it!=vertices_to_be_examined.end(); it++){
	    int temp_current_index = (*it);
	    std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(temp_current_index,*graph_steiner),*graph_steiner);
	    for (; neighbors.first!=neighbors.second; ++neighbors.first){
	      //std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
	      if (total_vertices_found.find(index(*neighbors.first))==total_vertices_found.end()){
		total_vertices_found.insert(index(*neighbors.first));
		vertices_saved_for_next.insert(index(*neighbors.first));
	      }
	    }
	  }
	  vertices_to_be_examined = vertices_saved_for_next;
	}
	
	WireCell::Point temp_p(wcp.x, wcp.y, wcp.z);
	std::vector<int> temp_results = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
	// std::cout << cloud.pts[wcp.index].index_u << " "  << temp_results.at(0) << " " << temp_results.at(1) << std::endl;
	int cur_time_slice = temp_results.at(0);
	int cur_wire_u = cloud.pts[wcp.index].index_u;
	int cur_wire_v = cloud.pts[wcp.index].index_v;
	int cur_wire_w = cloud.pts[wcp.index].index_w;
	
	SMGCSet nearby_mcells_set;
	std::vector<int> point_indices;
	std::vector<int> point_timeslices;
	
	for (auto it = total_vertices_found.begin(); it!=total_vertices_found.end(); it++){
	  SlimMergeGeomCell *mcell = cloud.pts[*it].mcell;
	  if (mcell!=0){
	    nearby_mcells_set.insert(mcell);
	  }else{
	    temp_p.x = cloud.pts[*it].x;
	    temp_p.y = cloud.pts[*it].y;
	    temp_p.z = cloud.pts[*it].z;
	    temp_results = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
	    
	    point_indices.push_back(*it);
	    point_timeslices.push_back(temp_results.at(0));
	  }
	}
	
	
	
	double dis_cut_u = dis_cut;
	double dis_cut_v = dis_cut;
	double dis_cut_w = dis_cut;
	
	double max_time_slice_u = 0;
	double max_time_slice_v = 0;
	double max_time_slice_w = 0;
	
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_u)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
		max_time_slice_u = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_v)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
		max_time_slice_v = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_w)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
		max_time_slice_w = fabs(this_time_slice-cur_time_slice);
	  }
	}
	for (size_t i=0; i!=point_indices.size(); i++){
	  int this_time_slice = point_timeslices.at(i);
	  if (fabs(cloud.pts[point_indices.at(i)].index_u - cur_wire_u)<=2){
	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
	      max_time_slice_u = fabs(this_time_slice-cur_time_slice);
	  }
	  if (fabs(cloud.pts[point_indices.at(i)].index_v - cur_wire_v)<=2){
	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
	      max_time_slice_v = fabs(this_time_slice-cur_time_slice);
	  }
	  if (fabs(cloud.pts[point_indices.at(i)].index_w - cur_wire_w)<=2){
	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
	      max_time_slice_w = fabs(this_time_slice-cur_time_slice);
	  }
	}
	
	if (max_time_slice_u * time_slice_width*1.2 < dis_cut_u)
	  dis_cut_u = max_time_slice_u * time_slice_width*1.2;
	if (max_time_slice_v * time_slice_width*1.2 < dis_cut_v)
	  dis_cut_v = max_time_slice_v * time_slice_width*1.2;
	if (max_time_slice_w * time_slice_width*1.2 < dis_cut_w)
	  dis_cut_w = max_time_slice_w * time_slice_width*1.2;
	
	// actual cut ...
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  
	  double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    float min_u_dis;
	    if (cur_wire_u < uwires.front()->index()){
	      min_u_dis = uwires.front()->index()-cur_wire_u;
	    }else if (cur_wire_u >= uwires.front()->index() &&
		      cur_wire_u <= uwires.back()->index()){
	      min_u_dis = 0;
	    }else{
	      min_u_dis = cur_wire_u-uwires.back()->index();
	    }
	    float min_v_dis;
	    if (cur_wire_v < vwires.front()->index()){
	      min_v_dis = vwires.front()->index()-cur_wire_v;
	    }else if (cur_wire_v >= vwires.front()->index() &&
		      cur_wire_v <= vwires.back()->index()){
	      min_v_dis = 0;
	    }else{
	      min_v_dis = cur_wire_v-vwires.back()->index();
	    }
	    float min_w_dis;
	    if (cur_wire_w < wwires.front()->index()){
	      min_w_dis = wwires.front()->index()-cur_wire_w;
	    }else if (cur_wire_w >= wwires.front()->index() &&
		      cur_wire_w <= wwires.back()->index()){
	      min_w_dis = 0;
	    }else{
	      min_w_dis = cur_wire_w-wwires.back()->index();
	    }
	    float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	    
	    if ( range_u > 0 && range_v >0 && range_w > 0){
	      float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
	      float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
	      float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
	      float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
	      float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
	      float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;
	      for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
		temp_2dut.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
		temp_2dvt.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
		temp_2dwt.insert(std::make_pair(j,this_time_slice));
	      }
	    } // cuts
	  } // cuts
	} // loop over all mcells
	
	for (size_t i=0; i!=point_indices.size(); i++){
	  int this_time_slice = point_timeslices.at(i);
	  int this_index_u = cloud.pts[point_indices.at(i)].index_u;
	  int this_index_v = cloud.pts[point_indices.at(i)].index_v;
	  int this_index_w = cloud.pts[point_indices.at(i)].index_w;
	  
	  double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
	    float min_u_dis;
	    if (cur_wire_u < this_index_u-1){
	      min_u_dis = this_index_u -1 - cur_wire_u;
	    }else if (cur_wire_u > this_index_u+1){
	      min_u_dis = cur_wire_u - this_index_u -1;
	    }else{
	      min_u_dis = 0;
	    }
	    float min_v_dis;
	    if (cur_wire_v < this_index_v-1){
	      min_v_dis = this_index_v -1 - cur_wire_v;
	    }else if (cur_wire_v > this_index_v+1){
	      min_v_dis = cur_wire_v - this_index_v -1;
	    }else{
	      min_v_dis = 0;
	    }
	    float min_w_dis;
	    if (cur_wire_w < this_index_w-1){
	      min_w_dis = this_index_w -1 - cur_wire_w;
	    }else if (cur_wire_w > this_index_w+1){
	      min_w_dis = cur_wire_w - this_index_w -1;
	    }else{
	      min_w_dis = 0;
	    }
	    float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	    
	    if ( range_u > 0 && range_v >0 && range_w > 0){
	      float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
	      float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
	      float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
	      float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
	      float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
	      float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;
	      for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
		temp_2dut.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
		temp_2dvt.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
		temp_2dwt.insert(std::make_pair(j,this_time_slice));
	      }
	    } // cuts
	  } // cuts
	  
	} // loop through points
	
      } // distance
    } // no steiner tree cloud
  }

  // just projection ...
  if (temp_2dut.size()==0 && temp_2dvt.size()==0 && temp_2dwt.size()==0){
    //std::cout << "haha " << std::endl;
    std::vector<int> temp_results = ct_point_cloud.convert_3Dpoint_time_ch(p);
    int cur_time_slice = temp_results.at(0);
    int cur_index_u = temp_results.at(1);
    int cur_index_v = temp_results.at(2);
    int cur_index_w = temp_results.at(3);

    for (int i=-time_cut; i!=time_cut+1;i++){
      for (int j=-time_cut; j!=time_cut+1; j++){
	if (abs(i)+abs(j) <= time_cut){
	  temp_2dut.insert(std::make_pair(cur_index_u+i, cur_time_slice+j));
	  temp_2dvt.insert(std::make_pair(cur_index_v+i, cur_time_slice+j));
	  temp_2dwt.insert(std::make_pair(cur_index_w+i, cur_time_slice+j));
	}
      }
    }
  }
  
  
}


void WireCellPID::PR3DCluster::form_map(WireCell::ToyCTPointCloud& ct_point_cloud, WireCell::PointVector& pts,
		  std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge,
		  std::map<int, std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set,
		  std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,
	       double end_point_factor, double mid_point_factor, int nlevel, double time_cut, double charge_cut){
  
  map_3D_2DU_set.clear();
  map_3D_2DV_set.clear();
  map_3D_2DW_set.clear();

  map_2DU_3D_set.clear();
  map_2DV_3D_set.clear();
  map_2DW_3D_set.clear();

 
  WireCell::PointVector saved_pts;
  int count = 0;
  
  // distance ...
  std::vector<double> distances;
  for (size_t i=0;i+1!=pts.size();i++){
    distances.push_back(sqrt(pow(pts.at(i+1).x-pts.at(i).x,2) +
			     pow(pts.at(i+1).y-pts.at(i).y,2) +
			     pow(pts.at(i+1).z-pts.at(i).z,2)));
  }

  

  // start to loop over the path ...
  for (size_t i=0;i!=pts.size();i++){
    double dis_cut;
    if (i==0){
      dis_cut = std::min(distances.at(i) * end_point_factor,4/3.*end_point_factor*units::cm);
    }else if (i+1==pts.size()){
      dis_cut = std::min(distances.back() * end_point_factor,4/3.*end_point_factor*units::cm);
    }else{
      dis_cut = std::min(std::max(distances.at(i-1)*mid_point_factor,distances.at(i)*mid_point_factor),4/3.*mid_point_factor*units::cm);
    }

    std::set<std::pair<int,int> > temp_2dut, temp_2dvt, temp_2dwt;
    form_point_association(pts.at(i), temp_2dut, temp_2dvt, temp_2dwt, ct_point_cloud, dis_cut, nlevel, time_cut);
    // examine ...
    std::vector<int> temp_flag = examine_point_association(temp_2dut, temp_2dvt, temp_2dwt, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,charge_cut);

    //    std::cout << temp_2dut.size() << " " << temp_2dvt.size() << " " << temp_2dwt.size() << " " << temp_flag.at(0) << " " << temp_flag.at(1) << " " << temp_flag.at(2) << std::endl;
    // just projection ...
    // fill the data ...
    if (temp_2dut.size() + temp_2dvt.size() + temp_2dwt.size() > 0){
      map_3D_2DU_set[count] = std::make_pair(temp_2dut,temp_flag.at(0));
      map_3D_2DV_set[count] = std::make_pair(temp_2dvt,temp_flag.at(1));
      map_3D_2DW_set[count] = std::make_pair(temp_2dwt,temp_flag.at(2));
      for (auto it = temp_2dut.begin(); it!=temp_2dut.end();it++){
	if (map_2DU_3D_set.find(*it)==map_2DU_3D_set.end()){
	  std::set<int>  temp_set;
	  temp_set.insert(count);
	  map_2DU_3D_set[*it] = temp_set;
	}else{
	  map_2DU_3D_set[*it].insert(count);
	}
      }
      for (auto it = temp_2dvt.begin(); it!=temp_2dvt.end();it++){
	if (map_2DV_3D_set.find(*it)==map_2DV_3D_set.end()){
	  std::set<int>  temp_set;
	  temp_set.insert(count);
	  map_2DV_3D_set[*it] = temp_set;
	}else{
	  map_2DV_3D_set[*it].insert(count);
	}
      }
      for (auto it = temp_2dwt.begin(); it!=temp_2dwt.end();it++){
	if (map_2DW_3D_set.find(*it)==map_2DW_3D_set.end()){
	  std::set<int>  temp_set;
	  temp_set.insert(count);
	  map_2DW_3D_set[*it] = temp_set;
	}else{
	  map_2DW_3D_set[*it].insert(count);
	}
      }

      
      saved_pts.push_back(pts.at(i));
      count ++;
    }
  }
  
  pts = saved_pts;
 }

void WireCellPID::PR3DCluster::prepare_data(WireCell::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge){
  
  std::vector<int> proj_channel;
  std::vector<int> proj_timeslice;
  std::vector<int> proj_charge;
  std::vector<int> proj_charge_err;
  std::vector<int> proj_flag;
  get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);

  int min_time = 1e9;
  int max_time = -1e9;
  int min_uch = 1e9;
  int max_uch = -1e9;
  int min_vch = 1e9;
  int max_vch = -1e9;
  int min_wch = 1e9;
  int max_wch = -1e9;
  
  // std::cout << proj_charge.size() << " " << proj_flag.size() << std::endl;
  for (size_t i=0;i!=proj_channel.size();i++){
    if (min_time > proj_timeslice.at(i)) min_time = proj_timeslice.at(i);
    if (max_time < proj_timeslice.at(i)) max_time = proj_timeslice.at(i);
	
    if (proj_channel.at(i)<2400){
      map_2D_ut_charge[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));

      if (min_uch > proj_channel.at(i)) min_uch = proj_channel.at(i);
      if (max_uch < proj_channel.at(i)) max_uch = proj_channel.at(i);
    }else if (proj_channel.at(i)<4800){
      map_2D_vt_charge[std::make_pair(proj_channel.at(i)-2400,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));

      if (min_vch > proj_channel.at(i)) min_vch = proj_channel.at(i);
      if (max_vch < proj_channel.at(i)) max_vch = proj_channel.at(i);
    }else{
      map_2D_wt_charge[std::make_pair(proj_channel.at(i)-4800,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));
      if (min_wch > proj_channel.at(i)) min_wch = proj_channel.at(i);
      if (max_wch < proj_channel.at(i)) max_wch = proj_channel.at(i);
    }
  }

  // flag 0 for dead or overlapping channels
  // flag 1 for good channels
  // flag 2 for isolated channels
  // flag 3 for additional channels ...
  
  // add the rest of live channels within range???
  std::map<std::pair<int,int>, std::pair<double,double> > map_u_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_uch, max_uch, 0);
  std::map<std::pair<int,int>, std::pair<double,double> > map_v_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_vch, max_vch, 1);
  std::map<std::pair<int,int>, std::pair<double,double> > map_w_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_wch, max_wch, 2);

  for (auto it = map_u_tcc.begin(); it!=map_u_tcc.end(); it++){
    if (map_2D_ut_charge.find(std::make_pair(it->first.second, it->first.first))==map_2D_ut_charge.end()){
      map_2D_ut_charge[std::make_pair(it->first.second, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }

  for (auto it = map_v_tcc.begin(); it!=map_v_tcc.end(); it++){
    if (map_2D_vt_charge.find(std::make_pair(it->first.second-2400, it->first.first))==map_2D_vt_charge.end()){
      map_2D_vt_charge[std::make_pair(it->first.second-2400, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }

  for (auto it = map_w_tcc.begin(); it!=map_w_tcc.end(); it++){
    if (map_2D_wt_charge.find(std::make_pair(it->first.second-4800, it->first.first))==map_2D_wt_charge.end()){
      map_2D_wt_charge[std::make_pair(it->first.second, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //      std::cout << it->first.first << " " << it->first.second << std::endl;
    }else{
      //      std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }
  
  
  //  for (auto it = map_w_tcc.begin(); it!=map_w_tcc.end(); it++){
  // std::cout << it->first.first << " " << it->first.second << std::endl;
  // }
  
}

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
