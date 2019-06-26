void WireCellPID::PR3DCluster::Create_steiner_tree(WireCell::GeomDataSource& gds){
  Create_graph();

  // find all the steiner terminal indices ...
  find_steiner_terminals(gds);
  
}

void WireCellPID::PR3DCluster::find_steiner_terminals(WireCell::GeomDataSource& gds){
  // reset ...
  steriner_terminal_indices.clear();
  
  // form all the maps ...
  form_cell_points_map();

  SMGCSelection temp_mcells;
  temp_mcells.push_back(mcells.at(0));
  find_peak_point_indices(temp_mcells, gds);
}

std::set<int> WireCellPID::PR3DCluster::find_peak_point_indices(SMGCSelection mcells, WireCell::GeomDataSource& gds){
  std::set<int> all_indices;
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::set<int>& indices = cell_point_indices_map[mcell];
    all_indices.insert(indices.begin(), indices.end());
  }
  
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();

  // form another set with the actual points according to their charge
  std::map<int,double> map_index_charge;
  std::set<std::pair<double, int>, std::greater<std::pair<double, int> > > candidates_set;
  for (auto it = all_indices.begin(); it!=all_indices.end(); it++){
    //std::cout << all_indices.size() << std::endl;
    WCPointCloud<double>::WCPoint& wcp = cloud.pts[(*it)];
    double charge = calc_charge_wcp(wcp,gds);
    map_index_charge[(*it)] = charge;
    if (charge > 4000){
      candidates_set.insert(std::make_pair(charge,*it));
    }
    //std::cout <<  << std::endl;
  }

  //std::cout << candidates_set.size() << std::endl;
  std::set<int> peak_indices;
  std::set<int> non_peak_indices;

  typedef boost::property_map<MCUGraph, boost::vertex_index_t>::type IndexMap;
  IndexMap index = get(boost::vertex_index,*graph);
  typedef boost::graph_traits<MCUGraph>::adjacency_iterator adjacency_iterator;
  
  for (auto it = candidates_set.begin(); it!=candidates_set.end(); it++){
    int current_index = it->second;
    double charge = it->first;
    // find the vertices with the point
    std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(current_index,*graph),*graph);
     	for (; neighbors.first!=neighbors.second; ++neighbors.first){
	  std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
	}
    // loop over the connected vertices (not in the dead or good list)
    // if charge smaller, push into dead list
    // if charge bigger, push current into dead list
    // if the current's charge is the biggest, push into good list 
  }
  
  
 
  
 


  return peak_indices;
}

double WireCellPID::PR3DCluster::calc_charge_wcp(WireCell::WCPointCloud<double>::WCPoint& wcp, WireCell::GeomDataSource& gds){
  double charge = 0;
  double ncharge = 0;
  SlimMergeGeomCell *mcell = wcp.mcell;
  
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),wcp.index_u);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),wcp.index_v);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),wcp.index_w);

  double charge_u = mcell->Get_Wire_Charge(uwire);
  double charge_v = mcell->Get_Wire_Charge(vwire);
  double charge_w = mcell->Get_Wire_Charge(wwire);

  charge += charge_u*charge_u; ncharge ++;
  charge += charge_v*charge_v; ncharge ++;
  charge += charge_w*charge_w; ncharge ++;
  //std::cout << charge_u << " " << charge_v << " " << charge_w << std::endl;

  // deal with bad planes ... 
  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
  for (size_t i=0;i!=bad_planes.size();i++){
    if (bad_planes.at(i)==WirePlaneType_t(0)){
      charge -= charge_u*charge_u; ncharge--;
    }else if (bad_planes.at(i)==WirePlaneType_t(1)){
      charge -= charge_v*charge_v; ncharge--;
    }else if (bad_planes.at(i)==WirePlaneType_t(2)){
      charge -= charge_w*charge_w; ncharge--;
    }
  }
  if (ncharge>0) {
    charge = sqrt(charge/ncharge);
  }else{
    charge = 0;
  }

  // get charge for each indices ...
  // how to average ??? 
  return charge;
}


void WireCellPID::PR3DCluster::form_cell_points_map(){
  cell_point_indices_map.clear();
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::vector<int>& wcps = point_cloud->get_mcell_indices(mcell);
    if (cell_point_indices_map.find(mcell)==cell_point_indices_map.end()){
      std::set<int> point_indices;
      cell_point_indices_map[mcell] = point_indices;
    }
    
    for (auto it1 = wcps.begin(); it1!=wcps.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp = cloud.pts[*it1]; 
      cell_point_indices_map[mcell].insert(wcp.index);
    }
  }

}
