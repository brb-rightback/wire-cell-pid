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

  std::cout << all_indices.size() << std::endl;
  WCPointCloud<double>::WCPoint& wcp = cloud.pts[(*all_indices.begin())]; 
  calc_charge_wcp(wcp,gds);
  //

  
  
  // form another set with the actual points according to their charge
  // find the vertices with the point
  // loop over the connected vertices (not in the dead or good list)
  // if charge smaller, push into dead list
  // if charge bigger, push current into dead list
  // if the current's charge is the biggest, push into good list 

  std::set<int> peak_indices;
  return peak_indices;
}

double WireCellPID::PR3DCluster::calc_charge_wcp(WireCell::WCPointCloud<double>::WCPoint& wcp, WireCell::GeomDataSource& gds){
  double charge = 0;
  SlimMergeGeomCell *mcell = wcp.mcell;
  
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),wcp.index_u);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),wcp.index_v);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),wcp.index_w);

  double charge_u = mcell->Get_Wire_Charge(uwire);
  double charge_v = mcell->Get_Wire_Charge(vwire);
  double charge_w = mcell->Get_Wire_Charge(wwire);

  std::cout << charge_u << " " << charge_v << " " << charge_w << std::endl;
  // deal with bad planes ... 
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
