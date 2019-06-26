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

  
}

std::set<int> WireCellPID::PR3DCluster::find_peak_point_indices(SMGCSelection mcells, WireCell::GeomDataSource& gds){
  std::set<int> all_indices;
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::set<int>& indices = cell_point_indices_map[mcell];
    all_indices.insert(indices.begin(), indices.end());
  }
  
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();

  //WCPointCloud<double>::WCPoint& wcp = cloud.pts[*it1]; 

  
  // form another set with the actual points according to their charge
  // find the vertices with the point
  // loop over the connected vertices (not in the dead or good list)
  // if charge smaller, push into dead list
  // if charge bigger, push current into dead list
  // if the current's charge is the biggest, push into good list 
  
}

double WireCellPID::PR3DCluster::calc_charge_wcp(WireCell::WCPointCloud<double>::WCPoint& wcp, WireCell::GeomDataSource& gds){
  double charge = 0;
  SlimMergeGeomCell *mcell = wcp.mcell;
  int index_u = wcp.index_u;
  int index_v = wcp.index_v;
  int index_w = wcp.index_w;
  
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
