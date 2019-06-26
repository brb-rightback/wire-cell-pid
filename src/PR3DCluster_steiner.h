void WireCellPID::PR3DCluster::Create_steiner_tree(){
  Create_graph();

  // find all the steiner terminal indices ...
  find_steiner_terminals();
  
}

void WireCellPID::PR3DCluster::find_steiner_terminals(){
  // reset ...
  steriner_terminal_indices.clear();
  
  // form all the maps ...
  form_cell_points_map();

  
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
