 
void WCPPID::NeutrinoID::determine_overall_main_vertex_DL(){
  // hack for now ... (need to be synced with main code ...)
  double dQdx_scale = 0.1;
  double dQdx_offset = -1000;
  
  // prepare the data for DL training
  std::vector<std::tuple<float, float, float, float> > vec_xyzq;

  // go through vertex first
  for ( auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    vec_xyzq.push_back(std::make_tuple(vtx->get_fit_pt().x/units::cm, vtx->get_fit_pt().y/units::cm, vtx->get_fit_pt().z/units::cm, vtx->get_dQ() * dQdx_scale + dQdx_offset));
  }
  // go through segment ...
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *seg = it->first;
    std::vector<WCP::Point>& pts = seg->get_point_vec();
    std::vector<double>& dQ_vec = seg->get_dQ_vec();
    for (size_t i=1; i+1<pts.size();i++){
      vec_xyzq.push_back(std::make_tuple(pts.at(i).x/units::cm, pts.at(i).y/units::cm, pts.at(i).z/units::cm, dQ_vec.at(i) * dQdx_scale + dQdx_offset));
    }
  }
  
  
  // prepare the candidate list ...
  std::vector<WCPPID::ProtoVertex* > cand_vertices;
  for ( auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    cand_vertices.push_back(vtx);
  }

  // find the vertex region (to be replaced by DL code ...)
  double x_reg = 0*units::cm;
  double y_reg = 0*units::cm;
  double z_reg = 0*units::cm;
  
  // find the main vertex
  double min_dis = 1e9;
  for (auto it = cand_vertices.begin(); it != cand_vertices.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    double dis = sqrt(pow(vtx->get_fit_pt().x - x_reg,2) + pow(vtx->get_fit_pt().y-y_reg,2) + pow(vtx->get_fit_pt().z - z_reg,2));
    if (dis < min_dis){
      min_dis = dis;
      main_vertex = vtx;
    }
  }
  
  // switch the main cluster to the one to main vertex
  WCPPID::PR3DCluster* tmp_new_main_cluster = 0;
  for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
    WCPPID::PR3DCluster *cluster = *it;
    if (cluster->get_cluster_id() == main_vertex->get_cluster_id()){
      tmp_new_main_cluster = cluster;
      break;
    }
  }
  if (tmp_new_main_cluster !=0) swap_main_cluster(tmp_new_main_cluster);
  
  
  // examine track connected to it
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end();it++){
    WCPPID::ProtoSegment *sg = *it;
    auto pair_results = calculate_num_daughter_showers(main_vertex, sg, false);
    if (pair_results.first==1 && sg->get_length() < 1.5*units::cm && sg->get_medium_dQ_dx()/(43e3/units::cm) > 1.6){
      TPCParams& mp = Singleton<TPCParams>::Instance();
      sg->set_particle_type(2212);
      sg->set_particle_mass(mp.get_mass_proton());
      sg->cal_4mom();
      //    std::cout << pair_results.first << " " << sg->get_length()/units::cm << " " << sg->get_medium_dQ_dx()/(43e3/units::cm) << std::endl;
      //    if (sg->get_length()<1*units::cm && pair_results.
    }
  }
  
  //clean up long muons
  {
    std::set<WCPPID::ProtoSegment* > tmp_segments;
    std::set<WCPPID::ProtoVertex* > tmp_vertices;
    for (auto it = segments_in_long_muon.begin(); it!= segments_in_long_muon.end(); it++){
      if ((*it)->get_cluster_id() != main_vertex->get_cluster_id()) tmp_segments.insert(*it);
    }
    for (auto it = vertices_in_long_muon.begin(); it!= vertices_in_long_muon.end(); it++){
      if ((*it)->get_cluster_id() != main_vertex->get_cluster_id()) tmp_vertices.insert(*it);
    }
    for (auto it = tmp_segments.begin(); it!=tmp_segments.end();it++){
      segments_in_long_muon.erase(*it);
    }
    for (auto it = tmp_vertices.begin(); it!=tmp_vertices.end(); it++){
      vertices_in_long_muon.erase(*it);
    }
  }
}