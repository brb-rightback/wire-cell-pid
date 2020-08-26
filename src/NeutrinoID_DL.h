#include "WCPPyUtil/read_numpy.h"
#include "WCPPyUtil/SCN_Vertex.h"

bool WCPPID::NeutrinoID::determine_overall_main_vertex_DL(){
  bool flag_change = false;
  
  // hack for now ... (need to be synced with main code ...)
  double dQdx_scale = 0.1;
  double dQdx_offset = -1000;
  
  // prepare the data for DL training
  std::vector< std::vector<float> > vec_xyzq;
  vec_xyzq.resize(4);

  // go through vertex first
  for ( auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    vec_xyzq[0].push_back(vtx->get_fit_pt().x/units::cm);
    vec_xyzq[1].push_back(vtx->get_fit_pt().y/units::cm);
    vec_xyzq[2].push_back(vtx->get_fit_pt().z/units::cm);
    vec_xyzq[3].push_back(vtx->get_dQ() * dQdx_scale + dQdx_offset);
  }
  // go through segment ...
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *seg = it->first;
    std::vector<WCP::Point>& pts = seg->get_point_vec();
    std::vector<double>& dQ_vec = seg->get_dQ_vec();
    for (size_t i=1; i+1<pts.size();i++){
      vec_xyzq[0].push_back(pts.at(i).x/units::cm);
      vec_xyzq[1].push_back(pts.at(i).y/units::cm);
      vec_xyzq[2].push_back(pts.at(i).z/units::cm);
      vec_xyzq[3].push_back(dQ_vec.at(i) * dQdx_scale + dQdx_offset);
    }
  }
  
   
  // prepare the candidate list ...
  std::vector<WCPPID::ProtoVertex* > cand_vertices;
  for ( auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    cand_vertices.push_back(vtx);
  }

  /*
  const std::string data_path = "../wire-cell-pydata/scn_vtx/nuecc-sample.npz";
    auto x = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "coords", 0, true);
    auto y = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "coords", 1);
    auto z = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "coords", 2);
    auto q = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "ft", 0);
  */

  const std::string module = "SCN_Vertex";
  const std::string function = "SCN_Vertex";
  const std::string weights =
      "./input_data_files/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth";

  auto dnn_vtx = WCPPyUtil::SCN_Vertex(module, function, weights, vec_xyzq, "float32", true);
  if (dnn_vtx.size()!=3) {
    throw std::runtime_error("vtx.size()!=3");
  }
  std::cout << "NeutrinoID_DL: DNN: " ;
  for(auto v : dnn_vtx) {
    std::cout << v << " ";
  }
  std::cout << " cm" << std::endl;
  
  double x_reg = dnn_vtx[0]*units::cm;
  double y_reg = dnn_vtx[1]*units::cm;
  double z_reg = dnn_vtx[2]*units::cm;
  
  // find the main vertex
  double min_dis = 1e9;
  WCPPID::ProtoVertex *min_vertex = 0;
  for (auto it = cand_vertices.begin(); it != cand_vertices.end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    double dis = sqrt(pow(vtx->get_fit_pt().x - x_reg,2) + pow(vtx->get_fit_pt().y-y_reg,2) + pow(vtx->get_fit_pt().z - z_reg,2));
    if (dis < min_dis){
      min_dis = dis;
      min_vertex = vtx;
    }
  }

  {
    auto tmp_point = min_vertex->get_fit_pt();
    std::cout << "NeutrinoID_DL: Hybrid: " << tmp_point.x/units::cm << ", " << tmp_point.y/units::cm << ", " << tmp_point.z/units::cm << " " << min_dis/units::cm << std::endl;
  }

  if (min_vertex !=0){
    bool flag_pass = true;
    int num_bad = 0;
    int num_tracks = 0;
    for (auto it = map_vertex_segments[min_vertex].begin(); it != map_vertex_segments[min_vertex].end(); it++){
      WCPPID::ProtoSegment *current_sg = *it;
      bool flag_start;
      if (current_sg->get_wcpt_vec().front().index == min_vertex->get_wcpt().index)
	flag_start = true;
      else if (current_sg->get_wcpt_vec().back().index == min_vertex->get_wcpt().index)
	flag_start = false;
      
      double length = current_sg->get_length();
      double direct_length = current_sg->get_direct_length();
      double medium_dQ_dx = current_sg->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);//0.85+0.95 *sqrt(25./ (length / units::cm));
      
      if (length > 15*units::cm && (!current_sg->is_dir_weak()) && (!flag_start) && medium_dQ_dx > dQ_dx_cut * 43e3/units::cm) num_bad ++;
      num_tracks ++;
    }
    if (num_bad >0 && num_bad == num_tracks) flag_pass = false;

    // 2 cm cut ..
    if (min_dis > dl_vtx_cut) flag_pass = false;

    //    std::cout <<"kaka: " << dl_vtx_cut/units::cm << std::endl;
    
    // now just do the change ....
    if (flag_pass){
      flag_change = true;
      main_vertex = min_vertex;
    }
  }

  if (flag_change){
    std::cout << "Change to DL vertex! " << std::endl;
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
    
    // hard examination ???
    bool flag_check = examine_direction(main_vertex);
    
    
    
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
  }else{
    std::cout << "Stay with the traditional vertex! " << std::endl;
  }

  return flag_change;
}
