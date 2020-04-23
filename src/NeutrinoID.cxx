#include "WCPPID/NeutrinoID.h"

#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include "WCPPaal/graph_metrics.h"
#include "WCPPaal/steiner_tree_greedy.h"

#include "WCPData/Line.h"


using namespace WCP;

#include "NeutrinoID_proto_vertex.h"
#include "NeutrinoID_improve_vertex.h"
#include "NeutrinoID_examine_structure.h"
#include "NeutrinoID_deghost.h"
#include "NeutrinoID_track_shower.h"
#include "NeutrinoID_energy_reco.h"
#include "NeutrinoID_shower_clustering.h"

WCPPID::NeutrinoID::NeutrinoID(WCPPID::PR3DCluster *main_cluster1, std::vector<WCPPID::PR3DCluster*>& other_clusters1, std::vector<WCPPID::PR3DCluster*>& all_clusters1, WCPPID::ToyFiducial* fid, WCPSst::GeomDataSource& gds, int nrebin, int frame_length, float unit_dis, ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, double flash_time, double offset_x)
  : acc_vertex_id(0)
  , acc_segment_id(0)
  , main_cluster(main_cluster1)
  , other_clusters(other_clusters1)
  , all_clusters(all_clusters1)
  , fid(fid)
  , ct_point_cloud(ct_point_cloud)
  , global_wc_map(global_wc_map)
  , flash_time(flash_time)
  , offset_x(offset_x)
  , type(0)
  , main_vertex(0)
  , main_cluster_initial_pair_vertices(std::make_pair((WCPPID::ProtoVertex*)0, (WCPPID::ProtoVertex*)0))
{
  bool flag_other_clusters = true;
  bool flag_main_cluster = true;

  PR3DCluster* max_length_cluster = 0;
  double max_length = 0;
  
  std::map<WCPPID::PR3DCluster*, double> map_cluster_length;
  std::set<WCPPID::PR3DCluster* > skip_clusters;
  // form id vs. cluster ...
  map_id_cluster[main_cluster->get_cluster_id()] = main_cluster;
  main_cluster->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  {
    std::pair<WCP::WCPointCloud<double>::WCPoint,WCP::WCPointCloud<double>::WCPoint> two_wcps = main_cluster->get_two_boundary_wcps();
    map_cluster_length[main_cluster] = sqrt(pow(two_wcps.first.x - two_wcps.second.x, 2) + pow(two_wcps.first.y - two_wcps.second.y, 2) + pow(two_wcps.first.z - two_wcps.second.z, 2));
  }
  for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
    map_id_cluster[(*it)->get_cluster_id()] = *it;
    (*it)->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
    std::pair<WCP::WCPointCloud<double>::WCPoint,WCP::WCPointCloud<double>::WCPoint> two_wcps = (*it)->get_two_boundary_wcps();
    double length = sqrt(pow(two_wcps.first.x - two_wcps.second.x, 2) + pow(two_wcps.first.y - two_wcps.second.y, 2) + pow(two_wcps.first.z - two_wcps.second.z, 2));
    map_cluster_length[*it] = length;
    if (length > max_length){
      max_length = length;
      max_length_cluster = *it;
    }
  }
  
  // std::cout << main_cluster->get_num_mcells() << " " << map_cluster_length[main_cluster]/units::cm << std::endl;
  // for  (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
  //  std::cout << (*it)->get_num_mcells() << " " << map_cluster_length[*it]/units::cm << std::endl;
  // }
  
  // hack switch 
  // {
  //   auto it = find(other_clusters.begin(), other_clusters.end(), max_length_cluster);
  //   other_clusters.erase(it);
  //   other_clusters.push_back(main_cluster);
  //   main_cluster = max_length_cluster;
  // }
  
  if (flag_main_cluster){
    // find the proto vertex ...
    find_proto_vertex(main_cluster, true, 2);
      
    // deal with shower ...
    clustering_points(main_cluster);
    separate_track_shower(main_cluster);
    determine_direction(main_cluster);
    
    shower_determing_in_main_cluster();
    determine_main_vertex(main_cluster);	
    
    if (max_length > map_cluster_length[main_cluster] * 0.8 ) check_switch_main_cluster(max_length_cluster, other_clusters, skip_clusters);        
    // fit the vertex in 3D 
    improve_vertex(main_cluster);
    clustering_points(main_cluster);
    examine_direction(main_vertex);
  }

  
  // loop over other clusters ...
  if (flag_other_clusters){
    for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
      if (skip_clusters.find(*it) != skip_clusters.end()) continue;
      // do not break track and find other tracks ...
      //      if ((*it)->get_cluster_id()!=50) continue;
      if (!find_proto_vertex(*it, false, 1)) init_point_segment(*it);
      
      clustering_points(*it);
      separate_track_shower(*it);
      determine_direction(*it);

      //      std::cout << (*it)->get_cluster_id() << " " << map_segment_vertices.size() << " " << map_vertex_segments.size() << std::endl;
      
    }
    //  deghost ...
    deghosting();
  }

  
  // for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
  //   WCPPID::ProtoVertex *vtx = it->first;
  //   if (vtx->get_cluster_id() == 50){
  //     std::cout << vtx->get_id() << " " << vtx->get_fit_pt() << " " << vtx->get_wcpt().x << " " << vtx->get_wcpt().y << " " << vtx->get_wcpt().z << std::endl;
  //   }
  // }
  
  
  if (flag_main_cluster){
    // overall
    separate_track_shower();
    // for charge based on calculation ...
    collect_2D_charges();
    // cluster E&M ...
    shower_clustering_with_nv();
  }
  
  //  std::cout << "Final Information: " << std::endl;
  //print_segs_info(main_vertex);

  
  // prepare output ...
  fill_fit_parameters();
}

void WCPPID::NeutrinoID::check_switch_main_cluster(WCPPID::PR3DCluster *max_length_cluster, WCPPID::PR3DClusterSelection& other_clusters, std::set<WCPPID::PR3DCluster*>& skip_clusters ){

  bool flag_switch = false;

  int n_showers = 0;
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end(); it++){
    if ((*it)->get_flag_shower()) n_showers ++;
  }
  if (n_showers == map_vertex_segments[main_vertex].size()) flag_switch = true;

  
  if (flag_switch){
    std::cout << "Switch Main Cluster! " << std::endl;
    auto it = find(other_clusters.begin(), other_clusters.end(), max_length_cluster);
    other_clusters.erase(it);
    skip_clusters.insert(main_cluster);
    other_clusters.push_back(main_cluster);
    main_cluster = max_length_cluster;
    
    // find the proto vertex ...
    find_proto_vertex(main_cluster, true, 2);    
    
    // deal with shower ...
    clustering_points(main_cluster);
    separate_track_shower(main_cluster);
    determine_direction(main_cluster);
    
    shower_determing_in_main_cluster();
    determine_main_vertex(main_cluster);	
  }
}
  

std::pair<int, int> WCPPID::NeutrinoID::count_num_tracks_showers(WCPPID::PR3DCluster* temp_cluster){
  int num[2] = {0,0};
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_flag_shower()) num[1] ++;
    else num[0] ++;
  }
  return std::make_pair(num[0], num[1]);
}



void WCPPID::NeutrinoID::fill_fit_parameters(){
  organize_vertices_segments();
  std::set<WCPPID::PR3DCluster*> clusters_set;
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    WCPPID::PR3DCluster *cluster = map_id_cluster[sg->get_cluster_id()];
    clusters_set.insert(cluster);
  }
  for (auto it = clusters_set.begin(); it!=clusters_set.end(); it++){
    (*it)->set_fit_parameters(proto_vertices, proto_segments, map_vertex_segments, map_segment_vertices);
  }

  PR3DClusterSelection all_clusters = other_clusters;
  all_clusters.push_back(main_cluster);
  for (auto it = all_clusters.begin(); it!=all_clusters.end(); it++){
    WCPPID::PR3DCluster* temp_cluster = *it;
    std::vector<int>& point_sub_cluster_ids = temp_cluster->get_point_sub_cluster_ids();
    std::vector<bool>& point_flag_showers = temp_cluster->get_point_flag_showers();
    if (point_flag_showers.size()==0)
      point_flag_showers.resize(point_sub_cluster_ids.size(), false);
  }
}

void WCPPID::NeutrinoID::clustering_points(WCPPID::PR3DCluster* temp_cluster){
  temp_cluster->clustering_points_master(map_vertex_segments, map_segment_vertices, *ct_point_cloud);

  std::map<int, WCPPID::ProtoSegment*> map_id_seg;
  std::map<WCPPID::ProtoSegment*, int> map_seg_id;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    map_id_seg[sg->get_id()] = sg;
    map_seg_id[sg] = sg->get_id();
    // std::cout << "A: " << sg->get_id() << std::endl;
    sg->reset_associate_points();
  }

  {
    // find the relevant point clouds ...
    WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud()->get_cloud();
    WCP::WC2DPointCloud<double>& cloud_u = temp_cluster->get_point_cloud()->get_cloud_u();
    WCP::WC2DPointCloud<double>& cloud_v = temp_cluster->get_point_cloud()->get_cloud_v();
    WCP::WC2DPointCloud<double>& cloud_w = temp_cluster->get_point_cloud()->get_cloud_w();
    
    std::vector<int>& point_sub_cluster_ids = temp_cluster->get_point_sub_cluster_ids();

    
    for (size_t i=0;i!=point_sub_cluster_ids.size();i++){
      //      std::cout << point_sub_cluster_ids.at(i) << std::endl;
      if (point_sub_cluster_ids.at(i) == -1) continue;
      if (map_id_seg.find(point_sub_cluster_ids.at(i))==map_id_seg.end()) continue;
      map_id_seg[point_sub_cluster_ids.at(i)]->add_associate_point(cloud.pts[i], cloud_u.pts[i], cloud_v.pts[i], cloud_w.pts[i]);
    }
  }

  {
    WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud_steiner()->get_cloud();
    std::vector<int>& point_steiner_sub_cluster_ids = temp_cluster->get_point_steiner_sub_cluster_ids();
    for (size_t i=0;i!=point_steiner_sub_cluster_ids.size();i++){
      if (point_steiner_sub_cluster_ids.at(i) == -1) continue;
      if (map_id_seg.find(point_steiner_sub_cluster_ids.at(i))==map_id_seg.end()) continue;
      map_id_seg[point_steiner_sub_cluster_ids.at(i)]->add_associate_point_steiner(cloud.pts[i]);
    }
  }
  

  // build kdtree
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    ToyPointCloud *pcloud_associate = sg->get_associated_pcloud();
    if (pcloud_associate !=0) pcloud_associate->build_kdtree_index();
    ToyPointCloud *pcloud_associate_steiner = sg->get_associated_pcloud_steiner();
    if (pcloud_associate_steiner !=0) pcloud_associate_steiner->build_kdtree_index();
  }
  
}



WCPPID::NeutrinoID::~NeutrinoID(){
  for (auto it = proto_vertices.begin(); it != proto_vertices.end(); it++){
    delete (*it);
  }
  for (auto it = proto_segments.begin(); it != proto_segments.end(); it++){
    delete (*it);
  }
  for (auto it = showers.begin(); it!= showers.end(); it++){
    delete (*it);
  }
}




  
int WCPPID::NeutrinoID::get_num_segments(WCPPID::ProtoVertex *pv){
  if (map_vertex_segments.find(pv) != map_vertex_segments.end()){
    return map_vertex_segments[pv].size();
  }else{
    return 0;
  }
}



void WCPPID::NeutrinoID::organize_vertices_segments(){
  WCPPID::ProtoVertexSelection temp_vertices = proto_vertices;
  WCPPID::ProtoSegmentSelection temp_segments = proto_segments;
  proto_vertices.clear();
  proto_segments.clear();
  WCPPID::ProtoVertexSet temp_vset;
  WCPPID::ProtoSegmentSet temp_sset;
  WCPPID::ProtoSegmentSet temp_delset;
  for (auto it = temp_vertices.begin(); it!=temp_vertices.end(); it++){
    if (map_vertex_segments.find(*it)!=map_vertex_segments.end() &&
	temp_vset.find(*it)==temp_vset.end()){
      proto_vertices.push_back(*it);
      temp_vset.insert(*it);
      //std::cout << "1: " << *it << std::endl;
      // std::cout <<"1: " << (*it)->get_wcpt().index << std::endl;
    }
  }
  for (auto it = temp_segments.begin(); it!=temp_segments.end(); it++){
    
    if (map_segment_vertices.find(*it)!=map_segment_vertices.end()){
      if (temp_sset.find(*it)==temp_sset.end()){
	proto_segments.push_back(*it);
	temp_sset.insert(*it);
	//std::cout << (*it)->get_wcpt_vec().front().index << " " << (*it)->get_wcpt_vec().back().index << " " << (*it)->get_point_vec().front() << " " << (*it)->get_wcpt_vec().front().x << " " << (*it)->get_wcpt_vec().front().y << " " << (*it)->get_wcpt_vec().front().z << " " << (*it)->get_point_vec().back() << " " << (*it)->get_wcpt_vec().back().x << " " << (*it)->get_wcpt_vec().back().y << " " << (*it)->get_wcpt_vec().back().z << std::endl;
	//std::cout << "2: " << *it << std::endl;
      }
    }else{
      temp_delset.insert(*it);
    }
  }

  for (auto it = temp_delset.begin(); it!=temp_delset.end();it++){
    delete *it;
  }
}


void WCPPID::NeutrinoID::process_main_cluster(){ 

  WCPPID::PR3DCluster *temp_cluster = main_cluster;
  if (temp_cluster->get_point_cloud_steiner()!=0){
    if (temp_cluster->get_point_cloud_steiner()->get_num_points() >= 2){
      std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2); 
      temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
      temp_cluster->cal_shortest_path(wcps.second,2);
    }
    if (temp_cluster->get_path_wcps().size()>=2){
      temp_cluster->collect_charge_trajectory(*ct_point_cloud);
      temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond);
    }
  }
  
}


void WCPPID::NeutrinoID::process_other_clusters(){
  for (auto it1 =other_clusters.begin(); it1!=other_clusters.end();it1++){
    WCPPID::PR3DCluster *temp_cluster = (*it1);
    if (temp_cluster->get_point_cloud_steiner()!=0){
      if (temp_cluster->get_point_cloud_steiner()->get_num_points() >= 2){
    	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2); 
    	temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
    	temp_cluster->cal_shortest_path(wcps.second,2);
      }
      if (temp_cluster->get_path_wcps().size()>=2){
    	temp_cluster->collect_charge_trajectory(*ct_point_cloud);
    	temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond);
      }
    }
  }
}

void WCPPID::NeutrinoID::fill_reco_tree(WCPPID::ProtoSegment* sg, WCRecoTree& rtree){
  rtree.mc_id[rtree.mc_Ntrack]  = sg->get_cluster_id()*1000 + sg->get_id();
  // particle code
  // e- 11  e+ -11
  // muon- 13  muon+ -13
  // gamma 22
  // pi+ 211, pi0 111, pi- -211
  // kaon+ 321, K- -321
  // p  2212
  // n 2112
  
  //    std::cout << sg->cal_kine_range()/units::MeV << " " << sg->cal_kine_dQdx()/units::MeV << std::endl;
  rtree.mc_pdg[rtree.mc_Ntrack] = sg->get_particle_type(); // all muons for now ...
  rtree.mc_process[rtree.mc_Ntrack] = 0;
  rtree.mc_mother[rtree.mc_Ntrack] = 0; 
  rtree.mc_dir_weak[rtree.mc_Ntrack] = sg->is_dir_weak();
  rtree.mc_kine_range[rtree.mc_Ntrack] = sg->cal_kine_range()/units::MeV;
  rtree.mc_kine_dQdx[rtree.mc_Ntrack] = sg->cal_kine_dQdx()/units::MeV;
  rtree.mc_kine_charge[rtree.mc_Ntrack] = cal_kine_charge(sg)/units::MeV;
  rtree.mc_length[rtree.mc_Ntrack] = sg->get_length()/units::cm;
  // std::cout << rtree.mc_id[rtree.mc_Ntrack] << " " << rtree.mc_dir_weak[rtree.mc_Ntrack] << " " << rtree.mc_kine_range[rtree.mc_Ntrack] << " " << rtree.mc_kine_dQdx[rtree.mc_Ntrack] << " " << rtree.mc_kine_charge[rtree.mc_Ntrack] << std::endl;
  
  sg->set_kine_charge( rtree.mc_kine_charge[rtree.mc_Ntrack] * units::MeV);
  //    std::cout << rtree.mc_kine_range[rtree.mc_Ntrack] << " " << rtree.mc_kine_dQdx[rtree.mc_Ntrack] << " " << rtree.mc_kine_charge[rtree.mc_Ntrack] << std::endl;
  
  if (sg->get_flag_dir()==0) return; // no direction not plot 

  else if (sg->get_flag_dir()==1){
    rtree.mc_startXYZT[rtree.mc_Ntrack][0] = sg->get_point_vec().front().x/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][1] = sg->get_point_vec().front().y/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][2] = sg->get_point_vec().front().z/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][3] = 0;
    
    rtree.mc_endXYZT[rtree.mc_Ntrack][0] = sg->get_point_vec().back().x/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][1] = sg->get_point_vec().back().y/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][2] = sg->get_point_vec().back().z/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][3] = 0;
  }else if (sg->get_flag_dir()==-1){
    rtree.mc_startXYZT[rtree.mc_Ntrack][0] = sg->get_point_vec().back().x/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][1] = sg->get_point_vec().back().y/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][2] = sg->get_point_vec().back().z/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][3] = 0;
    
    rtree.mc_endXYZT[rtree.mc_Ntrack][0] = sg->get_point_vec().front().x/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][1] = sg->get_point_vec().front().y/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][2] = sg->get_point_vec().front().z/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][3] = 0;
  }
  
  if (sg->get_particle_4mom(3)>0){
    //      std::cout << sg->get_particle_4mom(0)/units::MeV << " " << sg->get_particle_4mom(1)/units::MeV << " " << sg->get_particle_4mom(2)/units::MeV << " " << sg->get_particle_4mom(3)/units::MeV << " " << sg->get_particle_type() << " " << sg->get_particle_mass() << std::endl;
    rtree.mc_startMomentum[rtree.mc_Ntrack][0] = sg->get_particle_4mom(0)/units::GeV;
    rtree.mc_startMomentum[rtree.mc_Ntrack][1] = sg->get_particle_4mom(1)/units::GeV;
    rtree.mc_startMomentum[rtree.mc_Ntrack][2] = sg->get_particle_4mom(2)/units::GeV;
    rtree.mc_startMomentum[rtree.mc_Ntrack][3] = sg->get_particle_4mom(3)/units::GeV;
    
    rtree.mc_endMomentum[rtree.mc_Ntrack][0] = 0;
    rtree.mc_endMomentum[rtree.mc_Ntrack][1] = 0;
    rtree.mc_endMomentum[rtree.mc_Ntrack][2] = 0;
    rtree.mc_endMomentum[rtree.mc_Ntrack][3] = sg->get_particle_mass()/units::GeV;
  }else{
    rtree.mc_startMomentum[rtree.mc_Ntrack][0] = 0;
    rtree.mc_startMomentum[rtree.mc_Ntrack][1] = 0;
    rtree.mc_startMomentum[rtree.mc_Ntrack][2] = 0;
    rtree.mc_startMomentum[rtree.mc_Ntrack][3] = 0;
    
    rtree.mc_endMomentum[rtree.mc_Ntrack][0] = 0;
    rtree.mc_endMomentum[rtree.mc_Ntrack][1] = 0;
    rtree.mc_endMomentum[rtree.mc_Ntrack][2] = 0;
    rtree.mc_endMomentum[rtree.mc_Ntrack][3] = 0;
  }
  rtree.mc_Ntrack++;
  rtree.mc_daughters->resize(rtree.mc_Ntrack);
}

void WCPPID::NeutrinoID::fill_reco_tree(WCPPID::WCShower* shower, WCRecoTree& rtree){
  //  fill_reco_tree(shower->get_start_segment(), rtree);

  WCPPID::ProtoSegment *sg1 = shower->get_start_segment();
  
  rtree.mc_id[rtree.mc_Ntrack]  = sg1->get_cluster_id()*1000 + sg1->get_id();
  rtree.mc_pdg[rtree.mc_Ntrack] = shower->get_particle_type();
  rtree.mc_process[rtree.mc_Ntrack] = 0;
  
  rtree.mc_kine_range[rtree.mc_Ntrack] = shower->get_kine_range()/units::MeV;
  rtree.mc_kine_dQdx[rtree.mc_Ntrack] = shower->get_kine_dQdx()/units::MeV;
  rtree.mc_kine_charge[rtree.mc_Ntrack] = shower->get_kine_charge()/units::MeV;
  rtree.mc_length[rtree.mc_Ntrack] = shower->get_total_length()/units::cm;
  
  rtree.mc_startXYZT[rtree.mc_Ntrack][0] = shower->get_start_point().x/units::cm;
  rtree.mc_startXYZT[rtree.mc_Ntrack][1] = shower->get_start_point().y/units::cm;
  rtree.mc_startXYZT[rtree.mc_Ntrack][2] = shower->get_start_point().z/units::cm;
  rtree.mc_startXYZT[rtree.mc_Ntrack][3] = 0;
  
  rtree.mc_endXYZT[rtree.mc_Ntrack][0] = shower->get_end_point().x/units::cm;
  rtree.mc_endXYZT[rtree.mc_Ntrack][1] = shower->get_end_point().y/units::cm;
  rtree.mc_endXYZT[rtree.mc_Ntrack][2] = shower->get_end_point().z/units::cm;
  rtree.mc_endXYZT[rtree.mc_Ntrack][3] = 0;

  
  std::pair<ProtoVertex*, int> pair_start_vertex = shower->get_start_vertex();
  if (pair_start_vertex.second ==4 ){
    rtree.mc_dir_weak[rtree.mc_Ntrack] = 1;
    return;
  }else{
    rtree.mc_dir_weak[rtree.mc_Ntrack] = 0;
  }

  double kine_best = shower->get_kine_best();
  if (kine_best ==0 ) kine_best = shower->get_kine_charge();
  rtree.mc_startMomentum[rtree.mc_Ntrack][3] = (kine_best + sg1->get_particle_mass())/units::GeV;
  double momentum = sqrt(pow(kine_best + sg1->get_particle_mass(),2) - pow(sg1->get_particle_mass(),2));
  TVector3 init_dir = shower->get_init_dir();
  init_dir = init_dir.Unit();
  rtree.mc_startMomentum[rtree.mc_Ntrack][0] = momentum * init_dir.X()/units::GeV;
  rtree.mc_startMomentum[rtree.mc_Ntrack][1] = momentum * init_dir.Y()/units::GeV;
  rtree.mc_startMomentum[rtree.mc_Ntrack][2] = momentum * init_dir.Z()/units::GeV;
  

  rtree.mc_mother[rtree.mc_Ntrack] = 0;
  rtree.mc_Ntrack++;
  rtree.mc_daughters->resize(rtree.mc_Ntrack);
  
}

std::pair<int, int> WCPPID::NeutrinoID::fill_pi0_reco_tree(WCPPID::WCShower* shower, WCRecoTree& rtree){

  if (map_pio_id_saved_pair.find( map_shower_pio_id[shower]) == map_pio_id_saved_pair.end()){
    WCPPID::ProtoSegment *sg1 = shower->get_start_segment();
    rtree.mc_id[rtree.mc_Ntrack]  = sg1->get_cluster_id()*1000 + map_shower_pio_id[shower];
    rtree.mc_pdg[rtree.mc_Ntrack] = 111;
    
    rtree.mc_process[rtree.mc_Ntrack] = 0;  
    rtree.mc_kine_range[rtree.mc_Ntrack] = 0;
    rtree.mc_kine_dQdx[rtree.mc_Ntrack] = 0;
    rtree.mc_kine_charge[rtree.mc_Ntrack] = 0;
    rtree.mc_length[rtree.mc_Ntrack] = 0;
    rtree.mc_stopped[rtree.mc_Ntrack-1] = 0;
    
    rtree.mc_startXYZT[rtree.mc_Ntrack][0] = shower->get_start_vertex().first->get_fit_pt().x/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][1] = shower->get_start_vertex().first->get_fit_pt().y/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][2] = shower->get_start_vertex().first->get_fit_pt().z/units::cm;
    rtree.mc_startXYZT[rtree.mc_Ntrack][3] = 0;
    
    rtree.mc_endXYZT[rtree.mc_Ntrack][0] = shower->get_start_vertex().first->get_fit_pt().x/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][1] = shower->get_start_vertex().first->get_fit_pt().y/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][2] = shower->get_start_vertex().first->get_fit_pt().z/units::cm;
    rtree.mc_endXYZT[rtree.mc_Ntrack][3] = 0;
    
    std::pair<ProtoVertex*, int> pair_start_vertex = shower->get_start_vertex();
    
    rtree.mc_dir_weak[rtree.mc_Ntrack] = 0;
    
    rtree.mc_startMomentum[rtree.mc_Ntrack][0] = sqrt(pow(map_pio_id_mass[map_shower_pio_id[shower]]+135*units::MeV,2) - pow(135*units::MeV,2))/units::GeV;
    rtree.mc_startMomentum[rtree.mc_Ntrack][1] = 0;
    rtree.mc_startMomentum[rtree.mc_Ntrack][2] = 0;
    rtree.mc_startMomentum[rtree.mc_Ntrack][3] = (map_pio_id_mass[map_shower_pio_id[shower]]+135*units::MeV)/units::GeV;
    
    rtree.mc_mother[rtree.mc_Ntrack] = 0;
    rtree.mc_Ntrack++;
    rtree.mc_daughters->resize(rtree.mc_Ntrack);
    
    map_pio_id_saved_pair[ map_shower_pio_id[shower] ] = std::make_pair(rtree.mc_Ntrack-1, rtree.mc_id[rtree.mc_Ntrack-1]);
    return std::make_pair(rtree.mc_Ntrack-1, rtree.mc_id[rtree.mc_Ntrack-1]);
  }else{
    return map_pio_id_saved_pair[ map_shower_pio_id[shower] ];
  }
}



int WCPPID::NeutrinoID::fill_psuedo_reco_tree(WCPPID::WCShower* shower, WCRecoTree& rtree){
  
  WCPPID::ProtoSegment *sg1 = shower->get_start_segment();
  
  rtree.mc_id[rtree.mc_Ntrack]  = sg1->get_cluster_id()*1000 + acc_segment_id; acc_segment_id ++;
  if (fabs(shower->get_particle_type())==11 || fabs(shower->get_particle_type())==22){
    rtree.mc_pdg[rtree.mc_Ntrack] = 22;
  }else{
    rtree.mc_pdg[rtree.mc_Ntrack] = 2112;
  }
  rtree.mc_process[rtree.mc_Ntrack] = 0;
  
  rtree.mc_kine_range[rtree.mc_Ntrack] = 0;
  rtree.mc_kine_dQdx[rtree.mc_Ntrack] = 0;
  rtree.mc_kine_charge[rtree.mc_Ntrack] = 0;
  rtree.mc_length[rtree.mc_Ntrack] = 0;
  rtree.mc_stopped[rtree.mc_Ntrack-1] = 0;
  
  rtree.mc_startXYZT[rtree.mc_Ntrack][0] = shower->get_start_vertex().first->get_fit_pt().x/units::cm;
  rtree.mc_startXYZT[rtree.mc_Ntrack][1] = shower->get_start_vertex().first->get_fit_pt().y/units::cm;
  rtree.mc_startXYZT[rtree.mc_Ntrack][2] = shower->get_start_vertex().first->get_fit_pt().z/units::cm;
  rtree.mc_startXYZT[rtree.mc_Ntrack][3] = 0;
  
  rtree.mc_endXYZT[rtree.mc_Ntrack][0] = shower->get_start_point().x/units::cm;
  rtree.mc_endXYZT[rtree.mc_Ntrack][1] = shower->get_start_point().y/units::cm;
  rtree.mc_endXYZT[rtree.mc_Ntrack][2] = shower->get_start_point().z/units::cm;
  rtree.mc_endXYZT[rtree.mc_Ntrack][3] = 0;

  std::pair<ProtoVertex*, int> pair_start_vertex = shower->get_start_vertex();

  rtree.mc_dir_weak[rtree.mc_Ntrack] = 0;
  

  double kine_best = shower->get_kine_best();
  if (kine_best ==0 ) kine_best = shower->get_kine_charge();
  double momentum;
  if (fabs(shower->get_particle_type())==11 || fabs(shower->get_particle_type())==22){
    rtree.mc_startMomentum[rtree.mc_Ntrack][3] = (kine_best)/units::GeV;
    momentum = kine_best;
  }else{
    TPCParams& mp = Singleton<TPCParams>::Instance();
    double mass_neutron = mp.get_mass_neutron();
    rtree.mc_startMomentum[rtree.mc_Ntrack][3] = (kine_best + mass_neutron)/units::GeV;
    momentum = sqrt(pow(kine_best + mass_neutron,2) - pow(mass_neutron,2));
  }
  
  TVector3 init_dir(  rtree.mc_endXYZT[rtree.mc_Ntrack][0] - rtree.mc_startXYZT[rtree.mc_Ntrack][0], rtree.mc_endXYZT[rtree.mc_Ntrack][1] - rtree.mc_startXYZT[rtree.mc_Ntrack][1] , rtree.mc_endXYZT[rtree.mc_Ntrack][2] - rtree.mc_startXYZT[rtree.mc_Ntrack][2]);
  init_dir = init_dir.Unit();
  if (init_dir.Mag()==0) init_dir = shower->get_start_segment()->cal_dir_3vector();
  rtree.mc_startMomentum[rtree.mc_Ntrack][0] = momentum * init_dir.X()/units::GeV;
  rtree.mc_startMomentum[rtree.mc_Ntrack][1] = momentum * init_dir.Y()/units::GeV;
  rtree.mc_startMomentum[rtree.mc_Ntrack][2] = momentum * init_dir.Z()/units::GeV;
  

  rtree.mc_mother[rtree.mc_Ntrack] = 0;
  rtree.mc_Ntrack++;
  rtree.mc_daughters->resize(rtree.mc_Ntrack);

  return rtree.mc_id[rtree.mc_Ntrack-1];
}


void WCPPID::NeutrinoID::fill_reco_simple_tree(WCPPID::WCRecoTree& rtree){
  
  //start to fill  
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment* sg= it->first;
    // only save the main cluster
    if (sg->get_cluster_id() != main_cluster->get_cluster_id()) continue;
   
    fill_reco_tree(sg, rtree);
    std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> pair_vertices = find_vertices(sg);
    if (map_vertex_segments[pair_vertices.first].size()==1 || map_vertex_segments[pair_vertices.second].size()==1 ){
      rtree.mc_stopped[rtree.mc_Ntrack-1] = 1;
    }else{
      rtree.mc_stopped[rtree.mc_Ntrack-1] = 0;
    }
  }
  // std::cout << rtree.mc_Ntrack << std::endl;
}



void WCPPID::NeutrinoID::fill_particle_tree(WCPPID::WCRecoTree& rtree){
  if (main_vertex ==0 ) return;

  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment* sg= it->first;
    if (map_segment_in_shower.find(sg)!=map_segment_in_shower.end()) continue;
    fill_reco_tree(sg, rtree);
    std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> pair_vertices = find_vertices(sg);
    if (map_vertex_segments[pair_vertices.first].size()==1 || map_vertex_segments[pair_vertices.second].size()==1 ){
      rtree.mc_stopped[rtree.mc_Ntrack-1] = 1;
    }else{
      rtree.mc_stopped[rtree.mc_Ntrack-1] = 0;
    }
    //    std::cout << "kak " << sg->get_cluster_id() << " " << sg->get_id() << std::endl;
  }
  for (auto it = showers.begin(); it!=showers.end();it++){
    fill_reco_tree(*it, rtree);
    rtree.mc_stopped[rtree.mc_Ntrack-1] = 0;
    //std::cout << "gag " << *it << " " << (*it)->get_start_segment()->get_id() << std::endl;
  }

    // id vs. rtree id
  std::map<int, int> map_sgid_rtid;
  std::map<int, int> map_rtid_sgid;
  for (int i=0;i!=rtree.mc_Ntrack;i++){
    int sgid = rtree.mc_id[i];
    map_sgid_rtid[sgid] = i;
    map_rtid_sgid[i] = sgid;
  }
  
  // id vs. sg 
  std::map<int, WCPPID::ProtoSegment*> map_sgid_sg;
  std::map<WCPPID::ProtoSegment*, int> map_sg_sgid;
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment* sg= it->first;
    // only save the main cluster
    int sgid = sg->get_cluster_id()*1000 + sg->get_id();
    map_sg_sgid[sg] = sgid;
    map_sgid_sg[sgid] = sg;
  }
  
  // main_vertex figure out the daughters and mother ...
  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;

  for (auto it = showers.begin(); it!=showers.end();it++){
    (*it)->fill_sets(used_vertices, used_segments, false);
  }
  
  
  std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > segments_to_be_examined;
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    // parent are all zero now ...
    used_segments.insert(*it);
    WCPPID::ProtoVertex *other_vertex = find_other_vertex(*it, main_vertex);
    segments_to_be_examined.push_back(std::make_pair(other_vertex, *it));
  }
  used_vertices.insert(main_vertex);
  //std::cout << segments_to_be_examined.size() << " " << used_vertices.size() << " " << used_segments.size() << std::endl;
  
  while(segments_to_be_examined.size()>0){
    std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > temp_segments;
    for (auto it = segments_to_be_examined.begin(); it!= segments_to_be_examined.end(); it++){
      WCPPID::ProtoVertex *curr_vtx = it->first;
      WCPPID::ProtoSegment *prev_sg = it->second;
      if (used_vertices.find(curr_vtx)!=used_vertices.end()) continue;

      for (auto it1 = map_vertex_segments[curr_vtx].begin(); it1!=map_vertex_segments[curr_vtx].end(); it1++){
	WCPPID::ProtoSegment *curr_sg = *it1;
	if (used_segments.find(curr_sg)!=used_segments.end()) continue;
	used_segments.insert(curr_sg);
	// set mother ...
	rtree.mc_mother[map_sgid_rtid[map_sg_sgid[curr_sg]]] = map_sg_sgid[prev_sg];
	// set daughters ...
	rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(map_sg_sgid[curr_sg]);
	WCPPID::ProtoVertex *other_vertex = find_other_vertex(curr_sg, curr_vtx);
	if (used_vertices.find(other_vertex) == used_vertices.end())
	  temp_segments.push_back(std::make_pair(other_vertex, curr_sg));
      }
      used_vertices.insert(curr_vtx);
    }
    segments_to_be_examined = temp_segments;
  }

  // Now the showers
  for (auto it = showers.begin(); it!= showers.end(); it++){
    WCPPID::WCShower *shower = *it;
    WCPPID::ProtoSegment* curr_sg = shower->get_start_segment();
    std::pair<ProtoVertex*, int> pair_vertex = shower->get_start_vertex();

    //    std::cout << pair_vertex.second << std::endl;

    if (pi0_showers.find(shower)==pi0_showers.end()){ // not pi0
      if (pair_vertex.second == 1){ // direct connection
	if (pair_vertex.first == main_vertex){
	  rtree.mc_mother[ map_sgid_rtid[map_sg_sgid[curr_sg]] ] = 0;
	}else{
	  WCPPID::ProtoSegment* prev_sg = 0;
	  if (map_vertex_in_shower.find(pair_vertex.first) != map_vertex_in_shower.end()){
	    prev_sg = map_vertex_in_shower[pair_vertex.first]->get_start_segment();
	  }else{
	    prev_sg = find_incoming_segment(pair_vertex.first);
	  }
	  // set mother ...
	  rtree.mc_mother[map_sgid_rtid[map_sg_sgid[curr_sg]]] = map_sg_sgid[prev_sg];
	  // set daughters ...
	  rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(map_sg_sgid[curr_sg]);
	}
      }else if (pair_vertex.second == 2 || pair_vertex.second == 3){
	int psuedo_particle_id = fill_psuedo_reco_tree(shower, rtree);

	if (pair_vertex.first == main_vertex){
	  rtree.mc_mother[rtree.mc_Ntrack-1] = 0;
	  rtree.mc_daughters->at(rtree.mc_Ntrack-1).push_back(map_sg_sgid[curr_sg]);
	  rtree.mc_mother[ map_sgid_rtid[map_sg_sgid[curr_sg]] ] = psuedo_particle_id;
	}else{
	  WCPPID::ProtoSegment* prev_sg = 0;
	  if (map_vertex_in_shower.find(pair_vertex.first) != map_vertex_in_shower.end()){
	    prev_sg = map_vertex_in_shower[pair_vertex.first]->get_start_segment();
	  }else{
	    prev_sg = find_incoming_segment(pair_vertex.first);
	  }
	  

	  rtree.mc_mother[rtree.mc_Ntrack-1] = map_sg_sgid[prev_sg];
	  rtree.mc_daughters->at(rtree.mc_Ntrack-1).push_back(map_sg_sgid[curr_sg]);
	  rtree.mc_mother[ map_sgid_rtid[map_sg_sgid[curr_sg]] ] = psuedo_particle_id;
	  rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(psuedo_particle_id);
	}
      }
    }else{

       std::pair<ProtoVertex*, int> pair_vertex = shower->get_start_vertex();
       // pi0 ...
       // create a pi0 ...
       std::pair<int, int> pio_info_pair = fill_pi0_reco_tree(shower, rtree);
       
       int psuedo_particle_id = fill_psuedo_reco_tree(shower, rtree);
       
       
       if (pair_vertex.first == main_vertex){
	 // pio
	 rtree.mc_mother[pio_info_pair.first] = 0;
	 rtree.mc_daughters->at(pio_info_pair.first).push_back(psuedo_particle_id);
	 // gamma
	 rtree.mc_mother[rtree.mc_Ntrack-1] = pio_info_pair.second;
	 rtree.mc_daughters->at(rtree.mc_Ntrack-1).push_back(map_sg_sgid[curr_sg]);
	 // electron 
	 rtree.mc_mother[ map_sgid_rtid[map_sg_sgid[curr_sg]] ] = psuedo_particle_id;
       }else{
	 WCPPID::ProtoSegment* prev_sg = 0;
	 if (map_vertex_in_shower.find(pair_vertex.first) != map_vertex_in_shower.end()){
	   prev_sg = map_vertex_in_shower[pair_vertex.first]->get_start_segment();
	 }else{
	   prev_sg = find_incoming_segment(pair_vertex.first);
	 }
	 //	 std::cout << map_sg_sgid[prev_sg] << " " <<  psuedo_particle_id << " " <<  pio_info_pair.first << " " << pio_info_pair.second << std::endl;


	 // pio
	 rtree.mc_mother[pio_info_pair.first] = map_sg_sgid[prev_sg];
	 rtree.mc_daughters->at(pio_info_pair.first).push_back(psuedo_particle_id);
	 // gamma
	 rtree.mc_mother[rtree.mc_Ntrack-1] =  pio_info_pair.second;
	 rtree.mc_daughters->at(rtree.mc_Ntrack-1).push_back(map_sg_sgid[curr_sg]);
	 // electron
	 rtree.mc_mother[ map_sgid_rtid[map_sg_sgid[curr_sg]] ] = psuedo_particle_id;
	 if (find(rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).begin(), rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).end(), pio_info_pair.second) == rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).end())
	   rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(pio_info_pair.second);
	 //	 rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(psuedo_particle_id);
       }
       
      // create the second psuedo particle

      // fill the final daughter ...

      
    }
  }
}

WCPPID::ProtoSegment* WCPPID::NeutrinoID::find_incoming_segment(WCPPID::ProtoVertex *vtx){
  WCPPID::ProtoSegment* sg = 0;
  // find the first segment ...
  for ( auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
    WCPPID::ProtoSegment* current_sg = (*it);
    bool flag_start;
    if (current_sg->get_wcpt_vec().front().index == vtx->get_wcpt().index)
      flag_start = true;
    else if (current_sg->get_wcpt_vec().back().index == vtx->get_wcpt().index)
      flag_start = false;

    if (flag_start && current_sg->get_flag_dir()==-1
	|| (!flag_start) && current_sg->get_flag_dir()==1 ){
      sg = current_sg;
      break;
    }
    
  }
  
  return sg;
}


void WCPPID::NeutrinoID::fill_proto_main_tree(WCPPID::WCRecoTree& rtree){
  fill_reco_simple_tree(rtree);

  // id vs. rtree id
  std::map<int, int> map_sgid_rtid;
  std::map<int, int> map_rtid_sgid;
  for (int i=0;i!=rtree.mc_Ntrack;i++){
    int sgid = rtree.mc_id[i];
    map_sgid_rtid[sgid] = i;
    map_rtid_sgid[i] = sgid;
  }
  
  // id vs. sg 
  std::map<int, WCPPID::ProtoSegment*> map_sgid_sg;
  std::map<WCPPID::ProtoSegment*, int> map_sg_sgid;
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment* sg= it->first;
    // only save the main cluster
    if (sg->get_cluster_id() != main_cluster->get_cluster_id()) continue;
    int sgid = sg->get_cluster_id()*1000 + sg->get_id();

    map_sg_sgid[sg] = sgid;
    map_sgid_sg[sgid] = sg;
  }

  
  // main_vertex figure out the daughters and mother ...
  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;

  std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > segments_to_be_examined;
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    // parent are all zero now ...
    used_segments.insert(*it);
    WCPPID::ProtoVertex *other_vertex = find_other_vertex(*it, main_vertex);
    segments_to_be_examined.push_back(std::make_pair(other_vertex, *it));
  }
  used_vertices.insert(main_vertex);
  //std::cout << segments_to_be_examined.size() << " " << used_vertices.size() << " " << used_segments.size() << std::endl;
  
  while(segments_to_be_examined.size()>0){
    std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > temp_segments;
    for (auto it = segments_to_be_examined.begin(); it!= segments_to_be_examined.end(); it++){
      WCPPID::ProtoVertex *curr_vtx = it->first;
      WCPPID::ProtoSegment *prev_sg = it->second;
      if (used_vertices.find(curr_vtx)!=used_vertices.end()) continue;

      for (auto it1 = map_vertex_segments[curr_vtx].begin(); it1!=map_vertex_segments[curr_vtx].end(); it1++){
	WCPPID::ProtoSegment *curr_sg = *it1;
	if (used_segments.find(curr_sg)!=used_segments.end()) continue;
	used_segments.insert(curr_sg);
	// set mother ...
	rtree.mc_mother[map_sgid_rtid[map_sg_sgid[curr_sg]]] = map_sg_sgid[prev_sg];
	// set daughters ...
	rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(map_sg_sgid[curr_sg]);
	WCPPID::ProtoVertex *other_vertex = find_other_vertex(curr_sg, curr_vtx);
	if (used_vertices.find(other_vertex) == used_vertices.end())
	  temp_segments.push_back(std::make_pair(other_vertex, curr_sg));
      }
      used_vertices.insert(curr_vtx);
    }
    segments_to_be_examined = temp_segments;
  }
}


void WCPPID::NeutrinoID::fill_skeleton_info_magnify(int mother_cluster_id, WCPPID::WCPointTree& ptree, TTree *T, double dQdx_scale, double dQdx_offset, bool flag_skip_vertex){

  WCPPID::PR3DClusterSelection clusters_vec;
  clusters_vec.push_back(main_cluster);
  for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
    clusters_vec.push_back(*it);
  }
  
  ptree.reco_chi2 = 1;
  ptree.reco_mother_cluster_id = mother_cluster_id;
  

  for (auto it0 = clusters_vec.begin(); it0!=clusters_vec.end(); it0++){
    WCPPID::PR3DCluster* cluster = *it0;
    
    if (!flag_skip_vertex){
      // fill vertex information
      ptree.reco_flag_track_shower = 0;
      for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
	WCPPID::ProtoVertex *vtx = it->first;
	if (vtx->get_cluster_id() != cluster->get_cluster_id()) continue;
	if (map_vertex_in_shower.find(vtx)!=map_vertex_in_shower.end()) continue;
	ptree.reco_cluster_id = vtx->get_cluster_id();
	ptree.reco_proto_cluster_id = -1;
	ptree.reco_particle_id = -1;
	ptree.reco_ndf = vtx->get_cluster_id();
	ptree.reco_flag_vertex = 1;
      
	ptree.reco_x = vtx->get_fit_pt().x/units::cm;
	ptree.reco_y = vtx->get_fit_pt().y/units::cm;
	ptree.reco_z = vtx->get_fit_pt().z/units::cm;
	ptree.reco_dQ = vtx->get_dQ() * dQdx_scale + dQdx_offset;
	ptree.reco_dx = vtx->get_dx()/units::cm;
	ptree.reco_pu = vtx->get_pu();
	ptree.reco_pv = vtx->get_pv();
	ptree.reco_pw = vtx->get_pw();
	ptree.reco_pt = vtx->get_pt();    
	ptree.reco_reduced_chi2 = vtx->get_reduced_chi2();
	ptree.reco_rr = -1; // no residual range
	
	T->Fill();
      }
    }


    for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *seg = it->first;
      if (seg->get_cluster_id() != cluster->get_cluster_id()) continue;
      std::vector<WCP::Point>& pts = seg->get_point_vec();
      std::vector<double>& dQ_vec = seg->get_dQ_vec();
      std::vector<double>& dx_vec = seg->get_dx_vec();
      std::vector<double>& pu_vec = seg->get_pu_vec();
      std::vector<double>& pv_vec = seg->get_pv_vec();
      std::vector<double>& pw_vec = seg->get_pw_vec();
      std::vector<double>& pt_vec = seg->get_pt_vec();
      std::vector<double>& reduced_chi2_vec = seg->get_reduced_chi2_vec();
      
      ptree.reco_cluster_id = seg->get_cluster_id();
      ptree.reco_ndf = seg->get_cluster_id();
      ptree.reco_proto_cluster_id = seg->get_cluster_id()*1000 + seg->get_id();

      // calculate rr ...
      std::vector<double> rr_vec(pts.size(),0);
      {
	WCPPID::ProtoVertex *start_v = 0, *end_v = 0;
	for (auto it1 = map_segment_vertices[it->first].begin(); it1!=map_segment_vertices[it->first].end(); it1++){
	  if ((*it1)->get_wcpt().index == it->first->get_wcpt_vec().front().index) start_v = *it1;
	  if ((*it1)->get_wcpt().index == it->first->get_wcpt_vec().back().index) end_v = *it1;
	}
	int start_n = map_vertex_segments[start_v].size();
	int end_n = map_vertex_segments[end_v].size();

	std::vector<double> L(pts.size(),0);
	double acc_length = 0;
	for (size_t i=0;i+1<pts.size();i++){
	  acc_length += sqrt(pow(pts.at(i+1).x-pts.at(i).x,2) + pow(pts.at(i+1).y - pts.at(i).y ,2) + pow(pts.at(i+1).z - pts.at(i).z,2)); 
	  L.at(i+1) = acc_length;
	}
	if ((seg)->get_flag_dir()==1){// forward direction
	  for (size_t i=0;i!=pts.size();i++){
	    rr_vec.at(pts.size()-1-i) = L.back() - L.at(pts.size()-1-i);
	  }
	}else if ((seg)->get_flag_dir()==-1){ // reverse direction
	  rr_vec = L;
	}else{
	  rr_vec = L; // quick
	}
	if (start_n>1) rr_vec.front() = -1;
	if (end_n >1) rr_vec.back() = -1;
	
      }
      
      // hack for now ...
      //      ptree.reco_particle_id = seg->get_cluster_id()*1000 + seg->get_id();
      if (seg->get_flag_shower()){
	ptree.reco_flag_track_shower = 1;
      }else{
	ptree.reco_flag_track_shower = 0;
      }

      if (ptree.reco_flag_track_shower==1){
	ptree.reco_particle_id =1;
      }else{
	ptree.reco_particle_id =4;
      }
      
      ptree.reco_flag_vertex = 0;
      for (size_t i=0;i!=pts.size();i++){
	ptree.reco_x = pts.at(i).x/units::cm;
	ptree.reco_y = pts.at(i).y/units::cm;
	ptree.reco_z = pts.at(i).z/units::cm;
	ptree.reco_dQ = dQ_vec.at(i) * dQdx_scale + dQdx_offset;
	ptree.reco_dx = dx_vec.at(i)/units::cm;
	ptree.reco_pu = pu_vec.at(i);
	ptree.reco_pv = pv_vec.at(i);
	ptree.reco_pw = pw_vec.at(i);
	ptree.reco_pt = pt_vec.at(i);    
	ptree.reco_reduced_chi2 = reduced_chi2_vec.at(i);
	ptree.reco_rr = rr_vec.at(i)/units::cm; // no residual range
	T->Fill();
      }
    }
  }
  
}


void WCPPID::NeutrinoID::fill_skeleton_info(int mother_cluster_id, WCPPID::WCPointTree& ptree, TTree *T, double dQdx_scale, double dQdx_offset, bool flag_skip_vertex){
  
  ptree.reco_chi2 = 1;
  ptree.reco_mother_cluster_id = mother_cluster_id;
  
    
  if (!flag_skip_vertex){
    // fill vertex information
    ptree.reco_flag_track_shower = 0;
    for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
      WCPPID::ProtoVertex *vtx = it->first;
      ptree.reco_cluster_id = vtx->get_cluster_id();
      ptree.reco_proto_cluster_id = -1;
      ptree.reco_particle_id = -1;

      ptree.reco_ndf = vtx->get_cluster_id();
      ptree.reco_flag_vertex = 1;
      ptree.reco_x = vtx->get_fit_pt().x/units::cm;
      ptree.reco_y = vtx->get_fit_pt().y/units::cm;
      ptree.reco_z = vtx->get_fit_pt().z/units::cm;
      ptree.reco_dQ = vtx->get_dQ() * dQdx_scale + dQdx_offset;
      ptree.reco_dx = vtx->get_dx()/units::cm;
      ptree.reco_pu = vtx->get_pu();
      ptree.reco_pv = vtx->get_pv();
      ptree.reco_pw = vtx->get_pw();
      ptree.reco_pt = vtx->get_pt();    
      ptree.reco_reduced_chi2 = vtx->get_reduced_chi2();
            
      T->Fill();
    }
  }
  
  
  for (auto it = map_vertex_to_shower.begin(); it!=map_vertex_to_shower.end(); it++){
    WCPPID::ProtoVertex* vtx = it->first;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      WCPPID::WCShower *shower = (*it1);
    
      if (map_vertex_segments[vtx].find(shower->get_start_segment())!=map_vertex_segments[vtx].end()){
	ptree.reco_cluster_id = vtx->get_cluster_id();
	
	ptree.reco_proto_cluster_id = shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id();
	ptree.reco_particle_id = shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id();
	
	ptree.reco_ndf = vtx->get_cluster_id();
	ptree.reco_flag_vertex = 1;
	
	ptree.reco_x = vtx->get_fit_pt().x/units::cm;
	ptree.reco_y = vtx->get_fit_pt().y/units::cm;
	ptree.reco_z = vtx->get_fit_pt().z/units::cm;
	ptree.reco_dQ = vtx->get_dQ() * dQdx_scale + dQdx_offset;
	ptree.reco_dx = vtx->get_dx()/units::cm;
	ptree.reco_pu = vtx->get_pu();
	ptree.reco_pv = vtx->get_pv();
	ptree.reco_pw = vtx->get_pw();
	ptree.reco_pt = vtx->get_pt();    
	ptree.reco_reduced_chi2 = vtx->get_reduced_chi2();
      
	T->Fill();
      }
    }
  }
  
  for (auto it = map_vertex_in_shower.begin(); it != map_vertex_in_shower.end(); it++){
     WCPPID::ProtoVertex *vtx = it->first;
     WCPPID::WCShower *shower = it->second;

     ptree.reco_proto_cluster_id = shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id();
     if (shower->get_flag_shower()){
       ptree.reco_flag_track_shower = 1;
     }else{
       ptree.reco_flag_track_shower = 0;
     }
     ptree.reco_particle_id = shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id();
     ptree.reco_flag_vertex = 0;

     ptree.reco_x = vtx->get_fit_pt().x/units::cm;
     ptree.reco_y = vtx->get_fit_pt().y/units::cm;
     ptree.reco_z = vtx->get_fit_pt().z/units::cm;
     ptree.reco_dQ = vtx->get_dQ() * dQdx_scale + dQdx_offset;
     ptree.reco_dx = vtx->get_dx()/units::cm;
     ptree.reco_pu = vtx->get_pu();
     ptree.reco_pv = vtx->get_pv();
     ptree.reco_pw = vtx->get_pw();
     ptree.reco_pt = vtx->get_pt();    
     ptree.reco_reduced_chi2 = vtx->get_reduced_chi2();

     
     // std::cout << vtx->get_cluster_id() << " " << vtx->get_id() << " " << vtx->get_fit_pt() << " " << shower->get_start_segment()->get_cluster_id() << " " << shower->get_start_segment()->get_id() << std::endl;
     
     // std::cout <<  ptree.reco_proto_cluster_id << std::endl;
     // std::cout << shower->get_start_segment()->get_length()/units::cm << " " << vtx << " " << shower << std::endl;
     T->Fill();
  }
  
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *seg = it->first;

    std::vector<WCP::Point>& pts = seg->get_point_vec();
    std::vector<double>& dQ_vec = seg->get_dQ_vec();
    std::vector<double>& dx_vec = seg->get_dx_vec();
    std::vector<double>& pu_vec = seg->get_pu_vec();
    std::vector<double>& pv_vec = seg->get_pv_vec();
    std::vector<double>& pw_vec = seg->get_pw_vec();
    std::vector<double>& pt_vec = seg->get_pt_vec();
    std::vector<double>& reduced_chi2_vec = seg->get_reduced_chi2_vec();

    ptree.reco_cluster_id = seg->get_cluster_id();
    ptree.reco_ndf = seg->get_cluster_id();
    
    auto it1 = map_segment_in_shower.find(seg);
    if (it1!=map_segment_in_shower.end()){
      WCPPID::WCShower *shower = it1->second;
      
      ptree.reco_proto_cluster_id = shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id();
      if (shower->get_flag_shower()){
	ptree.reco_flag_track_shower = 1;
      }else{
	ptree.reco_flag_track_shower = 0;
      }
      ptree.reco_particle_id = shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id();
      ptree.reco_flag_vertex = 0;

      for (size_t i=1; i+1<pts.size();i++){
	ptree.reco_x = pts.at(i).x/units::cm;
	ptree.reco_y = pts.at(i).y/units::cm;
	ptree.reco_z = pts.at(i).z/units::cm;
	ptree.reco_dQ = dQ_vec.at(i) * dQdx_scale + dQdx_offset;
	ptree.reco_dx = dx_vec.at(i)/units::cm;
	ptree.reco_pu = pu_vec.at(i);
	ptree.reco_pv = pv_vec.at(i);
	ptree.reco_pw = pw_vec.at(i);
	ptree.reco_pt = pt_vec.at(i);    
	ptree.reco_reduced_chi2 = reduced_chi2_vec.at(i);
	T->Fill();
      }

          
    }else{
      
      ptree.reco_proto_cluster_id = seg->get_cluster_id()*1000 + seg->get_id();

      if (seg->get_flag_shower()){
	ptree.reco_flag_track_shower = 1;
      }else{
	ptree.reco_flag_track_shower = 0;
      }
      
      ptree.reco_particle_id = ptree.reco_particle_id = seg->get_cluster_id()*1000 + seg->get_id();
      
      ptree.reco_flag_vertex = 0;
      for (size_t i=0;i!=pts.size();i++){
	ptree.reco_x = pts.at(i).x/units::cm;
	ptree.reco_y = pts.at(i).y/units::cm;
	ptree.reco_z = pts.at(i).z/units::cm;
	ptree.reco_dQ = dQ_vec.at(i) * dQdx_scale + dQdx_offset;
	ptree.reco_dx = dx_vec.at(i)/units::cm;
	ptree.reco_pu = pu_vec.at(i);
	ptree.reco_pv = pv_vec.at(i);
	ptree.reco_pw = pw_vec.at(i);
	ptree.reco_pt = pt_vec.at(i);    
	ptree.reco_reduced_chi2 = reduced_chi2_vec.at(i);
	T->Fill();
      }
    }
    
  }
  
  
}


void WCPPID::NeutrinoID::fill_point_info(int mother_cluster_id, WCPPID::WCPointTree& ptree, TTree *T){
  
  ptree.reco_mother_cluster_id = mother_cluster_id;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *seg = it->first;

    auto it1 = map_segment_in_shower.find(seg);
    if (it1 != map_segment_in_shower.end()){
      ptree.reco_proto_cluster_id = (it1->second)->get_start_segment()->get_cluster_id() * 1000 + (it1->second)->get_start_segment()->get_id();
      if (it1->second->get_flag_shower()){
	ptree.reco_flag_track_shower = 1;
	ptree.reco_flag_track_shower_charge = 15000;
      }else{
	ptree.reco_flag_track_shower = 0;
	ptree.reco_flag_track_shower_charge = 0;
      }
    }else{
      ptree.reco_proto_cluster_id = seg->get_cluster_id() * 1000 + seg->get_id();
      if (seg->get_flag_shower()){
	ptree.reco_flag_track_shower = 1;
	ptree.reco_flag_track_shower_charge = 15000;
      }else{
	ptree.reco_flag_track_shower = 0;
	ptree.reco_flag_track_shower_charge = 0;
      }
    }
    ptree.reco_cluster_id = seg->get_cluster_id();
    
    ToyPointCloud *pcloud =  seg->get_associated_pcloud();
    if (pcloud != 0){
      auto cloud = pcloud->get_cloud();
      //    std::cout << cloud.pts.size() << std::endl;
      for (size_t i=0;i!=cloud.pts.size();i++){
	ptree.reco_x = cloud.pts.at(i).x/units::cm;
	ptree.reco_y = cloud.pts.at(i).y/units::cm;
	ptree.reco_z = cloud.pts.at(i).z/units::cm;
	T->Fill();
      }
    }
    
  }
}
    
