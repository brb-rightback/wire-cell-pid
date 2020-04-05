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
#include "NeutrinoID_deghost.h"
#include "NeutrinoID_track_shower.h"

WCPPID::NeutrinoID::NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, std::vector<WCPPID::PR3DCluster*>& all_clusters, WCPPID::ToyFiducial* fid, WCPSst::GeomDataSource& gds, int nrebin, int frame_length, float unit_dis, ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, double flash_time, double offset_x)
  : acc_vertex_id(0)
  , acc_segment_id(0)
  , main_cluster(main_cluster)
  , other_clusters(other_clusters)
  , all_clusters(all_clusters)
  , fid(fid)
  , ct_point_cloud(ct_point_cloud)
  , global_wc_map(global_wc_map)
  , flash_time(flash_time)
  , offset_x(offset_x)
  , type(0)
  , main_vertex(0)
{
  bool flag_other_clusters = false;
  bool flag_main_cluster = true;
  
  // form id vs. cluster ...
  map_id_cluster[main_cluster->get_cluster_id()] = main_cluster;
  //std::cout << main_cluster->get_cluster_id() << " " << main_cluster << std::endl;
  for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
    map_id_cluster[(*it)->get_cluster_id()] = *it;
    //    std::cout << (*it)->get_cluster_id() <<  " " << *it << std::endl;
  }
  //  std::cout << map_id_cluster.size() << std::endl;

  
  if (flag_main_cluster){
    // create Steiner-inspired graph
    main_cluster->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
    // find the proto vertex ...
    find_proto_vertex(main_cluster);
    // fit the vertex in 3D 
    improve_vertex(main_cluster);
  }

  // main_cluster->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  
  // for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
  //   std::cout << it->first->get_fit_pt() << std::endl;
  // }

  if (flag_other_clusters){
    //deal with the other clusters ...
    for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
      //   if ((*it)->get_cluster_id()>38) continue;
      (*it)->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
      //std::cout << map_vertex_segments.size() << " " << map_segment_vertices.size() << std::endl;
      // do not break track and find other tracks ...
      if (!find_proto_vertex(*it, false, 1)) init_point_segment(*it);
      //      std::cout << map_vertex_segments.size() << " " << map_segment_vertices.size() << std::endl;
    }
    //  deghost ...
    deghost_clusters();
  }
  
  // clustering points
  if (flag_main_cluster)
    clustering_points(main_cluster);
  if (flag_other_clusters){
    //deal with the other clusters ...
    for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
      clustering_points(*it);
    }
  }
  
  // track shower separation
  separate_track_shower();

  if (flag_main_cluster){
    determine_direction(main_cluster);
    determine_main_vertex(main_cluster);
  }
  if (flag_other_clusters){
    for (auto it = other_clusters.begin(); it != other_clusters.end(); it++){
      determine_direction(*it);
    }
  }
  
  // for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
  //   std::cout << it->first->get_fit_pt() << std::endl;
  // }
  // prepare output ...
  fill_fit_parameters();
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
  for (auto it = vertices.begin(); it!= vertices.end(); it++){
    delete (*it);
  }
  for (auto it = particles.begin(); it!=particles.end(); it++){
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


void WCPPID::NeutrinoID::fill_reco_simple_tree(WCPPID::WCRecoTree& rtree){
  rtree.mc_daughters->clear();

  rtree.mc_Ntrack = 0;
  //start to fill
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment* sg= it->first;
    // only save the main cluster
    if (sg->get_cluster_id() != main_cluster->get_cluster_id()) continue;
    rtree.mc_id[rtree.mc_Ntrack]  = sg->get_cluster_id()*1000 + sg->get_id();

    // particle code
    // e- 11  e+ -11
    // muon- 13  muon+ -13
    // gamma 22
    // pi+ 211, pi0 111, pi- -211
    // kaon+ 321, K- -321
    // p  2212
    // n 2112
    
    rtree.mc_pdg[rtree.mc_Ntrack] = sg->get_particle_type(); // all muons for now ...
    
    rtree.mc_process[rtree.mc_Ntrack] = 0;
    rtree.mc_mother[rtree.mc_Ntrack] = 0; 

    if (sg->get_flag_dir()==0) continue;
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
  }
  rtree.mc_daughters->resize(rtree.mc_Ntrack);

  // std::cout << rtree.mc_Ntrack << std::endl;
}



void WCPPID::NeutrinoID::fill_proto_tree(WCPPID::WCRecoTree& rtree){
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
    //std::cout << segments_to_be_examined.size() << " " << used_vertices.size() << " " << used_segments.size() << std::endl;
  }
  
}
