#include "WCPPID/NeutrinoID.h"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include "WCPPaal/graph_metrics.h"
#include "WCPPaal/steiner_tree_greedy.h"

#include "WCPData/Line.h"



using namespace WCP;

#include "NeutrinoID_proto_vertex.h"
#include "NeutrinoID_improve_vertex.h"

WCPPID::NeutrinoID::NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, WCPSst::GeomDataSource& gds, int nrebin, int frame_length, float unit_dis, ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, double flash_time)
  : acc_vertex_id(0)
  , acc_segment_id(0)
  , main_cluster(main_cluster)
  , other_clusters(other_clusters)
  , ct_point_cloud(ct_point_cloud)
  , global_wc_map(global_wc_map)
  , flash_time(flash_time)
  , type(0)
{
  bool flag_other_clusters = false;
  
  // create Steiner-inspired graph
  main_cluster->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  // find the proto vertex ...
  find_proto_vertex(main_cluster);
  // fit the vertex in 3D 
  improve_vertex(main_cluster);
  
  
  // clustering points
  clustering_points(main_cluster);
  
  // prepare output ...
  organize_vertices_segments();
  main_cluster->set_fit_parameters(proto_vertices, proto_segments);

  
  if (flag_other_clusters){
    //deal with the other clusters ...
    for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
      (*it)->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
      find_proto_vertex(*it, false, 1);
      clustering_points(*it);

      organize_vertices_segments();
      (*it)->set_fit_parameters(proto_vertices, proto_segments);
    }
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
      if (point_sub_cluster_ids.at(i) == -1) continue;
      map_id_seg[point_sub_cluster_ids.at(i)]->add_associate_point(cloud.pts[i], cloud_u.pts[i], cloud_v.pts[i], cloud_w.pts[i]);
    }
  }

  {
    WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud_steiner()->get_cloud();
    std::vector<int>& point_steiner_sub_cluster_ids = temp_cluster->get_point_steiner_sub_cluster_ids();
    for (size_t i=0;i!=point_steiner_sub_cluster_ids.size();i++){
      if (point_steiner_sub_cluster_ids.at(i) == -1) continue;
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
