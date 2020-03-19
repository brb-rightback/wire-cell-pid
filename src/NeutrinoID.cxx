#include "WCPPID/NeutrinoID.h"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include "WCPPaal/graph_metrics.h"
#include "WCPPaal/steiner_tree_greedy.h"

#include "WCPData/Line.h"



using namespace WCP;

#include "NeutrinoID_proto_vertex.h"


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
  // create Steiner-inspired graph
  main_cluster->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  find_proto_vertex(main_cluster);
  
  // deal with the other clusters ...
  // for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
  //   (*it)->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  //   find_proto_vertex(*it, false, 1);
  // }



  
  // deal with main cluster
  //  process_main_cluster();
  
  // deal with other clusters ...
  // process_other_clusters();
  

  
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


WCPPID::ProtoSegment* WCPPID::NeutrinoID::init_first_segment(WCPPID::PR3DCluster *temp_cluster){
  // do the first search of the trajectory ...
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2);
  temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
  temp_cluster->cal_shortest_path(wcps.second,2);
  
  WCPPID::ProtoVertex *v1=0;
  WCPPID::ProtoVertex *v2=0;
  WCPPID::ProtoSegment *sg1=0;
  
  if (temp_cluster->get_path_wcps().size()>1){
    v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.first, temp_cluster->get_cluster_id()); acc_vertex_id++;
    v2 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.second, temp_cluster->get_cluster_id()); acc_vertex_id++;
    sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
    add_proto_connection(v1,sg1,temp_cluster);
    add_proto_connection(v2,sg1,temp_cluster);

    temp_cluster->collect_charge_trajectory(*ct_point_cloud);
    temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    v1->set_fit(temp_cluster->get_fine_tracking_path().front(), temp_cluster->get_dQ().front(), temp_cluster->get_dx().front(), temp_cluster->get_pu().front(), temp_cluster->get_pv().front(), temp_cluster->get_pw().front(), temp_cluster->get_pt().front(), temp_cluster->get_reduced_chi2().front());
    v2->set_fit(temp_cluster->get_fine_tracking_path().back(), temp_cluster->get_dQ().back(), temp_cluster->get_dx().back(), temp_cluster->get_pu().back(), temp_cluster->get_pv().back(), temp_cluster->get_pw().back(), temp_cluster->get_pt().back(), temp_cluster->get_reduced_chi2().back());
    sg1->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2());
  }else{
     v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcps.first, temp_cluster->get_cluster_id()); acc_vertex_id++;
     sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
     add_proto_connection(v1,sg1,temp_cluster);
  }
  
  
  // std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg1].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
  //  del_proto_vertex(v1); // del_proto_vertex(v1);
  // std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg1].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
   // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v1].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;
   // del_proto_segment(sg1); // del_proto_segment(sg1);
   // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v1].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;

  // temporary fit ...
  
  // finish the temporary fit ...

  //  std::cout << temp_cluster->get_fine_tracking_path().size() << " " << temp_cluster->get_fine_tracking_path().at(40) << std::endl;

  
  // sg1->print_dis();
  // std::cout << v1->get_fit_init_dis()/units::cm << " " << v2->get_fit_init_dis()/units::cm << std::endl;
  return sg1;
}


void WCPPID::NeutrinoID::find_proto_vertex(WCPPID::PR3DCluster *temp_cluster, bool flag_break_track, int nrounds_find_other_tracks){
  
  if (temp_cluster->get_point_cloud_steiner()==0) return;
  if (temp_cluster->get_point_cloud_steiner()->get_num_points()<2) return;
  
  WCPPID::ProtoSegment* sg1 = init_first_segment(temp_cluster);


  // for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
  //   WCPPID::ProtoSegment *sg = it->first;
  //   std::cout << sg->get_point_vec().size() << " A " << sg->get_point_vec().front() << " " << sg->get_point_vec().back() << std::endl;
  // }

  if (sg1->get_wcpt_vec().size()>1){
    // break tracks ...
    if (flag_break_track){
      std::vector<WCPPID::ProtoSegment*> remaining_segments;
      remaining_segments.push_back(sg1);
      break_segments(remaining_segments, temp_cluster);
    }
    
    // find other segments ...
    for (size_t i=0;i!=nrounds_find_other_tracks;i++){
      find_other_segments(temp_cluster, flag_break_track);
    }
  }
  
  // prepare output ...
  organize_vertices_segments();
  temp_cluster->set_fit_parameters(proto_vertices, proto_segments);
  

  // for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
  //   WCPPID::ProtoSegment *sg = it->first;
  //   std::cout << sg->get_id() << " " << sg->get_point_vec().size() << " " << sg->get_point_vec().front() << " " << sg->get_point_vec().back() << std::endl;
  //   for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  //     std::cout << (*it1)->get_id() << " " << (*it1) << std::endl;
  //   }
  // }
  // for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
  //   std::cout << (it->first)->get_id() << std::endl;
  //   for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  //     std::cout << (*it1)->get_id() << " " << (*it1) << std::endl;
  //   }
  // }
  
  //  temp_cluster->set_fit_parameters(map_vertex_segments, map_segment_vertices);
  
  
  
  
  // practice other components ...
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
