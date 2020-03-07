#include "WCPPID/NeutrinoID.h"

using namespace WCP;

WCPPID::NeutrinoID::NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, WCPSst::GeomDataSource& gds, int nrebin, int frame_length, float unit_dis, ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, double flash_time)
  : main_cluster(main_cluster)
  , other_clusters(other_clusters)
  , ct_point_cloud(ct_point_cloud)
  , global_wc_map(global_wc_map)
  , flash_time(flash_time)
  , type(0)
{
  // create Steiner-inspired graph
  main_cluster->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
    (*it)->create_steiner_graph(*ct_point_cloud, gds, nrebin, frame_length, unit_dis);
  }


  find_proto_vertex(main_cluster);

  
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


void WCPPID::NeutrinoID::find_proto_vertex(WCPPID::PR3DCluster *temp_cluster){
  
  if (temp_cluster->get_point_cloud_steiner()==0) return;
  if (temp_cluster->get_point_cloud_steiner()->get_num_points()<2) return;
  
  
  // do the first search of the trajectory ...
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2);
  temp_cluster->dijkstra_shortest_paths(wcps.first,2); 
  temp_cluster->cal_shortest_path(wcps.second,2);
  WCPPID::ProtoVertex *v1 = new WCPPID::ProtoVertex(wcps.first);
  WCPPID::ProtoVertex *v2 = new WCPPID::ProtoVertex(wcps.second);
  WCPPID::ProtoSegment *sg1 = new WCPPID::ProtoSegment(temp_cluster->get_path_wcps());
  add_proto_connection(v1,sg1,temp_cluster);
  add_proto_connection(v2,sg1,temp_cluster);
  
  // std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg1].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
  //  del_proto_vertex(v1); // del_proto_vertex(v1);
  //std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg1].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
   // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v1].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;
   // del_proto_segment(sg1); // del_proto_segment(sg1);
   // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v1].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;


  // temporary fit ...
  temp_cluster->collect_charge_trajectory(*ct_point_cloud);
  temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
  v1->set_fit(temp_cluster->get_fine_tracking_path().front(), temp_cluster->get_dQ().front(), temp_cluster->get_dx().front(), temp_cluster->get_pu().front(), temp_cluster->get_pv().front(), temp_cluster->get_pw().front(), temp_cluster->get_pt().front(), temp_cluster->get_reduced_chi2().front());
  v2->set_fit(temp_cluster->get_fine_tracking_path().back(), temp_cluster->get_dQ().back(), temp_cluster->get_dx().back(), temp_cluster->get_pu().back(), temp_cluster->get_pv().back(), temp_cluster->get_pw().back(), temp_cluster->get_pt().back(), temp_cluster->get_reduced_chi2().back());
  sg1->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2());
  // finish the temporary fit ...
  //sg1->print_dis();
  //std::cout << v1->get_fit_init_dis()/units::cm << " " << v2->get_fit_init_dis()/units::cm << std::endl;


  
  std::tuple<Point, TVector3, bool> kink_tuple = sg1->search_kink(temp_cluster->get_fine_tracking_path().front());

  
  if (std::get<1>(kink_tuple).Mag()!=0 ){
    // find the extreme point ... PR3DCluster function
    WCP::WCPointCloud<double>::WCPoint break_wcp = temp_cluster->proto_extend_point(std::get<0>(kink_tuple), std::get<1>(kink_tuple), std::get<2>(kink_tuple));

    //std::cout << std::get<0>(kink_tuple) << " " << std::get<1>(kink_tuple).Mag() << " " << std::get<2>(kink_tuple) << " " << break_wcp.x/units::cm << " " << break_wcp.y/units::cm << " " << break_wcp.z/units::cm << " " << break_wcp.index << std::endl;
    
    // find the new path ... PR3DCluster function
    // break into two tracks and adding a ProtoVertex in this function ...
    std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list1;
    std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list2;
    bool flag_break = temp_cluster->proto_break_tracks(v1->get_wcpt(), break_wcp, v2->get_wcpt(), wcps_list1, wcps_list2);

    //    std::cout << wcps_list1.front().index << " " << wcps_list1.back().index << " " << wcps_list2.front().index << " " << wcps_list2.back().index << std::endl;

    if (flag_break){
      WCPPID::ProtoVertex *v3 = new WCPPID::ProtoVertex(break_wcp);
      WCPPID::ProtoSegment *sg2 = new WCPPID::ProtoSegment(wcps_list1);
      WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(wcps_list2);
      
      add_proto_connection(v1, sg2, temp_cluster);
      add_proto_connection(v3, sg2, temp_cluster);
      add_proto_connection(v3, sg3, temp_cluster);
      add_proto_connection(v2, sg3, temp_cluster);
  
      del_proto_segment(sg1);
      delete sg1;
      
      // temporary fit hack ...
      temp_cluster->set_path_wcps(wcps_list1);
      temp_cluster->collect_charge_trajectory(*ct_point_cloud);
      temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
      
      sg2->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2());
      v3->set_fit(temp_cluster->get_fine_tracking_path().back(), temp_cluster->get_dQ().back(), temp_cluster->get_dx().back(), temp_cluster->get_pu().back(), temp_cluster->get_pv().back(), temp_cluster->get_pw().back(), temp_cluster->get_pt().back(), temp_cluster->get_reduced_chi2().back());
      temp_cluster->set_path_wcps(wcps_list2);
      temp_cluster->collect_charge_trajectory(*ct_point_cloud);
      temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
      sg3->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2());
    }
    
    
      
    
    //    std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg2].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
    // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v3].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;
  }
  
  temp_cluster->set_fit_parameters(map_vertex_segments, map_segment_vertices);
  
  
 

  
  // practice other components ...
}

int WCPPID::NeutrinoID::get_num_segments(WCPPID::ProtoVertex *pv){
  if (map_vertex_segments.find(pv) != map_vertex_segments.end()){
    return map_vertex_segments[pv].size();
  }else{
    return 0;
  }
}

bool WCPPID::NeutrinoID::add_proto_connection(WCPPID::ProtoVertex *pv, WCPPID::ProtoSegment *ps, WCPPID::PR3DCluster* cluster){

  if (pv->get_wcpt().index != ps->get_wcpt_vec().front().index && pv->get_wcpt().index != ps->get_wcpt_vec().back().index){
    std::cout << "Error! Vertex and Segment does not match" << std::endl;
    return false;
  }

  
  proto_vertices.insert(pv);
  map_vertex_cluster[pv] = cluster;
  if (map_cluster_vertices.find(cluster)!=map_cluster_vertices.end()){
    map_cluster_vertices[cluster].insert(pv);
  }else{
    ProtoVertexSet temp_vertices;
    temp_vertices.insert(pv);
    map_cluster_vertices[cluster] = temp_vertices;
  }

  proto_segments.insert(ps);
  map_segment_cluster[ps] = cluster;
  if (map_cluster_segments.find(cluster)!=map_cluster_segments.end()){
    map_cluster_segments[cluster].insert(ps);
  }else{
    ProtoSegmentSet temp_segments;
    temp_segments.insert(ps);
    map_cluster_segments[cluster] = temp_segments;
  }

  if (map_vertex_segments.find(pv)!=map_vertex_segments.end()){
    map_vertex_segments[pv].insert(ps);
  }else{
    ProtoSegmentSet temp_segments;
    temp_segments.insert(ps);
    map_vertex_segments[pv] = temp_segments;
  }

  if (map_segment_vertices.find(ps)!=map_segment_vertices.end()){
    map_segment_vertices[ps].insert(pv);
  }else{
    ProtoVertexSet temp_vertices;
    map_segment_vertices[ps] = temp_vertices;
    map_segment_vertices[ps].insert(pv);
  }

  return true;
}


bool WCPPID::NeutrinoID::del_proto_vertex(WCPPID::ProtoVertex *pv){
  if (proto_vertices.find(pv)!=proto_vertices.end()){
    // overall
    proto_vertices.erase(proto_vertices.find(pv));
    // map between vertex and segment
    for (auto it = map_vertex_segments[pv].begin(); it!= map_vertex_segments[pv].end(); it++){ // loop very semeng *it
      map_segment_vertices[*it].erase(map_segment_vertices[*it].find(pv));
    }
    map_vertex_segments.erase(pv);
    // map between cluster and vertices
    map_cluster_vertices[map_vertex_cluster[pv]].erase(map_cluster_vertices[map_vertex_cluster[pv]].find(pv));
    map_vertex_cluster.erase(pv);
    
    return true;
  }else{
    return false;
  }
}

bool WCPPID::NeutrinoID::del_proto_segment(WCPPID::ProtoSegment *ps){
  if (proto_segments.find(ps)!=proto_segments.end()){
    // overall
    proto_segments.erase(proto_segments.find(ps));
    //map between segment and vertex
    for (auto it = map_segment_vertices[ps].begin(); it!=map_segment_vertices[ps].end(); it++){
      map_vertex_segments[*it].erase(map_vertex_segments[*it].find(ps));
    }
    map_segment_vertices.erase(ps);
    //map between cluster and segment
    map_cluster_segments[map_segment_cluster[ps]].erase(map_cluster_segments[map_segment_cluster[ps]].find(ps));
    map_segment_cluster.erase(ps);
    return true;
  }else{
    return false;
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
