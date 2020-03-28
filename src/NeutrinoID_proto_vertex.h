
void WCPPID::NeutrinoID::find_proto_vertex(WCPPID::PR3DCluster *temp_cluster, bool flag_break_track, int nrounds_find_other_tracks){
  
  if (temp_cluster->get_point_cloud_steiner()==0) return;
  if (temp_cluster->get_point_cloud_steiner()->get_num_points()<2) return;
  
  WCPPID::ProtoSegment* sg1 = init_first_segment(temp_cluster);

  // hack test
  //temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true);
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
    // examine the vertices ...
    examine_vertices(temp_cluster);
  }
  

  

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



WCPPID::ProtoSegment* WCPPID::NeutrinoID::init_first_segment(WCPPID::PR3DCluster *temp_cluster){
  // do the first search of the trajectory ...
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = temp_cluster->get_two_boundary_wcps(2);
  // good for the first track
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

    // fit dQ/dx and everything ...
    temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true);

    //std::cout << temp_cluster->get_fine_tracking_path().size() << " " << temp_cluster->get_dx().size() << std::endl;
    
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




void WCPPID::NeutrinoID::break_segments(std::vector<WCPPID::ProtoSegment*>& remaining_segments, WCPPID::PR3DCluster* temp_cluster){
  bool flag_print = false;
  
  int count = 0;
  
  while(remaining_segments.size()!=0 && count < 2){
    WCPPID::ProtoSegment* curr_sg = remaining_segments.back();
    remaining_segments.pop_back();
    WCPPID::ProtoVertex *start_v=0, *end_v=0;
    for (auto it = map_segment_vertices[curr_sg].begin(); it!=map_segment_vertices[curr_sg].end(); it++){
      if ((*it)->get_wcpt().index == curr_sg->get_wcpt_vec().front().index) start_v = *it;
      if ((*it)->get_wcpt().index == curr_sg->get_wcpt_vec().back().index) end_v = *it;
    }
    if (start_v==0 || end_v==0){
      std::cout << "Error in finding vertices for a segment" << std::endl; 
      //std::cout << "Vertex: " << start_v->get_wcpt().index << " " << end_v->get_wcpt().index << std::endl;
    }
    
    // initialize the start points
    WCP::WCPointCloud<double>::WCPoint break_wcp = start_v->get_wcpt();
    Point test_start_p = curr_sg->get_point_vec().front();
    
    while(sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2) +
	       pow(start_v->get_wcpt().y - break_wcp.y,2) +
	       pow(start_v->get_wcpt().z - break_wcp.z,2)) <= 1*units::cm &&
	  sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2) +
	       pow(end_v->get_wcpt().y - break_wcp.y,2) +
	       pow(end_v->get_wcpt().z - break_wcp.z,2)) > 1*units::cm){
      std::tuple<Point, TVector3, TVector3, bool> kink_tuple = curr_sg->search_kink(test_start_p);
      if (std::get<1>(kink_tuple).Mag()!=0 ){
	// find the extreme point ... PR3DCluster function
	break_wcp = temp_cluster->proto_extend_point(std::get<0>(kink_tuple), std::get<1>(kink_tuple), std::get<2>(kink_tuple), std::get<3>(kink_tuple));
	//	std::cout << "Break point: " << remaining_segments.size() << " " << std::get<0>(kink_tuple) << " " << std::get<1>(kink_tuple).Mag() << " " << std::get<3>(kink_tuple) << " " << break_wcp.x/units::cm << " " << break_wcp.y/units::cm << " " << break_wcp.z/units::cm << " " << break_wcp.index << std::endl;
        if (sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2) +
	       pow(start_v->get_wcpt().y - break_wcp.y,2) +
	       pow(start_v->get_wcpt().z - break_wcp.z,2)) <= 1*units::cm &&
	    sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2) +
		 pow(end_v->get_wcpt().y - break_wcp.y,2) +
		 pow(end_v->get_wcpt().z - break_wcp.z,2)) > 1*units::cm)
	  test_start_p = std::get<0>(kink_tuple);
      }else{
	break;
      }
    }

    if (sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2) +
	     pow(start_v->get_wcpt().y - break_wcp.y,2) +
	     pow(start_v->get_wcpt().z - break_wcp.z,2)) > 1*units::cm){
      // find the new path ... PR3DCluster function
      // break into two tracks and adding a ProtoVertex in this function ...
      std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list1;
      std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list2;
      bool flag_break = temp_cluster->proto_break_tracks(start_v->get_wcpt(), break_wcp, end_v->get_wcpt(), wcps_list1, wcps_list2);

      std::cout << "breaking point: " << flag_break << " " << wcps_list1.front().index << " " << wcps_list1.back().index << " " << wcps_list2.front().index << " " << wcps_list2.back().index << std::endl;
	  
      if (flag_break){
	WCPPID::ProtoVertex *v3 = new WCPPID::ProtoVertex(acc_vertex_id, break_wcp, temp_cluster->get_cluster_id()); acc_vertex_id ++;
	WCPPID::ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list1, temp_cluster->get_cluster_id()); acc_segment_id++;
	WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list2, temp_cluster->get_cluster_id()); acc_segment_id++;
	
	add_proto_connection(start_v, sg2, temp_cluster);
	add_proto_connection(v3, sg2, temp_cluster);
	add_proto_connection(v3, sg3, temp_cluster);
	add_proto_connection(end_v, sg3, temp_cluster);
	
	del_proto_segment(curr_sg);
	//delete curr_sg;

	flag_print = true;
	// fit dQ/dx here, do not exclude others
	temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, false);
	
	// temporary fit hack ...
	/* temp_cluster->set_path_wcps(wcps_list1); */
	/* temp_cluster->collect_charge_trajectory(*ct_point_cloud); */
	/* temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false); */
	/* sg2->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2()); */
	/* v3->set_fit(temp_cluster->get_fine_tracking_path().back(), temp_cluster->get_dQ().back(), temp_cluster->get_dx().back(), temp_cluster->get_pu().back(), temp_cluster->get_pv().back(), temp_cluster->get_pw().back(), temp_cluster->get_pt().back(), temp_cluster->get_reduced_chi2().back()); */
	/* temp_cluster->set_path_wcps(wcps_list2); */
	/* temp_cluster->collect_charge_trajectory(*ct_point_cloud); */
	/* temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false); */
	/* sg3->set_fit_vec(temp_cluster->get_fine_tracking_path(), temp_cluster->get_dQ(), temp_cluster->get_dx(), temp_cluster->get_pu(), temp_cluster->get_pv(), temp_cluster->get_pw(), temp_cluster->get_pt(), temp_cluster->get_reduced_chi2()); */
	//	count ++;
	remaining_segments.push_back(sg3); 
      }
    }
    
    //    std::cout << proto_vertices.size() << " " << map_vertex_segments.size() << " " << map_segment_vertices[sg2].size() << " " << map_cluster_vertices[main_cluster].size() << " " << map_vertex_cluster.size()  << std::endl;
    // std::cout << proto_segments.size() << " " << map_segment_vertices.size() << " " << map_vertex_segments[v3].size() << " " << map_cluster_segments[main_cluster].size() << " " << map_segment_cluster.size() << std::endl;
  }

  if (flag_print)
    std::cout << "Break tracks  -- # of Vertices: " << map_vertex_segments.size() << "; # of Segments: " << map_segment_vertices.size() << std::endl;
}




void WCPPID::NeutrinoID::find_other_segments(WCPPID::PR3DCluster* temp_cluster, bool flag_break_track, double search_range, double scaling_2d){
  

  const int N = temp_cluster->get_point_cloud_steiner()->get_num_points();
  std::vector<bool>& flag_tagged = temp_cluster->get_flag_tagged_steiner_graph();
  flag_tagged.resize(N, false);
  int num_tagged = 0;
  WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud_steiner()->get_cloud();
  for (size_t i=0;i!=N;i++){
    Point p(cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z);
    double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
    double min_3d_dis = 1e9;
    for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg1 = it->first;
      if (sg1->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
      std::pair<double, WCP::Point> closest_dis_point = sg1->get_closest_point(p);
      std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p);
      if (closest_dis_point.first < min_3d_dis) min_3d_dis = closest_dis_point.first;
      if (closest_dis_point.first < search_range){
	flag_tagged[i] = true;
	num_tagged ++;
	break;
      }
      if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis);
      if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis);
      if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis);
    }

    if (!flag_tagged[i]){
      if ((min_dis_u < scaling_2d * search_range || ct_point_cloud->get_closest_dead_chs(p, 0) ) &&
	  (min_dis_v < scaling_2d * search_range || ct_point_cloud->get_closest_dead_chs(p, 1) ) &&
	  (min_dis_w < scaling_2d * search_range || ct_point_cloud->get_closest_dead_chs(p, 2) ) )
	flag_tagged[i] = true;
	  // std::cout << min_dis_u/units::cm << " " << min_dis_v/units::cm << " " << min_dis_w/units::cm << std::endl;
    }
    
    //    if (cloud.pts[i].index_v>4440-2400 && cloud.pts[i].index_v < 4460-2400 && cloud.pts[i].index_w > 7895-4800 && cloud.pts[i].index_w < 7910-4800 && p.x < 980){
    //std::cout << p << " " << cloud.pts[i].index << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v+2400 << " " << cloud.pts[i].index_w+4800 << " " << min_dis_u/units::cm << " " << min_dis_v/units::cm << " " << min_dis_w/units::cm << " " <<  min_3d_dis/units::cm << " " << flag_tagged[i] << std::endl;
    //}
    
    
  }

  // Now figure out the Terminal Graph ...
  std::vector<int> terminals;
  std::map<int, int> map_oindex_tindex;
  std::vector<bool>& flag_steiner_terminal = temp_cluster->get_flag_steiner_terminal();
  for (size_t i = 0;i!=flag_steiner_terminal.size();i++){
    if (flag_steiner_terminal[i]){
      map_oindex_tindex[i] = terminals.size();
      terminals.push_back(i);
    }
  }

  using Vertex = typename boost::graph_traits<WCPPID::MCUGraph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<WCPPID::MCUGraph>::edge_descriptor;
  using Base = typename boost::property<edge_base_t, Edge>;
  using EdgeWeightMap = typename boost::property_map<WCPPID::MCUGraph, boost::edge_weight_t>::type;
  using Weight = typename boost::property_traits<EdgeWeightMap>::value_type;
  using WeightProperty =
    typename boost::property<boost::edge_weight_t, Weight, Base>;
  using TerminalGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, WeightProperty>;
  using EdgeTerminal =
    typename boost::graph_traits<TerminalGraph>::edge_descriptor;
  
  WCPPID::MCUGraph* graph_steiner = temp_cluster->get_graph_steiner();
  EdgeWeightMap edge_weight = boost::choose_pmap(boost::get_param(boost::no_named_parameters(), boost::edge_weight), *graph_steiner, boost::edge_weight);
  
  // distance array used in the dijkstra runs
  std::vector<Weight> distance(N);

  std::vector<Vertex> nearest_terminal(num_vertices(*graph_steiner));
  auto index = get(boost::vertex_index, *graph_steiner);
  auto nearest_terminal_map = boost::make_iterator_property_map(nearest_terminal.begin(), get(boost::vertex_index, *graph_steiner));
  for (auto terminal : terminals) {
    nearest_terminal_map[terminal] = terminal;
  }
  // compute voronoi diagram each vertex get nearest terminal and last edge on
  // path to nearest terminal
  auto distance_map = make_iterator_property_map(distance.begin(), index);
  std::vector<Edge> vpred(N);
  auto last_edge = boost::make_iterator_property_map(vpred.begin(), get(boost::vertex_index, *graph_steiner));
  boost::dijkstra_shortest_paths(*graph_steiner, terminals.begin(), terminals.end(), boost::dummy_property_map(),
				 distance_map, edge_weight, index, paal::utils::less(),
				 boost::closed_plus<Weight>(), std::numeric_limits<Weight>::max(), 0,
				 boost::make_dijkstra_visitor(paal::detail::make_nearest_recorder(
												  nearest_terminal_map, last_edge, boost::on_edge_relaxed{})));
  
  // computing distances between terminals
  // creating terminal_graph
  TerminalGraph terminal_graph(N);
  std::map<std::pair<int, int>, std::pair<Weight, Edge> > map_saved_edge;

  for (auto w : boost::as_array(edges(*graph_steiner))) {
    auto const &nearest_to_source = nearest_terminal_map[source(w, *graph_steiner)];
    auto const &nearest_to_target = nearest_terminal_map[target(w, *graph_steiner)];
    if (nearest_to_source != nearest_to_target) {
      Weight temp_weight = distance[source(w, *graph_steiner)] + distance[target(w, *graph_steiner)] + edge_weight[w];
      if (map_saved_edge.find(std::make_pair(nearest_to_source, nearest_to_target))!=map_saved_edge.end()){
	if (temp_weight < map_saved_edge[std::make_pair(nearest_to_source, nearest_to_target)].first)
	  map_saved_edge[std::make_pair(nearest_to_source, nearest_to_target)] = std::make_pair(temp_weight,w);
      }else if (map_saved_edge.find(std::make_pair(nearest_to_target, nearest_to_source))!=map_saved_edge.end()){
	if (temp_weight < map_saved_edge[std::make_pair(nearest_to_target, nearest_to_source)].first)
	  map_saved_edge[std::make_pair(nearest_to_target, nearest_to_source)] = std::make_pair(temp_weight,w);
      }else{
	map_saved_edge[std::make_pair(nearest_to_source, nearest_to_target)] = std::make_pair(temp_weight,w);
      }
    }
  }

  std::vector<Edge> terminal_edge;
  for (auto it = map_saved_edge.begin(); it!=map_saved_edge.end(); it++){
    std::pair<Edge, bool> p = add_edge(it->first.first, it->first.second,
				       WeightProperty(it->second.first, Base(it->second.second)),terminal_graph);
    //      terminal_edge.push_back(p.first);
  }
  // minimal spanning tree ...
  boost::kruskal_minimum_spanning_tree(terminal_graph,
				       std::back_inserter(terminal_edge));
  
  // new Graph ...
  std::map<int, std::set<int> > map_connection;
  TerminalGraph terminal_graph_cluster(terminals.size());
  for (auto it = terminal_edge.begin(); it!=terminal_edge.end(); it++){
    if (flag_tagged[source(*it, terminal_graph)] == flag_tagged[target(*it, terminal_graph)]){
      add_edge(map_oindex_tindex[source(*it, terminal_graph)],map_oindex_tindex[target(*it, terminal_graph)], terminal_graph_cluster);
    }else{
      if (map_connection.find(source(*it, terminal_graph))==map_connection.end()){
	std::set<int> temp_results;
	temp_results.insert(target(*it,terminal_graph));
	map_connection[source(*it, terminal_graph)] = temp_results;
      }else{
	map_connection[source(*it, terminal_graph)].insert(target(*it,terminal_graph));
      }
      if (map_connection.find(target(*it, terminal_graph))==map_connection.end()){
	std::set<int> temp_results;
	temp_results.insert(source(*it,terminal_graph));
	map_connection[target(*it, terminal_graph)] = temp_results;
      }else{
	map_connection[target(*it, terminal_graph)].insert(source(*it,terminal_graph));
      }
    }
    //      std::cout << source(*it, terminal_graph) << " " << target(*it, terminal_graph) << std::endl;
  }
  //std::cout << map_connection.size() << std::endl;
  //    std::cout << map_saved_edge.size() << " " << terminal_edge.size() << std::endl;

  std::vector<int> component(num_vertices(terminal_graph_cluster));
  const int num = connected_components(terminal_graph_cluster,&component[0]);
  //  std::cout << "jaja: " << num << std::endl;
  int ncounts[num];
  std::vector<std::vector<int> > sep_clusters(num);
  for (size_t i=0;i!=num;i++){
    ncounts[i] = 0;
    std::vector<int> temp;
    sep_clusters[i] = temp;
  }  
  {
    std::vector<int>::size_type i;
    for (i=0;i!=component.size(); ++i){
      ncounts[component[i]]++;
      sep_clusters[component[i]].push_back(terminals.at(i));
    }
  }

  //test
  //  for (size_t j=0;j!=sep_clusters.size();j++){
  // for(size_t i=0;i!=sep_clusters[j].size();i++){
  //  if (sep_clusters[j].at(i) > 530 && sep_clusters[j].at(i) < 570)
  //	std::cout << j << " " << i << " " << sep_clusters[j].at(i) << " " << flag_tagged[sep_clusters[j].at(i)] << std::endl;
  //}
  //}

  


  std::vector<WCPPID::Res_proto_segment> temp_segments(num);
  std::set<int> remaining_segments;
  for (int i=0;i!=num;i++){
    remaining_segments.insert(i);
  }

  for (size_t i=0;i!=num;i++){
    // inside the original track 
    if (flag_tagged[sep_clusters[i].front()] ) {
      remaining_segments.erase(i);
      continue;
    }
    temp_segments.at(i).group_num = i;
    temp_segments.at(i).number_points = ncounts[i];
    
    int special_A = -1;
    for (size_t j=0;j!=ncounts[i];j++){
      if (map_connection.find(sep_clusters[i].at(j))!=map_connection.end()){
	special_A = sep_clusters[i].at(j);
	break;
      }
    }

    //    if (special_A !=-1){
    //  std::cout << i << " " << special_A << " " << sep_clusters[i].size() << std::endl;
    //  for (auto it1 = map_connection[special_A].begin(); it1 != map_connection[special_A].end(); it1++){
    // 	std::cout << *it1 << std::endl;
    //  }
    // }
    
    int special_B = special_A;
    double min_dis = 0;
    
    int number_not_faked = 0;
    double max_dis_u = 0, max_dis_v = 0, max_dis_w = 0;
    for (size_t j=0;j!=ncounts[i];j++){
      double dis = sqrt(pow(cloud.pts[sep_clusters[i].at(j)].x - cloud.pts[special_A].x,2)
			+ pow(cloud.pts[sep_clusters[i].at(j)].y - cloud.pts[special_A].y,2)
			+ pow(cloud.pts[sep_clusters[i].at(j)].z - cloud.pts[special_A].z,2));
      if (dis > min_dis){
	min_dis = dis;
	special_B = sep_clusters[i].at(j); // furthest away from special_A
      }

      // also judge whether this track is fake ...
      WCP::Point p(cloud.pts[sep_clusters[i].at(j)].x, cloud.pts[sep_clusters[i].at(j)].y, cloud.pts[sep_clusters[i].at(j)].z);
      double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
      for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
	WCPPID::ProtoSegment *sg1 = it->first;
	if (sg1->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
	std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p);
	if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis);
	if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis);
	if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis);
      }
      int flag_num = 0;
      if (min_dis_u > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 0))) flag_num ++;
      if (min_dis_v > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 1))) flag_num ++;
      if (min_dis_w > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 2))) flag_num ++;
      if (min_dis_u > max_dis_u && (!ct_point_cloud->get_closest_dead_chs(p, 0))) max_dis_u = min_dis_u;
      if (min_dis_v > max_dis_v && (!ct_point_cloud->get_closest_dead_chs(p, 1))) max_dis_v = min_dis_v;
      if (min_dis_w > max_dis_w && (!ct_point_cloud->get_closest_dead_chs(p, 2))) max_dis_w = min_dis_w;
      if (flag_num>=2) number_not_faked ++;
      
    }

    
    double length = sqrt(pow(cloud.pts[special_A].x - cloud.pts[special_B].x,2) + pow(cloud.pts[special_A].y - cloud.pts[special_B].y,2) + pow(cloud.pts[special_A].z - cloud.pts[special_B].z,2)); 
    // check ...
    
    if ( length < 3*units::cm){
      int save_index = special_A;
      double save_dis = 1e9;
      for (auto it1 = map_connection[special_A].begin(); it1!=map_connection[special_A].end();it1++){
    	double temp_dis = sqrt(pow(cloud.pts[special_A].x - cloud.pts[*it1].x,2) + pow(cloud.pts[*it1].y - cloud.pts[*it1].y,2) + pow(cloud.pts[*it1].z - cloud.pts[*it1].z,2));
    	if (temp_dis < save_dis){
    	  save_index = *it1;
    	  save_dis = temp_dis;
    	}
      }
      special_A = save_index;
      length = sqrt(pow(cloud.pts[special_A].x - cloud.pts[special_B].x,2) + pow(cloud.pts[special_A].y - cloud.pts[special_B].y,2) + pow(cloud.pts[special_A].z - cloud.pts[special_B].z,2));
    }
    //    std::cout << special_A << " " << special_B << " " << length/units::cm << std::endl;
    
    temp_segments.at(i).special_A = special_A;
    temp_segments.at(i).special_B = special_B;
    temp_segments.at(i).length = length;
    temp_segments.at(i).number_not_faked = number_not_faked;
    temp_segments.at(i).max_dis_u = max_dis_u;
    temp_segments.at(i).max_dis_v = max_dis_v;
    temp_segments.at(i).max_dis_w = max_dis_w;

    //    std::cout << special_A << " " << special_B << " " << length/units::cm << std::endl;
    
    if (temp_segments.at(i).number_points ==1  //  only one point 
    	|| number_not_faked == 0 &&
	(length < 3.5*units::cm  // very short & fake
	 || (number_not_faked < 0.25 * temp_segments.at(i).number_points || number_not_faked < 0.4 * temp_segments.at(i).number_points && length < 7 * units::cm) && max_dis_u/units::cm < 3 && max_dis_v/units::cm < 3 && max_dis_w/units::cm < 3 && max_dis_u + max_dis_v + max_dis_w < 6*units::cm)  // many fake things and very close to each other ...
     	)
      remaining_segments.erase(i);
  }

  /* for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){  */
  /*   std::cout << "jaja1: " << *it << " " << temp_segments.at(*it).number_not_faked << " " << temp_segments.at(*it).number_points << " " <<  temp_segments.at(*it).max_dis_u/units::cm << " " << temp_segments.at(*it).max_dis_v/units::cm << " " << temp_segments.at(*it).max_dis_w/units::cm << " " << temp_segments.at(*it).length/units::cm << std::endl;  */
  /* }  */
  

  // plan to examine things ... 
  std::vector<int> saved_clusters;
  std::vector<std::pair<int, int> > saved_cluster_points(num);

  /* saved_clusters.push_back(8); */
  /* saved_cluster_points.at(8) = std::make_pair(temp_segments.at(8).special_A, temp_segments.at(8).special_B); */
  
  ProtoSegmentSelection temp_segments_1, temp_segments_2;
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    if (it->first->get_cluster_id() != temp_cluster->get_cluster_id()) continue;
    temp_segments_1.push_back(it->first);
  }

  while (remaining_segments.size()>0){
    //    std::cout << num << " " << remaining_segments.size() << std::endl;
    // find the maximal length cluster
    double max_length = 0; int max_length_cluster = -1;
    double max_number_not_faked = 0;
    for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){
      if (temp_segments.at(*it).number_not_faked > max_number_not_faked){
  	max_length_cluster = *it;
  	max_number_not_faked = temp_segments.at(*it).number_not_faked;
  	max_length = temp_segments.at(*it).length;
      }else if (temp_segments.at(*it).number_not_faked == max_number_not_faked){
  	if (temp_segments.at(*it).length > max_length){
  	  max_length_cluster = *it;
  	  max_number_not_faked = temp_segments.at(*it).number_not_faked;
  	  max_length = temp_segments.at(*it).length;
  	}
      }
    }
    // save things ...
    saved_clusters.push_back(max_length_cluster);
    saved_cluster_points.at(max_length_cluster) = std::make_pair(temp_segments.at(max_length_cluster).special_A, temp_segments.at(max_length_cluster).special_B);
    remaining_segments.erase(max_length_cluster);

    // form a new segment with this cluster ...
    // good for a new segement ...
    temp_cluster->dijkstra_shortest_paths(cloud.pts[temp_segments.at(max_length_cluster).special_A],2);
    temp_cluster->cal_shortest_path(cloud.pts[temp_segments.at(max_length_cluster).special_B],2);
    WCPPID::ProtoSegment *tmp_sg = new WCPPID::ProtoSegment(-1, temp_cluster->get_path_wcps() ,temp_cluster->get_cluster_id());
    temp_segments_1.push_back(tmp_sg);
    temp_segments_2.push_back(tmp_sg);

    // examine the rest of points ...
    for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){
      //reset
      temp_segments.at(*it).number_not_faked = 0;
      temp_segments.at(*it).max_dis_u = 0;
      temp_segments.at(*it).max_dis_v = 0;
      temp_segments.at(*it).max_dis_w = 0;

      for (size_t j=0;j!=ncounts[*it];j++){
  	// also judge whether this track is fake ...
  	WCP::Point p(cloud.pts[sep_clusters[*it].at(j)].x, cloud.pts[sep_clusters[*it].at(j)].y, cloud.pts[sep_clusters[*it].at(j)].z);
  	double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
  	for (auto it1 = temp_segments_1.begin(); it1 != temp_segments_1.end(); it1++){
  	  WCPPID::ProtoSegment *sg1 = *it1;
  	  std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p);
  	  if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis);
  	  if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis);
  	  if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis);
  	}
  	int flag_num = 0;
  	if (min_dis_u > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 0))) flag_num ++;
  	if (min_dis_v > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 1))) flag_num ++;
  	if (min_dis_w > scaling_2d * search_range   && (!ct_point_cloud->get_closest_dead_chs(p, 2))) flag_num ++;
  	if (flag_num>=2) temp_segments.at(*it).number_not_faked ++;
  	if (min_dis_u > temp_segments.at(*it).max_dis_u && (!ct_point_cloud->get_closest_dead_chs(p, 0))) temp_segments.at(*it).max_dis_u = min_dis_u;
  	if (min_dis_v > temp_segments.at(*it).max_dis_v && (!ct_point_cloud->get_closest_dead_chs(p, 1))) temp_segments.at(*it).max_dis_v = min_dis_v;
  	if (min_dis_w > temp_segments.at(*it).max_dis_w && (!ct_point_cloud->get_closest_dead_chs(p, 2))) temp_segments.at(*it).max_dis_w = min_dis_w;
      }

       if (temp_segments.at(*it).number_points ==1  //  only one point
       	|| temp_segments.at(*it).number_not_faked ==0 &&
       	   (temp_segments.at(*it).length < 3.5*units::cm  // very short & fake
       	    || (temp_segments.at(*it).number_not_faked < 0.25 * temp_segments.at(*it).number_points || temp_segments.at(*it).number_not_faked < 0.4 * temp_segments.at(*it).number_points && temp_segments.at(*it).length < 7 * units::cm) && temp_segments.at(*it).max_dis_u/units::cm < 3 && temp_segments.at(*it).max_dis_v/units::cm < 3 && temp_segments.at(*it).max_dis_w/units::cm < 3 && temp_segments.at(*it).max_dis_u + temp_segments.at(*it).max_dis_v + temp_segments.at(*it).max_dis_w < 6*units::cm)  // many fake things and very close to each other ...
       	)
       	 remaining_segments.erase(*it);
      
    }

    /* for (auto it = remaining_segments.begin(); it!=remaining_segments.end(); it++){ */
    /*   std::cout << "jaja1: " << *it << " " << temp_segments.at(*it).number_not_faked << " " << temp_segments.at(*it).number_points << " " <<  temp_segments.at(*it).max_dis_u/units::cm << " " << temp_segments.at(*it).max_dis_v/units::cm << " " << temp_segments.at(*it).max_dis_w/units::cm << " " << temp_segments.at(*it).length/units::cm << std::endl; */
    /* } */
    /* std::cout << std::endl; */
    //    break;
  }
  for (auto it = temp_segments_2.begin(); it != temp_segments_2.end(); it++){
    delete *it;
  }

  bool flag_final_fit = true;
  
  std::vector<WCPPID::ProtoSegment*> new_segments;
  // save segments ...
  for (auto it = saved_clusters.begin(); it!= saved_clusters.end(); it++){
    // within the new segment, OK ...
    temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).first],2); 
    temp_cluster->cal_shortest_path(cloud.pts[(saved_cluster_points.at(*it)).second],2);
    std::list<WCP::WCPointCloud<double>::WCPoint> temp_path_list = temp_cluster->get_path_wcps();
    temp_cluster->collect_charge_trajectory(*ct_point_cloud);
    // do not fit dQ/dx at beginning ...
    temp_cluster->do_tracking(*ct_point_cloud, global_wc_map, flash_time*units::microsecond, false, false);

    if (temp_cluster->get_fine_tracking_path().size() >1){      
      WCPPID::ProtoVertex *v1 = find_vertex_other_segment(temp_cluster, true, cloud.pts[(saved_cluster_points.at(*it)).first]);
      WCPPID::ProtoVertex *v2 = find_vertex_other_segment(temp_cluster, false, cloud.pts[(saved_cluster_points.at(*it)).second]);   
      if (v1->get_wcpt().index == cloud.pts[(saved_cluster_points.at(*it)).first].index && v2->get_wcpt().index == cloud.pts[(saved_cluster_points.at(*it)).second].index){
	// this is between v1 and v2 ...
      }else{
	
	//std::cout << temp_path_list.size() << std::endl;
	if (v1->get_wcpt().index != cloud.pts[(saved_cluster_points.at(*it)).first].index){
	  double tmp_length = sqrt(pow(v1->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).first].x,2) + pow(v1->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).first].y,2) + pow(v1->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).first].z,2));
	  int tmp_count = std::round(tmp_length/(0.6*units::cm));
	  for (int k=1;k<tmp_count;k++){
	    Point tmp_p;
	    tmp_p.x = cloud.pts[(saved_cluster_points.at(*it)).first].x + (v1->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).first].x)*1./tmp_count*k;
	    tmp_p.y = cloud.pts[(saved_cluster_points.at(*it)).first].y + (v1->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).first].y)*1./tmp_count*k;
	    tmp_p.z = cloud.pts[(saved_cluster_points.at(*it)).first].z + (v1->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).first].z)*1./tmp_count*k;
	    
	    WCP::WCPointCloud<double>::WCPoint& tmp_wcp = temp_cluster->get_point_cloud_steiner()->get_closest_wcpoint(tmp_p);
	    if (tmp_wcp.index != temp_path_list.front().index) temp_path_list.push_front(tmp_wcp);
	  }
	  if (v1->get_wcpt().index != temp_path_list.front().index) temp_path_list.push_front(v1->get_wcpt());

	  /*
	  // the starting point is already the first point
	  temp_cluster->cal_shortest_path(v1->get_wcpt(),2);
	  if (temp_cluster->get_path_wcps().size() >0 && temp_path_list.size() > 0){
	    WCP::WCPointCloud<double>::WCPoint curr_wcp;
	    bool flag_replace = false;
	    while(temp_cluster->get_path_wcps().front().index == temp_path_list.front().index){
	      curr_wcp =  temp_path_list.front();
	      flag_replace = true;
	      temp_cluster->get_path_wcps().pop_front();
	      temp_path_list.pop_front();
	      if (temp_cluster->get_path_wcps().size() ==0 || temp_path_list.size() == 0) break;
	    }
	    if (flag_replace) temp_path_list.push_front(curr_wcp);
	  }
	  for (auto it1 = temp_cluster->get_path_wcps().begin(); it1 != temp_cluster->get_path_wcps().end(); it1++){
	    temp_path_list.push_front(*it1);
	  }
	  // above needs to be replaced ...
	  */
	  
	}

	//for (auto it1 = temp_path_list.begin(); it1!=temp_path_list.end(); it1++){
	//   std::cout <<  " B: " << (*it1).index << std::endl;
	// }
	
	//std::cout << temp_path_list.size() << std::endl;
	
	if (v2->get_wcpt().index != cloud.pts[(saved_cluster_points.at(*it)).second].index){
	  double tmp_length = sqrt(pow(v2->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).second].x,2) + pow(v2->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).second].y,2) + pow(v2->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).second].z,2));
	  int tmp_count = std::round(tmp_length/(0.6*units::cm));
	  for (int k=1;k<tmp_count;k++){
	    Point tmp_p;
	    tmp_p.x = cloud.pts[(saved_cluster_points.at(*it)).second].x + (v2->get_wcpt().x - cloud.pts[(saved_cluster_points.at(*it)).second].x)*1./tmp_count*k;
	    tmp_p.y = cloud.pts[(saved_cluster_points.at(*it)).second].y + (v2->get_wcpt().y - cloud.pts[(saved_cluster_points.at(*it)).second].y)*1./tmp_count*k;
	    tmp_p.z = cloud.pts[(saved_cluster_points.at(*it)).second].z + (v2->get_wcpt().z - cloud.pts[(saved_cluster_points.at(*it)).second].z)*1./tmp_count*k;
	    
	    WCP::WCPointCloud<double>::WCPoint& tmp_wcp = temp_cluster->get_point_cloud_steiner()->get_closest_wcpoint(tmp_p);
	    if (tmp_wcp.index != temp_path_list.back().index) temp_path_list.push_back(tmp_wcp);
	  }
	  if (v2->get_wcpt().index != temp_path_list.back().index) temp_path_list.push_back(v2->get_wcpt());
	  
	  /*
	  // below needs to be replaced ...
	  temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).second],2); 
	  temp_cluster->cal_shortest_path(v2->get_wcpt(),2);

	  if (temp_cluster->get_path_wcps().size() >0 && temp_path_list.size() > 0){
	    WCP::WCPointCloud<double>::WCPoint curr_wcp;
	    bool flag_replace = false;
	    while(temp_cluster->get_path_wcps().front().index == temp_path_list.back().index){
	      curr_wcp = temp_path_list.back();
	      flag_replace = true;
	      temp_cluster->get_path_wcps().pop_front();
	      temp_path_list.pop_back();
	    }
	    if (flag_replace) temp_path_list.push_back(curr_wcp);
	  }
	  for (auto it1 = temp_cluster->get_path_wcps().begin(); it1 != temp_cluster->get_path_wcps().end(); it1++){
	    temp_path_list.push_back(*it1);
	  }
	  */
	  
	}

	/* std::cout << v1->get_wcpt().index << " " << temp_path_list.front().index << " " << temp_path_list.size() << " " << temp_path_list.back().index << " " << v2->get_wcpt().index << " " << temp_cluster->get_path_wcps().size() << " "<< temp_cluster->get_path_wcps().back().index << " " << temp_cluster->get_cluster_id() << " " << v2->get_cluster_id() << std::endl; */
	/* std::cout << cloud.pts[(saved_cluster_points.at(*it)).second].index << " " << cloud.pts[(saved_cluster_points.at(*it)).first].index << std::endl; */
	//temp_cluster->dijkstra_shortest_paths(v1->get_wcpt(),2);
	/* temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).first],2); */
	/* temp_cluster->cal_shortest_path(v2->get_wcpt(),2); */
	/* std::cout << " A: " << temp_cluster->get_path_wcps().size() << std::endl; */
	/* for (auto it1 = temp_cluster->get_path_wcps().begin(); it1!= temp_cluster->get_path_wcps().end();it1++){ */
	/*   std::cout << it1->index << std::endl; */
	/* } */
	/* temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).first],2); */
	/* temp_cluster->cal_shortest_path(cloud.pts[(saved_cluster_points.at(*it)).second],2); */
	/* std::cout << " A: " << temp_cluster->get_path_wcps().size() << std::endl; */
	/* for (auto it1 = temp_cluster->get_path_wcps().begin(); it1!= temp_cluster->get_path_wcps().end();it1++){ */
	/*   std::cout << it1->index << std::endl; */
	/* } */
	
	/* temp_cluster->dijkstra_shortest_paths(cloud.pts[(saved_cluster_points.at(*it)).second],2); */
	/* temp_cluster->cal_shortest_path(v2->get_wcpt(),2); */
	/* std::cout << " A: " << temp_cluster->get_path_wcps().size() << std::endl; */
	/* for (auto it1 = temp_cluster->get_path_wcps().begin(); it1!= temp_cluster->get_path_wcps().end();it1++){ */
	/*   std::cout << it1->index << std::endl; */
	/* } */
	
	//	std::cout << temp_path_list.size() << std::endl;
	temp_cluster->set_path_wcps(temp_path_list);
	
	
      }

      if (map_vertex_segments.find(v1)!=map_vertex_segments.end() || map_vertex_segments.find(v2) != map_vertex_segments.end()){
	WCPPID::ProtoSegment *sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
	add_proto_connection(v1,sg1,temp_cluster);
	add_proto_connection(v2,sg1,temp_cluster);
	if (sg1->get_length()/units::cm > 12*units::cm)	new_segments.push_back(sg1);

	// do dQ/dx fitting and exclude others ...
	temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
	flag_final_fit = false;
	// update output ...
	std::cout << "Other tracks -- # of Vertices: " << map_vertex_segments.size() << "; # of Segments: " << map_segment_vertices.size() << std::endl;
      }else{
	// hack
	/* WCPPID::ProtoSegment *sg1 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++; */
	/* add_proto_connection(v1,sg1,temp_cluster); */
	/* add_proto_connection(v2,sg1,temp_cluster); */
	/* temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true); */
	
	residual_segment_candidates.push_back(std::make_tuple(temp_cluster, v1->get_wcpt().index, v2->get_wcpt().index));
	std::cout << "Isolated residual segment found: " << sqrt(pow(v1->get_fit_pt().x - v2->get_fit_pt().x,2) + pow(v1->get_fit_pt().y - v2->get_fit_pt().y,2) + pow(v1->get_fit_pt().z - v2->get_fit_pt().z,2))/units::cm << " cm " << v1->get_fit_pt() << " " << v2->get_fit_pt() << std::endl;
	delete v1;
	delete v2;
      }

      //hack ...
      //break;
    }
  }
  
  if (flag_break_track)
    break_segments(new_segments, temp_cluster);

  //  if (flag_final_fit)
  // temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
}

WCPPID::ProtoVertex* WCPPID::NeutrinoID::find_vertex_other_segment(WCPPID::PR3DCluster* temp_cluster, bool flag_forward, WCP::WCPointCloud<double>::WCPoint& wcp){
  WCPPID::ProtoVertex *v1 = 0;
  
  std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, Point> check_results = check_end_point(temp_cluster->get_fine_tracking_path(), flag_forward);
  
 // hack ...
 // std::get<0>(check_results) = 0;
 //std::get<1>(check_results) = 0;

 
 if (std::get<0>(check_results) != 0){
   v1 = std::get<0>(check_results);
 }else if (std::get<1>(check_results) != 0){
   
   WCPPID::ProtoSegment *sg1 = std::get<1>(check_results);
   Point test_p = std::get<2>(check_results);
   WCP::WCPointCloud<double>::WCPoint break_wcp = sg1->get_closest_wcpt(test_p);
   WCPPID::ProtoVertex *start_v = 0, *end_v = 0;
   for (auto it2 = map_segment_vertices[sg1].begin(); it2!=map_segment_vertices[sg1].end(); it2++){
     if ((*it2)->get_wcpt().index == sg1->get_wcpt_vec().front().index) start_v = *it2;
     if ((*it2)->get_wcpt().index == sg1->get_wcpt_vec().back().index) end_v = *it2;
   }
   if (start_v ==0 || end_v ==0) std::cout << "Error in finding vertices for a sgement" << std::endl;
   
   //	if (break_wcp.index==start_v->get_wcpt().index){
   if (sqrt(pow(start_v->get_wcpt().x - break_wcp.x,2)+pow(start_v->get_wcpt().y - break_wcp.y,2) + pow(start_v->get_wcpt().z - break_wcp.z,2)) < 0.9*units::cm){
     v1 = start_v;
   }else if (sqrt(pow(end_v->get_wcpt().x - break_wcp.x,2)+pow(end_v->get_wcpt().y - break_wcp.y,2) + pow(end_v->get_wcpt().z - break_wcp.z,2)) < 0.9*units::cm){
     //}else if (break_wcp.index == end_v->get_wcpt().index){
     v1 = end_v;
   }else{
     std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list1;
     std::list<WCP::WCPointCloud<double>::WCPoint> wcps_list2;
     temp_cluster->proto_break_tracks(start_v->get_wcpt(), break_wcp, end_v->get_wcpt(), wcps_list1, wcps_list2, true);
	  
     v1 = new WCPPID::ProtoVertex(acc_vertex_id, break_wcp, temp_cluster->get_cluster_id()); acc_vertex_id++;
     WCPPID::ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list1, temp_cluster->get_cluster_id()); acc_segment_id++;
     WCPPID::ProtoSegment *sg3 = new WCPPID::ProtoSegment(acc_segment_id, wcps_list2, temp_cluster->get_cluster_id()); acc_segment_id++;
     
     //	  std::cout << start_v->get_wcpt().index << " " << break_wcp.index << " " << end_v->get_wcpt().index << " " << wcps_list1.size() << " " << wcps_list2.size() << std::endl;
     //std::cout << map_vertex_segments.size() << " " << map_segment_vertices.size() << std::endl;
     add_proto_connection(start_v, sg2, temp_cluster);
     add_proto_connection(v1, sg2, temp_cluster);
     add_proto_connection(v1, sg3, temp_cluster);
     add_proto_connection(end_v, sg3, temp_cluster);
     del_proto_segment(sg1);
   } 
 }else{
   v1 = new WCPPID::ProtoVertex(acc_vertex_id, wcp, temp_cluster->get_cluster_id()); acc_vertex_id++;
 }
 
 
 return v1;
}


std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, WCP::Point> WCPPID::NeutrinoID::check_end_point(WCP::PointVector& tracking_path, bool flag_front ){

  if (tracking_path.size()<2) std::cout << "vector size wrong!" << std::endl;
  int ncount = 5;
  Point test_p;
  if (flag_front){
    test_p = tracking_path.front();
  }else{
    test_p = tracking_path.back();
  }
  WCPPID::ProtoVertex *vtx=0;
  WCPPID::ProtoSegment *seg=0;
  
  Point dir_p(0,0,0);
  for (int i=0;i!=ncount;i++){
    if (flag_front){
      if (i+1<tracking_path.size()){
	dir_p.x = test_p.x - tracking_path.at(i+1).x;
	dir_p.y = test_p.y - tracking_path.at(i+1).y;
	dir_p.z = test_p.z - tracking_path.at(i+1).z;
      }else{
	dir_p.x = test_p.x - tracking_path.back().x;
	dir_p.y = test_p.y - tracking_path.back().y;
	dir_p.z = test_p.z - tracking_path.back().z;
      }
    }else{
      if (tracking_path.size() >= i+2){
	dir_p.x = test_p.x - tracking_path.at(tracking_path.size()-2-i).x;
	dir_p.y = test_p.y - tracking_path.at(tracking_path.size()-2-i).y;
	dir_p.z = test_p.z - tracking_path.at(tracking_path.size()-2-i).z;
      }else{
	dir_p.x = test_p.x - tracking_path.front().x;
	dir_p.x = test_p.y - tracking_path.front().y;
	dir_p.x = test_p.z - tracking_path.front().z;
      }
    }
    TVector3 temp_dir(dir_p.x, dir_p.y, dir_p.z);
    if (temp_dir.Mag() == 0 ) continue;
    
    temp_dir = temp_dir.Unit();
    Line temp_l(test_p, temp_dir);

    // check vertex
    for (auto it1 = map_vertex_segments.begin(); it1!= map_vertex_segments.end(); it1++){
      WCPPID::ProtoVertex *test_v = it1->first;
      Point p1(test_v->get_wcpt().x, test_v->get_wcpt().y, test_v->get_wcpt().z);
      Point p2 = test_v->get_fit_pt();
      
      double dis1 = sqrt(pow(test_p.x - p1.x,2)
			 + pow(test_p.y - p1.y,2)
			 + pow(test_p.z - p1.z,2));
      double dis2 = temp_l.closest_dis(p1);
      double dis3 = sqrt(pow(test_p.x - p2.x,2)
			 + pow(test_p.y - p2.y,2)
			 + pow(test_p.z - p2.z,2));
      double dis4 = temp_l.closest_dis(p2);
      
      //      std::cout << i << " " << test_p << " " << temp_dir.X() << " " << temp_dir.Y() << " " << temp_dir.Z() << " " << dis1/units::cm << " " << dis2/units::cm << " " << dis3/units::cm << " " << dis4/units::cm << " " << it1->second.size() << " "  << std::endl;
      
      if (std::min(dis1,dis2) < 0.9*units::cm && std::max(dis1,dis2) < 2.0*units::cm ||
	  std::min(dis3,dis4) < 0.9*units::cm && std::max(dis3,dis4) < 2.0*units::cm){
	vtx = test_v;
	break;
      }
    }

    if (vtx!=0) break;
  }
   
  // check segment 
  if (vtx == 0){
    WCPPID::ProtoSegment* min_sg=0;
    Point min_point(0,0,0);
    double min_dis = 1e9;
      
    for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
      std::pair<double, WCP::Point> results = it->first->get_closest_point(test_p);
      if (results.first < 2*units::cm){
	PointVector& points = it->first->get_point_vec();

	Point dir_p(0,0,0);
	for (int i=0;i!=ncount;i++){
	  if (flag_front){
	    if (i+1<tracking_path.size()){
	      dir_p.x = test_p.x - tracking_path.at(i+1).x;
	      dir_p.y = test_p.y - tracking_path.at(i+1).y;
	      dir_p.z = test_p.z - tracking_path.at(i+1).z;
	    }else{
	      dir_p.x = test_p.x - tracking_path.back().x;
	      dir_p.y = test_p.y - tracking_path.back().y;
	      dir_p.z = test_p.z - tracking_path.back().z;
	    }
	  }else{
	    if (tracking_path.size() >= i+2){
	      dir_p.x = test_p.x - tracking_path.at(tracking_path.size()-2-i).x;
	      dir_p.y = test_p.y - tracking_path.at(tracking_path.size()-2-i).y;
	      dir_p.z = test_p.z - tracking_path.at(tracking_path.size()-2-i).z;
	    }else{
	      dir_p.x = test_p.x - tracking_path.front().x;
	      dir_p.x = test_p.y - tracking_path.front().y;
	      dir_p.x = test_p.z - tracking_path.front().z;
	    }
	  }
	  TVector3 temp_dir(dir_p.x, dir_p.y, dir_p.z);
	  if (temp_dir.Mag()==0) continue;
	  temp_dir = temp_dir.Unit();
	  Line temp_l(test_p, temp_dir);
	
	  for (size_t i= 0; i!=points.size();i++){
	    double dis1 = temp_l.closest_dis(points.at(i) );
	    if (dis1 < 1.2*units::cm && dis1 < min_dis){
	      min_dis = dis1;
	      min_point = points.at(i);
	      min_sg = it->first;
	    }
	  }
	}
      }
    }
    
    if (min_sg!=0){
      seg = min_sg;
      test_p = min_point;
    }
  }
    
  
  // return ...
  return std::make_tuple(vtx, seg, test_p);
}



bool WCPPID::NeutrinoID::add_proto_connection(WCPPID::ProtoVertex *pv, WCPPID::ProtoSegment *ps, WCPPID::PR3DCluster* cluster){

  if (pv->get_wcpt().index != ps->get_wcpt_vec().front().index && pv->get_wcpt().index != ps->get_wcpt_vec().back().index){
    std::cout << "Error! Vertex and Segment does not match " << pv->get_wcpt().index << " " << ps->get_wcpt_vec().front().index << " " << ps->get_wcpt_vec().back().index << std::endl;
    return false;
  }

  if (map_vertex_cluster.find(pv)==map_vertex_cluster.end())
    proto_vertices.push_back(pv);

  map_vertex_cluster[pv] = cluster;
  if (map_cluster_vertices.find(cluster)!=map_cluster_vertices.end()){
    map_cluster_vertices[cluster].insert(pv);
  }else{
    ProtoVertexSet temp_vertices;
    temp_vertices.insert(pv);
    map_cluster_vertices[cluster] = temp_vertices;
  }

  if (map_segment_cluster.find(ps)==map_segment_cluster.end())
    proto_segments.push_back(ps);
  
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
  if (map_vertex_segments.find(pv)!=map_vertex_segments.end()){
    // overall
    //    proto_vertices.erase(proto_vertices.find(pv));
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
  if (map_segment_vertices.find(ps)!=map_segment_vertices.end()){
    // overall
    //    proto_segments.erase(proto_segments.find(ps));
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

void WCPPID::NeutrinoID::examine_vertices(WCPPID::PR3DCluster* temp_cluster){
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_offset = mp.get_time_offset();
  int nrebin = mp.get_nrebin();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w + 0.5;
  double offset_u = -first_u_dis/pitch_u + 0.5;
  double offset_v = -first_v_dis/pitch_v + 0.5;
  
  double first_t_dis = temp_cluster->get_point_cloud()->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - temp_cluster->get_point_cloud()->get_cloud().pts[0].x;
  double slope_xt = 1./time_slice_width;
  double offset_t =  first_t_dis/time_slice_width + 0.5;
  
  bool flag_continue = true;

  while (flag_continue){
    flag_continue = false;
    //std::map<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> map_vertex_replace;
    WCPPID::ProtoVertex *v1 = 0;
    WCPPID::ProtoVertex *v2 = 0;
    WCPPID::ProtoVertex *v3 = 0;
    
    for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
      WCPPID::ProtoVertex *vtx = it->first;
      if (it->second.size()==2){ // potential check
	for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	  WCPPID::ProtoSegment* sg = (*it1);
	  if (sg->get_length() > 4*units::cm) continue;
	  for (auto it2 = map_segment_vertices[sg].begin(); it2 != map_segment_vertices[sg].end(); it2++){
	    WCPPID::ProtoVertex *vtx1 = (*it2);
	    if (vtx1 == vtx) continue;
	    if (map_vertex_segments[vtx1].size() > 2){ // the other track needs to be larger than 2
	      if (examine_vertices(vtx, vtx1, offset_t, slope_xt, offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw)){
		v1 = vtx;
		v2 = vtx1;

		for (auto it3 = it->second.begin(); it3 != it->second.end(); it3++){
		  if (*it3 == sg) continue;
		  for (auto it4 = map_segment_vertices[*it3].begin(); it4!=map_segment_vertices[*it3].end();it4++){
		    if (*it4 == v1) continue;
		    v3 = *it4;
		  }
		}
		
		flag_continue = true;
		break;
	      }
	    }
	  }
	  if (flag_continue) break;
	}
	if (flag_continue) break;
      }
    }

    //    std::cout << map_vertex_replace.size() << std::endl;
    // replace ...
    if (v1!=0 && v2!=0){
      std::cout << "Merge vertices" << std::endl;
      
      ProtoSegment *sg = find_segment(v1,v2);
      ProtoSegment *sg1 = find_segment(v1,v3);

      //  std::cout << sg<< " " << sg1 << std::endl;
      temp_cluster->dijkstra_shortest_paths(v3->get_wcpt(),2); 
      temp_cluster->cal_shortest_path(v2->get_wcpt(),2);
      ProtoSegment *sg2 = new WCPPID::ProtoSegment(acc_segment_id, temp_cluster->get_path_wcps(), temp_cluster->get_cluster_id()); acc_segment_id++;
      add_proto_connection(v2, sg2, temp_cluster);
      add_proto_connection(v3, sg2, temp_cluster);
      
      del_proto_vertex(v1);
      del_proto_segment(sg);
      del_proto_segment(sg1);

      temp_cluster->do_multi_tracking(map_vertex_segments, map_segment_vertices, *ct_point_cloud, global_wc_map, flash_time*units::microsecond, true, true, true);
      
    }
    
    // hack for now
    // flag_continue = false;
  }
  

  
}

bool WCPPID::NeutrinoID::examine_vertices(WCPPID::ProtoVertex* v1, WCPPID::ProtoVertex *v2, double offset_t, double slope_xt, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw){
  // U, V, W ...
  Point& v1_p = v1->get_fit_pt();
  double v1_t = offset_t + slope_xt * v1_p.x;
  double v1_u = offset_u + slope_yu * v1_p.y + slope_zu * v1_p.z;
  double v1_v = offset_v + slope_yv * v1_p.y + slope_zv * v1_p.z;
  double v1_w = offset_w + slope_yw * v1_p.y + slope_zw * v1_p.z;
  
  Point& v2_p = v2->get_fit_pt();
  double v2_t = offset_t + slope_xt * v2_p.x;
  double v2_u = offset_u + slope_yu * v2_p.y + slope_zu * v2_p.z;
  double v2_v = offset_v + slope_yv * v2_p.y + slope_zv * v2_p.z;
  double v2_w = offset_w + slope_yw * v2_p.y + slope_zw * v2_p.z;

  // find the track segment in between
  WCPPID::ProtoSegment *sg = find_segment(v1,v2);
  if (sg == 0 || map_vertex_segments[v1].size()!=2) return false;
  
  
  int ncount_close = 0;
  int ncount_dead = 0;
  int ncount_line = 0;

  // one view must be close
  if (sqrt(pow(v1_u - v2_u, 2) + pow(v1_t - v2_t,2)) < 2.0){
    ncount_close ++;
  }else{
    // one view must be dead? or shorter distance
    bool flag_dead = true;
    for (size_t i=0;i!=sg->get_point_vec().size();i++){
      if (!ct_point_cloud->get_closest_dead_chs(sg->get_point_vec().at(i),0))
	flag_dead = false;
    }
    if (flag_dead) ncount_dead++;
    else{
      // third view must be like a line
      PointVector& pts_1 = sg->get_point_vec();
      WCPPID::ProtoSegment *sg1 = 0;
      for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	// sg is one track
	if (*it1 == sg) continue;
	sg1 = *it1;
      }

      PointVector& pts_2 = sg1->get_point_vec();
      TVector3 v1(v2_u-v1_u, v2_t-v1_t, 0);
      TVector3 v2(0,0,0);
      double min_dis = 1e9;
      
      for (size_t i=0;i!=pts_2.size();i++){
	double p_t = offset_t + slope_xt * pts_2.at(i).x;
	double p_u = offset_u + slope_yu * pts_2.at(i).y + slope_zu * pts_2.at(i).z;
	TVector3 v3(p_u-v1_u, p_t-v1_t,0);
	double dis = fabs(v3.Mag()-9);
	if (dis < min_dis){
	  min_dis = dis;
	  v2 = v3;
	}
      }

      if (180-v1.Angle(v2)/3.1415926*180. < 30) ncount_line++;
      // std::cout << v1.Mag() << " " << v2.Mag() << " " << v1.Angle(v2)/3.1415926*180. << std::endl;
    }
  }

  //  std::cout << ncount_close << " " << ncount_dead << " " << ncount_line << std::endl;
  
   // one view must be close
  if (sqrt(pow(v1_v - v2_v, 2) + pow(v1_t - v2_t,2)) < 2.0){
    ncount_close ++;
  }else{
    // one view must be dead? or shorter distance
    bool flag_dead = true;
    for (size_t i=0;i!=sg->get_point_vec().size();i++){
      if (!ct_point_cloud->get_closest_dead_chs(sg->get_point_vec().at(i),1))
	flag_dead = false;
    }
    if (flag_dead) ncount_dead++;
    else{
      // third view must be like a line
      PointVector& pts_1 = sg->get_point_vec();
      WCPPID::ProtoSegment *sg1 = 0;
      for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	// sg is one track
	if (*it1 == sg) continue;
	sg1 = *it1;
      }

      PointVector& pts_2 = sg1->get_point_vec();
      TVector3 v1(v2_v-v1_v, v2_t-v1_t, 0);
      TVector3 v2(0,0,0);
      double min_dis = 1e9;
      
      for (size_t i=0;i!=pts_2.size();i++){
	double p_t = offset_t + slope_xt * pts_2.at(i).x;
	double p_v = offset_v + slope_yv * pts_2.at(i).y + slope_zv * pts_2.at(i).z;
	TVector3 v3(p_v-v1_v, p_t-v1_t,0);
	double dis = fabs(v3.Mag()-9);
	if (dis < min_dis){
	  min_dis = dis;
	  v2 = v3;
	}
      }

      if (180-v1.Angle(v2)/3.1415926*180. < 30) ncount_line++;
      //      std::cout << v1.Mag() << " " << v2.Mag() << " " << v1.Angle(v2)/3.1415926*180. << std::endl;
    }
  }

  //std::cout << ncount_close << " " << ncount_dead << " " << ncount_line << std::endl;

   // one view must be close
  if (sqrt(pow(v1_w - v2_w, 2) + pow(v1_t - v2_t,2)) < 2.0){
    //    std::cout << "haha " << std::endl;
    ncount_close ++;
  }else{
    // one view must be dead? or shorter distance
    bool flag_dead = true;
    for (size_t i=0;i!=sg->get_point_vec().size();i++){
      if (!ct_point_cloud->get_closest_dead_chs(sg->get_point_vec().at(i),2))
	flag_dead = false;
    }

    
    
    if (flag_dead) ncount_dead++;
    else{
      // third view must be like a line
      PointVector& pts_1 = sg->get_point_vec();
      WCPPID::ProtoSegment *sg1 = 0;
      for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
	// sg is one track
	if (*it1 == sg) continue;
	sg1 = *it1;
      }

      PointVector& pts_2 = sg1->get_point_vec();
      TVector3 v1(v2_w-v1_w, v2_t-v1_t, 0);
      TVector3 v2(0,0,0);
      double min_dis = 1e9;
      
      for (size_t i=0;i!=pts_2.size();i++){
	double p_t = offset_t + slope_xt * pts_2.at(i).x;
	double p_w = offset_w + slope_yw * pts_2.at(i).y + slope_zw * pts_2.at(i).z;
	TVector3 v3(p_w-v1_w, p_t-v1_t,0);
	double dis = fabs(v3.Mag()-9);
	if (dis < min_dis){
	  min_dis = dis;
	  v2 = v3;
	}
      }

      if (180-v1.Angle(v2)/3.1415926*180. < 30) ncount_line++;

      //      std::cout << v1.Mag() << " " << v2.Mag() << " " << v1.Angle(v2)/3.1415926*180. << std::endl;
    }
  }

  
  
  if (ncount_close >=2 ||
      ncount_close ==1 && ncount_dead ==1 & ncount_line==1 ||
      ncount_close ==1 && ncount_dead == 2)
    return true;

  //  std::cout << ncount_close << " " << ncount_dead << " " << ncount_line << std::endl;
  //  std::cout <<  << " " << sqrt(pow(v1_u - v2_u, 2) + pow(v1_t - v2_t,2)) << " " << sqrt(pow(v1_v - v2_v, 2) + pow(v1_t - v2_t,2)) << std::endl;

  
  return false;
}



WCPPID::ProtoSegment* WCPPID::NeutrinoID::find_segment(WCPPID::ProtoVertex *v1, WCPPID::ProtoVertex *v2){

  WCPPID::ProtoSegment* sg = 0;
  for (auto it1 = map_vertex_segments[v1].begin(); it1 != map_vertex_segments[v1].end(); it1++){
    sg = *it1;
    if (map_segment_vertices[sg].find(v2) == map_segment_vertices[sg].end()){
      sg = 0;
    }else{
      break;
    }
  }
  return sg;
}
