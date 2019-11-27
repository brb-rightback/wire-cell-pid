
void WCPPID::PR3DCluster::search_other_tracks(double search_range){
  if (fit_tracks.size()==0 && fine_tracking_path.size() == 0) return;
    
  if (fit_tracks.size()==0){
    WCP::TrackInfo *track = new WCP::TrackInfo(fine_tracking_path, dQ, dx, pu, pv, pw, pt, reduced_chi2);
    fit_tracks.push_back(track);
  }
  
  // number of terminals ... 
  const int N = point_cloud_steiner->get_num_points();
  std::vector<bool> flag_tagged(N, false);
  int num_tagged = 0;
  WCP::WCPointCloud<double>& cloud = point_cloud_steiner->get_cloud();
  for (size_t i=0;i!=N;i++){
    Point p(cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z);
    double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
    for (size_t j=0;j!=fit_tracks.size();j++){
      std::pair<double, WCP::Point> closest_dis_point = fit_tracks.at(j)->get_closest_point(p);
      std::tuple<double, double, double> closest_2d_dis = fit_tracks.at(j)->get_closest_2d_dis(p);
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
      if (min_dis_u < 0.75 * search_range && min_dis_v < 0.75 * search_range && min_dis_w < 0.75*search_range)
	flag_tagged[i] = true;
	  // std::cout << min_dis_u/units::cm << " " << min_dis_v/units::cm << " " << min_dis_w/units::cm << std::endl;
    }
  }

  //  std::cout << num_tagged << " " << N << std::endl;


  // Now figure out the TerminalGraph ...
  std::vector<int> terminals;
  std::map<int, int> map_oindex_tindex;
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
    
    //    std::cout << N << " " << terminals.size() << " " << map_saved_edge.size() << " " << terminal_edge.size() << std::endl;
    
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
    std::cout << num <<  " " << component.size() << " " << map_connection.size() << std::endl;
}
