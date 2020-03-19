void WCPPID::PR3DCluster::clustering_points(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices, int choice){

  MCUGraph *temp_graph = 0;
  WCP::ToyPointCloud *temp_point_cloud = 0;
  std::vector<int> *temp_points_ids = 0;
  if (choice == 1){
    temp_graph = graph;
    temp_point_cloud = point_cloud;
    temp_points_ids = &point_sub_cluster_ids;
  }else if (choice == 2){
    temp_graph = graph_steiner;
    temp_point_cloud = point_cloud_steiner;
    temp_points_ids = &point_steiner_sub_cluster_ids;
  }
  temp_points_ids->resize(temp_point_cloud->get_num_points(),0);
  WCP::WCPointCloud<double>& cloud = temp_point_cloud->get_cloud();
  
  // define steiner terminal for vertices ...
  // for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
  //  if (it->first->get_cluster_id()!=cluster_id) continue;
  //  Point tmp_p = it->first->get_fit_pt();
  //  
  //}

  std::map<size_t, WCPPID::ProtoSegment*> map_pindex_segment;
  
  // define steiner terminal for segments ...
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    if (it->first->get_cluster_id()!=cluster_id) continue;
    std::vector<WCP::Point >& pts = it->first->get_point_vec();
    for (size_t i=1;i+1<pts.size();i++){
      std::vector<std::pair<size_t, double> > temp_results = temp_point_cloud->get_closest_index(pts.at(i),5);
      for (auto it1 = temp_results.begin(); it1!=temp_results.end(); it1++){
	if (map_pindex_segment.find(it1->first)==map_pindex_segment.end()){
	  map_pindex_segment[it1->first] = it->first;
	  break;
	}
      }
    }
  }

  // these are steiner terminals ...
  //  std::cout << map_pindex_segment.size() << std::endl;

  std::vector<int> terminals;
  for (auto it = map_pindex_segment.begin(); it!=map_pindex_segment.end(); it++){
    terminals.push_back(it->first);
  }
  
  const int N = temp_point_cloud->get_num_points();

  using Vertex = typename boost::graph_traits<WCPPID::MCUGraph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<WCPPID::MCUGraph>::edge_descriptor;
  using Base = typename boost::property<edge_base_t, Edge>;
  // using EdgeWeightMap = typename boost::choose_const_pmap(boost::get_param(boost::no_named_parameters(), boost::edge_weight), *graph, boost::edge_weight);
  using EdgeWeightMap = typename boost::property_map<WCPPID::MCUGraph, boost::edge_weight_t>::type;
  using Weight = typename boost::property_traits<EdgeWeightMap>::value_type;
  using WeightProperty =
        typename boost::property<boost::edge_weight_t, Weight, Base>;
  using TerminalGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, WeightProperty>;
  using EdgeTerminal =
    typename boost::graph_traits<TerminalGraph>::edge_descriptor;

  EdgeWeightMap edge_weight = boost::choose_pmap(boost::get_param(boost::no_named_parameters(), boost::edge_weight), *graph, boost::edge_weight);
  
  // distance array used in the dijkstra runs
  std::vector<Weight> distance(N);
  std::vector<Vertex> nearest_terminal(num_vertices(*temp_graph));
  auto index = get(boost::vertex_index, *temp_graph);

  auto nearest_terminal_map = boost::make_iterator_property_map(nearest_terminal.begin(), get(boost::vertex_index, *graph));
  for (auto terminal : terminals) {
    nearest_terminal_map[terminal] = terminal;
  }
  // compute voronoi diagram each vertex get nearest terminal and last edge on
  // path to nearest terminal
  auto distance_map = make_iterator_property_map(distance.begin(), index);
  std::vector<Edge> vpred(N);
  auto last_edge = boost::make_iterator_property_map(vpred.begin(), get(boost::vertex_index, *graph));
  boost::dijkstra_shortest_paths(*graph, terminals.begin(), terminals.end(), boost::dummy_property_map(),
  				 distance_map, edge_weight, index, paal::utils::less(),
  				 boost::closed_plus<Weight>(), std::numeric_limits<Weight>::max(), 0,
  				 boost::make_dijkstra_visitor(paal::detail::make_nearest_recorder(
  											    nearest_terminal_map, last_edge, boost::on_edge_relaxed{})));

  for (size_t i=0;i!=nearest_terminal.size();i++){
    int index = nearest_terminal.at(i);
    //if (map_pindex_segment.find(index)==map_pindex_segment.end()) std::cout << "Wrong!  " << std::endl;
    temp_points_ids->at(i) = map_pindex_segment[index]->get_id();
    //std::cout << i << " " << nearest_terminal.at(i) << std::endl;
    //std::cout << nearest_terminal.size() << " " << N << std::endl;
  }

  // now examine to remove ghost points ....


  
  //std::cout << nearest_terminal.size() << " " << temp_points_ids->size() << std::endl;
}
