void WCPPID::PR3DCluster::clustering_points_master(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud, double search_range, double scaling_2d){

  int ntrack = 0;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    if (it->first->get_cluster_id()==cluster_id) ntrack++;
  }
  
  if (ntrack > 0){
    /* // steiner-inspired graph */
    /* clustering_points(map_vertex_segments, map_segment_vertices, ct_point_cloud, 2); */
    
    /* //  for (size_t i=0;i!=point_steiner_sub_cluster_ids.size();i++){ */
    /* // std::cout << point_steiner_sub_cluster_ids.at(i) << std::endl; */
    /* //} */
    
    /* ToyPointCloud *pcloud = new ToyPointCloud(); */
    /* PointVector pts; */
    /* WCP::WCPointCloud<double>& cloud = point_cloud_steiner->get_cloud(); */
    /* for (size_t i=0;i!=cloud.pts.size();i++){ */
    /*   if (point_steiner_sub_cluster_ids.at(i)!=-1){ */
    /* 	Point p(cloud.pts[i].x,cloud.pts[i].y, cloud.pts[i].z); */
    /* 	pts.push_back(p); */
    /*   } */
    /* } */
    /* pcloud->AddPoints(pts); */
    /* pcloud->build_kdtree_index(); */
    
    // regular ...
    clustering_points(map_vertex_segments, map_segment_vertices, ct_point_cloud, 1);
    //delete pcloud;
  }

  
}

  void WCPPID::PR3DCluster::clustering_points(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud, int choice, WCP::ToyPointCloud* pcloud, double search_range, double scaling_2d){

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
  
  std::map<size_t, WCPPID::ProtoSegment*> map_pindex_segment;
  
  // define steiner terminal for segments ...
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){
    if (it->first->get_cluster_id()!=cluster_id) continue;
    std::vector<WCP::Point >& pts = it->first->get_point_vec();

    if (pts.size()>2){
      for (size_t i=1;i+1<pts.size();i++){
	std::vector<std::pair<size_t, double> > temp_results = temp_point_cloud->get_closest_index(pts.at(i),5);
	for (auto it1 = temp_results.begin(); it1!=temp_results.end(); it1++){
	  if (map_pindex_segment.find(it1->first)==map_pindex_segment.end()){
	    map_pindex_segment[it1->first] = it->first;
	    break;
	  }
	}
      }
    }else{
      Point p((pts.front().x+pts.back().x)/2., (pts.front().y+pts.back().y)/2., (pts.front().z+pts.back().z)/2.);
      std::vector<std::pair<size_t, double> > temp_results = temp_point_cloud->get_closest_index(p,5);
      for (auto it1 = temp_results.begin(); it1!=temp_results.end(); it1++){
	if (map_pindex_segment.find(it1->first)==map_pindex_segment.end()){
	  map_pindex_segment[it1->first] = it->first;
	  break;
	}
      }
    }
  }

 
  
  // these are steiner terminals ...
  if (map_pindex_segment.size()>0){
    std::vector<int> terminals;
    for (auto it = map_pindex_segment.begin(); it!=map_pindex_segment.end(); it++){
      // std::cout << it->first << " " << it->second->get_cluster_id() << std::endl;
      terminals.push_back(it->first);
    }

    
    const int N = temp_point_cloud->get_num_points();
    
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
    
    EdgeWeightMap edge_weight = boost::choose_pmap(boost::get_param(boost::no_named_parameters(), boost::edge_weight), *temp_graph, boost::edge_weight);
    
    // distance array used in the dijkstra runs
    std::vector<Weight> distance(N);
    std::vector<Vertex> nearest_terminal(num_vertices(*temp_graph));
    auto index = get(boost::vertex_index, *temp_graph);
    
    auto nearest_terminal_map = boost::make_iterator_property_map(nearest_terminal.begin(), get(boost::vertex_index, *temp_graph));
    for (auto terminal : terminals) {
      nearest_terminal_map[terminal] = terminal;
    }
    // compute voronoi diagram each vertex get nearest terminal and last edge on
    // path to nearest terminal
    auto distance_map = make_iterator_property_map(distance.begin(), index);
    std::vector<Edge> vpred(N);
    auto last_edge = boost::make_iterator_property_map(vpred.begin(), get(boost::vertex_index, *temp_graph));
    boost::dijkstra_shortest_paths(*temp_graph, terminals.begin(), terminals.end(), boost::dummy_property_map(),
				   distance_map, edge_weight, index, paal::utils::less(),
				   boost::closed_plus<Weight>(), std::numeric_limits<Weight>::max(), 0,
				   boost::make_dijkstra_visitor(paal::detail::make_nearest_recorder(
												    nearest_terminal_map, last_edge, boost::on_edge_relaxed{})));

    /* if (cluster_id==66)  { */
    /* /\*   typename boost::graph_traits<WCPPID::MCUGraph>::out_edge_iterator ei, ei_end; *\/ */
    /* /\*   for (boost::tie(ei, ei_end) = out_edges(1676, *temp_graph); ei != ei_end; ++ei){ *\/ */
    /* /\* 	auto source = boost::source(*ei, *temp_graph); *\/ */
    /* /\* 	auto target = boost::target(*ei, *temp_graph); *\/ */
    /* /\* 	std::cout << "There is an edge from " << source << " to " << target << std::endl;  *\/ */
    /* /\*   } *\/ */
    /*   std::vector<int> component(num_vertices(*temp_graph)); */
    /*   const int num = connected_components(*temp_graph,&component[0]); */
    /*   std::cout << num << std::endl; */
    /* }  */
    
    for (size_t i=0;i!=nearest_terminal.size();i++){
      int index = nearest_terminal.at(i);
      if (map_pindex_segment.find(index)==map_pindex_segment.end()) {
	//	std::cout << "Wrong, no terminal found!  " << index << " " << i << std::endl;
	temp_points_ids->at(i) = -1;
      }else{
	temp_points_ids->at(i) = map_pindex_segment[index]->get_id();
      }
      //std::cout << i << " " << nearest_terminal.at(i) << std::endl;
      //
    }
    //    std::cout << nearest_terminal.size() << " " << N << std::endl;
    
    WCPPID::ProtoSegmentSelection temp_segments;
    for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
      //if (it->first->get_cluster_id() != cluster_id) continue;
      // consider all segments ...
      temp_segments.push_back(it->first);
    }

    
    // now examine to remove ghost points ....
    for (size_t i=0;i!=N;i++){
      if (map_pindex_segment.find(nearest_terminal.at(i)) == map_pindex_segment.end()) continue;
      Point p(cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z);
      WCPPID::ProtoSegment *main_sg = map_pindex_segment[nearest_terminal.at(i)];

      std::pair<double, WCP::Point> closest_dis_point = main_sg->get_closest_point(p);
      std::tuple<double, double, double> closest_2d_dis = main_sg->get_closest_2d_dis(p);

      std::tuple<double, double, double> min_2d_dis = closest_2d_dis;
      
      // check against main_sg;
      bool flag_change = true;
      
      for (auto it = temp_segments.begin(); it!=temp_segments.end(); it++){
	if (main_sg == (*it)) continue;
	std::tuple<double, double, double> temp_2d_dis = (*it)->get_closest_2d_dis(p);
	if (std::get<0>(temp_2d_dis) < std::get<0>(min_2d_dis)) std::get<0>(min_2d_dis) = std::get<0>(temp_2d_dis);
	if (std::get<1>(temp_2d_dis) < std::get<1>(min_2d_dis)) std::get<1>(min_2d_dis) = std::get<1>(temp_2d_dis);
	if (std::get<2>(temp_2d_dis) < std::get<2>(min_2d_dis)) std::get<2>(min_2d_dis) = std::get<2>(temp_2d_dis);
      }
      
      if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) && std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis)) // all closest
	flag_change = false;
      else if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) ) //&& (std::get<2>(closest_2d_dis) < scaling_2d * search_range || closest_dis_point.first < search_range)) // 2 closest
	flag_change = false;
      else if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis) ) //&& (std::get<1>(closest_2d_dis) < scaling_2d * search_range || closest_dis_point.first < search_range))
	flag_change = false;
      else if (std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) && std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis) ) //&& (std::get<0>(closest_2d_dis) < scaling_2d * search_range || closest_dis_point.first < search_range))
	flag_change = false;
      else if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && (closest_dis_point.first < search_range ||std::get<1>(closest_2d_dis) < scaling_2d * search_range &&  std::get<2>(closest_2d_dis) < scaling_2d * search_range) )
	flag_change = false;
      else if (std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) && (closest_dis_point.first < search_range ||std::get<0>(closest_2d_dis) < scaling_2d * search_range &&  std::get<2>(closest_2d_dis) < scaling_2d * search_range) )
	flag_change = false;
      else if (std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis) && (closest_dis_point.first < search_range ||std::get<1>(closest_2d_dis) < scaling_2d * search_range &&  std::get<0>(closest_2d_dis) < scaling_2d * search_range) )
	flag_change = false;

      // deal with dead channels ...
      if (!flag_change){
	if (ct_point_cloud.get_closest_dead_chs(p,0) && std::get<0>(closest_2d_dis) > scaling_2d * search_range){
	  if (std::get<1>(closest_2d_dis) < scaling_2d * search_range ||  std::get<2>(closest_2d_dis) < scaling_2d * search_range)
	    flag_change = true;
	}else if (ct_point_cloud.get_closest_dead_chs(p,1) && std::get<1>(closest_2d_dis) > scaling_2d * search_range){
	  if (std::get<0>(closest_2d_dis) < scaling_2d * search_range ||  std::get<2>(closest_2d_dis) < scaling_2d * search_range)
	    flag_change = true;
	}else if (ct_point_cloud.get_closest_dead_chs(p,2) && std::get<2>(closest_2d_dis) > scaling_2d * search_range){
	  if (std::get<1>(closest_2d_dis) < scaling_2d * search_range ||  std::get<0>(closest_2d_dis) < scaling_2d * search_range)
	    flag_change = true;
	}
      }
      
      
      /* if (closest_dis_point.first < search_range){ */
      /* 	flag_change = false; */
      /* }else if (std::get<0>(closest_2d_dis) < scaling_2d * search_range && std::get<1>(closest_2d_dis) < scaling_2d * search_range && std::get<2>(closest_2d_dis) < scaling_2d * search_range){ */
      /* 	flag_change = false; */
	
      /* 	// need to be close to the steiner tree somehow ... */
      /* 	if(!flag_change  && pcloud != 0){ */
      /* 	  std::pair<int,double> test_u = pcloud->get_closest_2d_dis(p, 0); */
      /* 	  std::pair<int,double> test_v = pcloud->get_closest_2d_dis(p, 1); */
      /* 	  std::pair<int,double> test_w = pcloud->get_closest_2d_dis(p, 2); */
      /* 	  double dis = pcloud->get_closest_dis(p); */
      /* 	  int ncount = 0; */
      /* 	  if (test_u.second > search_range * scaling_2d) ncount++; */
      /* 	  if (test_v.second > search_range * scaling_2d) ncount++; */
      /* 	  if (test_w.second > search_range * scaling_2d) ncount++; */
      /* 	  if (dis > search_range) ncount ++; */
      /* 	  if (ncount > 0){ */
      /* 	    flag_change = true; */
      /* 	  } */
      /* 	} */
      /* }else{ */
      /* 	// check against all_sg;  These likely to be ghosts ... */
      /* 	double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9; */
      /* 	for (auto it = temp_segments.begin(); it!= temp_segments.end(); it++){ */
      /* 	  WCPPID::ProtoSegment *sg1 = *it; */
      /* 	  //std::pair<double, WCP::Point> closest_dis_point = sg1->get_closest_point(p); */
      /* 	  std::tuple<double, double, double> closest_2d_dis = sg1->get_closest_2d_dis(p); */
      /* 	  if (std::get<0>(closest_2d_dis) < min_dis_u) min_dis_u = std::get<0>(closest_2d_dis); */
      /* 	  if (std::get<1>(closest_2d_dis) < min_dis_v) min_dis_v = std::get<1>(closest_2d_dis); */
      /* 	  if (std::get<2>(closest_2d_dis) < min_dis_w) min_dis_w = std::get<2>(closest_2d_dis); */
      /* 	} */
      /* 	if ((min_dis_u < scaling_2d * search_range || ct_point_cloud.get_closest_dead_chs(p, 0) ) && */
      /* 	    (min_dis_v < scaling_2d * search_range || ct_point_cloud.get_closest_dead_chs(p, 1) ) && */
      /* 	    (min_dis_w < scaling_2d * search_range || ct_point_cloud.get_closest_dead_chs(p, 2) ) ) */
      /* 	  flag_change = true; */
      /* 	else */
      /* 	  flag_change = false; */

      /* 	// need to be close to the steiner tree somehow ... */
      /* 	if(!flag_change  && pcloud != 0){ */
      /* 	  std::pair<int,double> test_u = pcloud->get_closest_2d_dis(p, 0); */
      /* 	  std::pair<int,double> test_v = pcloud->get_closest_2d_dis(p, 1); */
      /* 	  std::pair<int,double> test_w = pcloud->get_closest_2d_dis(p, 2); */
      /* 	  double dis = pcloud->get_closest_dis(p); */
      /* 	  int ncount = 0; */
      /* 	  if (test_u.second > search_range * scaling_2d) ncount++; */
      /* 	  if (test_v.second > search_range * scaling_2d) ncount++; */
      /* 	  if (test_w.second > search_range * scaling_2d) ncount++; */
      /* 	  if (dis > search_range) ncount ++; */
      /* 	  if (ncount > 0){ */
      /* 	    flag_change = true; */
      /* 	  } */
      /* 	} */
      /* } */
      
     
      // change the point's clustering ... 
      if (flag_change) temp_points_ids->at(i) = -1;
    }
  }
  
  
  //std::cout << nearest_terminal.size() << " " << temp_points_ids->size() << std::endl;
}
