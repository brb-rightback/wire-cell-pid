void WCPPID::NeutrinoID::shower_clustering(){
  shower_clustering_in_main_cluster();

  
  update_shower_maps();  
}

void WCPPID::NeutrinoID::update_shower_maps(){

  map_vertex_to_shower.clear();
  map_vertex_in_shower.clear();
  map_segment_in_shower.clear();
  
  for (auto it = showers.begin(); it!=showers.end(); it++){
    WCPPID::WCShower* shower = *it;
    map_vertex_to_shower[shower->get_start_vertex().first] = shower;
    shower->fill_maps(map_vertex_in_shower, map_segment_in_shower);
  }
}


void WCPPID::NeutrinoID::shower_clustering_in_main_cluster(){
  // search from main vertex ...
  // search trees, if find an electron, then the rest are all added to it ...
  

  // main_vertex figure out the daughters and mother ... 
  std::set<WCPPID::ProtoSegment* > used_segments; 

  std::vector<std::pair<WCPPID::ProtoSegment*, WCPPID::ProtoVertex*> > segments_to_be_examined; // current_segment, daughter vertex
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){ 
    // parent are all zero now ...
    WCPPID::ProtoVertex *other_vertex = find_other_vertex(*it, main_vertex);
    segments_to_be_examined.push_back(std::make_pair(*it, other_vertex)); 
  } 
  
  while(segments_to_be_examined.size()>0){ 
    std::vector<std::pair<WCPPID::ProtoSegment*, WCPPID::ProtoVertex*> > temp_segments; 
    for (auto it = segments_to_be_examined.begin(); it!= segments_to_be_examined.end(); it++){ 
      WCPPID::ProtoSegment *curr_sg = it->first;
      WCPPID::ProtoVertex *daughter_vtx = it->second;
      used_segments.insert(curr_sg);
      
      if (curr_sg->get_flag_shower() ){
	WCPPID::ProtoVertex *parent_vtx = find_other_vertex(curr_sg, daughter_vtx);
	WCPPID::WCShower *shower = new WCPPID::WCShower();
	shower->set_start_vertex(parent_vtx, 1);
	shower->set_start_segment(curr_sg);
	showers.push_back(shower);
      }else{
	// keep searching its daughter
	for (auto it1 = map_vertex_segments[daughter_vtx].begin(); it1 != map_vertex_segments[daughter_vtx].end(); it1++){
	  WCPPID::ProtoSegment *next_sg = *it1;
	  if (used_segments.find(next_sg)!=used_segments.end()) continue;
	  WCPPID::ProtoVertex *other_vertex = find_other_vertex(next_sg, daughter_vtx);
	  temp_segments.push_back(std::make_pair(next_sg, other_vertex));
	}
      }
    } 
    segments_to_be_examined = temp_segments; 
  }

  // complete the shower construction
  for (size_t i=0;i!=showers.size();i++){
    showers.at(i)->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, used_segments); 
  }
  
  //  std::cout << showers.size() << " " << used_segments.size() << std::endl;
}
