void WCPPID::NeutrinoID::shower_clustering(){

  shower_clustering_in_main_cluster();
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


void WCPPID::WCShower::complete_structure_with_start_segment(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments){
  // fill the start segment ...
  std::vector<ProtoSegment* > new_segments;
  std::vector<ProtoVertex* > new_vertices;
  
  for (auto it = map_segment_vertices[start_segment].begin(); it!= map_segment_vertices[start_segment].end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (vtx == start_vertex) continue;
    map_vtx_segs[vtx].insert(start_segment);
    map_seg_vtxs[start_segment].insert(vtx);

    new_vertices.push_back(vtx);
  }

  while( new_vertices.size()>0 || new_segments.size()>0 ){
    if (new_vertices.size()>0){
      ProtoVertex *vtx = new_vertices.back();
      new_vertices.pop_back();
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
    	ProtoSegment *seg = *it;
    	if (used_segments.find(seg)!=used_segments.end()) continue;
    	map_vtx_segs[vtx].insert(seg);
    	map_seg_vtxs[seg].insert(vtx);
    	new_segments.push_back(seg);
	used_segments.insert(seg);
      }
    }

    if (new_segments.size()>0){
      ProtoSegment *seg = new_segments.back();
      new_segments.pop_back();
      for (auto it = map_segment_vertices[seg].begin(); it!= map_segment_vertices[seg].end(); it++){
  	ProtoVertex *vtx = *it;
  	if (map_vtx_segs.find(vtx)!= map_vtx_segs.end() || vtx == start_vertex) continue;
  	map_vtx_segs[vtx].insert(seg);
  	map_seg_vtxs[seg].insert(vtx);
  	new_vertices.push_back(vtx);
      }
    }
  }
  
  //std::cout << map_vtx_segs.size() << " " << map_seg_vtxs.size() << std::endl;
}
