void WCPPID::NeutrinoID::particle_clustering(){

  particle_clustering_in_main_cluster();
}


void WCPPID::NeutrinoID::particle_clustering_in_main_cluster(){
  // search from main vertex ...

  /* // main_vertex figure out the daughters and mother ... */
  /* std::set<WCPPID::ProtoVertex* > used_vertices; */
  /* std::set<WCPPID::ProtoSegment* > used_segments; */

  /* std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > segments_to_be_examined; */
  /* for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){ */
  /*   // parent are all zero now ... */
  /*   used_segments.insert(*it); */
  /*   WCPPID::ProtoVertex *other_vertex = find_other_vertex(*it, main_vertex); */
  /*   segments_to_be_examined.push_back(std::make_pair(other_vertex, *it)); */
  /* } */
  /* used_vertices.insert(main_vertex); */

  /*  while(segments_to_be_examined.size()>0){ */
  /*   std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > temp_segments; */
  /*   for (auto it = segments_to_be_examined.begin(); it!= segments_to_be_examined.end(); it++){ */
  /*     WCPPID::ProtoVertex *curr_vtx = it->first; */
  /*     WCPPID::ProtoSegment *prev_sg = it->second; */
  /*     if (used_vertices.find(curr_vtx)!=used_vertices.end()) continue; */

  /*     for (auto it1 = map_vertex_segments[curr_vtx].begin(); it1!=map_vertex_segments[curr_vtx].end(); it1++){ */
  /* 	WCPPID::ProtoSegment *curr_sg = *it1; */
  /* 	if (used_segments.find(curr_sg)!=used_segments.end()) continue; */
  /* 	used_segments.insert(curr_sg); */
  /* 	// set mother ... */
  /* 	rtree.mc_mother[map_sgid_rtid[map_sg_sgid[curr_sg]]] = map_sg_sgid[prev_sg]; */
  /* 	// set daughters ... */
  /* 	rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(map_sg_sgid[curr_sg]); */
  /* 	WCPPID::ProtoVertex *other_vertex = find_other_vertex(curr_sg, curr_vtx); */
  /* 	if (used_vertices.find(other_vertex) == used_vertices.end()) */
  /* 	  temp_segments.push_back(std::make_pair(other_vertex, curr_sg)); */
  /*     } */
  /*     used_vertices.insert(curr_vtx); */
  /*   } */
  /*   segments_to_be_examined = temp_segments; */
  /* } */
   
}
