void WCPPID::NeutrinoID::shower_clustering(){

  // connect to the main cluster ...
  shower_clustering_in_main_cluster();
  shower_clustering_from_main_cluster();
  
  // check remaining clusters ...
  shower_clustering_in_other_clusters();
  
  calculate_shower_kinematics();
}

void WCPPID::NeutrinoID::shower_clustering_from_main_cluster(){
  
}

void WCPPID::NeutrinoID::shower_clustering_from_vertices(){

}


void WCPPID::NeutrinoID::calculate_shower_kinematics(){
  for (size_t i=0;i!=showers.size();i++){
    WCPPID::WCShower *shower = showers.at(i);
    shower->calculate_kinematics();
    double kine_charge = cal_kine_charge(shower);
    shower->set_kine_charge(kine_charge);
    
    //  std::cout << shower->get_kine_range()/units::MeV << " " << shower->get_kine_dQdx()/units::MeV << " " << shower->get_kine_charge()/units::MeV << std::endl;
  }
}

void WCPPID::NeutrinoID::update_shower_maps(){

  map_vertex_to_shower.clear();
  map_vertex_in_shower.clear();
  map_segment_in_shower.clear();
  used_shower_clusters.clear();
  
  for (auto it = showers.begin(); it!=showers.end(); it++){
    WCPPID::WCShower* shower = *it;
    map_vertex_to_shower[shower->get_start_vertex().first] = shower;
    shower->fill_maps(map_vertex_in_shower, map_segment_in_shower);
  }
  for (auto it = map_segment_in_shower.begin(); it!=map_segment_in_shower.end();it++){
    used_shower_clusters.insert(it->first->get_cluster_id());
  }
  
}

std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> WCPPID::NeutrinoID::get_start_end_vertices(WCPPID::ProtoSegment* sg){
  WCPPID::ProtoVertex *start_v=0, *end_v=0;
  for (auto it = map_segment_vertices[sg].begin(); it!=map_segment_vertices[sg].end(); it++){
    if ((*it)->get_wcpt().index == sg->get_wcpt_vec().front().index) start_v = *it;
    if ((*it)->get_wcpt().index == sg->get_wcpt_vec().back().index) end_v = *it;
  }
  return std::make_pair(start_v, end_v);
}


// place holder ...
void WCPPID::NeutrinoID::shower_clustering_in_other_clusters(){
  for (auto it = other_clusters.begin(); it != other_clusters.end(); it++){
    WCPPID::PR3DCluster *cluster = *it;
    if (used_shower_clusters.find(cluster->get_cluster_id()) != used_shower_clusters.end()) continue;
    
    WCPPID::ProtoSegment *sg = 0;
    for (auto it1 = map_segment_vertices.begin(); it1 != map_segment_vertices.end(); it1++){
      if (it1->first->get_cluster_id() != cluster->get_cluster_id()) continue;
      sg = it1->first;
      break;
    }

    if (sg != 0){
      WCPPID::WCShower *shower = new WCPPID::WCShower();
      shower->set_start_vertex(main_vertex, 3);
      shower->set_start_segment(sg);
      
      if (sg->get_flag_dir()==0){
	auto tmp_vertices = get_start_end_vertices(sg);
	if (map_vertex_segments[tmp_vertices.first].size()==1 && map_vertex_segments[tmp_vertices.second].size()>1){
	  sg->set_flag_dir(1);
	}else if (map_vertex_segments[tmp_vertices.first].size()>1 && map_vertex_segments[tmp_vertices.second].size()==1){
	  sg->set_flag_dir(-1);
	}else{
	  // examine vertices 
	  double dis1 = sqrt(pow(main_vertex->get_fit_pt().x - sg->get_point_vec().front().x,2) + pow(main_vertex->get_fit_pt().y - sg->get_point_vec().front().y,2) + pow(main_vertex->get_fit_pt().z - sg->get_point_vec().front().z,2));
	  double dis2 = sqrt(pow(main_vertex->get_fit_pt().x - sg->get_point_vec().back().x,2) + pow(main_vertex->get_fit_pt().y - sg->get_point_vec().back().y,2) + pow(main_vertex->get_fit_pt().z - sg->get_point_vec().back().z,2));
	    
	  if (dis1 < dis2){
	    sg->set_flag_dir(1);
	  }else{
	    sg->set_flag_dir(-1);
	  }
	}
      }
      //      std::cout << sg->get_particle_type() << " " << std::endl;
      if (sg->get_particle_type()==0 || fabs(sg->get_particle_type())==13){
	sg->set_particle_type(11);
	TPCParams& mp = Singleton<TPCParams>::Instance();
	sg->set_particle_mass(mp.get_mass_electron());
	sg->cal_4mom();
      }

      
      std::set<WCPPID::ProtoSegment* > used_segments; 
      shower->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, used_segments); 
      showers.push_back(shower);
    }
    
  }
  update_shower_maps();  
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

  update_shower_maps();  
  //  std::cout << showers.size() << " " << used_segments.size() << std::endl;
}
