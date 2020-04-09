
struct cluster_point_info{
  WCPPID::PR3DCluster *cluster;
  double min_angle;
  double min_dis;
  Point min_point;
  WCPPID::ProtoVertex *min_vertex;
};

bool sortbydis(const cluster_point_info &a, const cluster_point_info &b){
  return (a.min_dis < b.min_dis);
}



void WCPPID::NeutrinoID::shower_clustering(){

  // connect to the main cluster ...
  shower_clustering_in_main_cluster();
  shower_clustering_from_main_cluster();
  shower_clustering_from_vertices();
  
  calculate_shower_kinematics();
  // check remaining clusters ...
  shower_clustering_in_other_clusters();
  
  calculate_shower_kinematics();
}

void WCPPID::NeutrinoID::shower_clustering_from_main_cluster(){
  std::map<WCPPID::WCShower*, TVector3> map_shower_dir;
  
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *seg = it->first;
    if (seg->get_cluster_id() != main_cluster->get_cluster_id()) continue;
    if (!seg->get_flag_shower_topology()) continue;
    WCPPID::WCShower *shower = map_segment_in_shower[seg];
    
    /* WCPPID::ProtoVertex *start_v=0, *end_v=0; */
    /* for (auto it = map_segment_vertices[seg].begin(); it!=map_segment_vertices[seg].end(); it++){ */
    /*   if ((*it)->get_wcpt().index == seg->get_wcpt_vec().front().index) start_v = *it; */
    /*   if ((*it)->get_wcpt().index == seg->get_wcpt_vec().back().index) end_v = *it; */
    /* } */
    /* if (seg->get_flag_dir()==1){ */
    /*   shower->set_start_vertex(start_v, 1); */
    /* }else if (seg->get_flag_dir()==-1){ */
    /*   shower->set_start_vertex(end_v, 1); */
    /* } */
    TVector3 dir_shower = seg->cal_dir_3vector(shower->get_start_vertex().first->get_fit_pt(), 15*units::cm);

    map_shower_dir[shower] = dir_shower;
  }
  
  
  // examine other segments first ... 
  for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){ 
    WCPPID::ProtoSegment *seg1 = it->first; 
    if (seg1->get_cluster_id() == main_cluster->get_cluster_id()) continue; 
    if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;

    // find the closest point 
    // judge the detection ...
    double min_dis = 1e9;
    WCPPID::WCShower* min_shower = 0;
    for (auto it1 = map_shower_dir.begin(); it1 != map_shower_dir.end(); it1++){
      WCPPID::WCShower *shower = it1->first;
      std::pair<double, WCP::Point> pair_dis_point = seg1->get_closest_point(shower->get_start_vertex().first->get_fit_pt()); 
      TVector3 v1(pair_dis_point.second.x - shower->get_start_vertex().first->get_fit_pt().x, pair_dis_point.second.y - shower->get_start_vertex().first->get_fit_pt().y, pair_dis_point.second.z - shower->get_start_vertex().first->get_fit_pt().z);
      double angle = it1->second.Angle(v1);
      if (angle/3.1415926*180.<15 && pair_dis_point.first < 80*units::cm ||
	  angle/3.1415926*180 < 7.5 && pair_dis_point.first < 120*units::cm ){
	double dis = pow(pair_dis_point.first*cos(angle),2)/pow(40*units::cm,2) + pow(pair_dis_point.first*sin(angle),2)/pow(5*units::cm,2);
	if (dis < min_dis){
	  min_dis = dis;
	  min_shower = shower;
	}
      }else{
	//std::cout << seg1->get_cluster_id()*1000 + seg1->get_id() << " " << pair_dis_point.first/units::cm << " " << angle/3.1415926*180. << " " << seg1->get_length()/units::cm << std::endl;
      }
    }
    if (min_shower!=0){
      min_shower->add_segment(seg1, map_segment_vertices);
    }
  }
  
  update_shower_maps();
}

/* void WCPPID::NeutrinoID::establish_cluster_segment_maps(){ */
/*   map_cluster_segments.clear(); */
/*   map_segment_cluster.clear(); */

/*   std::map<int, WCPPID::PR3DCluster*> map_cid_cluster; */
/*   std::map<WCPPID::ProtoSegment*, int> map_segment_cid; */

/*   map_cid_cluster[main_cluster->get_cluster_id()] = main_cluster; */
/*   for (size_t i=0; i!=other_clusters.size(); i++){ */
/*     map_cid_cluster[other_clusters.at(i)->get_cluster_id()] = other_clusters.at(i); */
/*   } */

/*   for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){ */
/*     WCPPID::ProtoSegment* seg = it->first; */
/*     map_segment_cid[seg] = seg->get_cluster_id(); */
/*   } */

  
  
/* } */

void WCPPID::NeutrinoID::shower_clustering_from_vertices(){
  // first map the sg vs. cluster existing ...
  //  std::cout << map_cluster_segments.size() << " " << other_clusters.size() << " " << map_segment_cluster.size() << " " << map_segment_vertices.size() << std::endl;
  std::map<PR3DCluster*, std::pair<WCP::Point, double> > map_cluster_center_point;
  
  for (auto it = other_clusters.begin(); it!=other_clusters.end(); it++){
    WCPPID::PR3DCluster *cluster = *it;
    auto it1 = map_cluster_segments.find(cluster);
    if (it1  == map_cluster_segments.end()) continue;
    double acc_length = 0;
    Point p(0,0,0);
    Int_t np = 0;
    for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
      WCPPID::ProtoSegment *seg = *it2;
      if (map_segment_in_shower.find(seg) != map_segment_in_shower.end()) continue;
      if (seg->get_flag_shower()==1){
	acc_length += seg->get_length();
	PointVector& pts = seg->get_point_vec();
	for (size_t i=0; i!= pts.size(); i++){
	  p.x += pts.at(i).x;
	  p.y += pts.at(i).y;
	  p.z += pts.at(i).z;
	  np ++;
	}
      }
    }
    if (acc_length > 1.0*units::cm){
      p.x /= np;
      p.y /= np;
      p.z /= np;
      map_cluster_center_point[cluster] = std::make_pair(p, acc_length);
    }
  }

  // list the main vertices ...
  std::vector<WCPPID::ProtoVertex*> main_cluster_vertices;
  for (auto it = map_vertex_segments.begin(); it != map_vertex_segments.end(); it++){
    if (it->first->get_cluster_id()==main_cluster->get_cluster_id()){
      //      if (it->first == main_vertex) std::cout << main_cluster_vertices.size() << std::endl;
      main_cluster_vertices.push_back(it->first);
    }
  }

 
  std::map<WCPPID::PR3DCluster*, cluster_point_info > map_cluster_pi;
  
  // start to check with each of clusters in main cluster ...
  for (auto it = map_cluster_center_point.begin(); it!= map_cluster_center_point.end(); it++){
    WCPPID::PR3DCluster *cluster = it->first;
    Point center_p = it->second.first;
    PointVector total_pts;
    for (auto it1 = map_cluster_segments[cluster].begin(); it1!= map_cluster_segments[cluster].end(); it1++){
      PointVector& pts = (*it1)->get_point_vec();
      total_pts.insert(total_pts.end(), pts.begin(), pts.end());
    }
    WCP::ToyPointCloud pcloud;
    pcloud.AddPoints(total_pts);
    pcloud.build_kdtree_index();

    cluster_point_info min_pi;
    min_pi.cluster = cluster;
    min_pi.min_angle = 1e9;
    min_pi.min_dis = 1e9;
    min_pi.min_vertex = 0;
    cluster_point_info main_pi;
    main_pi.cluster = cluster;
    main_pi.min_vertex = main_vertex;
    
    for (size_t i=0;i!=main_cluster_vertices.size();i++){
      std::pair<double, WCP::Point> result = pcloud.get_closest_point(main_cluster_vertices.at(i)->get_fit_pt());
      TVector3 v1(result.second.x - main_cluster_vertices.at(i)->get_fit_pt().x,
		  result.second.y - main_cluster_vertices.at(i)->get_fit_pt().y,
		  result.second.z - main_cluster_vertices.at(i)->get_fit_pt().z);
      TVector3 v2(center_p.x - result.second.x,
		  center_p.y - result.second.y,
		  center_p.z - result.second.z);
      double angle = v1.Angle(v2)/3.1415926*180.;
      if (angle < min_pi.min_angle){
	min_pi.min_angle = angle;
	min_pi.min_dis = result.first;
	min_pi.min_vertex = main_cluster_vertices.at(i);
	min_pi.min_point = result.second;
      }
      if (main_cluster_vertices.at(i) == main_vertex) {
	main_pi.min_angle = angle;
	main_pi.min_dis = result.first;
	main_pi.min_point = result.second;
      }
      //  std::cout << i << " " << v1.Angle(v2)/3.1415926*180. << " " << result.first/units::cm << std::endl;
    }
    if (main_pi.min_angle < min_pi.min_angle + 3 && min_pi.min_angle > 0.9 * main_pi.min_angle && main_pi.min_dis < min_pi.min_dis * 1.2){
      map_cluster_pi[cluster] = main_pi;
      //
    }else{
      map_cluster_pi[cluster] = min_pi;
    }
  }

  std::vector<cluster_point_info > vec_pi;
  for (auto it = map_cluster_pi.begin(); it!=map_cluster_pi.end(); it++){
    vec_pi.push_back(it->second);
  }
  std::sort(vec_pi.begin(), vec_pi.end(), sortbydis);
  //for (size_t i=0;i!=vec_pi.size();i++){
  // std::cout << vec_pi.at(i).min_dis << std::endl;
  //}

  for (size_t i=0;i!=vec_pi.size();i++){
    // find the segment and/or create the new vertices ...
    // create the new shower ...
    // cluster the other things ...
    
    // udate the map
    update_shower_maps();
  }

  
  
}


void WCPPID::NeutrinoID::calculate_shower_kinematics(){
  for (size_t i=0;i!=showers.size();i++){
    WCPPID::WCShower *shower = showers.at(i);
    if (!shower->get_flag_kinematics()){
      shower->calculate_kinematics();
      double kine_charge = cal_kine_charge(shower);
      shower->set_kine_charge(kine_charge);
      shower->set_flag_kinematics(true);
    }
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
  // vertices in main cluster as well as existing showers ...
  std::vector<WCPPID::ProtoVertex* > vertices;
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() == main_cluster->get_cluster_id() || map_vertex_in_shower.find(vtx) != map_vertex_in_shower.end()){
      vertices.push_back(vtx);
    }
  }
  //  std::cout << vertices.size() << std::endl;
  
  for (auto it = other_clusters.begin(); it != other_clusters.end(); it++){
    WCPPID::PR3DCluster *cluster = *it;
    if (used_shower_clusters.find(cluster->get_cluster_id()) != used_shower_clusters.end()) continue;

    // check against the main cluster first ...
    ToyPointCloud *pcloud = cluster->get_point_cloud();
    double min_dis = 1e9;
    WCPPID::ProtoVertex *min_vertex = 0;
    for (size_t i=0;i!=vertices.size();i++){
      double dis = pcloud->get_closest_dis(vertices.at(i)->get_fit_pt());
      if (dis < min_dis){
	min_dis = dis;
	min_vertex = vertices.at(i);
      }
    }
    
    
    WCPPID::ProtoSegment *sg = 0;
    for (auto it1 = map_segment_vertices.begin(); it1 != map_segment_vertices.end(); it1++){
      if (it1->first->get_cluster_id() != cluster->get_cluster_id()) continue;
      sg = it1->first;
      break;
    }

    int connection_type = 2;
    if (min_dis > 80*units::cm){
      connection_type = 3;
    }
    
    
    if (sg != 0){
      WCPPID::WCShower *shower = new WCPPID::WCShower();
      shower->set_start_vertex(min_vertex, connection_type);
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
