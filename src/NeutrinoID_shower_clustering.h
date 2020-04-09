
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
  // std::cout << showers.size() << std::endl;
  shower_clustering_from_main_cluster();
  //std::cout << showers.size() << std::endl;
  shower_clustering_from_vertices();
  //  std::cout << showers.size() << std::endl;  
  calculate_shower_kinematics();

  // check remaining clusters ...
  shower_clustering_in_other_clusters(true);
 
  calculate_shower_kinematics();

  id_pi0_with_vertex();
}

void WCPPID::NeutrinoID::id_pi0_with_vertex(){
  for (auto it = map_vertex_to_shower.begin(); it!= map_vertex_to_shower.end(); it++){
    std::vector<WCPPID::WCShower*> tmp_showers;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_start_vertex().second <3)
	tmp_showers.push_back(*it1);
    }
    
    if (tmp_showers.size()>1){
      std::map<std::pair<WCPPID::WCShower*, WCPPID::WCShower*>, double> map_shower_pair_mass;
      for (size_t i=0;i!= tmp_showers.size();i++){
	WCPPID::WCShower *shower_1 = tmp_showers.at(i);
	TVector3 dir1 = shower_1->get_init_dir();
	for (size_t j=i+1; j<tmp_showers.size();j++){
	  WCPPID::WCShower *shower_2 = tmp_showers.at(j);
	  TVector3 dir2 = shower_2->get_init_dir();
	  double angle = dir1.Angle(dir2);
	  double mass_pio = sqrt(4*shower_1->get_kine_charge()* shower_2->get_kine_charge()*pow(sin(angle/2.),2));
	  map_shower_pair_mass[std::make_pair(shower_1, shower_2)] = mass_pio;
	  //  std::cout << it->first << " " << i << " " << j << " " << shower_1->get_kine_charge()/units::MeV << " " << shower_2->get_kine_charge()/units::MeV << " " << dir1.Mag() << " " << dir2.Mag() << " " << angle/3.1415926*180. << " " << sqrt(4*shower_1->get_kine_charge()* shower_2->get_kine_charge()*pow(sin(angle/2.),2))/units::MeV<< std::endl;
	}
      }
      while(map_shower_pair_mass.size()>0){
	// find the one close to the pi0 mass ...
	double mass_diff = 1e9;
	double mass_save = 0;
	WCPPID::WCShower *shower_1 = 0;
	WCPPID::WCShower *shower_2 = 0;
	for (auto it = map_shower_pair_mass.begin(); it!= map_shower_pair_mass.end(); it++){
	  if (fabs(it->second - 135*units::MeV) < mass_diff){
	    mass_diff = fabs(it->second - 135*units::MeV);
	    mass_save = it->second;
	    shower_1 = it->first.first;
	    shower_2 = it->first.second;
	  }
	}

	if (mass_diff < 35*units::MeV){
	  pi0_showers.insert(shower_1);
	  pi0_showers.insert(shower_2);
	  int pio_id = acc_segment_id; acc_segment_id ++;
	  map_shower_pio_id[shower_1] = pio_id;
	  map_shower_pio_id[shower_2] = pio_id;
	  map_pio_id_mass[pio_id] = mass_save;
	  map_pio_id_showers[pio_id].push_back(shower_1);
	  map_pio_id_showers[pio_id].push_back(shower_2);
	  std::cout << "Pi0 found with mass: " << " " << mass_save/units::MeV << " MeV" << std::endl;
	}else{
	  break;
	}
	std::vector<std::pair<WCPPID::WCShower*, WCPPID::WCShower*> > to_be_removed;
	for (auto it = map_shower_pair_mass.begin(); it!= map_shower_pair_mass.end(); it++){
	  if (it->first.first == shower_1 || it->first.first == shower_2 || it->first.second == shower_1 || it->first.second == shower_2)
	    to_be_removed.push_back(it->first);
	}
	for (auto it = to_be_removed.begin(); it != to_be_removed.end(); it++){
	  map_shower_pair_mass.erase(*it);
	}
      }
      


      
    }
  }
}

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
      Point near_center = pcloud.get_center_point_radius(result.second, 2*units::cm);
      TVector3 v3(near_center.x - result.second.x,
		  near_center.y - result.second.y,
		  near_center.z - result.second.z);
      double angle = v1.Angle(v2)/3.1415926*180.;
      if (angle < 30)
	angle = std::min(v1.Angle(v2)/3.1415926*180., v1.Angle(v3)/3.1415926*180.);
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
      // std::cout << i << " " << v1.Angle(v2)/3.1415926*180. << " " << v1.Angle(v3)/3.1415926*180. << " " << result.first/units::cm << std::endl;
    }
    //    std::cout << main_pi.min_angle << " " << min_pi.min_angle << " " << " " << main_pi.min_dis/units::cm << " " << min_pi.min_dis/units::cm << std::endl;
    
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

  std::map<WCPPID::PR3DCluster*, WCPPID::ProtoVertex* > map_cluster_associated_vertex;
  for (size_t i=0;i!=vec_pi.size();i++){
    if (vec_pi.at(i).min_angle < 10)
      map_cluster_associated_vertex[vec_pi.at(i).cluster] = vec_pi.at(i).min_vertex;
  }

  /* for (auto it = map_segment_in_shower.begin(); it!= map_segment_in_shower.end(); it++){ */
  /*   std::cout << it->first->get_id() << " "; */
  /* } */
  /* std::cout << std::endl; */


  
  
  for (size_t i=0;i!=vec_pi.size();i++){
    //    if (i>=2) continue;
    
    // find the segment and/or create the new vertices ...
    WCPPID::PR3DCluster *cluster = vec_pi.at(i).cluster;
    WCPPID::ProtoVertex *vertex = vec_pi.at(i).min_vertex;
    WCP::Point point = vec_pi.at(i).min_point;
    WCPPID::ProtoSegment* sg1 = 0;
    double angle = vec_pi.at(i).min_angle;
    if (angle > 25) continue;
    
    for (auto it = map_cluster_segments[cluster].begin(); it!=map_cluster_segments[cluster].end();it++){
      WCPPID::ProtoSegment *sg = *it;
      if (map_segment_in_shower.find(sg) != map_segment_in_shower.end() ) continue;
      std::pair<double, WCP::Point> result = sg->get_closest_point(point);
      if (result.first < 0.01*units::cm){
	sg1 = sg;
	break;
      }
    }
    if (sg1 == 0 ) continue;

    // create the new shower ...
    WCPPID::WCShower *shower = new WCPPID::WCShower();
    //    std::cout << shower << " " << angle << std::endl;
    shower->set_start_vertex(vertex, 2); // second type ...
    showers.push_back(shower); 
    
    
    if (sqrt(pow(sg1->get_point_vec().front().x-point.x,2) + pow(sg1->get_point_vec().front().y-point.y,2) + pow(sg1->get_point_vec().front().z-point.z,2)) < 0.01*units::cm || sqrt(pow(sg1->get_point_vec().back().x-point.x,2) + pow(sg1->get_point_vec().back().y-point.y,2) + pow(sg1->get_point_vec().back().z-point.z,2)) < 0.01*units::cm){

      shower->set_start_segment(sg1, map_segment_vertices);
    }else{
      // break the segment ...
      //      std::cout << sg1 << " " << sg1->get_wcpt_vec().size() << " " << sg1->get_point_vec().size() << std::endl;
      // create two segments with filled information , and find the correct segments ...
      std::tuple<WCPPID::ProtoSegment*, WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> result1 = sg1->break_segment_at_point(point, acc_segment_id, acc_vertex_id);
      if (std::get<1>(result1) == 0){
	shower->set_start_segment(sg1, map_segment_vertices);
      }else{
	 std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> result2 = find_vertices(sg1);
	 TVector3 v3(point.x - vertex->get_fit_pt().x, point.y - vertex->get_fit_pt().y, point.z - vertex->get_fit_pt().z);
	 //std::cout << result2.first->get_wcpt().index << " " << result2.second->get_wcpt().index << " " << std::get<0>(result1)->get_wcpt_vec().front().index << " " << std::get<2>(result1)->get_wcpt_vec().back().index << std::endl;
	 if (result2.first->get_wcpt().index == std::get<0>(result1)->get_wcpt_vec().front().index){
	   add_proto_connection(result2.first, std::get<0>(result1), cluster);
	   add_proto_connection(result2.second, std::get<2>(result1), cluster);
	 }else{
	   add_proto_connection(result2.first, std::get<2>(result1), cluster);
	   add_proto_connection(result2.second, std::get<0>(result1), cluster);
	 }
	 add_proto_connection(std::get<1>(result1), std::get<0>(result1), cluster);
	 add_proto_connection(std::get<1>(result1), std::get<2>(result1), cluster);
	 del_proto_segment(sg1);

	 TVector3 v1 = std::get<0>(result1)->cal_dir_3vector(point, 5*units::cm);
	 TVector3 v2 = std::get<2>(result1)->cal_dir_3vector(point, 5*units::cm);
	
	 //std::cout << v3.Angle(v1) << " " << v3.Angle(v2) << std::endl;
	 if ( v3.Angle(v1) <v3.Angle(v2)) {
	   shower->set_start_segment(std::get<0>(result1), map_segment_vertices);
	 }else{
	   shower->set_start_segment(std::get<2>(result1), map_segment_vertices);
	 }
      }
    }


    
    
    if (shower->get_start_segment()->get_flag_dir()==0){
      // examine vertices 
      double dis1 = sqrt(pow(main_vertex->get_fit_pt().x - shower->get_start_segment()->get_point_vec().front().x,2) + pow(main_vertex->get_fit_pt().y - shower->get_start_segment()->get_point_vec().front().y,2) + pow(main_vertex->get_fit_pt().z - shower->get_start_segment()->get_point_vec().front().z,2));
      double dis2 = sqrt(pow(main_vertex->get_fit_pt().x - shower->get_start_segment()->get_point_vec().back().x,2) + pow(main_vertex->get_fit_pt().y - shower->get_start_segment()->get_point_vec().back().y,2) + pow(main_vertex->get_fit_pt().z - shower->get_start_segment()->get_point_vec().back().z,2));
      
      if (dis1 < dis2){
	shower->get_start_segment()->set_flag_dir(1);
      }else{
	shower->get_start_segment()->set_flag_dir(-1);
      }
    }
    
    //      std::cout << shower->get_start_segment()->get_particle_type() << " " << std::endl;
    if (shower->get_start_segment()->get_particle_type()==0 || fabs(shower->get_start_segment()->get_particle_type())==13){
      shower->get_start_segment()->set_particle_type(11);
      TPCParams& mp = Singleton<TPCParams>::Instance();
      shower->get_start_segment()->set_particle_mass(mp.get_mass_electron());
      shower->get_start_segment()->cal_4mom();
    }
    
    // cluster the other things ...
    
    TVector3 dir_shower = shower->get_start_segment()->cal_dir_3vector(point, 15*units::cm);
    for (auto it = map_segment_vertices.begin(); it!=map_segment_vertices.end(); it++){ 
      WCPPID::ProtoSegment *seg1 = it->first; 
      if (seg1->get_cluster_id() == main_cluster->get_cluster_id()) continue; 
      if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
      if (seg1 == shower->get_start_segment()) continue;
      auto it1 = map_cluster_associated_vertex.find(map_segment_cluster[seg1]);
      if (it1 != map_cluster_associated_vertex.end())
	if (it1->second != vertex) continue;
      // find the closest point 
      // judge the detection ...
      
      std::pair<double, WCP::Point> pair_dis_point = seg1->get_closest_point(shower->get_start_vertex().first->get_fit_pt()); 
      TVector3 v1(pair_dis_point.second.x - shower->get_start_vertex().first->get_fit_pt().x, pair_dis_point.second.y - shower->get_start_vertex().first->get_fit_pt().y, pair_dis_point.second.z - shower->get_start_vertex().first->get_fit_pt().z);
      double angle = dir_shower.Angle(v1);

      //      std::cout << seg1->get_cluster_id()*1000 + seg1->get_id() << " " << angle/3.1415926*180. << " " << pair_dis_point.first/units::cm << std::endl;
      
      if (angle/3.1415926*180.< 25 && pair_dis_point.first < 80*units::cm ||
	  angle/3.1415926*180 < 12.5 && pair_dis_point.first < 120*units::cm ){
	//	double dis = pow(pair_dis_point.first*cos(angle),2)/pow(40*units::cm,2) + pow(pair_dis_point.first*sin(angle),2)/pow(5*units::cm,2);
	//      min_dis = dis;
	shower->add_segment(seg1, map_segment_vertices);
      }
    }
    
    shower->update_particle_type();
    bool tmp_flag = (shower->get_start_vertex().first == main_vertex);
    
    std::cout << "Separated shower: " <<  shower->get_start_segment()->get_cluster_id()*1000 + shower->get_start_segment()->get_id() << " " << shower->get_start_segment()->get_particle_type() << " " << shower->get_num_segments() << " " << tmp_flag << " " << vec_pi.at(i).min_dis/units::cm << std::endl;
    
    // udate the map
    update_shower_maps();
  }

  std::cout << "With separated-cluster shower: " << showers.size() << std::endl;
    
  /* for (auto it = map_segment_in_shower.begin(); it!= map_segment_in_shower.end(); it++){ */
  /*   std::cout << it->first->get_id() << " "; */
  /* } */
  /* std::cout << std::endl; */
}

void WCPPID::NeutrinoID::update_shower_maps(){

  map_vertex_to_shower.clear();
  map_vertex_in_shower.clear();
  map_segment_in_shower.clear();
  used_shower_clusters.clear();
  
  for (auto it = showers.begin(); it!=showers.end(); it++){
    WCPPID::WCShower* shower = *it;
    map_vertex_to_shower[shower->get_start_vertex().first].insert(shower);
    shower->fill_maps(map_vertex_in_shower, map_segment_in_shower);
  }
  for (auto it = map_segment_in_shower.begin(); it!=map_segment_in_shower.end();it++){
    used_shower_clusters.insert(it->first->get_cluster_id());
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


std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> WCPPID::NeutrinoID::get_start_end_vertices(WCPPID::ProtoSegment* sg){
  WCPPID::ProtoVertex *start_v=0, *end_v=0;
  for (auto it = map_segment_vertices[sg].begin(); it!=map_segment_vertices[sg].end(); it++){
    if ((*it)->get_wcpt().index == sg->get_wcpt_vec().front().index) start_v = *it;
    if ((*it)->get_wcpt().index == sg->get_wcpt_vec().back().index) end_v = *it;
  }
  return std::make_pair(start_v, end_v);
}


// place holder ...
void WCPPID::NeutrinoID::shower_clustering_in_other_clusters(bool flag_save){
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
    double main_dis;
    for (size_t i=0;i!=vertices.size();i++){
      double dis = pcloud->get_closest_dis(vertices.at(i)->get_fit_pt());
      if (dis < min_dis){
	min_dis = dis;
	min_vertex = vertices.at(i);
      }
      if (vertices.at(i) == main_vertex)
	main_dis = dis;
    }
    
    // std::cout << min_dis/units::cm << " " << main_dis/units::cm << " " << min_dis/main_dis << std::endl;
    if (min_dis > 0.8 * main_dis){
      min_dis = main_dis;
      min_vertex = main_vertex;
    }
    
    WCPPID::ProtoSegment *sg = 0;
    for (auto it1 = map_segment_vertices.begin(); it1 != map_segment_vertices.end(); it1++){
      if (it1->first->get_cluster_id() != cluster->get_cluster_id()) continue;
      sg = it1->first;
      break;
    }



    int connection_type = 3;
    if (min_dis > 80*units::cm){
      connection_type = 4;
    }

    if (!flag_save)
      connection_type = 4;
    
    
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
	bool tmp_flag = parent_vtx == main_vertex;
	showers.push_back(shower);
	std::cout << "Main-cluster shower " << showers.size() << " : " << curr_sg->get_cluster_id()*1000 + curr_sg->get_id() << " " << curr_sg->get_particle_type() << " " << tmp_flag << " " << curr_sg->get_flag_shower_topology() << std::endl;
	

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
      if (angle/3.1415926*180.<25. && pair_dis_point.first < 80*units::cm ||
	  angle/3.1415926*180 < 12.5 && pair_dis_point.first < 120*units::cm ){
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
