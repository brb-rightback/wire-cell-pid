std::pair<bool, double> WCPPID::NeutrinoID::numu_tagger(){
  bool flag_long_muon = false;

  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);
  
  bool flag_numu_cc = false;

  bool flag_print = false;
  
  double dis_cut = 5*units::cm; // muon shorter than this is not useful ...
  // for electron, we cut on 60 MeV anyway ...
  double max_muon_length = 0;
  WCPPID::ProtoSegment* max_muon = 0;
  WCPPID::WCShower* max_long_muon = 0;
  
  // first round of check of everything connect to the vertex ...
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    double length = sg->get_length();
    double direct_length = sg->get_direct_length();
    double medium_dQ_dx = sg->get_medium_dQ_dx();

    /* if (sg->get_particle_type()==13){ */
    
    /* } */
    
    double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);//0.85+0.95 *sqrt(25./ (length / units::cm));

    // std::cout << sg->get_particle_type() << " Xin1 " << length/units::cm << " " << medium_dQ_dx/(43e3/units::cm) << " " << direct_length/units::cm << " " << dQ_dx_cut << std::endl; 
    
    if (abs(sg->get_particle_type())==13 && length > dis_cut && medium_dQ_dx < dQ_dx_cut * 43e3/units::cm && (length > 40*units::cm || length <= 40*units::cm && direct_length > 0.925 * length)){
      flag_numu_cc = true;
      
      if (length > max_muon_length) {
	max_muon_length = length;
	max_muon = sg;
      }
    }
  }



  //  if (!flag_numu_cc){
    // check long muon  in WCshowers ...
  for (auto it = showers.begin(); it != showers.end(); it++){
    WCPPID::WCShower *shower = *it;
    if (abs(shower->get_start_segment()->get_particle_type())==13 && main_cluster->get_cluster_id() == shower->get_start_segment()->get_cluster_id()){
      double length = shower->get_total_track_length();
      double total_length = shower->get_total_length();
      if (length > 18*units::cm) flag_numu_cc = true;
      
      if (length > max_muon_length) {
	max_muon_length = length;
	max_long_muon = shower;
	max_muon = 0;
      }
      //	std::cout << main_cluster->get_cluster_id() << " " << shower->get_start_segment()->get_cluster_id() << " " << length/units::cm << " " << total_length/units::cm << std::endl;
    }
  }
  // }
  
  
  
  //  if (!flag_numu_cc){
  double max_length = 0;
  double max_direct_length = 0;
  double max_medium_dQ_dx = 0;
  double acc_track_length = 0;
  double max_length_all = 0;
  WCPPID::ProtoSegment *tmp_max_muon = 0;
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    double length = sg->get_length();

    if (sg->get_cluster_id() != main_vertex->get_cluster_id()) continue;
    if ((!sg->get_flag_shower()) && max_length_all < length) max_length_all = length;
    if (abs(sg->get_particle_type())==211) continue;
    

    if (!sg->get_flag_shower()){
      acc_track_length += length;
    }
    if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end()) continue;
    
    
    double direct_length = sg->get_direct_length();
    double medium_dQ_dx = sg->get_medium_dQ_dx();
    double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
    
    //    std::cout << sg->get_id() << " " << sg->get_flag_shower() << " " << sg->get_flag_shower_topology() << " " << medium_dQ_dx/(43e3/units::cm) << " " << dQ_dx_cut << " " << length/units::cm << " " << direct_length/units::cm << std::endl;

    if (sg->get_flag_avoid_muon_check()) continue;
    if (sg->get_flag_shower() && (!sg->get_flag_shower_topology()) || (!sg->get_flag_shower()) || length > 50*units::cm){
      
      if (medium_dQ_dx < dQ_dx_cut * 1.05 * 43e3/units::cm  && medium_dQ_dx > 0.75 * 43e3/units::cm && (direct_length > 0.925*length || length > 50*units::cm)){
	if (length > max_length){
	  max_length = length;
	  max_direct_length = direct_length;
	  max_medium_dQ_dx = medium_dQ_dx/(43e3/units::cm);
	  tmp_max_muon  = sg;
	}
	
	// calculate the length of muon 
	if (length > 50*units::cm){
	  double tmp_length = length;
	  if (!sg->get_flag_shower_topology()){
	    auto pair_vertices = find_vertices(sg);
	    auto results1 = find_cont_muon_segment(sg, pair_vertices.first);
	    if (results1.first != 0){
	      tmp_length += results1.first->get_length();
	    }
	    auto results2 = find_cont_muon_segment(sg, pair_vertices.second);
	    if (results2.first != 0){
	      tmp_length += results2.first->get_length();
	    }
	  }
	  //std::cout << sg->get_particle_type() << " " << sg->get_flag_shower_topology() << " " << sg->get_id() << " " << sg->get_length()/units::cm << std::endl;
	  
	  if (tmp_length > max_muon_length) {
	    max_muon_length = tmp_length;
	    max_muon = sg;
	    max_long_muon = 0;
	  }
	}
      }
    }
  }
  if (max_length > 25*units::cm && acc_track_length >0 ||
      max_length > 30*units::cm && acc_track_length == 0){
    
    if (max_length > max_muon_length) {
      max_muon_length = max_length;
      max_muon = tmp_max_muon;
      max_long_muon = 0;
    }
    if (flag_print) std::cout << "Xin_A: " << std::endl;
    flag_numu_cc = true;
  }


  
  //  if (!flag_numu_cc) 
  //  std::cout <<"Xin " << max_length/units::cm << " " << max_direct_length/units::cm << " " << max_medium_dQ_dx << " " << acc_track_length/units::cm << std::endl; 

  if (max_muon !=0 || max_long_muon !=0){
    TVector3 muon_dir;
    int n_daughter_tracks = 0;
    int n_daughter_all = 0;
    
    if (max_muon != 0){
      auto pair_vertices = find_vertices(max_muon);
      double dis1 = sqrt(pow(pair_vertices.first->get_fit_pt().x - main_vertex->get_fit_pt().x,2) + pow(pair_vertices.first->get_fit_pt().y - main_vertex->get_fit_pt().y,2) + pow(pair_vertices.first->get_fit_pt().z - main_vertex->get_fit_pt().z,2));
      double dis2 = sqrt(pow(pair_vertices.second->get_fit_pt().x - main_vertex->get_fit_pt().x,2) + pow(pair_vertices.second->get_fit_pt().y - main_vertex->get_fit_pt().y,2) + pow(pair_vertices.second->get_fit_pt().z - main_vertex->get_fit_pt().z,2));
      if (dis1 < dis2){
	muon_dir.SetXYZ(pair_vertices.second->get_fit_pt().x - pair_vertices.first->get_fit_pt().x,
			pair_vertices.second->get_fit_pt().y - pair_vertices.first->get_fit_pt().y,
			pair_vertices.second->get_fit_pt().z - pair_vertices.first->get_fit_pt().z );
	n_daughter_tracks = calculate_num_daughter_tracks(pair_vertices.first, max_muon, false, 3*units::cm).first;
	n_daughter_all = calculate_num_daughter_tracks(pair_vertices.first, max_muon, true, 3*units::cm).first;
      }else{
	muon_dir.SetXYZ(pair_vertices.first->get_fit_pt().x - pair_vertices.second->get_fit_pt().x,
			pair_vertices.first->get_fit_pt().y - pair_vertices.second->get_fit_pt().y,
			pair_vertices.first->get_fit_pt().z - pair_vertices.second->get_fit_pt().z );
	n_daughter_tracks = calculate_num_daughter_tracks(pair_vertices.second, max_muon, false, 3*units::cm).first;
	n_daughter_all = calculate_num_daughter_tracks(pair_vertices.second, max_muon, true, 3*units::cm).first;
      }
      if (!max_muon->get_flag_shower()){
	n_daughter_tracks -= 1;
      }
      n_daughter_all -=1;
    }else if (max_long_muon !=0){
      muon_dir = max_long_muon->cal_dir_3vector(main_vertex->get_fit_pt(), 30*units::cm);
      auto pair_result = max_long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon);
      WCPPID::ProtoSegment *last_sg = pair_result.first;
      auto pair_vertices = find_vertices(last_sg);
      double dis1 = sqrt(pow(pair_vertices.first->get_fit_pt().x - main_vertex->get_fit_pt().x,2) + pow(pair_vertices.first->get_fit_pt().y - main_vertex->get_fit_pt().y,2) + pow(pair_vertices.first->get_fit_pt().z - main_vertex->get_fit_pt().z,2));
      double dis2 = sqrt(pow(pair_vertices.second->get_fit_pt().x - main_vertex->get_fit_pt().x,2) + pow(pair_vertices.second->get_fit_pt().y - main_vertex->get_fit_pt().y,2) + pow(pair_vertices.second->get_fit_pt().z - main_vertex->get_fit_pt().z,2));
      if (dis1 < dis2){
	n_daughter_tracks = calculate_num_daughter_tracks(pair_vertices.first, last_sg, false, 3*units::cm).first;
	n_daughter_all = calculate_num_daughter_tracks(pair_vertices.first, last_sg, true, 3*units::cm).first;
      }else{
	n_daughter_tracks = calculate_num_daughter_tracks(pair_vertices.second, last_sg, false, 3*units::cm).first;
	n_daughter_all = calculate_num_daughter_tracks(pair_vertices.second, last_sg, true, 3*units::cm).first;
      }
      //      std::cout << "last segment: " << last_sg->get_id() << std::endl;
      if (!last_sg->get_flag_shower()){
	n_daughter_tracks -= 1;
      }
      n_daughter_all -=1;
    }

    

    if (flag_numu_cc){
      double angle = muon_dir.Angle(dir_beam)/3.1415926*180.;

      //      std::cout <<  "Xin: " << max_length_all/units::cm << " " << max_muon_length/units::cm << " " << muon_dir.Angle(dir_beam)/3.1415926*180. << " " << n_daughter_tracks << " " << n_daughter_all - n_daughter_tracks << " " << max_muon_length/(80*units::cm) + angle/46 << " " << " " << (angle-120)/40. - max_muon_length/(80.*units::cm) << " " << (angle-120)/120 - (max_length_all-max_muon_length-120*units::cm)/(100.*units::cm)<< std::endl;
      
      if (n_daughter_tracks >1 || n_daughter_all - n_daughter_tracks > 2
	  || (max_length_all-max_muon_length > 100.*units::cm)){
	flag_numu_cc = false;
	if (flag_print) std::cout << "Xin_B: " << std::endl;
      }
    }
  }
  
  // std::cout << "Xin1: " << max_muon_length/units::cm << std::endl;
  
  if ((max_muon_length > 100*units::cm || max_length_all > 120*units::cm) && flag_numu_cc) flag_long_muon = true;
  
  if (flag_numu_cc){
    neutrino_type |= 1UL << 2; //numu
  }else{
    neutrino_type |= 1UL << 3; // nc ...
  }
  
  
  return std::make_pair(flag_long_muon, max_muon_length);
}
