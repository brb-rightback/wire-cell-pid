std::pair<bool, double> WCPPID::NeutrinoID::numu_tagger(){
  bool flag_long_muon = false;

  bool flag_numu_cc = false;
  
  double dis_cut = 18*units::cm; // muon shorter than this is not useful ...
  // for electron, we cut on 60 MeV anyway ...
  double max_muon_length = 0;
  
  // first round of check of everything connect to the vertex ...
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    double length = sg->get_length();
    double direct_length = sg->get_direct_length();
    double medium_dQ_dx = sg->get_medium_dQ_dx();

    /* if (sg->get_particle_type()==13){ */
    //    std::cout << sg->get_particle_type() << " Xin " << length/units::cm << " " << medium_dQ_dx/(43e3/units::cm) << " " << direct_length/units::cm << std::endl; 
    /* } */
    
    double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);//0.85+0.95 *sqrt(25./ (length / units::cm));
    if (abs(sg->get_particle_type())==13 && length > dis_cut && medium_dQ_dx < dQ_dx_cut * 43e3/units::cm && (length > 40*units::cm || length <= 40*units::cm && direct_length > 0.925 * length)){
      flag_numu_cc = true;
      
      if (length > max_muon_length) max_muon_length = length;
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
	
	if (length > max_muon_length) max_muon_length = length;
	//	std::cout << main_cluster->get_cluster_id() << " " << shower->get_start_segment()->get_cluster_id() << " " << length/units::cm << " " << total_length/units::cm << std::endl;
      }
    }
    // }

    

    //  if (!flag_numu_cc){
    double max_length = 0;
    double max_direct_length = 0;
    double max_medium_dQ_dx = 0;
    double acc_track_length = 0;
    for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      if (sg->get_cluster_id() != main_vertex->get_cluster_id()) continue;
      if (abs(sg->get_particle_type())==211) continue;
      
      double length = sg->get_length();
      if (!sg->get_flag_shower()){
	acc_track_length += length;
      }
      double direct_length = sg->get_direct_length();
      double medium_dQ_dx = sg->get_medium_dQ_dx();
      double dQ_dx_cut = 0.8866+0.9533 *pow(18*units::cm/length, 0.4234);
      
      // std::cout << sg->get_id() << " " << sg->get_flag_shower() << " " << sg->get_flag_shower_topology() << " " << medium_dQ_dx/(43e3/units::cm) << " " << dQ_dx_cut << " " << length/units::cm << " " << direct_length/units::cm << std::endl;
      
      if (sg->get_flag_shower() && (!sg->get_flag_shower_topology()) || (!sg->get_flag_shower()) || length > 50*units::cm){
	
	if (medium_dQ_dx < dQ_dx_cut * 1.05 * 43e3/units::cm  && medium_dQ_dx > 0.75 * 43e3/units::cm && direct_length > 0.925*length ){
	  if (length > max_length){
	    max_length = length;
	    max_direct_length = direct_length;
	    max_medium_dQ_dx = medium_dQ_dx/(43e3/units::cm);
	  }
	  
	  // calculate the length of muon 
	  if (length > 50*units::cm){
	    auto pair_vertices = find_vertices(sg);
	    double tmp_length = length;
	    auto results1 = find_cont_muon_segment(sg, pair_vertices.first);
	    if (results1.first != 0){
	      tmp_length += results1.first->get_length();
	    }
	    auto results2 = find_cont_muon_segment(sg, pair_vertices.second);
	    if (results2.first != 0){
	      tmp_length += results2.first->get_length();
	    }
	    if (tmp_length > max_muon_length) max_muon_length = tmp_length;
	  }
	}
      }
    }
    if (max_length > 25*units::cm && acc_track_length >0 ||
	max_length > 30*units::cm && acc_track_length == 0){

      if (max_length > max_muon_length) max_muon_length = max_length;
      
      flag_numu_cc = true;
    }

    /* if (!flag_numu_cc) */
    /*   std::cout <<"Xin " << max_length/units::cm << " " << max_direct_length/units::cm << " " << max_medium_dQ_dx << " " << acc_track_length/units::cm << std::endl; */
    //  }

    // std::cout << "Xin1: " << max_muon_length/units::cm << std::endl;

  if (max_muon_length > 100*units::cm) flag_long_muon = true;
  
  if (flag_numu_cc){
    neutrino_type |= 1UL << 2; //numu
  }else{
    neutrino_type |= 1UL << 3; // nc ...
  }
  
  
  return std::make_pair(flag_long_muon, max_muon_length);
}
