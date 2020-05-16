bool WCPPID::NeutrinoID::nue_tagger(){
  bool flag_nue = false;

  // check main_vertex ...
  {
    auto it = map_vertex_to_shower.find(main_vertex);
    if (it != map_vertex_to_shower.end()){
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	WCPPID::WCShower *shower = *it1;
	WCPPID::ProtoSegment *sg = shower->get_start_segment();

	if (bad_reconstruction(shower)) continue;
	
	// higher than 70 MeV, and connect to the vertex ...
	if (map_vertex_segments[main_vertex].find(sg) != map_vertex_segments[main_vertex].end() && shower->get_kine_charge() > 70*units::MeV){
	  //	  std::cout << shower->get_kine_charge()/units::MeV << std::endl;
	  if (!gap_identification(main_vertex, sg)){ // gap id
	    int mip_id = mip_identification(main_vertex,sg);
	    if (mip_id == 1){
	      flag_nue = true;
	    }else if (mip_id == 0){
	      if (!pi0_identification(main_vertex, sg)){
		flag_nue = true;
	      }
	    }
	  } // gap id
	} // 70 MeV, and connected
      } // loop over showers
    } // has a shower ...
  } // main vertex
  

  if (flag_nue){
    neutrino_type |= 1UL << 5; //numu
  }

  return flag_nue;
}

bool WCPPID::NeutrinoID::gap_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment* sg){
  bool flag_gap = false;

  return flag_gap;
}

int WCPPID::NeutrinoID::mip_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg){
  int mip_id = 1; 
  // 1 good, -1 bad, 0 not sure ...
  
  return mip_id;
}

bool WCPPID::NeutrinoID::pi0_identification(WCPPID::ProtoVertex* vertex, WCPPID::ProtoSegment *sg){
  bool flag_pi0 = false;

  return flag_pi0;
}

bool WCPPID::NeutrinoID::bad_reconstruction(WCPPID::WCShower* shower){
  bool flag_bad_shower = false;
  
  return flag_bad_shower;
}
