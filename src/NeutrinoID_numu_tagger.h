bool WCPPID::NeutrinoID::numu_tagger(){
  bool flag_long_muon = false;

  bool flag_numu_cc = false;
  
  double dis_cut = 18*units::cm; // muon shorter than this is not useful ...
  // for electron, we cut on 60 MeV anyway ...
  
  // first round of check of everything connect to the vertex ...
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    double length = sg->get_length();
    double direct_length = sg->get_direct_length();
    double medium_dQ_dx = sg->get_medium_dQ_dx();

    if (sg->get_particle_type()==13){
      std::cout << sg->get_particle_type() << " Xin " << length/units::cm << " " << medium_dQ_dx/(43e3/units::cm) << " " << direct_length/units::cm << std::endl;
    }
    
    double dQ_dx_cut = 0.85+0.95 *sqrt(25./ (length / units::cm));
    if (abs(sg->get_particle_type())==13 && length > dis_cut && medium_dQ_dx < dQ_dx_cut * 43e3/units::cm && (length > 40*units::cm || length <= 40*units::cm && direct_length > 0.925 * length)){
      flag_numu_cc = true;
    }
  }


  if (flag_numu_cc){
    neutrino_type |= 1UL << 2; //numu
  }else{
    neutrino_type |= 1UL << 3; // nc ...
  }
  
  
  return flag_long_muon;
}
