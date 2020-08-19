void WCPPID::NeutrinoID::fill_kine_tree(WCPPID::KineInfo& ktree){
  // neutrino vertex distance to wall ...
  Point nu_vtx = main_vertex->get_fit_pt();

  TPCParams& mp = Singleton<TPCParams>::Instance();
  Point corr_nu_vtx = mp.func_pos_SCE_correction(nu_vtx);
  ktree.kine_nu_x_corr = corr_nu_vtx.x/units::cm;
  ktree.kine_nu_y_corr = corr_nu_vtx.y/units::cm;
  ktree.kine_nu_z_corr = corr_nu_vtx.z/units::cm;

  //  std::cout << nu_vtx << " " << corr_nu_vtx << std::endl;
  
  // neutrino energy reconstruction ...
  ktree.kine_reco_Enu = 0;
  ktree.kine_reco_add_energy = 0;


  double ave_binding_energy = 8.6 * units::MeV;
  
  std::set<WCPPID::ProtoVertex* > used_vertices;
  std::set<WCPPID::ProtoSegment* > used_segments;
  std::set<WCPPID::WCShower* > used_showers;

  // put all shower related ones in ...
  for (auto it = showers.begin(); it!=showers.end();it++){
    (*it)->fill_sets(used_vertices, used_segments, false);
  }
  // map of sg w.r.t. shower
  std::map<WCPPID::ProtoSegment*, WCPPID::WCShower*> map_sg_shower;
  for (auto it = showers.begin(); it!= showers.end(); it++){
    WCPPID::WCShower *shower = *it;
    WCPPID::ProtoSegment* curr_sg = shower->get_start_segment();
    map_sg_shower[curr_sg] = shower;
  }

  
  std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > segments_to_be_examined;
  for (auto it = map_vertex_segments[main_vertex].begin(); it != map_vertex_segments[main_vertex].end(); it++){
    // parent are all zero now ...
    auto it1 = map_sg_shower.find(*it);

    if (it1 != map_sg_shower.end()){
      // shower ...
      WCPPID::WCShower* shower = it1->second;

      float kine_best = shower->get_kine_best();
      if (kine_best ==0 ) kine_best = shower->get_kine_charge();
      float kine_charge = shower->get_kine_charge();
      float kine_dQdx = shower->get_kine_dQdx();
      float kine_range = shower->get_kine_range();
      
      ktree.kine_energy_particle.push_back(kine_best/units::MeV);
      ktree.kine_particle_type.push_back(shower->get_start_segment()->get_particle_type());
      
      if (fabs(kine_best - kine_charge) < 0.001 * kine_best){
	ktree.kine_energy_info.push_back(2); // charge ...
      }else if (fabs(kine_best - kine_range) < 0.001*kine_best){
	ktree.kine_energy_info.push_back(1); // range ...
      }else{
	ktree.kine_energy_info.push_back(0); // dQ/dx ...
      }
      ktree.kine_energy_included.push_back(1);
      
      used_showers.insert(shower);

      // either muon or electron in the showers ...
      if (shower->get_start_segment()->get_particle_type()!=11)
	ktree.kine_reco_add_energy += shower->get_start_segment()->get_particle_mass()/units::MeV;
      
    }else{
      // track ...
      used_segments.insert(*it);
      WCPPID::ProtoVertex *other_vertex = find_other_vertex(*it, main_vertex);
      segments_to_be_examined.push_back(std::make_pair(other_vertex, *it));

      float kine_best = (*it)->get_kine_best();
      float kine_charge = cal_kine_charge(*it);
      float kine_dQdx = (*it)->cal_kine_dQdx();
      float kine_range = (*it)->cal_kine_range();
      
      ktree.kine_energy_particle.push_back(kine_best/units::MeV);
      ktree.kine_particle_type.push_back((*it)->get_particle_type());
      
      if (fabs(kine_best - kine_charge) < 0.001 * kine_best){
	ktree.kine_energy_info.push_back(2); // charge ...
      }else if (fabs(kine_best - kine_range) < 0.001*kine_best){
	ktree.kine_energy_info.push_back(1); // range ...
      }else{
	ktree.kine_energy_info.push_back(0); // dQ/dx ...
      }
      ktree.kine_energy_included.push_back(1);

      if ((*it)->get_particle_type()==2212){ // proton 
	ktree.kine_reco_add_energy += ave_binding_energy/units::MeV;
      }else if ((*it)->get_particle_type()!=11){ // pion, muon, electron
	ktree.kine_reco_add_energy += (*it)->get_particle_mass()/units::MeV;
      }
      
      
    }
  }
  used_vertices.insert(main_vertex);

  
  
  while(segments_to_be_examined.size()>0){
    std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*> > temp_segments;
    for (auto it = segments_to_be_examined.begin(); it!= segments_to_be_examined.end(); it++){
      WCPPID::ProtoVertex *curr_vtx = it->first;
      WCPPID::ProtoSegment *prev_sg = it->second;
      
      if (used_vertices.find(curr_vtx)!=used_vertices.end()) continue;

      bool flag_reduce = false;

      
      for (auto it1 = map_vertex_segments[curr_vtx].begin(); it1!=map_vertex_segments[curr_vtx].end(); it1++){
	WCPPID::ProtoSegment *curr_sg = *it1;
	if (curr_sg == prev_sg) continue;

	if (curr_sg->get_particle_type() == prev_sg->get_particle_type() ||
	    prev_sg->get_particle_type() == 211 && curr_sg->get_particle_type()==13 ||
	    prev_sg->get_particle_type() == 13 && curr_sg->get_particle_type()==211 
	    ) flag_reduce = true;
	
	auto it2 = map_sg_shower.find(curr_sg);
		
	if (it2 == map_sg_shower.end()){
	  // tracks ...
	  if (used_segments.find(curr_sg)!=used_segments.end()) continue;
	  used_segments.insert(curr_sg);
	  
	  /* // set mother ... */
	  /* rtree.mc_mother[map_sgid_rtid[map_sg_sgid[curr_sg]]] = map_sg_sgid[prev_sg]; */
	  /* // set daughters ... */
	  /* rtree.mc_daughters->at(map_sgid_rtid[map_sg_sgid[prev_sg]]).push_back(map_sg_sgid[curr_sg]); */

	  float kine_best = curr_sg->get_kine_best();
	  float kine_charge = cal_kine_charge(curr_sg);
	  float kine_dQdx = curr_sg->cal_kine_dQdx();
	  float kine_range = curr_sg->cal_kine_range();
	  
	  ktree.kine_energy_particle.push_back(kine_best/units::MeV);
	  ktree.kine_particle_type.push_back(curr_sg->get_particle_type());
	  
	  if (fabs(kine_best - kine_charge) < 0.001 * kine_best){
	    ktree.kine_energy_info.push_back(2); // charge ...
	  }else if (fabs(kine_best - kine_range) < 0.001*kine_best){
	    ktree.kine_energy_info.push_back(1); // range ...
	  }else{
	    ktree.kine_energy_info.push_back(0); // dQ/dx ...
	  }
	  ktree.kine_energy_included.push_back(1);
	  
	  if (curr_sg->get_particle_type()==2212){ // proton 
	    ktree.kine_reco_add_energy += ave_binding_energy/units::MeV;
	  }else if (curr_sg->get_particle_type()!=11){ // pion, muon
	    ktree.kine_reco_add_energy += curr_sg->get_particle_mass()/units::MeV;
	  }
	  
	  
	  WCPPID::ProtoVertex *other_vertex = find_other_vertex(curr_sg, curr_vtx);
	  if (used_vertices.find(other_vertex) == used_vertices.end())
	    temp_segments.push_back(std::make_pair(other_vertex, curr_sg));
	}else{
	  // shower ...
	  WCPPID::WCShower *shower = it2->second;
	  float kine_best = shower->get_kine_best();
	  if (kine_best ==0 ) kine_best = shower->get_kine_charge();
	  float kine_charge = shower->get_kine_charge();
	  float kine_dQdx = shower->get_kine_dQdx();
	  float kine_range = shower->get_kine_range();
	  
	  ktree.kine_energy_particle.push_back(kine_best/units::MeV);
	  ktree.kine_particle_type.push_back(shower->get_start_segment()->get_particle_type());
	  
	  if (fabs(kine_best - kine_charge) < 0.001 * kine_best){
	    ktree.kine_energy_info.push_back(2); // charge ...
	  }else if (fabs(kine_best - kine_range) < 0.001*kine_best){
	    ktree.kine_energy_info.push_back(1); // range ...
	  }else{
	    ktree.kine_energy_info.push_back(0); // dQ/dx ...
	  }
	  ktree.kine_energy_included.push_back(1);
	  
	  // either muon or electron in the showers ...
	  if (shower->get_start_segment()->get_particle_type()!=11)
	    ktree.kine_reco_add_energy += shower->get_start_segment()->get_particle_mass()/units::MeV;
	  used_showers.insert(shower);
	}
      }
      used_vertices.insert(curr_vtx);

      if (flag_reduce){
        if (prev_sg->get_particle_type()==2212){ // proton 
	  ktree.kine_reco_add_energy -= ave_binding_energy/units::MeV;
	}else if (prev_sg->get_particle_type()!=11){ // pion, muon, electron
	  ktree.kine_reco_add_energy -= prev_sg->get_particle_mass()/units::MeV;
	}
      }
      
      // calculate added energy ...
    }
    segments_to_be_examined = temp_segments;

    
  }

  // loop over remaining showers ...
  for (auto it = showers.begin(); it!= showers.end(); it++){
    WCPPID::WCShower *shower = *it;
    WCPPID::ProtoSegment* curr_sg = shower->get_start_segment();
    std::pair<ProtoVertex*, int> pair_vertex = shower->get_start_vertex();
    if (pair_vertex.second > 3) continue;
    if (used_showers.find(shower) != used_showers.end()) continue;
    
    float kine_best = shower->get_kine_best();
    if (kine_best ==0 ) kine_best = shower->get_kine_charge();
    float kine_charge = shower->get_kine_charge();
    float kine_dQdx = shower->get_kine_dQdx();
    float kine_range = shower->get_kine_range();
	  
    ktree.kine_energy_particle.push_back(kine_best/units::MeV);
    ktree.kine_particle_type.push_back(shower->get_start_segment()->get_particle_type());
    
    if (fabs(kine_best - kine_charge) < 0.001 * kine_best){
      ktree.kine_energy_info.push_back(2); // charge ...
    }else if (fabs(kine_best - kine_range) < 0.001*kine_best){
      ktree.kine_energy_info.push_back(1); // range ...
    }else{
      ktree.kine_energy_info.push_back(0); // dQ/dx ...
    }
    if (pair_vertex.second != 3){
      ktree.kine_energy_included.push_back(1);
    }else{
      ktree.kine_energy_included.push_back(pair_vertex.second);
    }
	  
    // either muon or electron in the showers ...
    if (shower->get_start_segment()->get_particle_type()==11){
      //ktree.kine_reco_add_energy += shower->get_start_segment()->get_particle_mass()/units::MeV;
    }else if (shower->get_start_segment()->get_particle_type()==2212 && shower->get_start_segment()->get_length() > 5.0*units::cm){
      //      std::cout << shower->get_start_segment()->get_length()/units::cm << std::endl;
      ktree.kine_reco_add_energy += ave_binding_energy/units::MeV;
    }

    used_showers.insert(shower);
  }

  ktree.kine_reco_Enu = 0;
  for (size_t i=0;i!=ktree.kine_energy_particle.size();i++){
    ktree.kine_reco_Enu += ktree.kine_energy_particle.at(i);

    // std::cout << ktree.kine_energy_particle.at(i) << " " << ktree.kine_particle_type.at(i) << " " << ktree.kine_energy_info.at(i) << " " << ktree.kine_energy_included.at(i) << std::endl;
  }
  
  ktree.kine_reco_Enu += ktree.kine_reco_add_energy;
  
  //std::cout << ktree.kine_reco_Enu << " " << ktree.kine_reco_add_energy << std::endl;
  
  // pio reconstruction ...
  
  ktree.kine_pio_mass = kine_pio_mass/units::MeV;
  ktree.kine_pio_flag = kine_pio_flag;
  ktree.kine_pio_vtx_dis = kine_pio_vtx_dis/units::cm;
  
  ktree.kine_pio_energy_1 = kine_pio_energy_1/units::MeV;
  ktree.kine_pio_theta_1 = kine_pio_theta_1/3.1415926*180.;
  ktree.kine_pio_phi_1 = kine_pio_phi_1/3.1415926*180.;
  ktree.kine_pio_dis_1 = kine_pio_dis_1/units::cm;
  
  ktree.kine_pio_energy_2 = kine_pio_energy_2/units::MeV;
  ktree.kine_pio_theta_2 = kine_pio_theta_2/3.1415926*180.;
  ktree.kine_pio_phi_2 = kine_pio_phi_2/3.1415926*180.;
  ktree.kine_pio_dis_2 = kine_pio_dis_2/units::cm;

  ktree.kine_pio_angle = kine_pio_angle/3.1415926*180.;

  // std::cout << ktree.kine_pio_flag << " " << ktree.kine_pio_mass << " " << ktree.kine_pio_energy_1 << " " << ktree.kine_pio_energy_2 << " " << ktree.kine_pio_angle << std::endl;
  
  
  
}
