bool WCPPID::NeutrinoID::ssm_tagger(){

  TVector3 dir_beam(0,0,1);
  TVector3 dir_drift(1,0,0);
  TVector3 dir_vertical(0,1,0);
  TVector3 target_dir(0.46, 0.05, 0.885);
  TVector3 absorber_dir(0.33, 0.75, -0.59);  

  std::cout<<"Starting SSM tagger"<<std::endl;
  int  Nsm = 0; //Number of Short Muons
  int  Nsm_wivtx = 0; //Number of Short Muons wi vertex activity
  bool backwards_muon = false;
  TVector3 init_dir(0,0,0);
  TVector3 init_dir_10(0,0,0);
  TVector3 init_dir_15(0,0,0);
  TVector3 init_dir_20(0,0,0);
  TVector3 init_dir_10_bp(0,0,0);
  TVector3 init_dir_20_bp(0,0,0);
  TVector3 mom(0,0,0);

  //only filled if there is one ssm
    //properties of the ssm
      //dq/dx info
  double dq_dx_fwd_1 = -999;//1st dq dx point wrt vertex
  double dq_dx_fwd_2 = -999;//2nd dq dx point wrt vertex
  double dq_dx_fwd_3 = -999;//3rd dq dx point wrt vertex
  double dq_dx_fwd_4 = -999;//4th dq dx point wrt vertex
  double dq_dx_fwd_5 = -999;//5th dq dx point wrt vertex
  double dq_dx_bck_1 = -999;//1st dq dx point wrt end
  double dq_dx_bck_2 = -999;//2nd dq dx point wrt end
  double dq_dx_bck_3 = -999;//3rd dq dx point wrt end
  double dq_dx_bck_4 = -999;//4th dq dx point wrt end
  double dq_dx_bck_5 = -999;//5th dq dx point wrt end
  double d_dq_dx_fwd_12 = -999;//1st d dq dx point wrt vertex
  double d_dq_dx_fwd_23 = -999;//2nd d dq dx point wrt vertex
  double d_dq_dx_fwd_34 = -999;//3rd d dq dx point wrt vertex
  double d_dq_dx_fwd_45 = -999;//4th d dq dx point wrt vertex
  double d_dq_dx_bck_12 = -999;//1st d dq dx point wrt end
  double d_dq_dx_bck_23 = -999;//2nd d dq dx point wrt end
  double d_dq_dx_bck_34 = -999;//3rd d dq dx point wrt end
  double d_dq_dx_bck_45 = -999;//4th d dq dx point wrt end
  double max_dq_dx_fwd_3 = -999;//max dq/dx in first 3 points wrt vertex
  double max_dq_dx_fwd_5 = -999;//max dq/dx in first 5 points wrt vertex
  double max_dq_dx_bck_3 = -999;//max dq/dx in first 3 points wrt end
  double max_dq_dx_bck_5 = -999;//max dq/dx in first 5 points wrt end
  double max_d_dq_dx_fwd_3 = -999;//max d dq/dx in first 3 points wrt vertex
  double max_d_dq_dx_fwd_5 = -999;//max d dq/dx in first 5 points wrt vertex
  double max_d_dq_dx_bck_3 = -999;//max d dq/dx in first 3 points wrt end
  double max_d_dq_dx_bck_5 = -999;//max d dq/dx in first 5 points wrt end
  double medium_dq_dx = -999;//medium dq/dx over the entire track
  double medium_dq_dx_bp = -999;//medium dq/dx over the entire track without the vertex acivity
      //angluar info
  double angle_to_z = -999;//angle to the z det cord
  double angle_to_target = -999;//angle to the target
  double angle_to_absorber = -999;//angle to the absorber
  double angle_to_vertical = -999;//angle to the vertical

  double angle_to_z_10 = -999;//angle to the z det cord
  double angle_to_target_10 = -999;//angle to the target
  double angle_to_absorber_10 = -999;//angle to the absorber
  double angle_to_vertical_10 = -999;//angle to the vertical

  double angle_to_z_15 = -999;//angle to the z det cord
  double angle_to_target_15 = -999;//angle to the target
  double angle_to_absorber_15 = -999;//angle to the absorber
  double angle_to_vertical_15 = -999;//angle to the vertical

  double angle_to_z_20 = -999;//angle to the z det cord
  double angle_to_target_20 = -999;//angle to the target
  double angle_to_absorber_20 = -999;//angle to the absorber
  double angle_to_vertical_20 = -999;//angle to the vertical
    //directional info
  double x_dir = -999;//x dir det cord
  double y_dir = -999;//y dir det cord
  double z_dir = -999;//z dir det cord
    //energy info
  double kine_energy = -999;//KE of the muon with range
  double kine_energy_reduced = -999;//KE of the muon with range, exluding high dqdx at vertex
      //general properties
  double vtx_activity = 0;//flag for identified vertex activity
  double pdg = -999;//initial pdg assigned to the track
  double dQ_dx_cut = -999;//cut normally used to choose muons
  double score_mu_fwd = -999;//ID score normally normally used to choose muons
  double score_p_fwd = -999;//ID score normally normally used to choose protons
  double score_e_fwd = -999;//ID score normally normally used to choose e
  double score_mu_bck = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck = -999;//ID score normally normally used to choose e but calculated backwards
  double score_mu_fwd_bp = -999;//ID score normally normally used to choose muons but calculated without the vertex activity
  double score_p_fwd_bp = -999;//ID score normally normally used to choose protons but calculated without the vertex activity
  double score_e_fwd_bp = -999;//ID score normally normally used to choose e but calculated without the vertex activity
      //track "straighness"
  double length = -999;//length of track direct caluclated from point to point
  double direct_length = -999;//length of track direct from end to end
  //double direct_length/length = -999;
  double max_dev = -999;//furthest distance from track's primary axis along the track
    //number of other particles
  double n_prim_tracks_1 = 0;//number of other primary tracks greater than 1 cm
  double n_prim_tracks_3 = 0;//number of other primary tracks greater than 3 cm
  double n_prim_tracks_5 = 0;//number of other primary tracks greater than 5 cm
  double n_prim_tracks_8 = 0;//number of other primary tracks greater than 8 cm
  double n_prim_tracks_11 = 0;//number of other primary tracks greater than 11 cm
  double n_prim_all_1 = 0;//number of other primary tracks and showers greater than 1 cm
  double n_prim_all_3 = 0;//number of other primary tracks and showers greater than 3 cm
  double n_prim_all_5 = 0;//number of other primary tracks and showers greater than 5 cm
  double n_prim_all_8 = 0;//number of other primary tracks and showers greater than 8 cm
  double n_prim_all_11 = 0;//number of other primary tracks and showers greater than 11 cm
  double n_daughter_tracks_1 = 0;//number of other daughter tracks greater than 1 cm
  double n_daughter_tracks_3 = 0;//number of other daughter tracks greater than 3 cm
  double n_daughter_tracks_5 = 0;//number of other daughter tracks greater than 5 cm
  double n_daughter_tracks_8 = 0;//number of other daughter tracks greater than 8 cm
  double n_daughter_tracks_11 = 0;//number of other daughter tracks greater than 11 cm
  double n_daughter_all_1 = 0;//number of other daughter tracks and showers greater than 1 cm
  double n_daughter_all_3 = 0;//number of other daughter tracks and showers greater than 3 cm
  double n_daughter_all_5 = 0;//number of other daughter tracks and showers greater than 5 cm
  double n_daughter_all_8 = 0;//number of other daughter tracks and showers greater than 8 cm
  double n_daughter_all_11 = 0;//number of other daughter tracks and showers greater than 11 cm
    //properties of leading other primary track
  double pdg_prim_track1 = -999;//initial pdg assigned to the track
  double score_mu_fwd_prim_track1 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_prim_track1 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_prim_track1 = -999;//ID score normally normally used to choose e
  double score_mu_bck_prim_track1 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_prim_track1 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_prim_track1 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_prim_track1 = -999;//length of track direct caluclated from point to point
  double direct_length_prim_track1 = -999;//length of track direct from end to end
  double max_dev_prim_track1 = -999;//furthest distance from track's primary axis along the track
  double kine_energy_range_prim_track1 = -999;//KE of the muon with range
  double kine_energy_range_mu_prim_track1 = -999;
  double kine_energy_range_p_prim_track1 = -999;
  double kine_energy_range_e_prim_track1 = -999;
  double kine_energy_cal_prim_track1 = -999;//KE of the track with dq/dx
  double medium_dq_dx_prim_track1 = -999;//medium dq/dx over the entire track
  double x_dir_prim_track1 = -999;//x dir det cord
  double y_dir_prim_track1 = -999;//y dir det cord
  double z_dir_prim_track1 = -999;//z dir det cord
  double add_daught_track_counts_1_prim_track1 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_prim_track1 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_prim_track1 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_prim_track1 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_prim_track1 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_prim_track1 = -999;//additional daughter tracks+showers longer than 11 cm
  //properties of sub-leading other primary track
  double pdg_prim_track2 = -999;//initial pdg assigned to the track
  double score_mu_fwd_prim_track2 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_prim_track2 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_prim_track2 = -999;//ID score normally normally used to choose e
  double score_mu_bck_prim_track2 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_prim_track2 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_prim_track2 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_prim_track2 = -999;//length of track direct caluclated from point to point
  double direct_length_prim_track2 = -999;//length of track direct from end to end
  double max_dev_prim_track2 = -999;//furthest distance from track's primary axis along the track
  double kine_energy_range_prim_track2 = -999;//KE of the track with range
  double kine_energy_range_mu_prim_track2 = -999;
  double kine_energy_range_p_prim_track2 = -999;
  double kine_energy_range_e_prim_track2 = -999;
  double kine_energy_cal_prim_track2 = -999;//KE of the track with dq/dx
  double medium_dq_dx_prim_track2 = -999;//medium dq/dx over the entire track
  double x_dir_prim_track2 = -999;//x dir det cord
  double y_dir_prim_track2 = -999;//y dir det cord
  double z_dir_prim_track2 = -999;//z dir det cord
  double add_daught_track_counts_1_prim_track2 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_prim_track2 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_prim_track2 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_prim_track2 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_prim_track2 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_prim_track2 = -999;//additional daughter tracks+showers longer than 11 cm
    //properties of leading daughter track
  double pdg_daught_track1 = -999;//initial pdg assigned to the track
  double score_mu_fwd_daught_track1 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_daught_track1 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_daught_track1 = -999;//ID score normally normally used to choose e
  double score_mu_bck_daught_track1 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_daught_track1 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_daught_track1 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_daught_track1 = -999;//length of track direct caluclated from point to point
  double direct_length_daught_track1 = -999;//length of track direct from end to end
  double max_dev_daught_track1 = -999;//furthest distance from track's daughtary axis along the track
  double kine_energy_range_daught_track1 = -999;//KE of the track with range
  double kine_energy_range_mu_daught_track1 = -999;
  double kine_energy_range_p_daught_track1 = -999;
  double kine_energy_range_e_daught_track1 = -999;
  double kine_energy_cal_daught_track1 = -999;//KE of the track with dq/dx
  double medium_dq_dx_daught_track1 = -999;//medium dq/dx over the entire track 
  double x_dir_daught_track1 = -999;//x dir det cord
  double y_dir_daught_track1 = -999;//y dir det cord
  double z_dir_daught_track1 = -999;//z dir det cord
  double add_daught_track_counts_1_daught_track1 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_daught_track1 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_daught_track1 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_daught_track1 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_daught_track1 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_daught_track1 = -999;//additional daughter tracks+showers longer than 11 cm
  //properties of sub-leading daughter track
  double pdg_daught_track2 = -999;//initial pdg assigned to the track
  double score_mu_fwd_daught_track2 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_daught_track2 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_daught_track2 = -999;//ID score normally normally used to choose e
  double score_mu_bck_daught_track2 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_daught_track2 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_daught_track2 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_daught_track2 = -999;//length of track direct caluclated from point to point
  double direct_length_daught_track2 = -999;//length of track direct from end to end
  double max_dev_daught_track2 = -999;//furthest distance from track's daughtary axis along the track
  double kine_energy_range_daught_track2 = -999;//KE of the muon with range
  double kine_energy_range_mu_daught_track2 = -999;
  double kine_energy_range_p_daught_track2 = -999;
  double kine_energy_range_e_daught_track2 = -999;
  double kine_energy_cal_daught_track2 = -999;//KE of the muon with dq/dx
  double medium_dq_dx_daught_track2 = -999;//medium dq/dx over the entire track
  double x_dir_daught_track2 = -999;//x dir det cord
  double y_dir_daught_track2 = -999;//y dir det cord
  double z_dir_daught_track2 = -999;//z dir det cord
  double add_daught_track_counts_1_daught_track2 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_daught_track2 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_daught_track2 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_daught_track2 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_daught_track2 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_daught_track2 = -999;//additional daughter tracks+showers longer than 11 cm
    //properties of leading other primary shower
  double pdg_prim_shw1 = -999;//initial pdg assigned to the shw
  double score_mu_fwd_prim_shw1 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_prim_shw1 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_prim_shw1 = -999;//ID score normally normally used to choose e
  double score_mu_bck_prim_shw1 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_prim_shw1 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_prim_shw1 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_prim_shw1 = -999;//length of shw direct caluclated from point to point
  double direct_length_prim_shw1 = -999;//length of shw direct from end to end
  double max_dev_prim_shw1 = -999;//furthest distance from shw's primary axis along the shw
  double kine_energy_range_prim_shw1 = -999;//KE of the muon with range
  double kine_energy_range_mu_prim_shw1 = -999;
  double kine_energy_range_p_prim_shw1 = -999;
  double kine_energy_range_e_prim_shw1 = -999;
  double kine_energy_cal_prim_shw1 = -999;//KE of the muon with dq/dx
  double medium_dq_dx_prim_shw1 = -999;//medium dq/dx over the entire shw
  double x_dir_prim_shw1 = -999;//x dir det cord
  double y_dir_prim_shw1 = -999;//y dir det cord
  double z_dir_prim_shw1 = -999;//z dir det cord
  double add_daught_track_counts_1_prim_shw1 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_prim_shw1 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_prim_shw1 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_prim_shw1 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_prim_shw1 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_prim_shw1 = -999;//additional daughter tracks+showers longer than 11 cm
    //properties of sub-leading other primary shower
  double pdg_prim_shw2 = -999;//initial pdg assigned to the shw
  double score_mu_fwd_prim_shw2 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_prim_shw2 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_prim_shw2 = -999;//ID score normally normally used to choose e
  double score_mu_bck_prim_shw2 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_prim_shw2 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_prim_shw2 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_prim_shw2 = -999;//length of shw direct caluclated from point to point
  double direct_length_prim_shw2 = -999;//length of shw direct from end to end
  double max_dev_prim_shw2 = -999;//furthest distance from shw's primary axis along the shw
  double kine_energy_range_prim_shw2 = -999;//KE of the shower with range
  double kine_energy_range_mu_prim_shw2 = -999;
  double kine_energy_range_p_prim_shw2 = -999;
  double kine_energy_range_e_prim_shw2 = -999;
  double kine_energy_cal_prim_shw2 = -999;//KE of the shower with dq/dx
  double medium_dq_dx_prim_shw2 = -999;//medium dq/dx over the entire shw
  double x_dir_prim_shw2 = -999;//x dir det cord
  double y_dir_prim_shw2 = -999;//y dir det cord
  double z_dir_prim_shw2 = -999;//z dir det cord
  double add_daught_track_counts_1_prim_shw2 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_prim_shw2 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_prim_shw2 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_prim_shw2 = -999;//additional daughter tracks+showers longer than 5 cm 
  double add_daught_track_counts_11_prim_shw2 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_prim_shw2 = -999;//additional daughter tracks+showers longer than 11 cm
   //properties of leading daughter shower
  double pdg_daught_shw1 = -999;//initial pdg assigned to the shw
  double score_mu_fwd_daught_shw1 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_daught_shw1 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_daught_shw1 = -999;//ID score normally normally used to choose e
  double score_mu_bck_daught_shw1 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_daught_shw1 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_daught_shw1 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_daught_shw1 = -999;//length of shw direct caluclated from point to point
  double direct_length_daught_shw1 = -999;//length of shw direct from end to end
  double max_dev_daught_shw1 = -999;//furthest distance from shw's daughtary axis along the shw
  double kine_energy_range_daught_shw1 = -999;//KE of the shower with range
  double kine_energy_range_mu_daught_shw1 = -999;
  double kine_energy_range_p_daught_shw1 = -999;
  double kine_energy_range_e_daught_shw1 = -999;
  double kine_energy_cal_daught_shw1 = -999;//KE of the shower with dq/dx
  double medium_dq_dx_daught_shw1 = -999;//medium dq/dx over the entire shw 
  double x_dir_daught_shw1 = -999;//x dir det cord
  double y_dir_daught_shw1 = -999;//y dir det cord
  double z_dir_daught_shw1 = -999;//z dir det cord
  double add_daught_track_counts_1_daught_shw1 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_daught_shw1 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_daught_shw1 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_daught_shw1 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_daught_shw1 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_daught_shw1 = -999;//additional daughter tracks+showers longer than 11 cm
    //properties of sub-leading daughter shower
  double pdg_daught_shw2 = -999;//initial pdg assigned to the shw
  double score_mu_fwd_daught_shw2 = -999;//ID score normally normally used to choose muons
  double score_p_fwd_daught_shw2 = -999;//ID score normally normally used to choose protons
  double score_e_fwd_daught_shw2 = -999;//ID score normally normally used to choose e
  double score_mu_bck_daught_shw2 = -999;//ID score normally normally used to choose muons but calculated backwards
  double score_p_bck_daught_shw2 = -999;//ID score normally normally used to choose protons but calculated backwards
  double score_e_bck_daught_shw2 = -999;//ID score normally normally used to choose e but calculated backwards
  double length_daught_shw2 = -999;//length of shw direct caluclated from point to point
  double direct_length_daught_shw2 = -999;//length of shw direct from end to end
  double max_dev_daught_shw2 = -999;//furthest distance from shw's daughtary axis along the shw
  double kine_energy_range_daught_shw2 = -999;//KE of the shower with range
  double kine_energy_range_mu_daught_shw2 = -999;
  double kine_energy_range_p_daught_shw2 = -999;
  double kine_energy_range_e_daught_shw2 = -999;
  double kine_energy_cal_daught_shw2 = -999;//KE of the shower with dq/dx
  double medium_dq_dx_daught_shw2 = -999;//medium dq/dx over the entire shw
  double x_dir_daught_shw2 = -999;//x dir det cord
  double y_dir_daught_shw2 = -999;//y dir det cord
  double z_dir_daught_shw2 = -999;//z dir det cord
  double add_daught_track_counts_1_daught_shw2 = -999;//additional daughter tracks longer than 1 cm
  double add_daught_all_counts_1_daught_shw2 = -999;//additional daughter tracks+showers longer than 1 cm
  double add_daught_track_counts_5_daught_shw2 = -999;//additional daughter tracks longer than 5 cm
  double add_daught_all_counts_5_daught_shw2 = -999;//additional daughter tracks+showers longer than 5 cm
  double add_daught_track_counts_11_daught_shw2 = -999;//additional daughter tracks longer than 11 cm
  double add_daught_all_counts_11_daught_shw2 = -999;//additional daughter tracks+showers longer than 11 cm
    //event level properties
  double nu_angle_to_z = -999;//angle to the z det cord
  double nu_angle_to_target = -999;
  double nu_angle_to_absorber = -999;
  double nu_angle_to_vertical = -999;
  double track_angle_to_z = -999;//angle to the z det cord
  double track_angle_to_target = -999;
  double track_angle_to_absorber = -999;
  double track_angle_to_vertical = -999;

  //saves all smm, if they have vtx activity, if they are backwards and their length
  std::map<WCPPID::ProtoSegment*, std::tuple<bool,bool,double> > all_ssm_sg;
  //saves the main ssm
  WCPPID::ProtoSegment* ssm_sg;
  //KDAR tagger, check everything connected to the vertex, make sure there is a ssm
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    double sg_length = sg->get_length()/units::cm;
    double sg_direct_length = sg->get_direct_length()/units::cm;
    std::cout<<"length check "<<sg_length<<"  dq/dx check "<<sg->get_medium_dQ_dx()/(43e3/units::cm)<<std::endl;
    bool sg_flag_backwards_muon = false;
    bool sg_flag_vtx_activity = false;
    //if( (length <= 31*units::cm && length >= 11*units::cm && abs(sg->get_particle_type())==13) || (direct_length > 0.9 * length && length <= 31*units::cm && length >= 11*units::cm && abs(sg->get_particle_type())==11)){
    //if( (sg_length <= 32 && sg_length >= 5 && abs(sg->get_particle_type())==13) || (sg_direct_length > 0.9 * sg_length && sg_length <= 32 && sg_length >= 5 && abs(sg->get_particle_type())==11)){
    if( (sg_length <= 34 && sg_length >= 1 && abs(sg->get_particle_type())==13) || (sg_direct_length > 0.9 * sg_length && sg_length <= 32 && sg_length >= 1 && abs(sg->get_particle_type())==11)){
    //veto stuff we still dont want
      if(sg_length<1) {std::cout<<"\n too short"<<std::endl; continue;}//prviously was 5
      //if(sg->get_medium_dQ_dx()/(43e3/units::cm)<1.1){ std::cout<<"\n bad dq/dx"<<std::endl; continue;}//make more conservative?
      if(sg->get_medium_dQ_dx()/(43e3/units::cm)<0.95){ std::cout<<"\n bad dq/dx"<<std::endl; continue;}//make more conservative?, prviously was 1
      Nsm+=1;
      if(abs(sg->get_particle_type())==13){ std::cout<<"Short Muon"<<std::endl;}
      else if(abs(sg->get_particle_type())==11){ std::cout<<"Short Straight Electron"<<std::endl;}

      //make the d dq/dx vector
      std::vector<double>& vec_dQ = sg->get_dQ_vec();
      std::vector<double>& vec_dx = sg->get_dx_vec();
      std::vector<double> vec_d_dq_dx;
      bool vtx_activity_fwd = false;
      bool vtx_activity_bck = false;
      double sg_max_d_dq_dx_fwd_5 = 0;
      double sg_max_d_dq_dx_bck_5 = 0;
      std::cout<<std::endl;
      double last_dq_dx = vec_dQ.at(0)/vec_dx.at(0)/(43e3/units::cm);
      for (int point = 1; point<vec_dQ.size(); point++){
        double dq_dx = vec_dQ.at(point)/vec_dx.at(point)/(43e3/units::cm);
        vec_d_dq_dx.push_back( dq_dx - last_dq_dx  );
        std::cout<<dq_dx - last_dq_dx<<", ";
        last_dq_dx = dq_dx;
      }
      std::cout<<std::endl;
      //check the first 5 points for a large change in dq/dx
      int end = 4;
      if(end>vec_d_dq_dx.size()) end=vec_d_dq_dx.size();
      std::cout<<"end "<<end<<std::endl;
      for (int point = 0; point<end; point++){
        if(abs(vec_d_dq_dx.at(point))>0.7 && abs(vec_d_dq_dx.at(point))>sg_max_d_dq_dx_fwd_5){
          vtx_activity_fwd = true;
          sg_max_d_dq_dx_fwd_5 = abs(vec_d_dq_dx.at(point));
        }
      }

      //check the last 5 points for a large change in dq/dx
      int start=vec_d_dq_dx.size()-4;
      if(start<0) start=0;
      std::cout<<"start: "<<start<<std::endl;
      for (int point = start; point<vec_d_dq_dx.size(); point++){
        if(abs(vec_d_dq_dx.at(point))>0.7 && abs(vec_d_dq_dx.at(point))>max_d_dq_dx_bck_5){
          vtx_activity_bck = true;
          max_d_dq_dx_bck_5 = abs(vec_d_dq_dx.at(point));
        }
      }

      //decide which direction has vertex activity 
      if(vtx_activity_bck && vtx_activity_fwd){
        if(sg_max_d_dq_dx_fwd_5<1.0) vtx_activity_fwd = false;
        if(sg_max_d_dq_dx_bck_5<1.0) vtx_activity_bck = false;
        if(!vtx_activity_bck && !vtx_activity_fwd){
          if(sg_max_d_dq_dx_fwd_5>max_d_dq_dx_bck_5){
            vtx_activity_fwd = true;
            vtx_activity_bck = false;
          }else{
            vtx_activity_fwd = false;
            vtx_activity_bck = true;
          }   
        }else if(vtx_activity_bck && vtx_activity_fwd){
          if(sg_max_d_dq_dx_fwd_5<1.3) vtx_activity_fwd = false;
          if(sg_max_d_dq_dx_bck_5<1.3) vtx_activity_bck = false;
          if(!vtx_activity_bck && !vtx_activity_fwd){
            if(sg_max_d_dq_dx_fwd_5>max_d_dq_dx_bck_5){
              vtx_activity_fwd = true;
              vtx_activity_bck = false;
            }else{
              vtx_activity_fwd = false;
              vtx_activity_bck = true;
            }
          }else if(vtx_activity_bck && vtx_activity_fwd){
            if(sg_max_d_dq_dx_fwd_5>max_d_dq_dx_bck_5){
              vtx_activity_fwd = true;
              vtx_activity_bck = false;
            }else{
              vtx_activity_fwd = false;
              vtx_activity_bck = true;
            }
          }
        }
      }

      //flag the muon as backwards if the direction with vertex activity does not match where the start of the muon is
      if(vtx_activity_bck){
        Nsm_wivtx+=1;
        if(sg->get_flag_dir()==1) sg_flag_backwards_muon=true;
        sg_flag_vtx_activity = true;
        ssm_sg = sg;
      }else if(vtx_activity_fwd){
        Nsm_wivtx+=1;
        if(sg->get_flag_dir()==-1) sg_flag_backwards_muon=true;
        sg_flag_vtx_activity = true;
      }
std::cout<<"sg_flag_backwards_muon "<<sg_flag_backwards_muon<<std::endl;
if(sg_flag_backwards_muon) std::cout<<"backwards_muon "<<std::endl;
      std::tuple<bool,bool,double> ssm_info = {sg_flag_vtx_activity,sg_flag_backwards_muon,sg_length};
      all_ssm_sg.insert(std::pair<WCPPID::ProtoSegment*, std::tuple<bool,bool,double>>(sg,ssm_info));
    }

  }//end loop over all segments connected to the vertex to determine if we have an ssm


  tagger_info.ssm_Nsm = Nsm;//number of short straight muons
  tagger_info.ssm_Nsm_wivtx = Nsm_wivtx;//number of short straight muons with vertex activity

  std::cout<<"tagger_info.ssm_Nsm "<<tagger_info.ssm_Nsm<<"  tagger_info.ssm_Nsm_wivtx "<<tagger_info.ssm_Nsm_wivtx<<std::endl;
  if( Nsm ==0 ){ std::cout<<"Exit Tagger"<<std::endl; return exit_ssm_tagger(); }

  //decision on which muon is the ssmi
  if(all_ssm_sg.begin()==all_ssm_sg.end()) std::cout<<"problem"<<std::endl;
  for (auto it = all_ssm_sg.begin(); it!=all_ssm_sg.end(); it++){
    std::tuple<bool,bool,double> ssm_info = it->second;
    std::cout<<"  get<0>(ssm_info) "<<get<0>(ssm_info)<<"  get<1>(ssm_info) "<<get<1>(ssm_info)<<"  get<2>(ssm_info) "<<get<2>(ssm_info)<<std::endl;
    if( get<0>(ssm_info)==true && ( get<2>(ssm_info) > length || vtx_activity==false ) ){      
      ssm_sg = it->first;
      length = get<2>(ssm_info);
      backwards_muon = get<1>(ssm_info);
      vtx_activity = get<0>(ssm_info);
    }
    else if( get<2>(ssm_info) > length && vtx_activity==false){
      ssm_sg = it->first;
      length = get<2>(ssm_info);
      backwards_muon = get<1>(ssm_info);
      vtx_activity = get<0>(ssm_info);
    } 
  }

  std::cout<<"abs(ssm_sg->get_particle_type() "<<ssm_sg->get_particle_type()<<std::endl;
  pdg = abs(ssm_sg->get_particle_type());
  direct_length = ssm_sg->get_direct_length()/units::cm;

  std::cout<<"Direction"<<std::endl;
  //get the direction and angle of the muon
  int dir = ssm_sg->get_flag_dir();
  std::cout<<"dir "<<dir<<std::endl;
  if(backwards_muon) dir=dir*-1;
  std::cout<<"dir "<<dir<<std::endl;
  init_dir = ssm_sg->cal_dir_3vector(dir,5);
  init_dir_10 = ssm_sg->cal_dir_3vector(dir,10);
  init_dir_15 = ssm_sg->cal_dir_3vector(dir,15);
  init_dir_20 = ssm_sg->cal_dir_3vector(dir,20);
  init_dir = init_dir.Unit();
  init_dir_10 = init_dir_10.Unit();
  init_dir_15 = init_dir_15.Unit();
  init_dir_20 = init_dir_20.Unit();
  x_dir = init_dir_10[0];
  y_dir = init_dir_10[1];
  z_dir = init_dir_10[2];
/*	
// old way, weird cases where I do not understand the output	
	std::cout<<"Counting Daughters and prims"<<std::endl;
	  //count our daughter and primary particles
  std::pair< std::pair<int, int>,std::pair<int, int> > pair_result_1 = count_daughters(ssm_sg,1.0);
  n_prim_tracks_1 = pair_result_1.first.first;
  n_prim_all_1 = pair_result_1.first.second;
  n_daughter_tracks_1 = pair_result_1.second.first;
  n_daughter_all_1 = pair_result_1.second.second;
  std::pair< std::pair<int, int>,std::pair<int, int> > pair_result_3 = count_daughters(ssm_sg,3.0);
  n_prim_tracks_3 = pair_result_3.first.first;
  n_prim_all_3 = pair_result_3.first.second;
  n_daughter_tracks_3 = pair_result_3.second.first;
  n_daughter_all_3 = pair_result_3.second.second;  
  std::pair< std::pair<int, int>,std::pair<int, int> > pair_result_5 = count_daughters(ssm_sg,5.0);
  n_prim_tracks_5 = pair_result_5.first.first;
  n_prim_all_5 = pair_result_5.first.second;
  n_daughter_tracks_5 = pair_result_5.second.first;
  n_daughter_all_5 = pair_result_5.second.second;
  std::pair< std::pair<int, int>,std::pair<int, int> > pair_result_8 = count_daughters(ssm_sg,8.0);
  n_prim_tracks_8 = pair_result_8.first.first;
  n_prim_all_8 = pair_result_8.first.second;
  n_daughter_tracks_8 = pair_result_8.second.first;
  n_daughter_all_8 = pair_result_8.second.second;
  std::pair< std::pair<int, int>,std::pair<int, int> > pair_result_11 = count_daughters(ssm_sg,11.0);
  n_prim_tracks_11 = pair_result_11.first.first;
  n_prim_all_11 = pair_result_11.first.second;
  n_daughter_tracks_11 = pair_result_11.second.first;
  n_daughter_all_11 = pair_result_11.second.second;
std::cout<<"backwards_muon "<<backwards_muon<<" dir "<<dir<<std::endl;
  if(backwards_muon){
  //if(dir==-1){
std::cout<<"backwards_muon "<<std::endl;
    n_prim_tracks_1 = pair_result_1.second.first;
    n_prim_all_1 = pair_result_1.second.second;
    n_daughter_tracks_1 = pair_result_1.first.first;
    n_daughter_all_1 = pair_result_1.first.second;
    n_prim_tracks_3 = pair_result_3.second.first;
    n_prim_all_3 = pair_result_3.second.second;
    n_daughter_tracks_3 = pair_result_3.first.first;
    n_daughter_all_3 = pair_result_3.first.second;
    n_prim_tracks_5 = pair_result_5.second.first;
    n_prim_all_5 = pair_result_5.second.second;
    n_daughter_tracks_5 = pair_result_5.first.first;
    n_daughter_all_5 = pair_result_5.first.second;
    n_prim_tracks_8 = pair_result_8.second.first;
    n_prim_all_8 = pair_result_8.second.second;
    n_daughter_tracks_8 = pair_result_8.first.first;
    n_daughter_all_8 = pair_result_8.first.second;
    n_prim_tracks_11 = pair_result_11.second.first;
    n_prim_all_11 = pair_result_11.second.second;
    n_daughter_tracks_11 = pair_result_11.first.first;
    n_daughter_all_11 = pair_result_11.first.second;
  }
*/
  max_dev = ssm_sg->get_max_deviation(int(0),int(ssm_sg->get_point_vec().size()))/units::cm;
/*
std::cout<<"Direction"<<std::endl;
  //get the direction and angle of the muon
  int dir = ssm_sg->get_flag_dir();
  if(backwards_muon) dir=dir*-1;
  init_dir = ssm_sg->cal_dir_3vector(dir,5);
  init_dir_10 = ssm_sg->cal_dir_3vector(dir,10);
  init_dir_15 = ssm_sg->cal_dir_3vector(dir,15);
  init_dir_20 = ssm_sg->cal_dir_3vector(dir,20);
  init_dir = init_dir.Unit();
  init_dir_10 = init_dir_10.Unit();
  init_dir_15 = init_dir_15.Unit();
  init_dir_20 = init_dir_20.Unit();
*/
    std::cout<<"Angle"<<std::endl;
  angle_to_z = acos( init_dir[0]*dir_beam[0] + init_dir[1]*dir_beam[1] + init_dir[2]*dir_beam[2] );
  angle_to_z_10 = acos( init_dir_10[0]*dir_beam[0] + init_dir_10[1]*dir_beam[1] + init_dir_10[2]*dir_beam[2] );
  angle_to_z_15 = acos( init_dir_15[0]*dir_beam[0] + init_dir_15[1]*dir_beam[1] + init_dir_15[2]*dir_beam[2] );
  angle_to_z_20 = acos( init_dir_20[0]*dir_beam[0] + init_dir_20[1]*dir_beam[1] + init_dir_20[2]*dir_beam[2] );

  angle_to_target = acos( init_dir[0]*target_dir[0] + init_dir[1]*target_dir[1] + init_dir[2]*target_dir[2] );
  angle_to_target_10 = acos( init_dir_10[0]*target_dir[0] + init_dir_10[1]*target_dir[1] + init_dir_10[2]*target_dir[2] );
  angle_to_target_15 = acos( init_dir_15[0]*target_dir[0] + init_dir_15[1]*target_dir[1] + init_dir_15[2]*target_dir[2] );
  angle_to_target_20 = acos( init_dir_20[0]*target_dir[0] + init_dir_20[1]*target_dir[1] + init_dir_20[2]*target_dir[2] );

  angle_to_absorber = acos( init_dir[0]*absorber_dir[0] + init_dir[1]*absorber_dir[1] + init_dir[2]*absorber_dir[2] );
  angle_to_absorber_10 = acos( init_dir_10[0]*absorber_dir[0] + init_dir_10[1]*absorber_dir[1] + init_dir_10[2]*absorber_dir[2] );
  angle_to_absorber_15 = acos( init_dir_15[0]*absorber_dir[0] + init_dir_15[1]*absorber_dir[1] + init_dir_15[2]*absorber_dir[2] );
  angle_to_absorber_20 = acos( init_dir_20[0]*absorber_dir[0] + init_dir_20[1]*absorber_dir[1] + init_dir_20[2]*absorber_dir[2] );

  angle_to_vertical = acos( init_dir[0]*dir_vertical[0] + init_dir[1]*dir_vertical[1] + init_dir[2]*dir_vertical[2] );
  angle_to_vertical_10 = acos( init_dir_10[0]*dir_vertical[0] + init_dir_10[1]*dir_vertical[1] + init_dir_10[2]*dir_vertical[2] );  
  angle_to_vertical_15 = acos( init_dir_15[0]*dir_vertical[0] + init_dir_15[1]*dir_vertical[1] + init_dir_15[2]*dir_vertical[2] );    
  angle_to_vertical_20 = acos( init_dir_20[0]*dir_vertical[0] + init_dir_20[1]*dir_vertical[1] + init_dir_20[2]*dir_vertical[2] );

std::cout<<"dqdx info"<<std::endl;
  //add dq/dx info
  int break_point_fwd = -1;
  int break_point_bck = -1;
  int break_point = 0;
  double reduced_muon_length = length;
  std::vector<double>& vec_dQ = ssm_sg->get_dQ_vec();
  std::vector<double>& vec_dx = ssm_sg->get_dx_vec();
  std::vector<double> vec_dq_dx;
  std::vector<double> vec_d_dq_dx;
  std::vector<double> vec_abs_d_dq_dx;
  double last_dq_dx = vec_dQ.at(0)/vec_dx.at(0)/(43e3/units::cm);
  for (int point = 0; point<vec_dQ.size(); point++){
    double dq_dx = vec_dQ.at(point)/vec_dx.at(point)/(43e3/units::cm);
    vec_dq_dx.push_back( dq_dx  );
    std::cout<<dq_dx<<", ";
    if (point>0) {vec_d_dq_dx.push_back( dq_dx - last_dq_dx  ); vec_abs_d_dq_dx.push_back( abs(dq_dx - last_dq_dx)  );}
    last_dq_dx = dq_dx;
  }
std::cout<<std::endl;
std::cout<<"fwd vtx"<<std::endl;
  //check the first 5 points for a large change in dq/dx
  int end = 4;
  if(vec_d_dq_dx.size()<end) end = vec_d_dq_dx.size();
  for (int point = 0; point<end; point++){
    if(fabs(vec_d_dq_dx.at(point))>0.7){
      vtx_activity = true;
      break_point_fwd = point+1;
    }
/*
    if(abs(vec_d_dq_dx.at(point))>max_d_dq_dx_fwd_5){
      max_d_dq_dx_fwd_5 = abs(vec_d_dq_dx.at(point));
    }
    if(abs(vec_dq_dx.at(point))>max_dq_dx_fwd_5){
      max_dq_dx_fwd_5 = abs(vec_dq_dx.at(point));
    }
    if(point>2) continue;
    if(abs(vec_d_dq_dx.at(point))>max_d_dq_dx_fwd_3){
      max_d_dq_dx_fwd_3 = abs(vec_d_dq_dx.at(point));
    }
    if(abs(vec_dq_dx.at(point))>max_dq_dx_fwd_3){
      max_dq_dx_fwd_3 = abs(vec_dq_dx.at(point));
    }
*/
  }
  if(dir==1 && vtx_activity == true){
    for (int point = 0; point<break_point_fwd; point++){
      reduced_muon_length-=abs(vec_dx.at(point)/units::cm);
      std::cout<<"fwd subtracting "<<abs(vec_dx.at(point)/units::cm)<<std::endl;
    }
    break_point = break_point_fwd;
  }
  max_dq_dx_fwd_3 = *std::max_element(vec_dq_dx.begin(), vec_dq_dx.begin()+std::min(int(vec_dq_dx.size()),3));
  max_dq_dx_fwd_5 = *std::max_element(vec_dq_dx.begin(), vec_dq_dx.begin()+std::min(int(vec_dq_dx.size()),5));
  max_d_dq_dx_fwd_3 = *std::max_element(vec_abs_d_dq_dx.begin(), vec_abs_d_dq_dx.begin()+std::min(int(vec_abs_d_dq_dx.size()),2));
  max_d_dq_dx_fwd_5 = *std::max_element(vec_abs_d_dq_dx.begin(), vec_abs_d_dq_dx.begin()+std::min(int(vec_abs_d_dq_dx.size()),4));

  //std::cout<<"max_dq_dx_fwd_5 old: "<<max_dq_dx_fwd_5<<"  new: "<<*std::max_element(vec_dq_dx.begin(), vec_dq_dx.begin()+std::min(int(vec_dq_dx.size()),5))<<std::endl;  
  //std::cout<<"max_dq_dx_fwd_3 old: "<<max_dq_dx_fwd_3<<"  new: "<<*std::max_element(vec_dq_dx.begin(), vec_dq_dx.begin()+std::min(int(vec_dq_dx.size()),3))<<std::endl;

  //std::cout<<"max_d_dq_dx_fwd_5 old: "<<max_d_dq_dx_fwd_5<<"  new: "<<*std::max_element(vec_abs_d_dq_dx.begin(), vec_abs_d_dq_dx.begin()+std::min(int(vec_abs_d_dq_dx.size()),4))<<std::endl;
  //std::cout<<"max_d_dq_dx_fwd_3 old: "<<max_d_dq_dx_fwd_3<<"  new: "<<*std::max_element(vec_abs_d_dq_dx.begin(), vec_abs_d_dq_dx.begin()+std::min(int(vec_abs_d_dq_dx.size()),2))<<std::endl;
  

  std::cout<<"bck vtx"<<std::endl;
  bool back_vtx_activity = false;
  //check the last 5 points for a large change in dq/dx
  int start = vec_d_dq_dx.size()-4;
  if(start<0) start=0;
  for (int point = start; point<vec_d_dq_dx.size(); point++){
    if(fabs(vec_d_dq_dx.at(point))>0.7){
      vtx_activity = true;
      if(back_vtx_activity==false) break_point_bck = point; back_vtx_activity = true;
    }
/* 
    if(abs(vec_d_dq_dx.at(point))>max_d_dq_dx_bck_5){
      max_d_dq_dx_bck_5 = abs(vec_d_dq_dx.at(point));
    }
    if(abs(vec_dq_dx.at(point))>max_dq_dx_bck_5){
      max_dq_dx_bck_5 = abs(vec_dq_dx.at(point));
    }
    if(point<vec_d_dq_dx.size()-3) continue;
    if(abs(vec_d_dq_dx.at(point))>max_d_dq_dx_bck_3){
      max_d_dq_dx_bck_3 = abs(vec_d_dq_dx.at(point));
    }
    if(abs(vec_dq_dx.at(point))>max_dq_dx_bck_3){
      max_dq_dx_bck_3 = abs(vec_dq_dx.at(point));
    }
  */
  }

  max_dq_dx_bck_3 = *std::max_element(vec_dq_dx.begin()+ std::max(int(vec_dq_dx.size()-3),0), vec_dq_dx.end());
  max_dq_dx_bck_5 = *std::max_element(vec_dq_dx.begin()+ std::max(int(vec_dq_dx.size()-5),0), vec_dq_dx.end());
  max_d_dq_dx_bck_3 = *std::max_element(vec_abs_d_dq_dx.begin()+ std::max(int(vec_abs_d_dq_dx.size()-2),0), vec_abs_d_dq_dx.end());
  max_d_dq_dx_bck_5 = *std::max_element(vec_abs_d_dq_dx.begin()+ std::max(int(vec_abs_d_dq_dx.size()-4),0), vec_abs_d_dq_dx.end());

  //start = vec_dq_dx.size()-5;
  //if(start<0) start=0;
  //std::cout<<"max_dq_dx_bck_5 old: "<<max_dq_dx_bck_5<<"  new: "<<*std::max_element(vec_dq_dx.begin()+start, vec_dq_dx.end())<<std::endl;
  //std::cout<<"max_dq_dx_bck_5 old: "<<max_dq_dx_bck_5<<"  new: "<<*std::max_element(vec_dq_dx.begin()+ std::max(int(vec_dq_dx.size()-5),0), vec_dq_dx.end())<<std::endl;
  //start = vec_dq_dx.size()-3;
  //if(start<0) start=0;
  //std::cout<<"max_dq_dx_bck_3 old: "<<max_dq_dx_bck_3<<"  new: "<<*std::max_element(vec_dq_dx.begin()+start, vec_dq_dx.end())<<std::endl;
  //std::cout<<"max_dq_dx_bck_3 old: "<<max_dq_dx_bck_3<<"  new: "<<*std::max_element(vec_dq_dx.begin()+ std::max(int(vec_dq_dx.size()-3),0), vec_dq_dx.end())<<std::endl;

  //start = vec_d_dq_dx.size()-4;
  //if(start<0) start=0;
  //std::cout<<"max_d_dq_dx_bck_5 old: "<<max_d_dq_dx_bck_5<<"  new: "<<*std::max_element(vec_abs_d_dq_dx.begin()+start, vec_abs_d_dq_dx.end())<<std::endl;
  //std::cout<<"max_d_dq_dx_bck_5 old: "<<max_d_dq_dx_bck_5<<"  new: "<<*std::max_element(vec_abs_d_dq_dx.begin()+ std::max(int(vec_abs_d_dq_dx.size()-4),0), vec_abs_d_dq_dx.end())<<std::endl;
  //start = vec_d_dq_dx.size()-2;
  //if(start<0) start=0;
  //std::cout<<"max_d_dq_dx_bck_3 old: "<<max_d_dq_dx_bck_3<<"  new: "<<*std::max_element(vec_abs_d_dq_dx.begin()+start, vec_abs_d_dq_dx.end())<<std::endl;
  //std::cout<<"max_d_dq_dx_bck_3 old: "<<max_d_dq_dx_bck_3<<"  new: "<<*std::max_element(vec_abs_d_dq_dx.begin()+ std::max(int(vec_abs_d_dq_dx.size()-2),0), vec_abs_d_dq_dx.end())<<std::endl;
  if(dir==-1 && vtx_activity == true){
    for (int point = break_point_bck; point<vec_dx.size(); point++){
      reduced_muon_length-=abs(vec_dx.at(point)/units::cm);
      std::cout<<"bck subtracting "<<abs(vec_dx.at(point)/units::cm)<<std::endl;
    }
    break_point = break_point_bck;
  }
/*
  //for short stuff, put the empty values on the front for backwards, and the back for forwards
  while(vec_dq_dx.size()<5){
	  if(dir==-1) {vec_dq_dx.insert(vec_dq_dx.begin(), 0);}
          else {vec_dq_dx.push_back(0);}
  }
  while(vec_d_dq_dx.size()<4){
	  if(dir==-1) {vec_d_dq_dx.insert(vec_d_dq_dx.begin(), 0);}
	  else{ vec_d_dq_dx.push_back(0);}
  }
*/
  if(vec_dq_dx.size()<1){dq_dx_fwd_1=0;}
  else {dq_dx_fwd_1 = vec_dq_dx.at(0);}
  if(vec_dq_dx.size()<2){dq_dx_fwd_2=0;}
  else {dq_dx_fwd_2 = vec_dq_dx.at(1);}
  if(vec_dq_dx.size()<3){dq_dx_fwd_3=0;}
  else {dq_dx_fwd_3 = vec_dq_dx.at(2);}
  if(vec_dq_dx.size()<4){dq_dx_fwd_4=0;}
  else {dq_dx_fwd_4 = vec_dq_dx.at(3);}
  if(vec_dq_dx.size()<5){dq_dx_fwd_5=0;}
  else {dq_dx_fwd_5 = vec_dq_dx.at(4);}

  if(vec_dq_dx.size()<1){dq_dx_bck_1=0;}
  else {dq_dx_bck_1 = vec_dq_dx.at( vec_dq_dx.size()-1 );}
  if(vec_dq_dx.size()<2){dq_dx_bck_2=0;}
  else {dq_dx_bck_2 = vec_dq_dx.at( vec_dq_dx.size()-2 );}
  if(vec_dq_dx.size()<3){dq_dx_bck_3=0;}
  else {dq_dx_bck_3 = vec_dq_dx.at( vec_dq_dx.size()-3 );}
  if(vec_dq_dx.size()<4){dq_dx_bck_4=0;}
  else {dq_dx_bck_4 = vec_dq_dx.at( vec_dq_dx.size()-4 );}
  if(vec_dq_dx.size()<5){dq_dx_bck_5=0;}
  else {dq_dx_bck_5 = vec_dq_dx.at( vec_dq_dx.size()-5 );}

  if(vec_d_dq_dx.size()<1){d_dq_dx_fwd_12=0;}
  else {d_dq_dx_fwd_12 = vec_d_dq_dx.at(0);}
  if(vec_d_dq_dx.size()<2){d_dq_dx_fwd_23=0;}
  else {d_dq_dx_fwd_23 = vec_d_dq_dx.at(1);}
  if(vec_d_dq_dx.size()<3){d_dq_dx_fwd_34=0;}
  else {d_dq_dx_fwd_34 = vec_d_dq_dx.at(2);}
  if(vec_d_dq_dx.size()<4){d_dq_dx_fwd_45=0;}
  else {d_dq_dx_fwd_45 = vec_d_dq_dx.at(3);}

  if(vec_d_dq_dx.size()<1){d_dq_dx_bck_12=0;}
  else {d_dq_dx_bck_12 = vec_d_dq_dx.at( vec_d_dq_dx.size()-1 );}
  if(vec_d_dq_dx.size()<2){d_dq_dx_bck_23=0;}
  else {d_dq_dx_bck_23 = vec_d_dq_dx.at( vec_d_dq_dx.size()-2 );}
  if(vec_d_dq_dx.size()<3){d_dq_dx_bck_34=0;}
  else {d_dq_dx_bck_34 = vec_d_dq_dx.at( vec_d_dq_dx.size()-3 );}
  if(vec_d_dq_dx.size()<4){d_dq_dx_bck_45=0;}
  else {d_dq_dx_bck_45 = vec_d_dq_dx.at( vec_d_dq_dx.size()-4 );}

  if(dir==-1){
    std::swap(dq_dx_fwd_1,dq_dx_bck_1);
    std::swap(dq_dx_fwd_2,dq_dx_bck_2);
    std::swap(dq_dx_fwd_3,dq_dx_bck_3);
    std::swap(dq_dx_fwd_4,dq_dx_bck_4);
    std::swap(dq_dx_fwd_5,dq_dx_bck_5);
    std::swap(d_dq_dx_fwd_12,d_dq_dx_bck_12);
    std::swap(d_dq_dx_fwd_23,d_dq_dx_bck_23);
    std::swap(d_dq_dx_fwd_34,d_dq_dx_bck_34);
    std::swap(d_dq_dx_fwd_45,d_dq_dx_bck_45);
    std::swap(max_dq_dx_fwd_3,max_dq_dx_bck_3);
    std::swap(max_dq_dx_fwd_5,max_dq_dx_bck_5);  
    std::swap(max_d_dq_dx_fwd_3,max_d_dq_dx_bck_5);
    std::swap(max_d_dq_dx_fwd_5,max_d_dq_dx_bck_5);
  }
  std::cout<<"Scores"<<std::endl;
  std::vector<double> scores = get_scores(ssm_sg);
  score_mu_fwd = scores.at(0);
  score_p_fwd = scores.at(1);
  score_e_fwd= scores.at(2);
  score_mu_bck = scores.at(3);
  score_p_bck = scores.at(4);
  score_e_bck = scores.at(5);
  std::cout<<score_mu_bck<<" "<<score_p_bck<<std::endl;
  if(dir==-1){
    std::swap(score_mu_fwd,score_mu_bck);
    std::swap(score_p_fwd,score_p_bck);
    std::swap(score_e_fwd,score_e_bck);
  }
  std::cout<<score_mu_bck<<" "<<score_p_bck<<std::endl;
  score_mu_fwd_bp = score_mu_fwd;
  score_p_fwd_bp = score_p_fwd;
  score_e_fwd_bp = score_e_fwd;
  std::cout<<"BP Scores"<<std::endl;
  if(vtx_activity == true){
    std::vector<double> bp_scores = get_scores(ssm_sg,break_point,dir);
    score_mu_fwd_bp = bp_scores.at(1);
    score_p_fwd_bp = bp_scores.at(2);
    score_e_fwd_bp = bp_scores.at(3);
  }

  std::cout<<"dqdx"<<std::endl;
  dQ_dx_cut = 0.8866+0.9533 *pow(18/length/units::cm, 0.4234);
  medium_dq_dx = ssm_sg->get_medium_dQ_dx()/(43e3/units::cm);
  medium_dq_dx_bp = medium_dq_dx;
  if(dir==1 && vtx_activity == true){medium_dq_dx_bp = ssm_sg->get_medium_dQ_dx(break_point,ssm_sg->get_dQ_vec().size())/(43e3/units::cm);}
  else if(vtx_activity == true){medium_dq_dx_bp = ssm_sg->get_medium_dQ_dx(0,ssm_sg->get_dQ_vec().size()-break_point)/(43e3/units::cm);}

  //catch and fix cases where we have a large d(dq/dx) for all points
  //if(reduced_muon_length<=0){
  if(break_point==vec_d_dq_dx.size()){
    std::cout<<"false break point"<<std::endl;
    reduced_muon_length = length;
    score_mu_fwd_bp = score_mu_fwd;
    score_p_fwd_bp = score_p_fwd;
    score_e_fwd_bp = score_e_fwd;
    medium_dq_dx_bp = medium_dq_dx;
  }

  //get the energy of the ssm
  TPCParams& mp = Singleton<TPCParams>::Instance();
  TGraph *g_range = 0;
  g_range = mp.get_muon_r2ke();
  kine_energy = g_range->Eval(length) * units::MeV;
  kine_energy_reduced = g_range->Eval(reduced_muon_length) * units::MeV;




  WCPPID::ProtoSegment *prim_track1_sg;
  WCPPID::ProtoSegment *prim_track2_sg;
  WCPPID::ProtoSegment *daught_track1_sg;
  WCPPID::ProtoSegment *daught_track2_sg;
  WCPPID::ProtoSegment *prim_shw1_sg;
  WCPPID::ProtoSegment *prim_shw2_sg;
  WCPPID::ProtoSegment *daught_shw1_sg;
  WCPPID::ProtoSegment *daught_shw2_sg;

  TVector3 dir_prim_track1(0,0,0);
  TVector3 dir_prim_track2(0,0,0);
  TVector3 dir_daught_track1(0,0,0);
  TVector3 dir_daught_track2(0,0,0);
  TVector3 dir_prim_shw1(0,0,0);
  TVector3 dir_prim_shw2(0,0,0);
  TVector3 dir_daught_shw1(0,0,0);
  TVector3 dir_daught_shw2(0,0,0);

  TVector3 mom_prim_track1(0,0,0);
  TVector3 mom_prim_track2(0,0,0);
  TVector3 mom_daught_track1(0,0,0);
  TVector3 mom_daught_track2(0,0,0);
  TVector3 mom_prim_shw1(0,0,0);
  TVector3 mom_prim_shw2(0,0,0);
  TVector3 mom_daught_shw1(0,0,0);
  TVector3 mom_daught_shw2(0,0,0);

std::cout<<"Other particles"<<std::endl;
  //add another loop to fill the prim and duaghter stuff
  for (auto it = map_vertex_segments[main_vertex].begin(); it!= map_vertex_segments[main_vertex].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    if(sg==ssm_sg) continue;
    double sg_length = sg->get_length()/units::cm;
    if(sg_length>1){
      n_prim_all_1 +=1;
      if (!sg->get_flag_shower()){
        n_prim_tracks_1 +=1;
      }
    }
    if(sg_length>3){
      n_prim_all_3 +=1;
      if (!sg->get_flag_shower()){
        n_prim_tracks_3 +=1;
      }
    }
    if(sg_length>5){
      n_prim_all_5 +=1;
      if (!sg->get_flag_shower()){
        n_prim_tracks_5 +=1;
      }
    }
    if(sg_length>8){
      n_prim_all_8 +=1;
      if (!sg->get_flag_shower()){
        n_prim_tracks_8 +=1;
      }
    }
    if(sg_length>11){
      n_prim_all_11 +=1;
      if (!sg->get_flag_shower()){
        n_prim_tracks_11 +=1;
      }
    }
    std::cout<<"sg_length "<<"sg->get_flag_shower() "<<sg->get_flag_shower()<<std::endl;
    double sg_direct_length = sg->get_direct_length()/units::cm;
    double sg_pdg = abs(sg->get_particle_type()); 
    double sg_medium_dq_dx =  sg->get_medium_dQ_dx()/(43e3/units::cm);
    double sg_max_dev = sg->get_max_deviation(int(0),int(sg->get_point_vec().size()))/units::cm;
    double sg_kine_energy_cal = sg->cal_kine_dQdx();
    double sg_kine_energy_range = sg->cal_kine_range();
    std::vector<double> sg_kine_energy_range_pdg = calc_kine_range_multi_pdg(sg_length);
    TVector3 sg_dir = sg->cal_dir_3vector();
    //track with vertex at the muon butt vertex
    if (!sg->get_flag_shower()){
      if(sg_length<length_prim_track2) continue;
      int sg_add_daught_track_counts_1 = calculate_num_daughter_tracks(main_vertex, sg, false, 1*units::cm).first-1;
      if(sg_length<1) sg_add_daught_track_counts_1+=1;//didn't subtract, was too short
      if(sg_add_daught_track_counts_1<0) sg_add_daught_track_counts_1=0;
      int sg_add_daught_all_counts_1 = calculate_num_daughter_tracks(main_vertex, sg, true, 1*units::cm).first-1;
      if(sg_length<1) sg_add_daught_all_counts_1+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_1<0) sg_add_daught_all_counts_1=0;
      int sg_add_daught_track_counts_5 = calculate_num_daughter_tracks(main_vertex, sg, false, 5*units::cm).first-1;
      if(sg_length<5) sg_add_daught_track_counts_5+=1;//didn't subtract, was too short
      if(sg_add_daught_track_counts_5<0) sg_add_daught_track_counts_5=0;
      int sg_add_daught_all_counts_5 = calculate_num_daughter_tracks(main_vertex, sg, true, 5*units::cm).first-1;
      if(sg_length<5) sg_add_daught_all_counts_5+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_5<0) sg_add_daught_all_counts_5=0;
      int sg_add_daught_track_counts_11 = calculate_num_daughter_tracks(main_vertex, sg, false, 11*units::cm).first-1;
      if(sg_length<11) sg_add_daught_track_counts_11+=1;//didn't subtract, was too short
      if(sg_add_daught_track_counts_11<0) sg_add_daught_track_counts_11=0;
      int sg_add_daught_all_counts_11 = calculate_num_daughter_tracks(main_vertex, sg, true, 11*units::cm).first-1;
      if(sg_length<11) sg_add_daught_all_counts_11+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_11<0) sg_add_daught_all_counts_11=0;
      if(sg_length<length_prim_track1){//this is the second longest track
        length_prim_track2 = sg_length;
        direct_length_prim_track2 = sg_direct_length;
        pdg_prim_track2 = sg_pdg;
        medium_dq_dx_prim_track2 = sg_medium_dq_dx;
        max_dev_prim_track2 = sg_max_dev;
        prim_track2_sg = sg;
        kine_energy_cal_prim_track2 = sg_kine_energy_cal;
        kine_energy_range_prim_track2 = sg_kine_energy_range;
	kine_energy_range_mu_prim_track2 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_prim_track2 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_prim_track2 = sg_kine_energy_range_pdg.at(2);
        dir_prim_track2 = sg_dir;
	add_daught_track_counts_1_prim_track2 = sg_add_daught_track_counts_1;
	add_daught_all_counts_1_prim_track2 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_prim_track2 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_prim_track2 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_prim_track2 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_prim_track2 = sg_add_daught_all_counts_11;
      }
      else{//longest track
        length_prim_track2 = length_prim_track1;
        direct_length_prim_track2 = direct_length_prim_track1;
        pdg_prim_track2 = pdg_prim_track1;
        medium_dq_dx_prim_track2 = medium_dq_dx_prim_track1;
        max_dev_prim_track2 = max_dev_prim_track1;
        prim_track2_sg = prim_track1_sg;
        kine_energy_cal_prim_track2 = kine_energy_cal_prim_track1;
        kine_energy_range_prim_track2 = kine_energy_range_prim_track1;
	kine_energy_range_mu_prim_track2 = kine_energy_range_mu_prim_track1;
	kine_energy_range_p_prim_track2 = kine_energy_range_p_prim_track1;
	kine_energy_range_e_prim_track2 = kine_energy_range_e_prim_track1;
        dir_prim_track2 = dir_prim_track1;
        add_daught_track_counts_1_prim_track2 = add_daught_track_counts_1_prim_track1;
        add_daught_all_counts_1_prim_track2 = add_daught_all_counts_1_prim_track1;
        add_daught_track_counts_5_prim_track2 = add_daught_track_counts_5_prim_track1;
        add_daught_all_counts_5_prim_track2 = add_daught_all_counts_5_prim_track1;
        add_daught_track_counts_11_prim_track2 = add_daught_track_counts_11_prim_track1;
        add_daught_all_counts_11_prim_track2 = add_daught_all_counts_11_prim_track1;
	length_prim_track1 = sg_length;
        direct_length_prim_track1 = sg_direct_length;
        pdg_prim_track1 = sg_pdg;
        medium_dq_dx_prim_track1 = sg_medium_dq_dx;
        max_dev_prim_track1 = sg_max_dev;
        prim_track1_sg = sg;
        kine_energy_cal_prim_track1 = sg_kine_energy_cal;
        kine_energy_range_prim_track1 = sg_kine_energy_range;
        kine_energy_range_mu_prim_track1 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_prim_track1 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_prim_track1 = sg_kine_energy_range_pdg.at(2);
	dir_prim_track1 = sg_dir;
	add_daught_track_counts_1_prim_track1 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_prim_track1 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_prim_track1 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_prim_track1 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_prim_track1 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_prim_track1 = sg_add_daught_all_counts_11;
      }
    }
    else{//shower at the muon butt vertex
      if(sg_length<length_prim_shw2) continue;
      int sg_add_daught_track_counts_1 = calculate_num_daughter_tracks(main_vertex, sg, false, 1*units::cm).first;//not a track, don't subtract
      if(sg_add_daught_track_counts_1<0) sg_add_daught_track_counts_1=0;
      int sg_add_daught_all_counts_1 = calculate_num_daughter_tracks(main_vertex, sg, true, 1*units::cm).first-1;
      if(sg_length<1) sg_add_daught_all_counts_1+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_1<0) sg_add_daught_all_counts_1=0;
      int sg_add_daught_track_counts_5 = calculate_num_daughter_tracks(main_vertex, sg, false, 5*units::cm).first;//not a track, don't subtract
      if(sg_add_daught_track_counts_5<0) sg_add_daught_track_counts_5=0;
      int sg_add_daught_all_counts_5 = calculate_num_daughter_tracks(main_vertex, sg, true, 5*units::cm).first-1;
      if(sg_length<5) sg_add_daught_all_counts_5+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_5<0) sg_add_daught_all_counts_5=0;
      int sg_add_daught_track_counts_11 = calculate_num_daughter_tracks(main_vertex, sg, false, 11*units::cm).first;//not a track, don't subtract
      if(sg_add_daught_track_counts_11<0) sg_add_daught_track_counts_11=0;
      int sg_add_daught_all_counts_11 = calculate_num_daughter_tracks(main_vertex, sg, true, 11*units::cm).first-1;
      if(sg_length<11) sg_add_daught_all_counts_11+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_11<0) sg_add_daught_all_counts_11=0;
      if(sg_length<length_prim_shw1){//this is the second longest shw
        length_prim_shw2 = sg_length;
        direct_length_prim_shw2 = sg_direct_length;
        pdg_prim_shw2 = sg_pdg; 
	medium_dq_dx_prim_shw2 = sg_medium_dq_dx;
        max_dev_prim_shw2 = sg_max_dev;
        prim_shw2_sg = sg;
        kine_energy_cal_prim_shw2 = sg_kine_energy_cal;
        kine_energy_range_prim_shw2 = sg_kine_energy_range;
        kine_energy_range_mu_prim_shw2 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_prim_shw2 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_prim_shw2 = sg_kine_energy_range_pdg.at(2);
	dir_prim_shw2 = sg_dir;
	add_daught_track_counts_1_prim_shw2 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_prim_shw2 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_prim_shw2 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_prim_shw2 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_prim_shw2 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_prim_shw2 = sg_add_daught_all_counts_11;
      }
      else{//longest shw
        length_prim_shw2 = length_prim_shw1;
        direct_length_prim_shw2 = direct_length_prim_shw1;
        pdg_prim_shw2 = pdg_prim_shw1;
        medium_dq_dx_prim_shw2 = medium_dq_dx_prim_shw1;
        max_dev_prim_shw2 = max_dev_prim_shw1;
        prim_shw2_sg = prim_shw1_sg;
        kine_energy_cal_prim_shw2 = kine_energy_cal_prim_shw1;
        kine_energy_range_prim_shw2 = kine_energy_range_prim_shw1;
        kine_energy_range_mu_prim_shw2 = kine_energy_range_mu_prim_shw1;
        kine_energy_range_p_prim_shw2 = kine_energy_range_p_prim_shw1;
        kine_energy_range_e_prim_shw2 = kine_energy_range_e_prim_shw1;
	dir_prim_shw2 = dir_prim_shw1;
	add_daught_track_counts_1_prim_shw2 = add_daught_track_counts_1_prim_shw1;
        add_daught_all_counts_1_prim_shw2 = add_daught_all_counts_1_prim_shw1;
        add_daught_track_counts_5_prim_shw2 = add_daught_track_counts_5_prim_shw1;
        add_daught_all_counts_5_prim_shw2 = add_daught_all_counts_5_prim_shw1;
        add_daught_track_counts_11_prim_shw2 = add_daught_track_counts_11_prim_shw1;
        add_daught_all_counts_11_prim_shw2 = add_daught_all_counts_11_prim_shw1;
        length_prim_shw1 = sg_length;
        direct_length_prim_shw1 = sg_direct_length;
        pdg_prim_shw1 = sg_pdg;
        medium_dq_dx_prim_shw1 = sg_medium_dq_dx;
        max_dev_prim_shw1 = sg_max_dev;
        prim_shw1_sg = sg;
        kine_energy_cal_prim_shw1 = sg_kine_energy_cal;
        kine_energy_range_prim_shw1 = sg_kine_energy_range;
        kine_energy_range_mu_prim_shw1 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_prim_shw1 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_prim_shw1 = sg_kine_energy_range_pdg.at(2);
	dir_prim_shw1 = sg_dir;
	add_daught_track_counts_1_prim_shw1 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_prim_shw1 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_prim_shw1 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_prim_shw1 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_prim_shw1 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_prim_shw1 = sg_add_daught_all_counts_11;
      }
    }
  }
std::cout<<"Other particles at the butt"<<std::endl;
  //loop the vertex corresponding to the muon butt
  WCPPID::ProtoVertex* second_vtx = find_other_vertex(ssm_sg,main_vertex);
  for (auto it = map_vertex_segments[second_vtx].begin(); it!= map_vertex_segments[second_vtx].end(); it++){
    WCPPID::ProtoSegment *sg = *it;
    if(sg==ssm_sg) continue;
    double sg_length = sg->get_length()/units::cm;
    if(sg_length>1){
      n_daughter_all_1 +=1;
      if (!sg->get_flag_shower()){
        n_daughter_tracks_1 +=1;
      }
    }
    if(sg_length>3){
      n_daughter_all_3 +=1;
      if (!sg->get_flag_shower()){
        n_daughter_tracks_3 +=1;
      }
    }
    if(sg_length>5){
      n_daughter_all_5 +=1;
      if (!sg->get_flag_shower()){
        n_daughter_tracks_5 +=1;
      }
    }
    if(sg_length>8){
      n_daughter_all_8 +=1;
      if (!sg->get_flag_shower()){
        n_daughter_tracks_8 +=1;
      }
    }
    if(sg_length>11){
      n_daughter_all_11 +=1;
      if (!sg->get_flag_shower()){
        n_daughter_tracks_11 +=1;
      }
    }
std::cout<<"sg_length "<<"sg->get_flag_shower() "<<sg->get_flag_shower()<<std::endl;
    double sg_direct_length = sg->get_direct_length()/units::cm;
    double sg_pdg = abs(sg->get_particle_type());
    double sg_medium_dq_dx =  sg->get_medium_dQ_dx()/(43e3/units::cm);
    double sg_max_dev = sg->get_max_deviation(int(0),int(sg->get_point_vec().size()))/units::cm;
    double sg_kine_energy_cal = sg->cal_kine_dQdx();
    double sg_kine_energy_range = sg->cal_kine_range();
    std::vector<double> sg_kine_energy_range_pdg = calc_kine_range_multi_pdg(sg_length);
    TVector3 sg_dir = sg->cal_dir_3vector();
    //track with vertex at the muon butt vertex
    if (!sg->get_flag_shower()){
      if(sg_length<length_daught_track2) continue;
      int sg_add_daught_track_counts_1 = calculate_num_daughter_tracks(second_vtx, sg, false, 1*units::cm).first-1;
      if(sg_length<1) sg_add_daught_track_counts_1+=1;//didn't subtract, was too short
      if(sg_add_daught_track_counts_1<0) sg_add_daught_track_counts_1=0;
      int sg_add_daught_all_counts_1 = calculate_num_daughter_tracks(second_vtx, sg, true, 1*units::cm).first-1;
      if(sg_length<1) sg_add_daught_all_counts_1+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_1<0) sg_add_daught_all_counts_1=0;
      int sg_add_daught_track_counts_5 = calculate_num_daughter_tracks(second_vtx, sg, false, 5*units::cm).first-1;
      if(sg_length<5) sg_add_daught_track_counts_5+=1;//didn't subtract, was too short
      if(sg_add_daught_track_counts_5<0) sg_add_daught_track_counts_5=0;
      int sg_add_daught_all_counts_5 = calculate_num_daughter_tracks(second_vtx, sg, true, 5*units::cm).first-1;
      if(sg_length<5) sg_add_daught_all_counts_5+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_5<0) sg_add_daught_all_counts_5=0;
      int sg_add_daught_track_counts_11 = calculate_num_daughter_tracks(second_vtx, sg, false, 11*units::cm).first-1;
      if(sg_length<11) sg_add_daught_track_counts_11+=1;//didn't subtract, was too short
      if(sg_add_daught_track_counts_11<0) sg_add_daught_track_counts_11=0;
      int sg_add_daught_all_counts_11 = calculate_num_daughter_tracks(second_vtx, sg, true, 11*units::cm).first-1;
      if(sg_length<11) sg_add_daught_all_counts_11+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_11<0) sg_add_daught_all_counts_11=0;
      if(sg_length<length_daught_track1){//this is the second longest track
        length_daught_track2 = sg_length;
        direct_length_daught_track2 = sg_direct_length;
        pdg_daught_track2 = sg_pdg;
        medium_dq_dx_daught_track2 = sg_medium_dq_dx;
        max_dev_daught_track2 = sg_max_dev;
        daught_track2_sg = sg;
        kine_energy_cal_daught_track2 = sg_kine_energy_cal;
	kine_energy_range_daught_track2 = sg_kine_energy_range;
	kine_energy_range_mu_daught_track2 = sg_kine_energy_range_pdg.at(0);
	kine_energy_range_p_daught_track2 = sg_kine_energy_range_pdg.at(1);
	kine_energy_range_e_daught_track2 = sg_kine_energy_range_pdg.at(2);
	dir_daught_track2 = sg_dir;
	add_daught_track_counts_1_daught_track2 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_daught_track2 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_daught_track2 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_daught_track2 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_daught_track2 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_daught_track2 = sg_add_daught_all_counts_11;
      }
      else{//longest track
        length_daught_track2 = length_daught_track1;
        direct_length_daught_track2 = direct_length_daught_track1;
        pdg_daught_track2 = pdg_daught_track1;
        medium_dq_dx_daught_track2 = medium_dq_dx_daught_track1;
        max_dev_daught_track2 = max_dev_daught_track1;
        daught_track2_sg = daught_track1_sg;
	kine_energy_cal_daught_track2 = kine_energy_cal_daught_track1;
	kine_energy_range_daught_track2 = kine_energy_range_daught_track1;
        kine_energy_range_mu_daught_track2 = kine_energy_range_mu_daught_track1;
        kine_energy_range_p_daught_track2 = kine_energy_range_p_daught_track1;
        kine_energy_range_e_daught_track2 = kine_energy_range_e_daught_track1;
	dir_daught_track2 = dir_daught_track1;
	add_daught_track_counts_1_daught_track2 = add_daught_track_counts_1_daught_track1;
        add_daught_all_counts_1_daught_track2 = add_daught_all_counts_1_daught_track1;
        add_daught_track_counts_5_daught_track2 = add_daught_track_counts_5_daught_track1;
        add_daught_all_counts_5_daught_track2 = add_daught_all_counts_5_daught_track1;
        add_daught_track_counts_11_daught_track2 = add_daught_track_counts_11_daught_track1;
        add_daught_all_counts_11_daught_track2 = add_daught_all_counts_11_daught_track1;
        length_daught_track1 = sg_length;
        direct_length_daught_track1 = sg_direct_length;
        pdg_daught_track1 = sg_pdg;
        medium_dq_dx_daught_track1 = sg_medium_dq_dx;
        max_dev_daught_track1 = sg_max_dev;  
        daught_track1_sg = sg;
	kine_energy_cal_daught_track1 = sg_kine_energy_cal;
	kine_energy_range_daught_track1 = sg_kine_energy_range;
        kine_energy_range_mu_daught_track1 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_daught_track1 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_daught_track1 = sg_kine_energy_range_pdg.at(2);
	dir_daught_track1 = sg_dir;
	add_daught_track_counts_1_daught_track1 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_daught_track1 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_daught_track1 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_daught_track1 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_daught_track1 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_daught_track1 = sg_add_daught_all_counts_11;
      }
    }
    else{//shower at the muon butt vertex
      if(sg_length<length_daught_shw2) continue;
      int sg_add_daught_track_counts_1 = calculate_num_daughter_tracks(second_vtx, sg, false, 1*units::cm).first;//not a track, don't subtract
      if(sg_add_daught_track_counts_1<0) sg_add_daught_track_counts_1=0;
      int sg_add_daught_all_counts_1 = calculate_num_daughter_tracks(second_vtx, sg, true, 1*units::cm).first-1;
      if(sg_length<1) sg_add_daught_all_counts_1+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_1<0) sg_add_daught_all_counts_1=0;
      int sg_add_daught_track_counts_5 = calculate_num_daughter_tracks(second_vtx, sg, false, 5*units::cm).first;//not a track, don't subtract
      if(sg_add_daught_track_counts_5<0) sg_add_daught_track_counts_5=0;
      int sg_add_daught_all_counts_5 = calculate_num_daughter_tracks(second_vtx, sg, true, 5*units::cm).first-1;
      if(sg_length<5) sg_add_daught_all_counts_5+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_5<0) sg_add_daught_all_counts_5=0;
      int sg_add_daught_track_counts_11 = calculate_num_daughter_tracks(second_vtx, sg, false, 11*units::cm).first;//not a track, don't subtract
      if(sg_add_daught_track_counts_11<0) sg_add_daught_track_counts_11=0;
      int sg_add_daught_all_counts_11 = calculate_num_daughter_tracks(second_vtx, sg, true, 11*units::cm).first-1;
      if(sg_length<11) sg_add_daught_all_counts_11+=1;//didn't subtract, was too short
      if(sg_add_daught_all_counts_11<0) sg_add_daught_all_counts_11=0;
      if(sg_length<length_daught_shw1){//this is the second longest shw
        length_daught_shw2 = sg_length;
        direct_length_daught_shw2 = sg_direct_length;
        pdg_daught_shw2 = sg_pdg;
        medium_dq_dx_daught_shw2 = sg_medium_dq_dx;
        max_dev_daught_shw2 = sg_max_dev;
        daught_shw2_sg = sg;
        kine_energy_cal_daught_shw2 = sg_kine_energy_cal;
        kine_energy_range_daught_shw2 = sg_kine_energy_range;
	kine_energy_range_mu_daught_shw2 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_daught_shw2 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_daught_shw2 = sg_kine_energy_range_pdg.at(2);
        dir_daught_shw2 = sg_dir;
	add_daught_track_counts_1_daught_shw2 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_daught_shw2 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_daught_shw2 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_daught_shw2 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_daught_shw2 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_daught_shw2 = sg_add_daught_all_counts_11;
      }
      else{//longest shw
        length_daught_shw2 = length_daught_shw1;
        direct_length_daught_shw2 = direct_length_daught_shw1;
        pdg_daught_shw2 = pdg_daught_shw1;
        medium_dq_dx_daught_shw2 = medium_dq_dx_daught_shw1;
        max_dev_daught_shw2 = max_dev_daught_shw1;
        daught_shw2_sg = daught_shw1_sg;
        kine_energy_cal_daught_shw2 = kine_energy_cal_daught_shw1;
        kine_energy_range_daught_shw2 = kine_energy_range_daught_shw1;
	kine_energy_range_mu_daught_shw2 = kine_energy_range_mu_daught_shw1;
        kine_energy_range_p_daught_shw2 = kine_energy_range_p_daught_shw1;
        kine_energy_range_e_daught_shw2 = kine_energy_range_e_daught_shw1;
        dir_daught_shw2 = dir_daught_shw1;
	add_daught_track_counts_1_daught_shw2 = add_daught_track_counts_1_daught_shw1;
        add_daught_all_counts_1_daught_shw2 = add_daught_all_counts_1_daught_shw1;
        add_daught_track_counts_5_daught_shw2 = add_daught_track_counts_5_daught_shw1;
        add_daught_all_counts_5_daught_shw2 = add_daught_all_counts_5_daught_shw1;
        add_daught_track_counts_11_daught_shw2 = add_daught_track_counts_11_daught_shw1;
        add_daught_all_counts_11_daught_shw2 = add_daught_all_counts_11_daught_shw1;
        length_daught_shw1 = sg_length;
        direct_length_daught_shw1 = sg_direct_length;
        pdg_daught_shw1 = sg_pdg;
        medium_dq_dx_daught_shw1 = sg_medium_dq_dx;
        max_dev_daught_shw1 = sg_max_dev;        
        daught_shw1_sg = sg;
        kine_energy_cal_daught_shw1 = sg_kine_energy_cal;
        kine_energy_range_daught_shw1 = sg_kine_energy_range;
        kine_energy_range_mu_daught_shw1 = sg_kine_energy_range_pdg.at(0);
        kine_energy_range_p_daught_shw1 = sg_kine_energy_range_pdg.at(1);
        kine_energy_range_e_daught_shw1 = sg_kine_energy_range_pdg.at(2);
	dir_daught_shw1 = sg_dir;
	add_daught_track_counts_1_daught_shw1 = sg_add_daught_track_counts_1;
        add_daught_all_counts_1_daught_shw1 = sg_add_daught_all_counts_1;
        add_daught_track_counts_5_daught_shw1 = sg_add_daught_track_counts_5;
        add_daught_all_counts_5_daught_shw1 = sg_add_daught_all_counts_5;
        add_daught_track_counts_11_daught_shw1 = sg_add_daught_track_counts_11;
        add_daught_all_counts_11_daught_shw1 = sg_add_daught_all_counts_11;
      }
    }
  }
    std::cout<<"Other particles scores"<<std::endl;
  if(length_prim_track1>0){
    std::vector<double> scores_prim_track1 = get_scores(prim_track1_sg);
    score_mu_fwd_prim_track1 = scores_prim_track1.at(0);
    score_p_fwd_prim_track1 = scores_prim_track1.at(1);
    score_e_fwd_prim_track1 = scores_prim_track1.at(2);
    score_mu_bck_prim_track1 = scores_prim_track1.at(3);
    score_p_bck_prim_track1 = scores_prim_track1.at(4);
    score_e_bck_prim_track1 = scores_prim_track1.at(5);
  }
  if(length_prim_track2>0){
    std::vector<double> scores_prim_track2 = get_scores(prim_track2_sg);
    score_mu_fwd_prim_track2 = scores_prim_track2.at(0);
    score_p_fwd_prim_track2 = scores_prim_track2.at(1);
    score_e_fwd_prim_track2 = scores_prim_track2.at(2);
    score_mu_bck_prim_track2 = scores_prim_track2.at(3);
    score_p_bck_prim_track2 = scores_prim_track2.at(4);
    score_e_bck_prim_track2 = scores_prim_track2.at(5);
  }
  if(length_daught_track1>0){
    std::vector<double> scores_daught_track1 = get_scores(daught_track1_sg);
    score_mu_fwd_daught_track1 = scores_daught_track1.at(0);
    score_p_fwd_daught_track1 = scores_daught_track1.at(1);
    score_e_fwd_daught_track1 = scores_daught_track1.at(2);
    score_mu_bck_daught_track1 = scores_daught_track1.at(3);
    score_p_bck_daught_track1 = scores_daught_track1.at(4);
    score_e_bck_daught_track1 = scores_daught_track1.at(5);
  }
  std::cout<<"daught track2"<<std::endl;
  if(length_daught_track2>0){
  std::vector<double> scores_daught_track2 = get_scores(daught_track2_sg);
  score_mu_fwd_daught_track2 = scores_daught_track2.at(0);
  score_p_fwd_daught_track2 = scores_daught_track2.at(1);
  score_e_fwd_daught_track2 = scores_daught_track2.at(2);
  score_mu_bck_daught_track2 = scores_daught_track2.at(3);
  score_p_bck_daught_track2 = scores_daught_track2.at(4);
  score_e_bck_daught_track2 = scores_daught_track2.at(5);
  }
  if(length_prim_shw1>0){
    std::vector<double> scores_prim_shw1 = get_scores(prim_shw1_sg);
    score_mu_fwd_prim_shw1 = scores_prim_shw1.at(0);
    score_p_fwd_prim_shw1 = scores_prim_shw1.at(1);
    score_e_fwd_prim_shw1 = scores_prim_shw1.at(2);
    score_mu_bck_prim_shw1 = scores_prim_shw1.at(3);
    score_p_bck_prim_shw1 = scores_prim_shw1.at(4);
    score_e_bck_prim_shw1 = scores_prim_shw1.at(5);
  }
  if(length_prim_shw2>0){
    std::vector<double> scores_prim_shw2 = get_scores(prim_shw2_sg);
    score_mu_fwd_prim_shw2 = scores_prim_shw2.at(0);
    score_p_fwd_prim_shw2 = scores_prim_shw2.at(1);
    score_e_fwd_prim_shw2 = scores_prim_shw2.at(2);
    score_mu_bck_prim_shw2 = scores_prim_shw2.at(3);
    score_p_bck_prim_shw2 = scores_prim_shw2.at(4);
    score_e_bck_prim_shw2 = scores_prim_shw2.at(5);
  }
  if(length_daught_shw1>0){
    std::vector<double> scores_daught_shw1 = get_scores(daught_shw1_sg);
    score_mu_fwd_daught_shw1 = scores_daught_shw1.at(0);
    score_p_fwd_daught_shw1 = scores_daught_shw1.at(1);
    score_e_fwd_daught_shw1 = scores_daught_shw1.at(2);
    score_mu_bck_daught_shw1 = scores_daught_shw1.at(3);
    score_p_bck_daught_shw1 = scores_daught_shw1.at(4);
    score_e_bck_daught_shw1 = scores_daught_shw1.at(5);
  }
  if(length_daught_shw2>0){
    std::vector<double> scores_daught_shw2 = get_scores(daught_shw2_sg);
    score_mu_fwd_daught_shw2 = scores_daught_shw2.at(0);
    score_p_fwd_daught_shw2 = scores_daught_shw2.at(1);
    score_e_fwd_daught_shw2 = scores_daught_shw2.at(2);
    score_mu_bck_daught_shw2 = scores_daught_shw2.at(3);
    score_p_bck_daught_shw2 = scores_daught_shw2.at(4);
    score_e_bck_daught_shw2 = scores_daught_shw2.at(5);
  }

    //figure our condition needed to swap, for now just daught and prim if muon is backwards
    //will be more complex if we look for off vertex stuff, here fwd/bckg may als switch if the muon in non-prim
std::cout<<"backwards_muon "<<backwards_muon<<"dir "<<dir<<std::endl;
  if(backwards_muon){
  //if(dir==-1){
  std::cout<<"backwards_muon "<<std::endl;
     std::swap(pdg_prim_track1, pdg_daught_track1);
     std::swap(length_prim_track1, length_daught_track1);
     std::swap(direct_length_prim_track1, direct_length_daught_track1);
     std::swap(max_dev_prim_track1, max_dev_daught_track1);
     std::swap(kine_energy_range_prim_track1, kine_energy_range_daught_track1);
     std::swap(kine_energy_range_mu_prim_track1, kine_energy_range_mu_daught_track1);
     std::swap(kine_energy_range_p_prim_track1, kine_energy_range_p_daught_track1);
     std::swap(kine_energy_range_e_prim_track1, kine_energy_range_e_daught_track1);
     std::swap(kine_energy_cal_prim_track1, kine_energy_cal_daught_track1);
     std::swap(medium_dq_dx_prim_track1, medium_dq_dx_daught_track1);
     std::swap(dir_prim_track1, dir_daught_track1);
  std::cout<<"1"<<std::endl;
    std::swap(n_prim_all_1, n_daughter_all_1);
    std::swap(n_prim_tracks_1, n_daughter_tracks_1);
    std::swap(n_prim_all_3, n_daughter_all_3);
    std::swap(n_prim_tracks_3, n_daughter_tracks_3);
    std::swap(n_prim_all_5, n_daughter_all_5);
    std::swap(n_prim_tracks_5, n_daughter_tracks_5);
    std::swap(n_prim_all_8, n_daughter_all_8);
    std::swap(n_prim_tracks_8, n_daughter_tracks_8);
    std::swap(n_prim_all_11, n_daughter_all_11);
    std::swap(n_prim_tracks_11, n_daughter_tracks_11);
  std::cout<<"2"<<std::endl;
    std::swap(add_daught_track_counts_1_prim_track1, add_daught_track_counts_1_daught_track1);
    std::swap(add_daught_all_counts_1_prim_track1, add_daught_all_counts_1_daught_track1);
    std::swap(add_daught_track_counts_1_prim_track2, add_daught_track_counts_1_daught_track2);
    std::swap(add_daught_all_counts_1_prim_track2, add_daught_all_counts_1_daught_track2);
    std::swap(add_daught_track_counts_1_prim_shw1, add_daught_track_counts_1_daught_shw1);
    std::swap(add_daught_all_counts_1_prim_shw1, add_daught_all_counts_1_daught_shw1);
    std::swap(add_daught_track_counts_1_prim_shw2, add_daught_track_counts_1_daught_shw2);
    std::swap(add_daught_all_counts_1_prim_shw2, add_daught_all_counts_1_daught_shw2);
  std::cout<<"3"<<std::endl;
    std::swap(add_daught_track_counts_5_prim_track1, add_daught_track_counts_5_daught_track1);
    std::swap(add_daught_all_counts_5_prim_track1, add_daught_all_counts_5_daught_track1);
    std::swap(add_daught_track_counts_5_prim_track2, add_daught_track_counts_5_daught_track2);
    std::swap(add_daught_all_counts_5_prim_track2, add_daught_all_counts_5_daught_track2);
    std::swap(add_daught_track_counts_5_prim_shw1, add_daught_track_counts_5_daught_shw1);
    std::swap(add_daught_all_counts_5_prim_shw1, add_daught_all_counts_5_daught_shw1);
    std::swap(add_daught_track_counts_5_prim_shw2, add_daught_track_counts_5_daught_shw2);
    std::swap(add_daught_all_counts_5_prim_shw2, add_daught_all_counts_5_daught_shw2);
  std::cout<<"4"<<std::endl;
    std::swap(add_daught_track_counts_11_prim_track1, add_daught_track_counts_11_daught_track1);
    std::swap(add_daught_all_counts_11_prim_track1, add_daught_all_counts_11_daught_track1);
    std::swap(add_daught_track_counts_11_prim_track2, add_daught_track_counts_11_daught_track2);
    std::swap(add_daught_all_counts_11_prim_track2, add_daught_all_counts_11_daught_track2);
    std::swap(add_daught_track_counts_11_prim_shw1, add_daught_track_counts_11_daught_shw1);
    std::swap(add_daught_all_counts_11_prim_shw1, add_daught_all_counts_11_daught_shw1);
    std::swap(add_daught_track_counts_11_prim_shw2, add_daught_track_counts_11_daught_shw2);
    std::swap(add_daught_all_counts_11_prim_shw2, add_daught_all_counts_11_daught_shw2);
  std::cout<<"5"<<std::endl;
     std::swap(pdg_prim_track2, pdg_daught_track2);
     std::swap(length_prim_track2, length_daught_track2);
     std::swap(direct_length_prim_track2, direct_length_daught_track2);
     std::swap(max_dev_prim_track2, max_dev_daught_track2);
     std::swap(kine_energy_range_prim_track2, kine_energy_range_daught_track2);
     std::swap(kine_energy_range_mu_prim_track2, kine_energy_range_mu_daught_track2);
     std::swap(kine_energy_range_p_prim_track2, kine_energy_range_p_daught_track2);
     std::swap(kine_energy_range_e_prim_track2, kine_energy_range_e_daught_track2);
     std::swap(kine_energy_cal_prim_track2, kine_energy_cal_daught_track2);
     std::swap(medium_dq_dx_prim_track2, medium_dq_dx_daught_track2);
     std::swap(dir_prim_track2, dir_daught_track2);
  std::cout<<"6"<<std::endl;
     std::swap(pdg_prim_shw1, pdg_daught_shw1);
     std::swap(length_prim_shw1, length_daught_shw1);
     std::swap(direct_length_prim_shw1, direct_length_daught_shw1);
     std::swap(max_dev_prim_shw1, max_dev_daught_shw1);
     std::swap(kine_energy_range_prim_shw1, kine_energy_range_daught_shw1);
     std::swap(kine_energy_range_mu_prim_shw1, kine_energy_range_mu_daught_shw1);
     std::swap(kine_energy_range_p_prim_shw1, kine_energy_range_p_daught_shw1);
     std::swap(kine_energy_range_e_prim_shw1, kine_energy_range_e_daught_shw1);
     std::swap(kine_energy_cal_prim_shw1, kine_energy_cal_daught_shw1);
     std::swap(medium_dq_dx_prim_shw1, medium_dq_dx_daught_shw1);
     std::swap(dir_prim_shw1, dir_daught_shw1);
  std::cout<<"7"<<std::endl;
     std::swap(pdg_prim_shw2, pdg_daught_shw2);
     std::swap(length_prim_shw2, length_daught_shw2);
     std::swap(direct_length_prim_shw2, direct_length_daught_shw2);
     std::swap(max_dev_prim_shw2, max_dev_daught_shw2);
     std::swap(kine_energy_range_prim_shw2, kine_energy_range_daught_shw2);
     std::swap(kine_energy_range_mu_prim_shw2, kine_energy_range_mu_daught_shw2);
     std::swap(kine_energy_range_p_prim_shw2, kine_energy_range_p_daught_shw2);
     std::swap(kine_energy_range_e_prim_shw2, kine_energy_range_e_daught_shw2);
     std::swap(kine_energy_cal_prim_shw2, kine_energy_cal_daught_shw2);
     std::swap(medium_dq_dx_prim_shw2, medium_dq_dx_daught_shw2);
     std::swap(dir_prim_shw2, dir_daught_shw2);
  std::cout<<"8"<<std::endl;
    std::swap(score_mu_fwd_prim_track1, score_mu_fwd_daught_track1);  
    std::swap(score_p_fwd_prim_track1, score_p_fwd_daught_track1); 
    std::swap(score_e_fwd_prim_track1, score_e_fwd_daught_track1);
    std::swap(score_mu_bck_prim_track1, score_mu_bck_daught_track1); 
    std::swap(score_p_bck_prim_track1, score_p_bck_daught_track1);
    std::swap(score_e_bck_prim_track1, score_e_bck_daught_track1);
  std::cout<<"9"<<std::endl;
    std::swap(score_mu_fwd_prim_track2, score_mu_fwd_daught_track2); 
    std::swap(score_p_fwd_prim_track2, score_p_fwd_daught_track2); 
    std::swap(score_e_fwd_prim_track2, score_e_fwd_daught_track2);
    std::swap(score_mu_bck_prim_track2, score_mu_bck_daught_track2); 
    std::swap(score_p_bck_prim_track2, score_p_bck_daught_track2);
    std::swap(score_e_bck_prim_track2, score_e_bck_daught_track2);
  std::cout<<"10"<<std::endl;
    std::swap(score_mu_fwd_prim_shw1, score_mu_fwd_daught_shw1);
    std::swap(score_p_fwd_prim_shw1, score_p_fwd_daught_shw1);
    std::swap(score_e_fwd_prim_shw1, score_e_fwd_daught_shw1);
    std::swap(score_mu_bck_prim_shw1, score_mu_bck_daught_shw1);
    std::swap(score_p_bck_prim_shw1, score_p_bck_daught_shw1); 
    std::swap(score_e_bck_prim_shw1, score_e_bck_daught_shw1);
  std::cout<<"11"<<std::endl;
    std::swap(score_mu_fwd_prim_shw2, score_mu_fwd_daught_shw2);
    std::swap(score_p_fwd_prim_shw2, score_p_fwd_daught_shw2);
    std::swap(score_e_fwd_prim_shw2, score_e_fwd_daught_shw2);
    std::swap(score_mu_bck_prim_shw2, score_mu_bck_daught_shw2);
    std::swap(score_p_bck_prim_shw2, score_p_bck_daught_shw2);
    std::swap(score_e_bck_prim_shw2, score_e_bck_daught_shw2);
  std::cout<<"12"<<std::endl;
  std::swap(prim_track1_sg,daught_track1_sg);
  std::swap(prim_track2_sg,daught_track2_sg);
  std::swap(prim_shw1_sg,daught_shw1_sg);
  std::swap(prim_shw2_sg,daught_shw2_sg);
  std::cout<<"13"<<std::endl;
  }

  if(length_prim_track1>0){
  std::cout<<"length_prim_track1 "<<length_prim_track1<<std::endl;
  double mass_prim_track1 = prim_track1_sg->get_particle_mass();
  std::cout<<"mass_prim_track1 "<<mass_prim_track1<<std::endl;
  double mom_mag_prim_track1 = sqrt(pow(kine_energy_range_prim_track1,2) + 2*kine_energy_range_prim_track1*mass_prim_track1);
  std::cout<<"mom_mag_prim_track1 "<<mom_mag_prim_track1<<std::endl;
  x_dir_prim_track1 = dir_prim_track1[0];
  y_dir_prim_track1 = dir_prim_track1[1];
  z_dir_prim_track1 = dir_prim_track1[2];
  mom_prim_track1[0] = dir_prim_track1[0]*mom_mag_prim_track1;
  mom_prim_track1[1] = dir_prim_track1[1]*mom_mag_prim_track1;
  mom_prim_track1[2] = dir_prim_track1[2]*mom_mag_prim_track1;
  }
  if(length_prim_track2>0){
  double mass_prim_track2 = prim_track2_sg->get_particle_mass();
  std::cout<<"mass_prim_track2 "<<mass_prim_track2<<std::endl;
  double mom_mag_prim_track2 = sqrt(pow(kine_energy_range_prim_track2,2) + 2*kine_energy_range_prim_track2*mass_prim_track2);
  x_dir_prim_track2 = dir_prim_track2[0];
  y_dir_prim_track2 = dir_prim_track2[1];
  z_dir_prim_track2 = dir_prim_track2[2];
  mom_prim_track2[0] = dir_prim_track2[0]*mom_mag_prim_track2;
  mom_prim_track2[1] = dir_prim_track2[1]*mom_mag_prim_track2;
  mom_prim_track2[2] = dir_prim_track2[2]*mom_mag_prim_track2;
  }
  if(length_daught_track1>0){
  double mass_daught_track1 = daught_track1_sg->get_particle_mass();
  std::cout<<"mass_daught_track1 "<<mass_daught_track1<<std::endl;
  double mom_mag_daught_track1 = sqrt(pow(kine_energy_range_daught_track1,2) + 2*kine_energy_range_daught_track1*mass_daught_track1);
  x_dir_daught_track1 = dir_prim_track1[0];
  y_dir_daught_track1 = dir_prim_track1[1];
  z_dir_daught_track1 = dir_prim_track1[2];
  mom_daught_track1[0] = dir_daught_track1[0]*mom_mag_daught_track1;
  mom_daught_track1[1] = dir_daught_track1[1]*mom_mag_daught_track1;
  mom_daught_track1[2] = dir_daught_track1[2]*mom_mag_daught_track1;
  }
  if(length_daught_track2>0){
  double mass_daught_track2 = daught_track2_sg->get_particle_mass();
  std::cout<<"mass_daught_track2 "<<mass_daught_track2<<std::endl;
  double mom_mag_daught_track2 = sqrt(pow(kine_energy_range_daught_track2,2) + 2*kine_energy_range_daught_track2*mass_daught_track2);
  x_dir_daught_track2 = dir_daught_track2[0];
  y_dir_daught_track2 = dir_daught_track2[1];
  z_dir_daught_track2 = dir_daught_track2[2];
  mom_daught_track2[0] = dir_daught_track2[0]*mom_mag_daught_track2;
  mom_daught_track2[1] = dir_daught_track2[1]*mom_mag_daught_track2;
  mom_daught_track2[2] = dir_daught_track2[2]*mom_mag_daught_track2;
  }
  if(length_prim_shw1>0){
  double mass_prim_shw1 = prim_shw1_sg->get_particle_mass();
  std::cout<<"mass_prim_shw1 "<<mass_prim_shw1<<std::endl; 
 double mom_mag_prim_shw1 = sqrt(pow(kine_energy_range_prim_shw1,2) + 2*kine_energy_range_prim_shw1*mass_prim_shw1);
  x_dir_prim_shw1 = dir_prim_shw1[0];
  y_dir_prim_shw1 = dir_prim_shw1[1];
  z_dir_prim_shw1 = dir_prim_shw1[2];
  mom_prim_shw1[0] = dir_prim_shw1[0]*mom_mag_prim_shw1;
  mom_prim_shw1[1] = dir_prim_shw1[1]*mom_mag_prim_shw1;
  mom_prim_shw1[2] = dir_prim_shw1[2]*mom_mag_prim_shw1;
  }
  if(length_prim_shw2>0){
  double mass_prim_shw2 = prim_shw2_sg->get_particle_mass();
  std::cout<<"mass_prim_shw2 "<<mass_prim_shw2<<std::endl;
  double mom_mag_prim_shw2 = sqrt(pow(kine_energy_range_prim_shw2,2) + 2*kine_energy_range_prim_shw2*mass_prim_shw2);
  x_dir_prim_shw2 = dir_prim_shw2[0];
  y_dir_prim_shw2 = dir_prim_shw2[1];
  z_dir_prim_shw2 = dir_prim_shw2[2];
  mom_prim_shw2[0] = dir_prim_shw2[0]*mom_mag_prim_shw2;
  mom_prim_shw2[1] = dir_prim_shw2[1]*mom_mag_prim_shw2;
  mom_prim_shw2[2] = dir_prim_shw2[2]*mom_mag_prim_shw2;
  }
  if(length_daught_shw1>0){
  double mass_daught_shw1 = daught_shw1_sg->get_particle_mass();
  std::cout<<"mass_daught_shw1 "<<mass_daught_shw1<<std::endl;
  double mom_mag_daught_shw1 = sqrt(pow(kine_energy_range_daught_shw1,2) + 2*kine_energy_range_daught_shw1*mass_daught_shw1);
  x_dir_daught_shw1 = dir_daught_shw1[0];
  y_dir_daught_shw1 = dir_daught_shw1[1];
  z_dir_daught_shw1 = dir_daught_shw1[2];
  mom_daught_shw1[0] = dir_daught_shw1[0]*mom_mag_daught_shw1;
  mom_daught_shw1[1] = dir_daught_shw1[1]*mom_mag_daught_shw1;
  mom_daught_shw1[2] = dir_daught_shw1[2]*mom_mag_daught_shw1;
  }
  if(length_daught_shw2>0){
  std::cout<<"length_daught_shw2 "<<length_daught_shw2<<std::endl;
  double mass_daught_shw2 = daught_shw2_sg->get_particle_mass();
  std::cout<<"mass_daught_shw2 "<<mass_daught_shw2<<std::endl;
  double mom_mag_daught_shw2 = sqrt(pow(kine_energy_range_daught_shw2,2) + 2*kine_energy_range_daught_shw2*mass_daught_shw2);
  std::cout<<"mom_mag_daught_shw2 "<<mom_mag_daught_shw2<<std::endl;
  x_dir_daught_shw2 = dir_daught_shw2[0];
  y_dir_daught_shw2 = dir_daught_shw2[1];
  z_dir_daught_shw2 = dir_daught_shw2[2];
  mom_daught_shw2[0] = dir_daught_shw2[0]*mom_mag_daught_shw2;
  mom_daught_shw2[1] = dir_daught_shw2[1]*mom_mag_daught_shw2;
  mom_daught_shw2[2] = dir_daught_shw2[2]*mom_mag_daught_shw2;
  }

  std::cout<<"14"<<std::endl;
  std::cout<<"nu dir"<<std::endl;
  //need to think about how to account for momentum here
  double mom_mag = sqrt(pow(kine_energy,2) + 2*kine_energy*105.66);
  std::cout<<"mom_mag "<<mom_mag<<std::endl;
  mom[0] = init_dir_10[0]*mom_mag;
  mom[1] = init_dir_10[1]*mom_mag;
  mom[2] = init_dir_10[2]*mom_mag;
  std::cout<<"mom "<<mom[0]<<", "<<mom[1]<<", "<<mom[2]<<std::endl;
  TVector3 nu_dir = mom+mom_prim_track1+mom_prim_track2+mom_prim_shw1+mom_prim_shw2;
  nu_dir = nu_dir.Unit();
  std::cout<<"mom "<<nu_dir[0]<<", "<<nu_dir[1]<<", "<<nu_dir[2]<<std::endl;
  nu_angle_to_z = acos( nu_dir[0]*dir_beam[0] + nu_dir[1]*dir_beam[1] + nu_dir[2]*dir_beam[2] );
  nu_angle_to_target = acos( nu_dir[0]*target_dir[0] + nu_dir[1]*target_dir[1] + nu_dir[2]*target_dir[2] );
  nu_angle_to_absorber = acos( nu_dir[0]*absorber_dir[0] + nu_dir[1]*absorber_dir[1] + nu_dir[2]*absorber_dir[2] );
  nu_angle_to_vertical = acos( nu_dir[0]*dir_vertical[0] + nu_dir[1]*dir_vertical[1] + nu_dir[2]*dir_vertical[2] );

  //need to think about how to account for momentum here
  TVector3 track_dir = mom+mom_prim_track1+mom_prim_track2;
  track_dir = track_dir.Unit();
  track_angle_to_z = acos( track_dir[0]*dir_beam[0] + track_dir[1]*dir_beam[1] + track_dir[2]*dir_beam[2] );
  track_angle_to_target = acos( track_dir[0]*target_dir[0] + track_dir[1]*target_dir[1] + track_dir[2]*target_dir[2] );
  track_angle_to_absorber = acos( track_dir[0]*absorber_dir[0] + track_dir[1]*absorber_dir[1] + track_dir[2]*absorber_dir[2] );
  track_angle_to_vertical = acos( track_dir[0]*dir_vertical[0] + track_dir[1]*dir_vertical[1] + track_dir[2]*dir_vertical[2] );

  //reset if we have too long of adaughter and or prim
  //if( length_prim_track1>11 || kine_energy_cal_daught_track1>70 || length_prim_track2>5 || kine_energy_cal_daught_track1>60 ||kine_energy_cal_daught_shw1>70 || kine_energy_cal_prim_shw1>100 || kine_energy_cal_daught_shw2>60 || kine_energy_cal_prim_shw2>60){ tagger_info.ssm_flag_st_kdar = exit_ssm_tagger();}

    bool flag_st_kdar = false;
    if(n_prim_tracks_1==0 && n_prim_all_3==0 && n_daughter_tracks_5==0 && n_daughter_all_5<2 && Nsm_wivtx==1 && !(kine_pio_mass>70 && kine_pio_mass<200)){
    //if(n_prim_tracks_1==0 && n_prim_tracks_3==0 && n_daughter_tracks_5==0 && n_daughter_tracks_5<2 && Nsm_wivtx==1 && kine_pio_flag==0){
      std::cout<<"Pass KDAR Tagger"<<std::endl;
      flag_st_kdar = true;
    }
    else if(Nsm_wivtx!=1){std::cout<<"Fail: no short muons with vertex activity"<<std::endl;}
    //else if(kine_pio_flag!=0){std::cout<<"Fail: pi0 spotted"<<std::endl;}
    else if(kine_pio_mass>70 && kine_pio_mass<200){std::cout<<"Fail: pi0 spotted"<<std::endl;}
    else if(n_daughter_all_5>=2){std::cout<<"Fail: too many daughter showers"<<std::endl;}
    else if(n_daughter_tracks_5!=0){std::cout<<"Fail: long daughter track"<<std::endl;}
    else if(n_prim_tracks_1!=0){std::cout<<"Fail: other primary track"<<std::endl;}
    else if(n_prim_all_3!=0){std::cout<<"Fail: other primary shower"<<std::endl;}
    std::cout<<std::endl;


  tagger_info.ssm_flag_st_kdar = flag_st_kdar; //pass cutbased

  //tagger_info.ssm_Nsm = Nsm;//number of short straight muons
  //tagger_info.ssm_Nsm_wivtx = Nsm_wivtx;//number of short straight muons with vertex activity

  //only filled if there is one ssm
    //properties of the ssm
      //dq/dx info
  tagger_info.ssm_dq_dx_fwd_1 = dq_dx_fwd_1;//1st dq dx point wrt vertex
  tagger_info.ssm_dq_dx_fwd_2 = dq_dx_fwd_2;//2nd dq dx point wrt vertex
  tagger_info.ssm_dq_dx_fwd_3 = dq_dx_fwd_3;//3rd dq dx point wrt vertex
  tagger_info.ssm_dq_dx_fwd_4 = dq_dx_fwd_4;//4th dq dx point wrt vertex
  tagger_info.ssm_dq_dx_fwd_5 = dq_dx_fwd_5;//5th dq dx point wrt vertex
  tagger_info.ssm_dq_dx_bck_1 = dq_dx_bck_1;//1st dq dx point wrt end
  tagger_info.ssm_dq_dx_bck_2 = dq_dx_bck_2;//2nd dq dx point wrt end
  tagger_info.ssm_dq_dx_bck_3 = dq_dx_bck_3;//3rd dq dx point wrt end
  tagger_info.ssm_dq_dx_bck_4 = dq_dx_bck_4;//4th dq dx point wrt end
  tagger_info.ssm_dq_dx_bck_5 = dq_dx_bck_5;//5th dq dx point wrt end
  tagger_info.ssm_d_dq_dx_fwd_12 = d_dq_dx_fwd_12;//1st d dq dx point wrt vertex
  tagger_info.ssm_d_dq_dx_fwd_23 = d_dq_dx_fwd_23;//2nd d dq dx point wrt vertex
  tagger_info.ssm_d_dq_dx_fwd_34 = d_dq_dx_fwd_34;//3rd d dq dx point wrt vertex
  tagger_info.ssm_d_dq_dx_fwd_45 = d_dq_dx_fwd_45;//4th d dq dx point wrt vertex
  tagger_info.ssm_d_dq_dx_bck_12 = d_dq_dx_bck_12;//1st d dq dx point wrt end
  tagger_info.ssm_d_dq_dx_bck_23 = d_dq_dx_bck_23;//2nd d dq dx point wrt end
  tagger_info.ssm_d_dq_dx_bck_34 = d_dq_dx_bck_34;//3rd d dq dx point wrt end
  tagger_info.ssm_d_dq_dx_bck_45 = d_dq_dx_bck_45;//4th d dq dx point wrt end
  tagger_info.ssm_max_dq_dx_fwd_3 = max_dq_dx_fwd_3;//max dq/dx in first 3 points wrt vertex
  tagger_info.ssm_max_dq_dx_fwd_5 = max_dq_dx_fwd_5;//max dq/dx in first 5 points wrt vertex
  tagger_info.ssm_max_dq_dx_bck_3 = max_dq_dx_bck_3;//max dq/dx in first 3 points wrt end
  tagger_info.ssm_max_dq_dx_bck_5 = max_dq_dx_bck_5;//max dq/dx in first 5 points wrt end
  tagger_info.ssm_max_d_dq_dx_fwd_3 = max_d_dq_dx_fwd_3;//max d dq/dx in first 3 points wrt vertex
  tagger_info.ssm_max_d_dq_dx_fwd_5 = max_d_dq_dx_fwd_5;//max d dq/dx in first 5 points wrt vertex
  tagger_info.ssm_max_d_dq_dx_bck_3 = max_d_dq_dx_bck_3;//max d dq/dx in first 3 points wrt end
  tagger_info.ssm_max_d_dq_dx_bck_5 = max_d_dq_dx_bck_5;//max d dq/dx in first 5 points wrt end
  tagger_info.ssm_medium_dq_dx = medium_dq_dx;//medium dq/dx over the entire track
  tagger_info.ssm_medium_dq_dx_bp = medium_dq_dx_bp;//medium dq/dx over the entire track without the vertex activity
      //angluar info
  tagger_info.ssm_angle_to_z = angle_to_z_10;//angle to the z det cord
  tagger_info.ssm_angle_to_target = angle_to_target_10;//angle to the target
  tagger_info.ssm_angle_to_absorber = angle_to_absorber_10;//angle to the absorber
  tagger_info.ssm_angle_to_vertical = angle_to_vertical_10;//angle to the vertical
      //direction info
  tagger_info.ssm_x_dir = x_dir;//x det cord direction
  tagger_info.ssm_y_dir = y_dir;//y det cord direction
  tagger_info.ssm_z_dir = z_dir;//z det cord direction
      //energy info
  tagger_info.ssm_kine_energy = kine_energy;//KE of the muon with range
  tagger_info.ssm_kine_energy_reduced = kine_energy_reduced;//KE of the muon with range, exluding high dqdx at vertex
      //general properties
  tagger_info.ssm_vtx_activity = vtx_activity;//flag for identified vertex activity
  tagger_info.ssm_pdg = pdg;//initial pdg assigned to the track
  tagger_info.ssm_dQ_dx_cut = dQ_dx_cut;//cut normally used to choose muons
  tagger_info.ssm_score_mu_fwd = score_mu_fwd;//ID score normally normally used to choose muons
  tagger_info.ssm_score_p_fwd = score_p_fwd;//ID score normally normally used to choose protons
  tagger_info.ssm_score_e_fwd = score_e_fwd;//ID score normally normally used to choose e
  tagger_info.ssm_score_mu_bck = score_mu_bck;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_score_p_bck = score_p_bck;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_score_e_bck = score_e_bck;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_score_mu_fwd_bp = score_mu_fwd_bp;//ID score normally normally used to choose muons but calculated without the vertex activity
  tagger_info.ssm_score_p_fwd_bp = score_p_fwd_bp;//ID score normally normally used to choose protons but calculated without the vertex activity
  tagger_info.ssm_score_e_fwd_bp = score_e_fwd_bp;//ID score normally normally used to choose e but calculated without the vertex activity
      //track "straighness"
  tagger_info.ssm_length = length;//length of track direct caluclated from point to point
  tagger_info.ssm_direct_length = direct_length;//length of track direct from end to end
  double length_ratio = -999;
  if(direct_length>0 && length>0) length_ratio = direct_length/length;
  tagger_info.ssm_length_ratio = length_ratio;
  tagger_info.ssm_max_dev = max_dev;//furthest distance from track's primary axis along the track
    //number of other particles
  tagger_info.ssm_n_prim_tracks_1 = n_prim_tracks_1;//number of other primary tracks greater than 1 cm
  tagger_info.ssm_n_prim_tracks_3 = n_prim_tracks_3;//number of other primary tracks greater than 3 cm
  tagger_info.ssm_n_prim_tracks_5 = n_prim_tracks_5;//number of other primary tracks greater than 5 cm
  tagger_info.ssm_n_prim_tracks_8 = n_prim_tracks_8;//number of other primary tracks greater than 8 cm
  tagger_info.ssm_n_prim_tracks_11 = n_prim_tracks_11;//number of other primary tracks greater than 11 cm
  tagger_info.ssm_n_all_tracks_1 = n_prim_all_1;//number of other primary tracks and showers greater than 1 cm
  tagger_info.ssm_n_all_tracks_3 = n_prim_all_3;//number of other primary tracks and showers greater than 3 cm
  tagger_info.ssm_n_all_tracks_5 = n_prim_all_5;//number of other primary tracks and showers greater than 5 cm
  tagger_info.ssm_n_all_tracks_8 = n_prim_all_8;//number of other primary tracks and showers greater than 8 cm
  tagger_info.ssm_n_all_tracks_11 = n_prim_all_11;//number of other primary tracks and showers greater than 11 cm
  tagger_info.ssm_n_daughter_tracks_1 = n_daughter_tracks_1;//number of other daughter tracks greater than 1 cm
  tagger_info.ssm_n_daughter_tracks_3 = n_daughter_tracks_3;//number of other daughter tracks greater than 3 cm
  tagger_info.ssm_n_daughter_tracks_5 = n_daughter_tracks_5;//number of other daughter tracks greater than 5 cm
  tagger_info.ssm_n_daughter_tracks_8 = n_daughter_tracks_8;//number of other daughter tracks greater than 8 cm
  tagger_info.ssm_n_daughter_tracks_11 = n_daughter_tracks_11;//number of other daughter tracks greater than 11 cm
  tagger_info.ssm_n_daughter_all_1 = n_daughter_all_1;//number of other daughter tracks and showers greater than 1 cm
  tagger_info.ssm_n_daughter_all_3 = n_daughter_all_3;//number of other daughter tracks and showers greater than 3 cm
  tagger_info.ssm_n_daughter_all_5 = n_daughter_all_5;//number of other daughter tracks and showers greater than 5 cm
  tagger_info.ssm_n_daughter_all_8 = n_daughter_all_8;//number of other daughter tracks and showers greater than 8 cm
  tagger_info.ssm_n_daughter_all_11 = n_daughter_all_11;//number of other daughter tracks and showers greater than 11 cm
    //properties of leading other primary track
  tagger_info.ssm_prim_track1_pdg = pdg_prim_track1;//initial pdg assigned to the track
  tagger_info.ssm_prim_track1_score_mu_fwd = score_mu_fwd_prim_track1;//ID score normally normally used to choose muons
  tagger_info.ssm_prim_track1_score_p_fwd = score_p_fwd_prim_track1;//ID score normally normally used to choose protons
  tagger_info.ssm_prim_track1_score_e_fwd = score_e_fwd_prim_track1;//ID score normally normally used to choose e
  tagger_info.ssm_prim_track1_score_mu_bck = score_mu_bck_prim_track1;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_prim_track1_score_p_bck = score_p_bck_prim_track1;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_prim_track1_score_e_bck = score_e_bck_prim_track1;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_prim_track1_length = length_prim_track1;//length of track direct caluclated from point to point
  tagger_info.ssm_prim_track1_direct_length = direct_length_prim_track1;//length of track direct from end to end
  double length_ratio_prim_track1 = -999;
  if(direct_length_prim_track1>0 && length_prim_track1>0) length_ratio_prim_track1 = direct_length_prim_track1/length_prim_track1;
  tagger_info.ssm_prim_track1_length_ratio = length_ratio_prim_track1;
  tagger_info.ssm_prim_track1_max_dev = max_dev_prim_track1;//furthest distance from track's primary axis along the track
  tagger_info.ssm_prim_track1_kine_energy_range = kine_energy_range_prim_track1;//KE of the muon with range
  tagger_info.ssm_prim_track1_kine_energy_range_mu = kine_energy_range_mu_prim_track1;
  tagger_info.ssm_prim_track1_kine_energy_range_p = kine_energy_range_p_prim_track1;
  tagger_info.ssm_prim_track1_kine_energy_range_e = kine_energy_range_e_prim_track1;
  tagger_info.ssm_prim_track1_kine_energy_cal = kine_energy_cal_prim_track1;//KE of the track with dq/dx
  tagger_info.ssm_prim_track1_medium_dq_dx = medium_dq_dx_prim_track1;//medium dq/dx over the entire track
  tagger_info.ssm_prim_track1_x_dir = x_dir_prim_track1;//x det cord direction
  tagger_info.ssm_prim_track1_y_dir = y_dir_prim_track1;//y det cord direction
  tagger_info.ssm_prim_track1_z_dir = z_dir_prim_track1;//z det cord direction
  tagger_info.ssm_prim_track1_add_daught_track_counts_1 = add_daught_track_counts_1_prim_track1;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_prim_track1_add_daught_all_counts_1 = add_daught_all_counts_1_prim_track1;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_prim_track1_add_daught_track_counts_5 = add_daught_track_counts_5_prim_track1;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_prim_track1_add_daught_all_counts_5 = add_daught_all_counts_5_prim_track1;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_prim_track1_add_daught_track_counts_11 = add_daught_track_counts_11_prim_track1;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_prim_track1_add_daught_all_counts_11 = add_daught_all_counts_11_prim_track1;//additional daughter tracks+showers longer than 11 cm
  //properties of sub-leading other primary track
  tagger_info.ssm_prim_track2_pdg = pdg_prim_track2;//initial pdg assigned to the track
  tagger_info.ssm_prim_track2_score_mu_fwd = score_mu_fwd_prim_track2;//ID score normally normally used to choose muons
  tagger_info.ssm_prim_track2_score_p_fwd = score_p_fwd_prim_track2;//ID score normally normally used to choose protons
  tagger_info.ssm_prim_track2_score_e_fwd = score_e_fwd_prim_track2;//ID score normally normally used to choose e
  tagger_info.ssm_prim_track2_score_mu_bck = score_mu_bck_prim_track2;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_prim_track2_score_p_bck = score_p_bck_prim_track2;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_prim_track2_score_e_bck = score_e_bck_prim_track2;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_prim_track2_length = length_prim_track2;//length of track direct caluclated from point to point
  tagger_info.ssm_prim_track2_direct_length = direct_length_prim_track2;//length of track direct from end to end
  double length_ratio_prim_track2 = -999;
  if(direct_length_prim_track2>0 && length_prim_track2>0) length_ratio_prim_track2 = direct_length_prim_track2/length_prim_track2;
  tagger_info.ssm_prim_track2_length_ratio = length_ratio_prim_track2;
  tagger_info.ssm_prim_track2_max_dev = max_dev_prim_track2;//furthest distance from track's primary axis along the track
  tagger_info.ssm_prim_track2_kine_energy_range = kine_energy_range_prim_track2;//KE of the track with range
  tagger_info.ssm_prim_track2_kine_energy_range_mu = kine_energy_range_mu_prim_track2;
  tagger_info.ssm_prim_track2_kine_energy_range_p = kine_energy_range_p_prim_track2;
  tagger_info.ssm_prim_track2_kine_energy_range_e = kine_energy_range_e_prim_track2;
  tagger_info.ssm_prim_track2_kine_energy_cal = kine_energy_cal_prim_track2;//KE of the track with dq/dx
  tagger_info.ssm_prim_track2_medium_dq_dx = medium_dq_dx_prim_track2;//medium dq/dx over the entire track
  tagger_info.ssm_prim_track2_x_dir = x_dir_prim_track2;//x det cord direction
  tagger_info.ssm_prim_track2_y_dir = y_dir_prim_track2;//y det cord direction
  tagger_info.ssm_prim_track2_z_dir = z_dir_prim_track2;//z det cord direction
  tagger_info.ssm_prim_track2_add_daught_track_counts_1 = add_daught_track_counts_1_prim_track2;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_prim_track2_add_daught_all_counts_1 = add_daught_all_counts_1_prim_track2;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_prim_track2_add_daught_track_counts_5 = add_daught_track_counts_5_prim_track2;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_prim_track2_add_daught_all_counts_5 = add_daught_all_counts_5_prim_track2;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_prim_track2_add_daught_track_counts_11 = add_daught_track_counts_11_prim_track2;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_prim_track2_add_daught_all_counts_11 = add_daught_all_counts_11_prim_track2;//additional daughter tracks+showers longer than 11 cm
  //properties of leading daughter track
  tagger_info.ssm_daught_track1_pdg = pdg_daught_track1;//initial pdg assigned to the track
  tagger_info.ssm_daught_track1_score_mu_fwd = score_mu_fwd_daught_track1;//ID score normally normally used to choose muons
  tagger_info.ssm_daught_track1_score_p_fwd = score_p_fwd_daught_track1;//ID score normally normally used to choose protons
  tagger_info.ssm_daught_track1_score_e_fwd = score_e_fwd_daught_track1;//ID score normally normally used to choose e
  tagger_info.ssm_daught_track1_score_mu_bck = score_mu_bck_daught_track1;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_daught_track1_score_p_bck = score_p_bck_daught_track1;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_daught_track1_score_e_bck = score_e_bck_daught_track1;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_daught_track1_length = length_daught_track1;//length of track direct caluclated from point to point
  tagger_info.ssm_daught_track1_direct_length = direct_length_daught_track1;//length of track direct from end to end
  double length_ratio_daught_track1 = -999;
  if(direct_length_daught_track1>0 && length_daught_track1>0) length_ratio_daught_track1 = direct_length_daught_track1/length_daught_track1;
  tagger_info.ssm_daught_track1_length_ratio = length_ratio_daught_track1;
  tagger_info.ssm_daught_track1_max_dev = max_dev_daught_track1;//furthest distance from track's daughtary axis along the track
  tagger_info.ssm_daught_track1_kine_energy_range = kine_energy_range_daught_track1;//KE of the track with range
  tagger_info.ssm_daught_track1_kine_energy_range_mu = kine_energy_range_mu_daught_track1;
  tagger_info.ssm_daught_track1_kine_energy_range_p = kine_energy_range_p_daught_track1;
  tagger_info.ssm_daught_track1_kine_energy_range_e = kine_energy_range_e_daught_track1;
  tagger_info.ssm_daught_track1_kine_energy_cal = kine_energy_cal_daught_track1;//KE of the track with dq/dx
  tagger_info.ssm_daught_track1_medium_dq_dx = medium_dq_dx_daught_track1;//medium dq/dx over the entire track 
  tagger_info.ssm_daught_track1_x_dir = x_dir_daught_track1;//x det cord direction
  tagger_info.ssm_daught_track1_y_dir = y_dir_daught_track1;//y det cord direction
  tagger_info.ssm_daught_track1_z_dir = z_dir_daught_track1;//z det cord direction
  tagger_info.ssm_daught_track1_add_daught_track_counts_1 = add_daught_track_counts_1_daught_track1;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_daught_track1_add_daught_all_counts_1 = add_daught_all_counts_1_daught_track1;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_daught_track1_add_daught_track_counts_5 = add_daught_track_counts_5_daught_track1;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_daught_track1_add_daught_all_counts_5 = add_daught_all_counts_5_daught_track1;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_daught_track1_add_daught_track_counts_11 = add_daught_track_counts_11_daught_track1;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_daught_track1_add_daught_all_counts_11 = add_daught_all_counts_11_daught_track1;//additional daughter tracks+showers longer than 11 cm
  //properties of sub-leading daughter track
  tagger_info.ssm_daught_track2_pdg = pdg_daught_track2;//initial pdg assigned to the track
  tagger_info.ssm_daught_track2_score_mu_fwd = score_mu_fwd_daught_track2;//ID score normally normally used to choose muons
  tagger_info.ssm_daught_track2_score_p_fwd = score_p_fwd_daught_track2;//ID score normally normally used to choose protons
  tagger_info.ssm_daught_track2_score_e_fwd = score_e_fwd_daught_track2;//ID score normally normally used to choose e
  tagger_info.ssm_daught_track2_score_mu_bck = score_mu_bck_daught_track2;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_daught_track2_score_p_bck = score_p_bck_daught_track2;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_daught_track2_score_e_bck = score_e_bck_daught_track2;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_daught_track2_length = length_daught_track2;//length of track direct caluclated from point to point
  tagger_info.ssm_daught_track2_direct_length = direct_length_daught_track2;//length of track direct from end to end
  double length_ratio_daught_track2 = -999;
  if(direct_length_daught_track2>0 && length_daught_track2>0) length_ratio_daught_track2 = direct_length_daught_track2/length_daught_track2;
  tagger_info.ssm_daught_track2_length_ratio = length_ratio_daught_track2;
  tagger_info.ssm_daught_track2_max_dev = max_dev_daught_track2;//furthest distance from track's daughtary axis along the track
  tagger_info.ssm_daught_track2_kine_energy_range = kine_energy_range_daught_track2;//KE of the muon with range
  tagger_info.ssm_daught_track2_kine_energy_range_mu = kine_energy_range_mu_daught_track2;
  tagger_info.ssm_daught_track2_kine_energy_range_p = kine_energy_range_p_daught_track2;
  tagger_info.ssm_daught_track2_kine_energy_range_e = kine_energy_range_e_daught_track2;
  tagger_info.ssm_daught_track2_kine_energy_cal = kine_energy_cal_daught_track2;//KE of the muon with dq/dx
  tagger_info.ssm_daught_track2_medium_dq_dx = medium_dq_dx_daught_track2;//medium dq/dx over the entire track
  tagger_info.ssm_daught_track2_x_dir = x_dir_daught_track2;//x det cord direction
  tagger_info.ssm_daught_track2_y_dir = y_dir_daught_track2;//y det cord direction
  tagger_info.ssm_daught_track2_z_dir = z_dir_daught_track2;//z det cord direction
  tagger_info.ssm_daught_track2_add_daught_track_counts_1 = add_daught_track_counts_1_daught_track2;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_daught_track2_add_daught_all_counts_1 = add_daught_all_counts_1_daught_track2;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_daught_track2_add_daught_track_counts_5 = add_daught_track_counts_5_daught_track2;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_daught_track2_add_daught_all_counts_5 = add_daught_all_counts_5_daught_track2;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_daught_track2_add_daught_track_counts_11 = add_daught_track_counts_11_daught_track2;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_daught_track2_add_daught_all_counts_11 = add_daught_all_counts_11_daught_track2;//additional daughter tracks+showers longer than 11 cm
  //properties of leading other primary shower
  tagger_info.ssm_prim_shw1_pdg = pdg_prim_shw1;//initial pdg assigned to the shw
  tagger_info.ssm_prim_shw1_score_mu_fwd = score_mu_fwd_prim_shw1;//ID score normally normally used to choose muons
  tagger_info.ssm_prim_shw1_score_p_fwd = score_p_fwd_prim_shw1;//ID score normally normally used to choose protons
  tagger_info.ssm_prim_shw1_score_e_fwd = score_e_fwd_prim_shw1;//ID score normally normally used to choose e
  tagger_info.ssm_prim_shw1_score_mu_bck = score_mu_bck_prim_shw1;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_prim_shw1_score_p_bck = score_p_bck_prim_shw1;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_prim_shw1_score_e_bck = score_e_bck_prim_shw1;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_prim_shw1_length = length_prim_shw1;//length of shw direct caluclated from point to point
  tagger_info.ssm_prim_shw1_direct_length = direct_length_prim_shw1;//length of shw direct from end to end
  double length_ratio_prim_shw1 = -999;
  if(direct_length_prim_shw1>0 && length_prim_shw1>0) length_ratio_prim_shw1 = direct_length_prim_shw1/length_prim_shw1;
  tagger_info.ssm_prim_shw1_length_ratio = length_ratio_prim_shw1;
  tagger_info.ssm_prim_shw1_max_dev = max_dev_prim_shw1;//furthest distance from shw's primary axis along the shw
  tagger_info.ssm_prim_shw1_kine_energy_range = kine_energy_range_prim_shw1;//KE of the muon with range
  tagger_info.ssm_prim_shw1_kine_energy_range_mu = kine_energy_range_mu_prim_shw1;
  tagger_info.ssm_prim_shw1_kine_energy_range_p = kine_energy_range_p_prim_shw1;
  tagger_info.ssm_prim_shw1_kine_energy_range_e = kine_energy_range_e_prim_shw1;
  tagger_info.ssm_prim_shw1_kine_energy_cal = kine_energy_cal_prim_shw1;//KE of the muon with dq/dx
  tagger_info.ssm_prim_shw1_medium_dq_dx = medium_dq_dx_prim_shw1;//medium dq/dx over the entire shw
  tagger_info.ssm_prim_shw1_x_dir = x_dir_prim_shw1;//x det cord direction
  tagger_info.ssm_prim_shw1_y_dir = y_dir_prim_shw1;//y det cord direction
  tagger_info.ssm_prim_shw1_z_dir = z_dir_prim_shw1;//z det cord direction
  tagger_info.ssm_prim_shw1_add_daught_track_counts_1 = add_daught_track_counts_1_prim_shw1;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_prim_shw1_add_daught_all_counts_1 = add_daught_all_counts_1_prim_shw1;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_prim_shw1_add_daught_track_counts_5 = add_daught_track_counts_5_prim_shw1;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_prim_shw1_add_daught_all_counts_5 = add_daught_all_counts_5_prim_shw1;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_prim_shw1_add_daught_track_counts_11 = add_daught_track_counts_11_prim_shw1;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_prim_shw1_add_daught_all_counts_11 = add_daught_all_counts_11_prim_shw1;//additional daughter tracks+showers longer than 11 cm
  //properties of sub-leading other primary shower
  tagger_info.ssm_prim_shw2_pdg = pdg_prim_shw2;//initial pdg assigned to the shw
  tagger_info.ssm_prim_shw2_score_mu_fwd = score_mu_fwd_prim_shw2;//ID score normally normally used to choose muons
  tagger_info.ssm_prim_shw2_score_p_fwd = score_p_fwd_prim_shw2;//ID score normally normally used to choose protons
  tagger_info.ssm_prim_shw2_score_e_fwd = score_e_fwd_prim_shw2;//ID score normally normally used to choose e
  tagger_info.ssm_prim_shw2_score_mu_bck = score_mu_bck_prim_shw2;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_prim_shw2_score_p_bck = score_p_bck_prim_shw2;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_prim_shw2_score_e_bck = score_e_bck_prim_shw2;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_prim_shw2_length = length_prim_shw2;//length of shw direct caluclated from point to point
  tagger_info.ssm_prim_shw2_direct_length = direct_length_prim_shw2;//length of shw direct from end to end
  double length_ratio_prim_shw2 = -999;
  if(direct_length_prim_shw2>0 && length_prim_shw2>0) length_ratio_prim_shw2 = direct_length_prim_shw2/length_prim_shw2;
  tagger_info.ssm_prim_shw2_length_ratio = length_ratio_prim_shw2;
  tagger_info.ssm_prim_shw2_max_dev = max_dev_prim_shw2;//furthest distance from shw's primary axis along the shw
  tagger_info.ssm_prim_shw2_kine_energy_range = kine_energy_range_prim_shw2;//KE of the shower with range
  tagger_info.ssm_prim_shw2_kine_energy_range_mu = kine_energy_range_mu_prim_shw2;
  tagger_info.ssm_prim_shw2_kine_energy_range_p = kine_energy_range_p_prim_shw2;
  tagger_info.ssm_prim_shw2_kine_energy_range_e = kine_energy_range_e_prim_shw2;
  tagger_info.ssm_prim_shw2_kine_energy_cal = kine_energy_cal_prim_shw2;//KE of the shower with dq/dx
  tagger_info.ssm_prim_shw2_medium_dq_dx = medium_dq_dx_prim_shw2;//medium dq/dx over the entire shw
  tagger_info.ssm_prim_shw2_x_dir = x_dir_prim_shw2;//x det cord direction
  tagger_info.ssm_prim_shw2_y_dir = y_dir_prim_shw2;//y det cord direction
  tagger_info.ssm_prim_shw2_z_dir = z_dir_prim_shw2;//z det cord direction
  tagger_info.ssm_prim_shw2_add_daught_track_counts_1 = add_daught_track_counts_1_prim_shw2;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_prim_shw2_add_daught_all_counts_1 = add_daught_all_counts_1_prim_shw2;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_prim_shw2_add_daught_track_counts_5 = add_daught_track_counts_5_prim_shw2;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_prim_shw2_add_daught_all_counts_5 = add_daught_all_counts_5_prim_shw2;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_prim_shw2_add_daught_track_counts_11 = add_daught_track_counts_11_prim_shw2;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_prim_shw2_add_daught_all_counts_11 = add_daught_all_counts_11_prim_shw2;//additional daughter tracks+showers longer than 11 cm
  //properties of leading daughter shower
  tagger_info.ssm_daught_shw1_pdg = pdg_daught_shw1;//initial pdg assigned to the shw
  tagger_info.ssm_daught_shw1_score_mu_fwd = score_mu_fwd_daught_shw1;//ID score normally normally used to choose muons
  tagger_info.ssm_daught_shw1_score_p_fwd = score_p_fwd_daught_shw1;//ID score normally normally used to choose protons
  tagger_info.ssm_daught_shw1_score_e_fwd = score_e_fwd_daught_shw1;//ID score normally normally used to choose e
  tagger_info.ssm_daught_shw1_score_mu_bck = score_mu_bck_daught_shw1;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_daught_shw1_score_p_bck = score_p_bck_daught_shw1;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_daught_shw1_score_e_bck = score_e_bck_daught_shw1;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_daught_shw1_length = length_daught_shw1;//length of shw direct caluclated from point to point
  tagger_info.ssm_daught_shw1_direct_length = direct_length_daught_shw1;//length of shw direct from end to end
  double length_ratio_daught_shw1 = -999;
  if(direct_length_daught_shw1>0 && length_daught_shw1>0) length_ratio_daught_shw1 = direct_length_daught_shw1/length_daught_shw1;
  tagger_info.ssm_daught_shw1_length_ratio = length_ratio_daught_shw1;
  tagger_info.ssm_daught_shw1_max_dev = max_dev_daught_shw1;//furthest distance from shw's daughtary axis along the shw
  tagger_info.ssm_daught_shw1_kine_energy_range = kine_energy_range_daught_shw1;//KE of the shower with range
  tagger_info.ssm_daught_shw1_kine_energy_range_mu = kine_energy_range_mu_daught_shw1;
  tagger_info.ssm_daught_shw1_kine_energy_range_p = kine_energy_range_p_daught_shw1;
  tagger_info.ssm_daught_shw1_kine_energy_range_e = kine_energy_range_e_daught_shw1;
  tagger_info.ssm_daught_shw1_kine_energy_cal = kine_energy_cal_daught_shw1;//KE of the shower with dq/dx
  tagger_info.ssm_daught_shw1_medium_dq_dx = medium_dq_dx_daught_shw1;//medium dq/dx over the entire shw 
  tagger_info.ssm_daught_shw1_x_dir = x_dir_daught_shw1;//x det cord direction
  tagger_info.ssm_daught_shw1_y_dir = y_dir_daught_shw1;//y det cord direction
  tagger_info.ssm_daught_shw1_z_dir = z_dir_daught_shw1;//z det cord direction
  tagger_info.ssm_daught_shw1_add_daught_track_counts_1 = add_daught_track_counts_1_daught_shw1;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_daught_shw1_add_daught_all_counts_1 = add_daught_all_counts_1_daught_shw1;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_daught_shw1_add_daught_track_counts_5 = add_daught_track_counts_5_daught_shw1;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_daught_shw1_add_daught_all_counts_5 = add_daught_all_counts_5_daught_shw1;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_daught_shw1_add_daught_track_counts_11 = add_daught_track_counts_11_daught_shw1;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_daught_shw1_add_daught_all_counts_11 = add_daught_all_counts_11_daught_shw1;//additional daughter tracks+showers longer than 11 cm
  //properties of sub-leading daughter shower
  tagger_info.ssm_daught_shw2_pdg = pdg_daught_shw2;//initial pdg assigned to the shw
  tagger_info.ssm_daught_shw2_score_mu_fwd = score_mu_fwd_daught_shw2;//ID score normally normally used to choose muons
  tagger_info.ssm_daught_shw2_score_p_fwd = score_p_fwd_daught_shw2;//ID score normally normally used to choose protons
  tagger_info.ssm_daught_shw2_score_e_fwd = score_e_fwd_daught_shw2;//ID score normally normally used to choose e
  tagger_info.ssm_daught_shw2_score_mu_bck = score_mu_bck_daught_shw2;//ID score normally normally used to choose muons but calculated backwards
  tagger_info.ssm_daught_shw2_score_p_bck = score_p_bck_daught_shw2;//ID score normally normally used to choose protons but calculated backwards
  tagger_info.ssm_daught_shw2_score_e_bck = score_e_bck_daught_shw2;//ID score normally normally used to choose e but calculated backwards
  tagger_info.ssm_daught_shw2_length = length_daught_shw2;//length of shw direct caluclated from point to point
  tagger_info.ssm_daught_shw2_direct_length = direct_length_daught_shw2;//length of shw direct from end to end
  double length_ratio_daught_shw2 = -999;
  if(direct_length_daught_shw2>0 && length_daught_shw2>0) length_ratio_daught_shw2 = direct_length_daught_shw2/length_daught_shw2;
  tagger_info.ssm_daught_shw2_length_ratio = length_ratio_daught_shw2;
  tagger_info.ssm_daught_shw2_max_dev = max_dev_daught_shw2;//furthest distance from shw's daughtary axis along the shw
  tagger_info.ssm_daught_shw2_kine_energy_range = kine_energy_range_daught_shw2;//KE of the shower with range
  tagger_info.ssm_daught_shw2_kine_energy_range_mu = kine_energy_range_mu_daught_shw2;
  tagger_info.ssm_daught_shw2_kine_energy_range_p = kine_energy_range_p_daught_shw2;
  tagger_info.ssm_daught_shw2_kine_energy_range_e = kine_energy_range_e_daught_shw2;
  tagger_info.ssm_daught_shw2_kine_energy_cal = kine_energy_cal_daught_shw2;//KE of the shower with dq/dx
  tagger_info.ssm_daught_shw2_medium_dq_dx = medium_dq_dx_daught_shw2;//medium dq/dx over the entire shw
  tagger_info.ssm_daught_shw2_x_dir = x_dir_daught_shw2;//x det cord direction
  tagger_info.ssm_daught_shw2_y_dir = y_dir_daught_shw2;//y det cord direction
  tagger_info.ssm_daught_shw2_z_dir = z_dir_daught_shw2;//z det cord direction
  tagger_info.ssm_daught_shw2_add_daught_track_counts_1 = add_daught_track_counts_1_daught_shw2;//additional daughter tracks longer than 1 cm
  tagger_info.ssm_daught_shw2_add_daught_all_counts_1 = add_daught_all_counts_1_daught_shw2;//additional daughter tracks+showers longer than 1 cm
  tagger_info.ssm_daught_shw2_add_daught_track_counts_5 = add_daught_track_counts_5_daught_shw2;//additional daughter tracks longer than 5 cm
  tagger_info.ssm_daught_shw2_add_daught_all_counts_5 = add_daught_all_counts_5_daught_shw2;//additional daughter tracks+showers longer than 5 cm
  tagger_info.ssm_daught_shw2_add_daught_track_counts_11 = add_daught_track_counts_11_daught_shw2;//additional daughter tracks longer than 11 cm
  tagger_info.ssm_daught_shw2_add_daught_all_counts_11 = add_daught_all_counts_11_daught_shw2;//additional daughter tracks+showers longer than 11 cm
  //event level properties
  tagger_info.ssm_nu_angle_z = nu_angle_to_z;
  tagger_info.ssm_nu_angle_target = nu_angle_to_target;
  tagger_info.ssm_nu_angle_absorber = nu_angle_to_absorber;
  tagger_info.ssm_nu_angle_vertical = nu_angle_to_vertical;
  tagger_info.ssm_track_angle_z = track_angle_to_z;
  tagger_info.ssm_track_angle_target = track_angle_to_target;
  tagger_info.ssm_track_angle_absorber = track_angle_to_absorber;
  tagger_info.ssm_track_angle_vertical = track_angle_to_vertical;

  print_ssm_tagger();

  return true;
};

std::vector<double> WCPPID::NeutrinoID::get_scores(WCPPID::ProtoSegment *sg){

  std::vector<WCP::Point >& fit_pt_vec =  sg->get_point_vec();
  std::vector<double>& vec_dQ = sg->get_dQ_vec();
  std::vector<double>& vec_dx = sg->get_dx_vec();

  int size = fit_pt_vec.size();
  std::cout<<fit_pt_vec.size()<<" fit_pt_vec.size()"<<std::endl;

  std::vector<double> L(size,0);
  std::vector<double> dQ_dx(size,0);
  std::vector<double> rL(size,0);
  std::vector<double> rdQ_dx(size,0);

  std::cout<<std::endl;
  double dis = 0;
  for (int i = 0; i<size; i++){
    double dq_dx = vec_dQ.at(i)/(vec_dx.at(i)/units::cm+1e-9);
    dQ_dx.at(i) = dq_dx;
    rdQ_dx.at(size-1-i) = dq_dx;
    L.at(i) = dis;
    if (i+1 < size)
      dis += sqrt(pow(fit_pt_vec.at(i+1).x-fit_pt_vec.at(i).x,2) + pow(fit_pt_vec.at(i+1).y-fit_pt_vec.at(i).y,2) + pow(fit_pt_vec.at(i+1).z - fit_pt_vec.at(i).z,2));
    std::cout<<dis<<"' ";
  }
  std::cout<<std::endl;
  // get reverse  vector for L
  for (size_t i=0;i!=L.size();i++){
    rL.at(i) = L.back() - L.at(L.size()-1-i);
  std::cout<<rL.at(i)<<"' ";
  }
  std::cout<<std::endl;
  std::vector<double> result_forward = sg->do_track_comp(L, dQ_dx, 15*units::cm);
  std::vector<double> result_backward = sg->do_track_comp(rL, rdQ_dx, 15*units::cm);

  std::vector<double> result;

  result.push_back( result_forward.at(1) );
  result.push_back( result_forward.at(2) );
  result.push_back( result_forward.at(3) );
  result.push_back( result_backward.at(1) );
  result.push_back( result_backward.at(2) );
  result.push_back( result_backward.at(3) );

  return result;

}

std::vector<double> WCPPID::NeutrinoID::get_scores(WCPPID::ProtoSegment *sg, int break_point, bool dir){

  std::vector<WCP::Point >& fit_pt_vec =  sg->get_point_vec();
  std::vector<double>& vec_dQ = sg->get_dQ_vec();
  std::vector<double>& vec_dx = sg->get_dx_vec();

  int size = fit_pt_vec.size();

  std::vector<double> L(size,0);
  std::vector<double> dQ_dx(size,0);
  std::vector<double> rL(size,0);
  std::vector<double> rdQ_dx(size,0);

  //if(break_point>=size) break_point==0;

  std::cout<<std::endl;
  //if(dir==1){
  double dis = 0;
    for (int i = break_point; i<size; i++){
      double dq_dx = vec_dQ.at(i)/(vec_dx.at(i)/units::cm+1e-9);
      dQ_dx.at(i) = dq_dx;
      rdQ_dx.at(size-1-i) = dq_dx;
      L.at(i) = dis;
      if (i+1 < size)
        dis += sqrt(pow(fit_pt_vec.at(i+1).x-fit_pt_vec.at(i).x,2) + pow(fit_pt_vec.at(i+1).y-fit_pt_vec.at(i).y,2) + pow(fit_pt_vec.at(i+1).z - fit_pt_vec.at(i).z,2));
    std::cout<<dis<<"' ";
    }
    std::cout<<std::endl;
    // get reverse  vector for L
    for (size_t i=break_point;i<L.size();i++){
      //rL.at(i) = L.back() - L.at(L.size()-1-i);
        rL.at(i) = L.at(break_point) - L.at(L.size()-1-i);
    }
  //}

  std::vector<double> result = sg->do_track_comp(L, dQ_dx, 15*units::cm);
  if(dir==-1){result = sg->do_track_comp(rL, rdQ_dx, 15*units::cm);}

  return result;

}


std::vector<double> WCPPID::NeutrinoID::calc_kine_range_multi_pdg(double length){

  std::vector<double> result;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  TGraph *g_range = 0;

  g_range = mp.get_muon_r2ke();
  double kine_energy_mu = g_range->Eval(length) * units::MeV;
  result.push_back(kine_energy_mu);

  g_range = mp.get_proton_r2ke();
  double kine_energy_p = g_range->Eval(length) * units::MeV;
  result.push_back(kine_energy_p);

  g_range = mp.get_electron_r2ke();
  double kine_energy_e = g_range->Eval(length) * units::MeV;
  result.push_back(kine_energy_e);

  return result;
}



bool WCPPID::NeutrinoID::exit_ssm_tagger(){
  tagger_info.ssm_dq_dx_fwd_1 = -999;
  tagger_info.ssm_dq_dx_fwd_2 = -999;
  tagger_info.ssm_dq_dx_fwd_3 = -999;
  tagger_info.ssm_dq_dx_fwd_4 = -999;
  tagger_info.ssm_dq_dx_fwd_5 = -999;
  tagger_info.ssm_dq_dx_bck_1 = -999;
  tagger_info.ssm_dq_dx_bck_2 = -999;
  tagger_info.ssm_dq_dx_bck_3 = -999;
  tagger_info.ssm_dq_dx_bck_4 = -999;
  tagger_info.ssm_dq_dx_bck_5 = -999;
  tagger_info.ssm_d_dq_dx_fwd_12 = -999;
  tagger_info.ssm_d_dq_dx_fwd_23 = -999;
  tagger_info.ssm_d_dq_dx_fwd_34 = -999;
  tagger_info.ssm_d_dq_dx_fwd_45 = -999;
  tagger_info.ssm_d_dq_dx_bck_12 = -999;
  tagger_info.ssm_d_dq_dx_bck_23 = -999;
  tagger_info.ssm_d_dq_dx_bck_34 = -999;
  tagger_info.ssm_d_dq_dx_bck_45 = -999;
  tagger_info.ssm_max_dq_dx_fwd_3 = -999;
  tagger_info.ssm_max_dq_dx_fwd_5 = -999;
  tagger_info.ssm_max_dq_dx_bck_3 = -999;
  tagger_info.ssm_max_dq_dx_bck_5 = -999;
  tagger_info.ssm_max_d_dq_dx_fwd_3 = -999;
  tagger_info.ssm_max_d_dq_dx_fwd_5 = -999;
  tagger_info.ssm_max_d_dq_dx_bck_3 = -999;
  tagger_info.ssm_max_d_dq_dx_bck_5 = -999;
  tagger_info.ssm_medium_dq_dx = -999;
  tagger_info.ssm_medium_dq_dx_bp = -999;
      //angluar info
  tagger_info.ssm_angle_to_z = -999;
  tagger_info.ssm_angle_to_target = -999;
  tagger_info.ssm_angle_to_absorber = -999;
  tagger_info.ssm_angle_to_vertical = -999;
      //directional info
  tagger_info.ssm_x_dir = -999;
  tagger_info.ssm_y_dir = -999;
  tagger_info.ssm_z_dir = -999;
      //energy info
  tagger_info.ssm_kine_energy = -999;
  tagger_info.ssm_kine_energy_reduced = -999;
      //general properties
  tagger_info.ssm_vtx_activity = -999;
  tagger_info.ssm_pdg = -999;
  tagger_info.ssm_dQ_dx_cut = -999;
  tagger_info.ssm_score_mu_fwd = -999;
  tagger_info.ssm_score_p_fwd = -999;
  tagger_info.ssm_score_e_fwd = -999;
  tagger_info.ssm_score_mu_bck = -999;
  tagger_info.ssm_score_p_bck = -999;
  tagger_info.ssm_score_e_bck = -999;
  tagger_info.ssm_score_mu_fwd_bp = -999;
  tagger_info.ssm_score_p_fwd_bp = -999;
  tagger_info.ssm_score_e_fwd_bp = -999;
      //track "straighness"
  tagger_info.ssm_length = -999;
  tagger_info.ssm_direct_length = -999;
  tagger_info.ssm_length_ratio = -999;
  tagger_info.ssm_max_dev = -999;
    //number of other particles
  tagger_info.ssm_n_prim_tracks_1 = -999;
  tagger_info.ssm_n_prim_tracks_3 = -999;
  tagger_info.ssm_n_prim_tracks_5 = -999;
  tagger_info.ssm_n_prim_tracks_8 = -999;
  tagger_info.ssm_n_prim_tracks_11 = -999;
  tagger_info.ssm_n_all_tracks_1 = -999;
  tagger_info.ssm_n_all_tracks_3 = -999;
  tagger_info.ssm_n_all_tracks_5 = -999;
  tagger_info.ssm_n_all_tracks_8 = -999;
  tagger_info.ssm_n_all_tracks_11 = -999;
  tagger_info.ssm_n_daughter_tracks_1 = -999;
  tagger_info.ssm_n_daughter_tracks_3 = -999;
  tagger_info.ssm_n_daughter_tracks_5 = -999;
  tagger_info.ssm_n_daughter_tracks_8 = -999;
  tagger_info.ssm_n_daughter_tracks_11 = -999;
  tagger_info.ssm_n_daughter_all_1 = -999;
  tagger_info.ssm_n_daughter_all_3 = -999;
  tagger_info.ssm_n_daughter_all_5 = -999;
  tagger_info.ssm_n_daughter_all_8 = -999;
  tagger_info.ssm_n_daughter_all_11 = -999;
    //properties of leading other primary track
  tagger_info.ssm_prim_track1_pdg = -999;
  tagger_info.ssm_prim_track1_score_mu_fwd = -999;
  tagger_info.ssm_prim_track1_score_p_fwd = -999;
  tagger_info.ssm_prim_track1_score_e_fwd = -999;
  tagger_info.ssm_prim_track1_score_mu_bck = -999;
  tagger_info.ssm_prim_track1_score_p_bck = -999;
  tagger_info.ssm_prim_track1_score_e_bck = -999;
  tagger_info.ssm_prim_track1_length = -999;
  tagger_info.ssm_prim_track1_direct_length = -999;
  tagger_info.ssm_prim_track1_length_ratio = -999;
  tagger_info.ssm_prim_track1_max_dev = -999;
  tagger_info.ssm_prim_track1_kine_energy_range = -999;
  tagger_info.ssm_prim_track1_kine_energy_range_mu = -999;
  tagger_info.ssm_prim_track1_kine_energy_range_p = -999;
  tagger_info.ssm_prim_track1_kine_energy_range_e = -999;
  tagger_info.ssm_prim_track1_kine_energy_cal = -999;
  tagger_info.ssm_prim_track1_medium_dq_dx = -999;
  tagger_info.ssm_prim_track1_x_dir = -999;
  tagger_info.ssm_prim_track1_y_dir = -999;
  tagger_info.ssm_prim_track1_z_dir = -999;
  tagger_info.ssm_prim_track1_add_daught_track_counts_1 = -999;
  tagger_info.ssm_prim_track1_add_daught_all_counts_1 = -999;
  tagger_info.ssm_prim_track1_add_daught_track_counts_5 = -999;
  tagger_info.ssm_prim_track1_add_daught_all_counts_5 = -999;
  tagger_info.ssm_prim_track1_add_daught_track_counts_11 = -999;
  tagger_info.ssm_prim_track1_add_daught_all_counts_11 = -999;
  //properties of sub-leading other primary track
  tagger_info.ssm_prim_track2_pdg = -999;
  tagger_info.ssm_prim_track2_score_mu_fwd = -999;
  tagger_info.ssm_prim_track2_score_p_fwd = -999;
  tagger_info.ssm_prim_track2_score_e_fwd = -999;
  tagger_info.ssm_prim_track2_score_mu_bck = -999;
  tagger_info.ssm_prim_track2_score_p_bck = -999;
  tagger_info.ssm_prim_track2_score_e_bck = -999;
  tagger_info.ssm_prim_track2_length = -999;
  tagger_info.ssm_prim_track2_direct_length = -999;
  tagger_info.ssm_prim_track2_length_ratio = -999;
  tagger_info.ssm_prim_track2_max_dev = -999;
  tagger_info.ssm_prim_track2_kine_energy_range = -999;
  tagger_info.ssm_prim_track2_kine_energy_range_mu = -999;
  tagger_info.ssm_prim_track2_kine_energy_range_p = -999;
  tagger_info.ssm_prim_track2_kine_energy_range_e = -999;
  tagger_info.ssm_prim_track2_kine_energy_cal = -999;
  tagger_info.ssm_prim_track2_medium_dq_dx = -999;
  tagger_info.ssm_prim_track2_x_dir = -999;
  tagger_info.ssm_prim_track2_y_dir = -999;
  tagger_info.ssm_prim_track2_z_dir = -999;
  tagger_info.ssm_prim_track2_add_daught_track_counts_1 = -999;
  tagger_info.ssm_prim_track2_add_daught_all_counts_1 = -999;
  tagger_info.ssm_prim_track2_add_daught_track_counts_5 = -999;
  tagger_info.ssm_prim_track2_add_daught_all_counts_5 = -999;
  tagger_info.ssm_prim_track2_add_daught_track_counts_11 = -999;
  tagger_info.ssm_prim_track2_add_daught_all_counts_11 = -999;
  //properties of leading daughter track
  tagger_info.ssm_daught_track1_pdg = -999;
  tagger_info.ssm_daught_track1_score_mu_fwd = -999;
  tagger_info.ssm_daught_track1_score_p_fwd = -999;
  tagger_info.ssm_daught_track1_score_e_fwd = -999;
  tagger_info.ssm_daught_track1_score_mu_bck = -999;
  tagger_info.ssm_daught_track1_score_p_bck = -999;
  tagger_info.ssm_daught_track1_score_e_bck = -999;
  tagger_info.ssm_daught_track1_length = -999;
  tagger_info.ssm_daught_track1_direct_length = -999;
  tagger_info.ssm_daught_track1_length_ratio = -999;
  tagger_info.ssm_daught_track1_max_dev = -999;
  tagger_info.ssm_daught_track1_kine_energy_range = -999;
  tagger_info.ssm_daught_track1_kine_energy_range_mu = -999;
  tagger_info.ssm_daught_track1_kine_energy_range_p = -999;
  tagger_info.ssm_daught_track1_kine_energy_range_e = -999;
  tagger_info.ssm_daught_track1_kine_energy_cal = -999;
  tagger_info.ssm_daught_track1_medium_dq_dx = -999;
  tagger_info.ssm_daught_track1_x_dir = -999;
  tagger_info.ssm_daught_track1_y_dir = -999;
  tagger_info.ssm_daught_track1_z_dir = -999;
  tagger_info.ssm_daught_track1_add_daught_track_counts_1 = -999;
  tagger_info.ssm_daught_track1_add_daught_all_counts_1 = -999;
  tagger_info.ssm_daught_track1_add_daught_track_counts_5 = -999;
  tagger_info.ssm_daught_track1_add_daught_all_counts_5 = -999;
  tagger_info.ssm_daught_track1_add_daught_track_counts_11 = -999;
  tagger_info.ssm_daught_track1_add_daught_all_counts_11 = -999;
  //properties of sub-leading daughter track
  tagger_info.ssm_daught_track2_pdg = -999;
  tagger_info.ssm_daught_track2_score_mu_fwd = -999;
  tagger_info.ssm_daught_track2_score_p_fwd = -999;
  tagger_info.ssm_daught_track2_score_e_fwd = -999;
  tagger_info.ssm_daught_track2_score_mu_bck = -999;
  tagger_info.ssm_daught_track2_score_p_bck = -999;
  tagger_info.ssm_daught_track2_score_e_bck = -999;
  tagger_info.ssm_daught_track2_length = -999;
  tagger_info.ssm_daught_track2_direct_length = -999;
  tagger_info.ssm_daught_track2_length_ratio = -999;
  tagger_info.ssm_daught_track2_max_dev = -999;
  tagger_info.ssm_daught_track2_kine_energy_range = -999;
  tagger_info.ssm_daught_track2_kine_energy_range_mu = -999;
  tagger_info.ssm_daught_track2_kine_energy_range_p = -999;
  tagger_info.ssm_daught_track2_kine_energy_range_e = -999;
  tagger_info.ssm_daught_track2_kine_energy_cal = -999;
  tagger_info.ssm_daught_track2_medium_dq_dx = -999;
  tagger_info.ssm_daught_track2_x_dir = -999;
  tagger_info.ssm_daught_track2_y_dir = -999;
  tagger_info.ssm_daught_track2_z_dir = -999;
  tagger_info.ssm_daught_track2_add_daught_track_counts_1 = -999;
  tagger_info.ssm_daught_track2_add_daught_all_counts_1 = -999;
  tagger_info.ssm_daught_track2_add_daught_track_counts_5 = -999;
  tagger_info.ssm_daught_track2_add_daught_all_counts_5 = -999;
  tagger_info.ssm_daught_track2_add_daught_track_counts_11 = -999;
  tagger_info.ssm_daught_track2_add_daught_all_counts_11 = -999;
  //properties of leading other primary shower
  tagger_info.ssm_prim_shw1_pdg = -999;
  tagger_info.ssm_prim_shw1_score_mu_fwd = -999;
  tagger_info.ssm_prim_shw1_score_p_fwd = -999;
  tagger_info.ssm_prim_shw1_score_e_fwd = -999;
  tagger_info.ssm_prim_shw1_score_mu_bck = -999;
  tagger_info.ssm_prim_shw1_score_p_bck = -999;
  tagger_info.ssm_prim_shw1_score_e_bck = -999;
  tagger_info.ssm_prim_shw1_length = -999;
  tagger_info.ssm_prim_shw1_direct_length = -999;
  tagger_info.ssm_prim_shw1_length_ratio = -999;
  tagger_info.ssm_prim_shw1_max_dev = -999;
  tagger_info.ssm_prim_shw1_kine_energy_range = -999;
  tagger_info.ssm_prim_shw1_kine_energy_range_mu = -999;
  tagger_info.ssm_prim_shw1_kine_energy_range_p = -999;
  tagger_info.ssm_prim_shw1_kine_energy_range_e = -999;
  tagger_info.ssm_prim_shw1_kine_energy_cal = -999;
  tagger_info.ssm_prim_shw1_medium_dq_dx = -999;
  tagger_info.ssm_prim_shw1_x_dir = -999;
  tagger_info.ssm_prim_shw1_y_dir = -999;
  tagger_info.ssm_prim_shw1_z_dir = -999;
  tagger_info.ssm_prim_shw1_add_daught_track_counts_1 = -999;
  tagger_info.ssm_prim_shw1_add_daught_all_counts_1 = -999;
  tagger_info.ssm_prim_shw1_add_daught_track_counts_5 = -999;
  tagger_info.ssm_prim_shw1_add_daught_all_counts_5 = -999;
  tagger_info.ssm_prim_shw1_add_daught_track_counts_11 = -999;
  tagger_info.ssm_prim_shw1_add_daught_all_counts_11 = -999;
    //properties of sub-leading other primary shower
  tagger_info.ssm_prim_shw2_pdg = -999;
  tagger_info.ssm_prim_shw2_score_mu_fwd = -999;
  tagger_info.ssm_prim_shw2_score_p_fwd = -999;
  tagger_info.ssm_prim_shw2_score_e_fwd = -999;
  tagger_info.ssm_prim_shw2_score_mu_bck = -999;
  tagger_info.ssm_prim_shw2_score_p_bck = -999;
  tagger_info.ssm_prim_shw2_score_e_bck = -999;
  tagger_info.ssm_prim_shw2_length = -999;
  tagger_info.ssm_prim_shw2_direct_length = -999;
  tagger_info.ssm_prim_shw2_length_ratio = -999;
  tagger_info.ssm_prim_shw2_max_dev = -999;
  tagger_info.ssm_prim_shw2_kine_energy_range = -999;
  tagger_info.ssm_prim_shw2_kine_energy_range_mu = -999;
  tagger_info.ssm_prim_shw2_kine_energy_range_p = -999;
  tagger_info.ssm_prim_shw2_kine_energy_range_e = -999;
  tagger_info.ssm_prim_shw2_kine_energy_cal = -999;
  tagger_info.ssm_prim_shw2_medium_dq_dx = -999;
  tagger_info.ssm_prim_shw2_x_dir = -999;
  tagger_info.ssm_prim_shw2_y_dir = -999;
  tagger_info.ssm_prim_shw2_z_dir = -999;
  tagger_info.ssm_prim_shw2_add_daught_track_counts_1 = -999;
  tagger_info.ssm_prim_shw2_add_daught_all_counts_1 = -999;
  tagger_info.ssm_prim_shw2_add_daught_track_counts_5 = -999;
  tagger_info.ssm_prim_shw2_add_daught_all_counts_5 = -999;
  tagger_info.ssm_prim_shw2_add_daught_track_counts_11 = -999;
  tagger_info.ssm_prim_shw2_add_daught_all_counts_11 = -999;
  //properties of leading daughter shower
  tagger_info.ssm_daught_shw1_pdg = -999;
  tagger_info.ssm_daught_shw1_score_mu_fwd = -999;
  tagger_info.ssm_daught_shw1_score_p_fwd = -999;
  tagger_info.ssm_daught_shw1_score_e_fwd = -999;
  tagger_info.ssm_daught_shw1_score_mu_bck = -999;
  tagger_info.ssm_daught_shw1_score_p_bck = -999;
  tagger_info.ssm_daught_shw1_score_e_bck = -999;
  tagger_info.ssm_daught_shw1_length = -999;
  tagger_info.ssm_daught_shw1_direct_length = -999;
  tagger_info.ssm_daught_shw1_length_ratio = -999;
  tagger_info.ssm_daught_shw1_max_dev = -999;
  tagger_info.ssm_daught_shw1_kine_energy_range = -999;
  tagger_info.ssm_daught_shw1_kine_energy_range_mu = -999;
  tagger_info.ssm_daught_shw1_kine_energy_range_p = -999;
  tagger_info.ssm_daught_shw1_kine_energy_range_e = -999;
  tagger_info.ssm_daught_shw1_kine_energy_cal = -999;
  tagger_info.ssm_daught_shw1_medium_dq_dx = -999;
  tagger_info.ssm_daught_shw1_x_dir = -999;
  tagger_info.ssm_daught_shw1_y_dir = -999;
  tagger_info.ssm_daught_shw1_z_dir = -999;
  tagger_info.ssm_daught_shw1_add_daught_track_counts_1 = -999;
  tagger_info.ssm_daught_shw1_add_daught_all_counts_1 = -999;
  tagger_info.ssm_daught_shw1_add_daught_track_counts_5 = -999;
  tagger_info.ssm_daught_shw1_add_daught_all_counts_5 = -999;
  tagger_info.ssm_daught_shw1_add_daught_track_counts_11 = -999;
  tagger_info.ssm_daught_shw1_add_daught_all_counts_11 = -999;
  //properties of sub-leading daughter shower
  tagger_info.ssm_daught_shw2_pdg = -999;
  tagger_info.ssm_daught_shw2_score_mu_fwd = -999;
  tagger_info.ssm_daught_shw2_score_p_fwd = -999;
  tagger_info.ssm_daught_shw2_score_e_fwd = -999;
  tagger_info.ssm_daught_shw2_score_mu_bck = -999;
  tagger_info.ssm_daught_shw2_score_p_bck = -999;
  tagger_info.ssm_daught_shw2_score_e_bck = -999;
  tagger_info.ssm_daught_shw2_length = -999;
  tagger_info.ssm_daught_shw2_direct_length = -999;
  tagger_info.ssm_daught_shw2_length_ratio = -999;
  tagger_info.ssm_daught_shw2_max_dev = -999;
  tagger_info.ssm_daught_shw2_kine_energy_range = -999;
  tagger_info.ssm_daught_shw2_kine_energy_range_mu = -999;
  tagger_info.ssm_daught_shw2_kine_energy_range_p = -999;
  tagger_info.ssm_daught_shw2_kine_energy_range_e = -999;
  tagger_info.ssm_daught_shw2_kine_energy_cal = -999;
  tagger_info.ssm_daught_shw2_medium_dq_dx = -999;
  tagger_info.ssm_daught_shw2_x_dir = -999;
  tagger_info.ssm_daught_shw2_y_dir = -999;
  tagger_info.ssm_daught_shw2_z_dir = -999;
  tagger_info.ssm_daught_shw2_add_daught_track_counts_1 = -999;
  tagger_info.ssm_daught_shw2_add_daught_all_counts_1 = -999;
  tagger_info.ssm_daught_shw2_add_daught_track_counts_5 = -999;
  tagger_info.ssm_daught_shw2_add_daught_all_counts_5 = -999;
  tagger_info.ssm_daught_shw2_add_daught_track_counts_11 = -999;
  tagger_info.ssm_daught_shw2_add_daught_all_counts_11 = -999;
  //event level properties
  tagger_info.ssm_nu_angle_z = -999;
  tagger_info.ssm_nu_angle_target = -999;
  tagger_info.ssm_nu_angle_absorber = -999;
  tagger_info.ssm_nu_angle_vertical = -999;
  tagger_info.ssm_track_angle_z = -999;
  tagger_info.ssm_track_angle_target = -999;
  tagger_info.ssm_track_angle_absorber = -999;
  tagger_info.ssm_track_angle_vertical = -999;
	return false;
}

void WCPPID::NeutrinoID::print_ssm_tagger(){
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<"tagger_info.ssm_flag_st_kdar: "<<tagger_info.ssm_flag_st_kdar<<std::endl; //pass cutbased

  std::cout<<"tagger_info.ssm_Nsm: "<<tagger_info.ssm_Nsm<<std::endl;//number of short straight muons
  std::cout<<"tagger_info.ssm_Nsm_wivtx: "<<tagger_info.ssm_Nsm_wivtx<<std::endl;//number of short straight muons with vertex activity

  std::cout<<"tagger_info.ssm_dq_dx_fwd_1: "<<tagger_info.ssm_dq_dx_fwd_1<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_fwd_2: "<<tagger_info.ssm_dq_dx_fwd_2<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_fwd_3: "<<tagger_info.ssm_dq_dx_fwd_3<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_fwd_4: "<<tagger_info.ssm_dq_dx_fwd_4<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_fwd_5: "<<tagger_info.ssm_dq_dx_fwd_5<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_bck_1: "<<tagger_info.ssm_dq_dx_bck_1<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_bck_2: "<<tagger_info.ssm_dq_dx_bck_2<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_bck_3: "<<tagger_info.ssm_dq_dx_bck_3<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_bck_4: "<<tagger_info.ssm_dq_dx_bck_4<<std::endl;
  std::cout<<"tagger_info.ssm_dq_dx_bck_5: "<<tagger_info.ssm_dq_dx_bck_5<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_fwd_12: "<<tagger_info.ssm_d_dq_dx_fwd_12<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_fwd_23: "<<tagger_info.ssm_d_dq_dx_fwd_23<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_fwd_34: "<<tagger_info.ssm_d_dq_dx_fwd_34<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_fwd_45: "<<tagger_info.ssm_d_dq_dx_fwd_45<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_bck_12: "<<tagger_info.ssm_d_dq_dx_bck_12<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_bck_23: "<<tagger_info.ssm_d_dq_dx_bck_23<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_bck_34: "<<tagger_info.ssm_d_dq_dx_bck_34<<std::endl;
  std::cout<<"tagger_info.ssm_d_dq_dx_bck_45: "<<tagger_info.ssm_d_dq_dx_bck_45<<std::endl;
  std::cout<<"tagger_info.ssm_max_dq_dx_fwd_3: "<<tagger_info.ssm_max_dq_dx_fwd_3<<std::endl;
  std::cout<<"tagger_info.ssm_max_dq_dx_fwd_5: "<<tagger_info.ssm_max_dq_dx_fwd_5<<std::endl;
  std::cout<<"tagger_info.ssm_max_dq_dx_bck_3: "<<tagger_info.ssm_max_dq_dx_bck_3<<std::endl;
  std::cout<<"tagger_info.ssm_max_dq_dx_bck_5: "<<tagger_info.ssm_max_dq_dx_bck_5<<std::endl;
  std::cout<<"tagger_info.ssm_max_d_dq_dx_fwd_3: "<<tagger_info.ssm_max_d_dq_dx_fwd_3<<std::endl;
  std::cout<<"tagger_info.ssm_max_d_dq_dx_fwd_5: "<<tagger_info.ssm_max_d_dq_dx_fwd_5<<std::endl;
  std::cout<<"tagger_info.ssm_max_d_dq_dx_bck_3: "<<tagger_info.ssm_max_d_dq_dx_bck_3<<std::endl;
  std::cout<<"tagger_info.ssm_max_d_dq_dx_bck_5: "<<tagger_info.ssm_max_d_dq_dx_bck_5<<std::endl;
  std::cout<<"tagger_info.ssm_medium_dq_dx: "<<tagger_info.ssm_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_medium_dq_dx_bp: "<<tagger_info.ssm_medium_dq_dx_bp<<std::endl;
      //angluar info
  std::cout<<"tagger_info.ssm_angle_to_z: "<<tagger_info.ssm_angle_to_z<<std::endl;
  std::cout<<"tagger_info.ssm_angle_to_target: "<<tagger_info.ssm_angle_to_target<<std::endl;
  std::cout<<"tagger_info.ssm_angle_to_absorber: "<<tagger_info.ssm_angle_to_absorber<<std::endl;
  std::cout<<"tagger_info.ssm_angle_to_vertical: "<<tagger_info.ssm_angle_to_vertical<<std::endl;
     //directional info
  std::cout<<"tagger_info.ssm_x_dir: "<<tagger_info.ssm_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_y_dir: "<<tagger_info.ssm_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_z_dir: "<<tagger_info.ssm_z_dir<<std::endl; 
     //energy info
  std::cout<<"tagger_info.ssm_kine_energy: "<<tagger_info.ssm_kine_energy<<std::endl;
  std::cout<<"tagger_info.ssm_kine_energy_reduced: "<<tagger_info.ssm_kine_energy_reduced<<std::endl;
      //general properties
  std::cout<<"tagger_info.ssm_vtx_activity: "<<tagger_info.ssm_vtx_activity<<std::endl;
  std::cout<<"tagger_info.ssm_pdg: "<<tagger_info.ssm_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_dQ_dx_cut: "<<tagger_info.ssm_dQ_dx_cut<<std::endl;
  std::cout<<"tagger_info.ssm_score_mu_fwd: "<<tagger_info.ssm_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_score_p_fwd: "<<tagger_info.ssm_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_score_e_fwd: "<<tagger_info.ssm_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_score_mu_bck: "<<tagger_info.ssm_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_score_p_bck: "<<tagger_info.ssm_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_score_e_bck: "<<tagger_info.ssm_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_score_mu_fwd_bp: "<<tagger_info.ssm_score_mu_fwd_bp<<std::endl;
  std::cout<<"tagger_info.ssm_score_p_fwd_bp: "<<tagger_info.ssm_score_p_fwd_bp<<std::endl;
  std::cout<<"tagger_info.ssm_score_e_fwd_bp: "<<tagger_info.ssm_score_e_fwd_bp<<std::endl;
      //track "straighness"
  std::cout<<"tagger_info.ssm_length: "<<tagger_info.ssm_length<<std::endl;
  std::cout<<"tagger_info.ssm_direct_length: "<<tagger_info.ssm_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_length_ratio: "<<tagger_info.ssm_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_max_dev: "<<tagger_info.ssm_max_dev<<std::endl;
    //number of other particles
  std::cout<<"tagger_info.ssm_n_prim_tracks_1: "<<tagger_info.ssm_n_prim_tracks_1<<std::endl;
  std::cout<<"tagger_info.ssm_n_prim_tracks_3: "<<tagger_info.ssm_n_prim_tracks_3<<std::endl;
  std::cout<<"tagger_info.ssm_n_prim_tracks_5: "<<tagger_info.ssm_n_prim_tracks_5<<std::endl;
  std::cout<<"tagger_info.ssm_n_prim_tracks_8: "<<tagger_info.ssm_n_prim_tracks_8<<std::endl;
  std::cout<<"tagger_info.ssm_n_prim_tracks_11: "<<tagger_info.ssm_n_prim_tracks_11<<std::endl;
  std::cout<<"tagger_info.ssm_n_all_tracks_1: "<<tagger_info.ssm_n_all_tracks_1<<std::endl;
  std::cout<<"tagger_info.ssm_n_all_tracks_3: "<<tagger_info.ssm_n_all_tracks_3<<std::endl;
  std::cout<<"tagger_info.ssm_n_all_tracks_5: "<<tagger_info.ssm_n_all_tracks_5<<std::endl;
  std::cout<<"tagger_info.ssm_n_all_tracks_8: "<<tagger_info.ssm_n_all_tracks_8<<std::endl;
  std::cout<<"tagger_info.ssm_n_all_tracks_11: "<<tagger_info.ssm_n_all_tracks_11<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_tracks_1: "<<tagger_info.ssm_n_daughter_tracks_1<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_tracks_3: "<<tagger_info.ssm_n_daughter_tracks_3<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_tracks_5: "<<tagger_info.ssm_n_daughter_tracks_5<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_tracks_8: "<<tagger_info.ssm_n_daughter_tracks_8<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_tracks_11: "<<tagger_info.ssm_n_daughter_tracks_11<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_all_1: "<<tagger_info.ssm_n_daughter_all_1<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_all_3: "<<tagger_info.ssm_n_daughter_all_3<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_all_5: "<<tagger_info.ssm_n_daughter_all_5<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_all_8: "<<tagger_info.ssm_n_daughter_all_8<<std::endl;
  std::cout<<"tagger_info.ssm_n_daughter_all_11: "<<tagger_info.ssm_n_daughter_all_11<<std::endl;
    //properties of leading other primary track
  std::cout<<"tagger_info.ssm_prim_track1_pdg: "<<tagger_info.ssm_prim_track1_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_score_mu_fwd: "<<tagger_info.ssm_prim_track1_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_score_p_fwd: "<<tagger_info.ssm_prim_track1_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_score_e_fwd: "<<tagger_info.ssm_prim_track1_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_score_mu_bck: "<<tagger_info.ssm_prim_track1_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_score_p_bck: "<<tagger_info.ssm_prim_track1_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_score_e_bck: "<<tagger_info.ssm_prim_track1_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_length: "<<tagger_info.ssm_prim_track1_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_direct_length: "<<tagger_info.ssm_prim_track1_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_length_ratio: "<<tagger_info.ssm_prim_track1_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_max_dev: "<<tagger_info.ssm_prim_track1_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_kine_energy_range: "<<tagger_info.ssm_prim_track1_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_kine_energy_range_mu: "<<tagger_info.ssm_prim_track1_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_kine_energy_range_p: "<<tagger_info.ssm_prim_track1_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_kine_energy_range_e: "<<tagger_info.ssm_prim_track1_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_kine_energy_cal: "<<tagger_info.ssm_prim_track1_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_medium_dq_dx: "<<tagger_info.ssm_prim_track1_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_x_dir: "<<tagger_info.ssm_prim_track1_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_y_dir: "<<tagger_info.ssm_prim_track1_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_z_dir: "<<tagger_info.ssm_prim_track1_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_add_daught_track_counts_1: "<<tagger_info.ssm_prim_track1_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_add_daught_all_counts_1: "<<tagger_info.ssm_prim_track1_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_add_daught_track_counts_5: "<<tagger_info.ssm_prim_track1_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_add_daught_all_counts_5: "<<tagger_info.ssm_prim_track1_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_add_daught_track_counts_11: "<<tagger_info.ssm_prim_track1_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track1_add_daught_all_counts_11: "<<tagger_info.ssm_prim_track1_add_daught_all_counts_11<<std::endl;
  //properties of sub-leading other primary track
  std::cout<<"tagger_info.ssm_prim_track2_pdg: "<<tagger_info.ssm_prim_track2_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_score_mu_fwd: "<<tagger_info.ssm_prim_track2_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_score_p_fwd: "<<tagger_info.ssm_prim_track2_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_score_e_fwd: "<<tagger_info.ssm_prim_track2_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_score_mu_bck: "<<tagger_info.ssm_prim_track2_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_score_p_bck: "<<tagger_info.ssm_prim_track2_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_score_e_bck: "<<tagger_info.ssm_prim_track2_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_length: "<<tagger_info.ssm_prim_track2_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_direct_length: "<<tagger_info.ssm_prim_track2_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_length_ratio: "<<tagger_info.ssm_prim_track2_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_max_dev: "<<tagger_info.ssm_prim_track2_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_kine_energy_range: "<<tagger_info.ssm_prim_track2_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_kine_energy_range_mu: "<<tagger_info.ssm_prim_track2_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_kine_energy_range_p: "<<tagger_info.ssm_prim_track2_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_kine_energy_range_e: "<<tagger_info.ssm_prim_track2_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_kine_energy_cal: "<<tagger_info.ssm_prim_track2_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_medium_dq_dx: "<<tagger_info.ssm_prim_track2_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_x_dir: "<<tagger_info.ssm_prim_track2_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_y_dir: "<<tagger_info.ssm_prim_track2_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_z_dir: "<<tagger_info.ssm_prim_track2_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_add_daught_track_counts_1: "<<tagger_info.ssm_prim_track2_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_add_daught_all_counts_1: "<<tagger_info.ssm_prim_track2_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_add_daught_track_counts_5: "<<tagger_info.ssm_prim_track2_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_add_daught_all_counts_5: "<<tagger_info.ssm_prim_track2_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_add_daught_track_counts_11: "<<tagger_info.ssm_prim_track2_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_prim_track2_add_daught_all_counts_11: "<<tagger_info.ssm_prim_track2_add_daught_all_counts_11<<std::endl;
  //properties of leading daughter track
  std::cout<<"tagger_info.ssm_daught_track1_pdg: "<<tagger_info.ssm_daught_track1_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_score_mu_fwd: "<<tagger_info.ssm_daught_track1_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_score_p_fwd: "<<tagger_info.ssm_daught_track1_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_score_e_fwd: "<<tagger_info.ssm_daught_track1_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_score_mu_bck: "<<tagger_info.ssm_daught_track1_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_score_p_bck: "<<tagger_info.ssm_daught_track1_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_score_e_bck: "<<tagger_info.ssm_daught_track1_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_length: "<<tagger_info.ssm_daught_track1_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_direct_length: "<<tagger_info.ssm_daught_track1_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_length_ratio: "<<tagger_info.ssm_daught_track1_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_max_dev: "<<tagger_info.ssm_daught_track1_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_kine_energy_range: "<<tagger_info.ssm_daught_track1_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_kine_energy_range_mu: "<<tagger_info.ssm_daught_track1_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_kine_energy_range_p: "<<tagger_info.ssm_daught_track1_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_kine_energy_range_e: "<<tagger_info.ssm_daught_track1_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_kine_energy_cal: "<<tagger_info.ssm_daught_track1_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_medium_dq_dx: "<<tagger_info.ssm_daught_track1_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_x_dir: "<<tagger_info.ssm_daught_track1_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_y_dir: "<<tagger_info.ssm_daught_track1_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_z_dir: "<<tagger_info.ssm_daught_track1_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_add_daught_track_counts_1: "<<tagger_info.ssm_daught_track1_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_add_daught_all_counts_1: "<<tagger_info.ssm_daught_track1_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_add_daught_track_counts_5: "<<tagger_info.ssm_daught_track1_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_add_daught_all_counts_5: "<<tagger_info.ssm_daught_track1_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_add_daught_track_counts_11: "<<tagger_info.ssm_daught_track1_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track1_add_daught_all_counts_11: "<<tagger_info.ssm_daught_track1_add_daught_all_counts_11<<std::endl;
  //properties of sub-leading daughter track
  std::cout<<"tagger_info.ssm_daught_track2_pdg: "<<tagger_info.ssm_daught_track2_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_score_mu_fwd: "<<tagger_info.ssm_daught_track2_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_score_p_fwd: "<<tagger_info.ssm_daught_track2_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_score_e_fwd: "<<tagger_info.ssm_daught_track2_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_score_mu_bck: "<<tagger_info.ssm_daught_track2_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_score_p_bck: "<<tagger_info.ssm_daught_track2_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_score_e_bck: "<<tagger_info.ssm_daught_track2_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_length: "<<tagger_info.ssm_daught_track2_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_direct_length: "<<tagger_info.ssm_daught_track2_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_length_ratio: "<<tagger_info.ssm_daught_track2_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_max_dev: "<<tagger_info.ssm_daught_track2_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_kine_energy_range: "<<tagger_info.ssm_daught_track2_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_kine_energy_range_mu: "<<tagger_info.ssm_daught_track2_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_kine_energy_range_p: "<<tagger_info.ssm_daught_track2_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_kine_energy_range_e: "<<tagger_info.ssm_daught_track2_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_kine_energy_cal: "<<tagger_info.ssm_daught_track2_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_medium_dq_dx: "<<tagger_info.ssm_daught_track2_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_x_dir: "<<tagger_info.ssm_daught_track2_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_y_dir: "<<tagger_info.ssm_daught_track2_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_z_dir: "<<tagger_info.ssm_daught_track2_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_add_daught_track_counts_1: "<<tagger_info.ssm_daught_track2_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_add_daught_all_counts_1: "<<tagger_info.ssm_daught_track2_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_add_daught_track_counts_5: "<<tagger_info.ssm_daught_track2_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_add_daught_all_counts_5: "<<tagger_info.ssm_daught_track2_add_daught_all_counts_5<<std::endl; 
  std::cout<<"tagger_info.ssm_daught_track2_add_daught_track_counts_11: "<<tagger_info.ssm_daught_track2_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_daught_track2_add_daught_all_counts_11: "<<tagger_info.ssm_daught_track2_add_daught_all_counts_11<<std::endl;
  //properties of leading other primary shower
  std::cout<<"tagger_info.ssm_prim_shw1_pdg: "<<tagger_info.ssm_prim_shw1_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_score_mu_fwd: "<<tagger_info.ssm_prim_shw1_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_score_p_fwd: "<<tagger_info.ssm_prim_shw1_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_score_e_fwd: "<<tagger_info.ssm_prim_shw1_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_score_mu_bck: "<<tagger_info.ssm_prim_shw1_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_score_p_bck: "<<tagger_info.ssm_prim_shw1_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_score_e_bck: "<<tagger_info.ssm_prim_shw1_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_length: "<<tagger_info.ssm_prim_shw1_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_direct_length: "<<tagger_info.ssm_prim_shw1_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_length_ratio: "<<tagger_info.ssm_prim_shw1_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_max_dev: "<<tagger_info.ssm_prim_shw1_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_kine_energy_range: "<<tagger_info.ssm_prim_shw1_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_kine_energy_range_mu: "<<tagger_info.ssm_prim_shw1_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_kine_energy_range_p: "<<tagger_info.ssm_prim_shw1_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_kine_energy_range_e: "<<tagger_info.ssm_prim_shw1_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_kine_energy_cal: "<<tagger_info.ssm_prim_shw1_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_medium_dq_dx: "<<tagger_info.ssm_prim_shw1_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_x_dir: "<<tagger_info.ssm_prim_shw1_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_y_dir: "<<tagger_info.ssm_prim_shw1_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_z_dir: "<<tagger_info.ssm_prim_shw1_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_add_daught_track_counts_1: "<<tagger_info.ssm_prim_shw1_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_add_daught_all_counts_1: "<<tagger_info.ssm_prim_shw1_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_add_daught_track_counts_5: "<<tagger_info.ssm_prim_shw1_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_add_daught_all_counts_5: "<<tagger_info.ssm_prim_shw1_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_add_daught_track_counts_11: "<<tagger_info.ssm_prim_shw1_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw1_add_daught_all_counts_11: "<<tagger_info.ssm_prim_shw1_add_daught_all_counts_11<<std::endl;
  //properties of sub-leading other primary shower
  std::cout<<"tagger_info.ssm_prim_shw2_pdg: "<<tagger_info.ssm_prim_shw2_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_score_mu_fwd: "<<tagger_info.ssm_prim_shw2_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_score_p_fwd: "<<tagger_info.ssm_prim_shw2_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_score_e_fwd: "<<tagger_info.ssm_prim_shw2_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_score_mu_bck: "<<tagger_info.ssm_prim_shw2_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_score_p_bck: "<<tagger_info.ssm_prim_shw2_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_score_e_bck: "<<tagger_info.ssm_prim_shw2_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_length: "<<tagger_info.ssm_prim_shw2_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_direct_length: "<<tagger_info.ssm_prim_shw2_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_length_ratio: "<<tagger_info.ssm_prim_shw2_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_max_dev: "<<tagger_info.ssm_prim_shw2_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_kine_energy_range: "<<tagger_info.ssm_prim_shw2_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_kine_energy_range_mu: "<<tagger_info.ssm_prim_shw2_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_kine_energy_range_p: "<<tagger_info.ssm_prim_shw2_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_kine_energy_range_e: "<<tagger_info.ssm_prim_shw2_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_kine_energy_cal: "<<tagger_info.ssm_prim_shw2_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_medium_dq_dx: "<<tagger_info.ssm_prim_shw2_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_x_dir: "<<tagger_info.ssm_prim_shw2_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_y_dir: "<<tagger_info.ssm_prim_shw2_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_z_dir: "<<tagger_info.ssm_prim_shw2_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_add_daught_track_counts_1: "<<tagger_info.ssm_prim_shw2_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_add_daught_all_counts_1: "<<tagger_info.ssm_prim_shw2_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_add_daught_track_counts_5: "<<tagger_info.ssm_prim_shw2_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_add_daught_all_counts_5: "<<tagger_info.ssm_prim_shw2_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_add_daught_track_counts_11: "<<tagger_info.ssm_prim_shw2_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_prim_shw2_add_daught_all_counts_11: "<<tagger_info.ssm_prim_shw2_add_daught_all_counts_11<<std::endl;
  //properties of leading daughter shower
  std::cout<<"tagger_info.ssm_daught_shw1_pdg: "<<tagger_info.ssm_daught_shw1_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_score_mu_fwd: "<<tagger_info.ssm_daught_shw1_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_score_p_fwd: "<<tagger_info.ssm_daught_shw1_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_score_e_fwd: "<<tagger_info.ssm_daught_shw1_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_score_mu_bck: "<<tagger_info.ssm_daught_shw1_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_score_p_bck: "<<tagger_info.ssm_daught_shw1_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_score_e_bck: "<<tagger_info.ssm_daught_shw1_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_length: "<<tagger_info.ssm_daught_shw1_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_direct_length: "<<tagger_info.ssm_daught_shw1_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_length_ratio: "<<tagger_info.ssm_daught_shw1_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_max_dev: "<<tagger_info.ssm_daught_shw1_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_kine_energy_range: "<<tagger_info.ssm_daught_shw1_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_kine_energy_range_mu: "<<tagger_info.ssm_daught_shw1_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_kine_energy_range_p: "<<tagger_info.ssm_daught_shw1_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_kine_energy_range_e: "<<tagger_info.ssm_daught_shw1_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_kine_energy_cal: "<<tagger_info.ssm_daught_shw1_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_medium_dq_dx: "<<tagger_info.ssm_daught_shw1_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_x_dir: "<<tagger_info.ssm_daught_shw1_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_y_dir: "<<tagger_info.ssm_daught_shw1_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_z_dir: "<<tagger_info.ssm_daught_shw1_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_add_daught_track_counts_1: "<<tagger_info.ssm_daught_shw1_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_add_daught_all_counts_1: "<<tagger_info.ssm_daught_shw1_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_add_daught_track_counts_5: "<<tagger_info.ssm_daught_shw1_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_add_daught_all_counts_5: "<<tagger_info.ssm_daught_shw1_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_add_daught_track_counts_11: "<<tagger_info.ssm_daught_shw1_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw1_add_daught_all_counts_11: "<<tagger_info.ssm_daught_shw1_add_daught_all_counts_11<<std::endl;
  //properties of sub-leading daughter shower
  std::cout<<"tagger_info.ssm_daught_shw2_pdg: "<<tagger_info.ssm_daught_shw2_pdg<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_score_mu_fwd: "<<tagger_info.ssm_daught_shw2_score_mu_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_score_p_fwd: "<<tagger_info.ssm_daught_shw2_score_p_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_score_e_fwd: "<<tagger_info.ssm_daught_shw2_score_e_fwd<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_score_mu_bck: "<<tagger_info.ssm_daught_shw2_score_mu_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_score_p_bck: "<<tagger_info.ssm_daught_shw2_score_p_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_score_e_bck: "<<tagger_info.ssm_daught_shw2_score_e_bck<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_length: "<<tagger_info.ssm_daught_shw2_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_direct_length: "<<tagger_info.ssm_daught_shw2_direct_length<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_length_ratio: "<<tagger_info.ssm_daught_shw2_length_ratio<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_max_dev: "<<tagger_info.ssm_daught_shw2_max_dev<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_kine_energy_range: "<<tagger_info.ssm_daught_shw2_kine_energy_range<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_kine_energy_range_mu: "<<tagger_info.ssm_daught_shw2_kine_energy_range_mu<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_kine_energy_range_p: "<<tagger_info.ssm_daught_shw2_kine_energy_range_p<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_kine_energy_range_e: "<<tagger_info.ssm_daught_shw2_kine_energy_range_e<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_kine_energy_cal: "<<tagger_info.ssm_daught_shw2_kine_energy_cal<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_medium_dq_dx: "<<tagger_info.ssm_daught_shw2_medium_dq_dx<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_x_dir: "<<tagger_info.ssm_daught_shw2_x_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_y_dir: "<<tagger_info.ssm_daught_shw2_y_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_z_dir: "<<tagger_info.ssm_daught_shw2_z_dir<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_add_daught_track_counts_1: "<<tagger_info.ssm_daught_shw2_add_daught_track_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_add_daught_all_counts_1: "<<tagger_info.ssm_daught_shw2_add_daught_all_counts_1<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_add_daught_track_counts_5: "<<tagger_info.ssm_daught_shw2_add_daught_track_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_add_daught_all_counts_5: "<<tagger_info.ssm_daught_shw2_add_daught_all_counts_5<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_add_daught_track_counts_11: "<<tagger_info.ssm_daught_shw2_add_daught_track_counts_11<<std::endl;
  std::cout<<"tagger_info.ssm_daught_shw2_add_daught_all_counts_11: "<<tagger_info.ssm_daught_shw2_add_daught_all_counts_11<<std::endl;

  //event level properties
  std::cout<<"tagger_info.ssm_nu_angle_z: "<<tagger_info.ssm_nu_angle_z<<std::endl;
  std::cout<<"tagger_info.ssm_nu_angle_target: "<<tagger_info.ssm_nu_angle_target<<std::endl;
  std::cout<<"tagger_info.ssm_nu_angle_absorber: "<<tagger_info.ssm_nu_angle_absorber<<std::endl;
  std::cout<<"tagger_info.ssm_nu_angle_vertical: "<<tagger_info.ssm_nu_angle_vertical<<std::endl;
  std::cout<<"tagger_info.ssm_track_angle_z: "<<tagger_info.ssm_track_angle_z<<std::endl;
  std::cout<<"tagger_info.ssm_track_angle_target: "<<tagger_info.ssm_track_angle_target<<std::endl;
  std::cout<<"tagger_info.ssm_track_angle_absorber: "<<tagger_info.ssm_track_angle_absorber<<std::endl;
  std::cout<<"tagger_info.ssm_track_angle_vertical: "<<tagger_info.ssm_track_angle_vertical<<std::endl;


  std::cout<<std::endl;
  std::cout<<std::endl;

        
}
