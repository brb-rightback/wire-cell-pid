float WCPPID::NeutrinoID::cal_numu_bdts_xgboost(){
  float val = -10;
  
  tagger_info.numu_1_score = cal_numu_1_bdt(-0.4);
  tagger_info.numu_2_score = cal_numu_2_bdt(-0.1);
  tagger_info.numu_3_score = cal_numu_3_bdt(-0.2);
  tagger_info.cosmict_10_score = cal_cosmict_10_bdt(0.7);

  TMVA::Reader reader;
  
  reader.AddVariable("numu_cc_flag_3", &tagger_info.numu_cc_flag_3);
  reader.AddVariable("numu_cc_3_particle_type", &tagger_info.numu_cc_3_particle_type);
  reader.AddVariable("numu_cc_3_max_length", &tagger_info.numu_cc_3_max_length);
  reader.AddVariable("numu_cc_3_track_length",&tagger_info.numu_cc_3_acc_track_length);
  reader.AddVariable("numu_cc_3_max_length_all",&tagger_info.numu_cc_3_max_length_all);
  reader.AddVariable("numu_cc_3_max_muon_length",&tagger_info.numu_cc_3_max_muon_length);
  reader.AddVariable("numu_cc_3_n_daughter_tracks",&tagger_info.numu_cc_3_n_daughter_tracks);
  reader.AddVariable("numu_cc_3_n_daughter_all",&tagger_info.numu_cc_3_n_daughter_all);
  reader.AddVariable("cosmict_flag_2", &tagger_info.cosmict_flag_2);
  reader.AddVariable("cosmict_2_filled", &tagger_info.cosmict_2_filled);
  reader.AddVariable("cosmict_2_particle_type",&tagger_info.cosmict_2_particle_type);
  reader.AddVariable("cosmict_2_n_muon_tracks",&tagger_info.cosmict_2_n_muon_tracks);
  reader.AddVariable("cosmict_2_total_shower_length",&tagger_info.cosmict_2_total_shower_length);
  reader.AddVariable("cosmict_2_flag_inside",&tagger_info.cosmict_2_flag_inside);
  reader.AddVariable("cosmict_2_angle_beam",&tagger_info.cosmict_2_angle_beam);
  reader.AddVariable("cosmict_2_flag_dir_weak", &tagger_info.cosmict_2_flag_dir_weak);
  reader.AddVariable("cosmict_2_dQ_dx_end", &tagger_info.cosmict_2_dQ_dx_end);
  reader.AddVariable("cosmict_2_dQ_dx_front", &tagger_info.cosmict_2_dQ_dx_front);
  reader.AddVariable("cosmict_2_theta", &tagger_info.cosmict_2_theta);
  reader.AddVariable("cosmict_2_phi", &tagger_info.cosmict_2_phi);
  reader.AddVariable("cosmict_2_valid_tracks", &tagger_info.cosmict_2_valid_tracks);
  reader.AddVariable("cosmict_flag_4", &tagger_info.cosmict_flag_4);
  reader.AddVariable("cosmict_4_filled", &tagger_info.cosmict_4_filled);
  reader.AddVariable("cosmict_4_flag_inside", &tagger_info.cosmict_4_flag_inside);
  reader.AddVariable("cosmict_4_angle_beam", &tagger_info.cosmict_4_angle_beam);
  reader.AddVariable("cosmict_4_connected_showers", &tagger_info.cosmict_4_connected_showers);
  reader.AddVariable("cosmict_flag_3", &tagger_info.cosmict_flag_3);
  reader.AddVariable("cosmict_3_filled", &tagger_info.cosmict_3_filled);
  reader.AddVariable("cosmict_3_flag_inside", &tagger_info.cosmict_3_flag_inside);
  reader.AddVariable("cosmict_3_angle_beam", &tagger_info.cosmict_3_angle_beam);
  reader.AddVariable("cosmict_3_flag_dir_weak", &tagger_info.cosmict_3_flag_dir_weak);
  reader.AddVariable("cosmict_3_dQ_dx_end", &tagger_info.cosmict_3_dQ_dx_end);
  reader.AddVariable("cosmict_3_dQ_dx_front", &tagger_info.cosmict_3_dQ_dx_front);
  reader.AddVariable("cosmict_3_theta", &tagger_info.cosmict_3_theta);
  reader.AddVariable("cosmict_3_phi", &tagger_info.cosmict_3_phi);
  reader.AddVariable("cosmict_3_valid_tracks", &tagger_info.cosmict_3_valid_tracks);
  reader.AddVariable("cosmict_flag_5", &tagger_info.cosmict_flag_5);
  reader.AddVariable("cosmict_5_filled", &tagger_info.cosmict_5_filled);
  reader.AddVariable("cosmict_5_flag_inside", &tagger_info.cosmict_5_flag_inside);
  reader.AddVariable("cosmict_5_angle_beam", &tagger_info.cosmict_5_angle_beam);
  reader.AddVariable("cosmict_5_connected_showers", &tagger_info.cosmict_5_connected_showers);
  reader.AddVariable("cosmict_flag_6", &tagger_info.cosmict_flag_6);
  reader.AddVariable("cosmict_6_filled", &tagger_info.cosmict_6_filled);
  reader.AddVariable("cosmict_6_flag_dir_weak", &tagger_info.cosmict_6_flag_dir_weak);
  reader.AddVariable("cosmict_6_flag_inside", &tagger_info.cosmict_6_flag_inside);
  reader.AddVariable("cosmict_6_angle", &tagger_info.cosmict_6_angle);
  reader.AddVariable("cosmict_flag_7", &tagger_info.cosmict_flag_7);
  reader.AddVariable("cosmict_7_filled", &tagger_info.cosmict_7_filled);
  reader.AddVariable("cosmict_7_flag_sec", &tagger_info.cosmict_7_flag_sec);
  reader.AddVariable("cosmict_7_n_muon_tracks", &tagger_info.cosmict_7_n_muon_tracks);
  reader.AddVariable("cosmict_7_total_shower_length", &tagger_info.cosmict_7_total_shower_length);
  reader.AddVariable("cosmict_7_flag_inside", &tagger_info.cosmict_7_flag_inside);
  reader.AddVariable("cosmict_7_angle_beam", &tagger_info.cosmict_7_angle_beam);
  reader.AddVariable("cosmict_7_flag_dir_weak", &tagger_info.cosmict_7_flag_dir_weak);
  reader.AddVariable("cosmict_7_dQ_dx_end", &tagger_info.cosmict_7_dQ_dx_end);
  reader.AddVariable("cosmict_7_dQ_dx_front", &tagger_info.cosmict_7_dQ_dx_front);
  reader.AddVariable("cosmict_7_theta", &tagger_info.cosmict_7_theta);
  reader.AddVariable("cosmict_7_phi", &tagger_info.cosmict_7_phi);
  reader.AddVariable("cosmict_flag_8", &tagger_info.cosmict_flag_8);
  reader.AddVariable("cosmict_8_filled", &tagger_info.cosmict_8_filled);
  reader.AddVariable("cosmict_8_flag_out", &tagger_info.cosmict_8_flag_out);
  reader.AddVariable("cosmict_8_muon_length", &tagger_info.cosmict_8_muon_length);
  reader.AddVariable("cosmict_8_acc_length", &tagger_info.cosmict_8_acc_length);
  reader.AddVariable("cosmict_flag_9", &tagger_info.cosmict_flag_9);
  reader.AddVariable("cosmic_flag", &tagger_info.cosmic_flag);
  reader.AddVariable("cosmic_filled", &tagger_info.cosmic_filled);
  reader.AddVariable("cosmict_flag", &tagger_info.cosmict_flag);
  reader.AddVariable("numu_cc_flag", &tagger_info.numu_cc_flag);
  reader.AddVariable("cosmict_flag_1", &tagger_info.cosmict_flag_1);
  reader.AddVariable("kine_reco_Enu",&kine_info.kine_reco_Enu);
  reader.AddVariable("match_isFC",&match_isFC);
  reader.AddVariable("cosmict_10_score", &tagger_info.cosmict_10_score);
  reader.AddVariable("numu_1_score", &tagger_info.numu_1_score);
  reader.AddVariable("numu_2_score", &tagger_info.numu_2_score);

  reader.BookMVA( "MyBDT", "input_data_files/numu_scalars_scores_0923.xml");
  double val1 = reader.EvaluateMVA("MyBDT");
  
  val = TMath::Log10( (1+val1)/(1-val1) );
  
  return val;
}

float WCPPID::NeutrinoID::cal_numu_bdts(){
  float val = -10;

  tagger_info.numu_1_score = cal_numu_1_bdt(-0.4);
  tagger_info.numu_2_score = cal_numu_2_bdt(-0.1);
  tagger_info.numu_3_score = cal_numu_3_bdt(-0.2);
  tagger_info.cosmict_2_4_score = cal_cosmict_2_4_bdt(0.3);
  tagger_info.cosmict_3_5_score = cal_cosmict_3_5_bdt(0.5);
  tagger_info.cosmict_6_score = cal_cosmict_6_bdt(0.15);
  tagger_info.cosmict_7_score = cal_cosmict_7_bdt(0.3);
  tagger_info.cosmict_8_score = cal_cosmict_8_bdt(0.15);
  tagger_info.cosmict_10_score = cal_cosmict_10_bdt(0.7);

  TMVA::Reader reader;
  reader.AddVariable("cosmict_flag_1",&tagger_info.cosmict_flag_1);
  reader.AddVariable("cosmict_flag_9",&tagger_info.cosmict_flag_9);
  reader.AddVariable("cosmic_flag",&tagger_info.cosmic_flag);
  reader.AddVariable("cosmic_filled",&tagger_info.cosmic_filled);
  reader.AddVariable("cosmict_2_4_score",&tagger_info.cosmict_2_4_score);
  reader.AddVariable("cosmict_3_5_score",&tagger_info.cosmict_3_5_score);
  reader.AddVariable("cosmict_6_score",&tagger_info.cosmict_6_score);
  reader.AddVariable("cosmict_7_score",&tagger_info.cosmict_7_score);
  reader.AddVariable("cosmict_8_score",&tagger_info.cosmict_8_score);
  reader.AddVariable("cosmict_10_score",&tagger_info.cosmict_10_score);
  reader.AddVariable("numu_1_score",&tagger_info.numu_1_score);
  reader.AddVariable("numu_2_score",&tagger_info.numu_2_score);
  reader.AddVariable("numu_3_score",&tagger_info.numu_3_score);
  reader.BookMVA( "MyBDT", "input_data_files/weights/numucc.weights.xml");

  val = reader.EvaluateMVA("MyBDT");
  
  return val;
}

float WCPPID::NeutrinoID::cal_cosmict_2_4_bdt(float default_val){
  double val = default_val;
  TMVA::Reader reader_cosmict_2_4;
  reader_cosmict_2_4.AddVariable( "cosmict_2_particle_type", &tagger_info.cosmict_2_particle_type);
  reader_cosmict_2_4.AddVariable( "cosmict_2_n_muon_tracks", &tagger_info.cosmict_2_n_muon_tracks);
  reader_cosmict_2_4.AddVariable( "cosmict_2_total_shower_length", &tagger_info.cosmict_2_total_shower_length);
  reader_cosmict_2_4.AddVariable( "cosmict_2_flag_inside", &tagger_info.cosmict_2_flag_inside);
  reader_cosmict_2_4.AddVariable( "cosmict_2_angle_beam", &tagger_info.cosmict_2_angle_beam);
  reader_cosmict_2_4.AddVariable( "cosmict_2_flag_dir_weak", &tagger_info.cosmict_2_flag_dir_weak);
  reader_cosmict_2_4.AddVariable( "cosmict_2_dQ_dx_end", &tagger_info.cosmict_2_dQ_dx_end);
  reader_cosmict_2_4.AddVariable( "cosmict_2_dQ_dx_front", &tagger_info.cosmict_2_dQ_dx_front);
  reader_cosmict_2_4.AddVariable( "cosmict_2_theta", &tagger_info.cosmict_2_theta);
  reader_cosmict_2_4.AddVariable( "cosmict_2_phi", &tagger_info.cosmict_2_phi);
  reader_cosmict_2_4.AddVariable( "cosmict_2_valid_tracks", &tagger_info.cosmict_2_valid_tracks);
  reader_cosmict_2_4.AddVariable( "cosmict_4_flag_inside", &tagger_info.cosmict_4_flag_inside);
  reader_cosmict_2_4.AddVariable( "cosmict_4_connected_showers",  &tagger_info.cosmict_4_connected_showers);
  reader_cosmict_2_4.BookMVA( "MyBDT", "input_data_files/weights/cos_tagger_2_4.weights.xml");


  
  if (tagger_info.cosmict_2_filled==1){
    if (std::isnan(tagger_info.cosmict_2_angle_beam)) tagger_info.cosmict_2_angle_beam = 0;
    if (std::isnan(tagger_info.cosmict_2_theta)) tagger_info.cosmict_2_theta = 0;
    if (std::isnan(tagger_info.cosmict_2_phi)) tagger_info.cosmict_2_phi = 0;
    
    
    val = reader_cosmict_2_4.EvaluateMVA("MyBDT");
  }
  
  return val;
}


float WCPPID::NeutrinoID::cal_cosmict_3_5_bdt(float default_val){
  float val = default_val;

  TMVA::Reader reader_cosmict_3_5;
  reader_cosmict_3_5.AddVariable( "cosmict_3_flag_inside", &tagger_info.cosmict_3_flag_inside);
  reader_cosmict_3_5.AddVariable( "cosmict_3_angle_beam", &tagger_info.cosmict_3_angle_beam);
  reader_cosmict_3_5.AddVariable( "cosmict_3_flag_dir_weak", &tagger_info.cosmict_3_flag_dir_weak);
  reader_cosmict_3_5.AddVariable( "cosmict_3_dQ_dx_end", &tagger_info.cosmict_3_dQ_dx_end);
  reader_cosmict_3_5.AddVariable( "cosmict_3_dQ_dx_front", &tagger_info.cosmict_3_dQ_dx_front);
  reader_cosmict_3_5.AddVariable( "cosmict_3_theta", &tagger_info.cosmict_3_theta);
  reader_cosmict_3_5.AddVariable( "cosmict_3_phi", &tagger_info.cosmict_3_phi);
  reader_cosmict_3_5.AddVariable( "cosmict_3_valid_tracks", &tagger_info.cosmict_3_valid_tracks);
  reader_cosmict_3_5.AddVariable( "cosmict_5_connected_showers", &tagger_info.cosmict_5_connected_showers);
  reader_cosmict_3_5.BookMVA( "MyBDT", "input_data_files/weights/cos_tagger_3_5.weights.xml");
  
  if (tagger_info.cosmict_3_filled==1){
    val = reader_cosmict_3_5.EvaluateMVA("MyBDT");
  }

  return val;
}
float WCPPID::NeutrinoID::cal_cosmict_6_bdt(float default_val){
  float val = default_val;

  TMVA::Reader reader_cosmict_6;
  reader_cosmict_6.AddVariable( "cosmict_6_flag_dir_weak", &tagger_info.cosmict_6_flag_dir_weak);
  reader_cosmict_6.AddVariable( "cosmict_6_flag_inside", &tagger_info.cosmict_6_flag_inside);
  reader_cosmict_6.AddVariable( "cosmict_6_angle", &tagger_info.cosmict_6_angle);
  reader_cosmict_6.BookMVA( "MyBDT", "input_data_files/weights/cos_tagger_6.weights.xml");
  
  if (tagger_info.cosmict_6_filled==1){
    val = reader_cosmict_6.EvaluateMVA("MyBDT");
  }

  return val;
}
float WCPPID::NeutrinoID::cal_cosmict_7_bdt(float default_val){
  float val = default_val;

   TMVA::Reader reader_cosmict_7;
  reader_cosmict_7.AddVariable( "cosmict_7_flag_sec", &tagger_info.cosmict_7_flag_sec);
  reader_cosmict_7.AddVariable( "cosmict_7_n_muon_tracks", &tagger_info.cosmict_7_n_muon_tracks);
  reader_cosmict_7.AddVariable( "cosmict_7_total_shower_length", &tagger_info.cosmict_7_total_shower_length);
  reader_cosmict_7.AddVariable( "cosmict_7_flag_inside", &tagger_info.cosmict_7_flag_inside);
  reader_cosmict_7.AddVariable( "cosmict_7_angle_beam", &tagger_info.cosmict_7_angle_beam);
  reader_cosmict_7.AddVariable( "cosmict_7_flag_dir_weak", &tagger_info.cosmict_7_flag_dir_weak);
  reader_cosmict_7.AddVariable( "cosmict_7_dQ_dx_end", &tagger_info.cosmict_7_dQ_dx_end);
  reader_cosmict_7.AddVariable( "cosmict_7_dQ_dx_front", &tagger_info.cosmict_7_dQ_dx_front);
  reader_cosmict_7.AddVariable( "cosmict_7_theta", &tagger_info.cosmict_7_theta);
  reader_cosmict_7.AddVariable( "cosmict_7_phi", &tagger_info.cosmict_7_phi);
  reader_cosmict_7.BookMVA( "MyBDT", "input_data_files/weights/cos_tagger_7.weights.xml");
  
  if (tagger_info.cosmict_7_filled==1){
    val = reader_cosmict_7.EvaluateMVA("MyBDT");
  }
  return val;
}
float WCPPID::NeutrinoID::cal_cosmict_8_bdt(float default_val){
  float val = default_val;

  TMVA::Reader reader_cosmict_8;
  reader_cosmict_8.AddVariable( "cosmict_8_flag_out", &tagger_info.cosmict_8_flag_out);
  reader_cosmict_8.AddVariable( "cosmict_8_muon_length", &tagger_info.cosmict_8_muon_length);
  reader_cosmict_8.AddVariable( "cosmict_8_acc_length", &tagger_info.cosmict_8_acc_length);
  reader_cosmict_8.BookMVA( "MyBDT", "input_data_files/weights/cos_tagger_8.weights.xml");
  
  if (tagger_info.cosmict_8_filled==1){
    val = reader_cosmict_8.EvaluateMVA("MyBDT");
  }
  return val;
}
float WCPPID::NeutrinoID::cal_cosmict_10_bdt(float default_val){
  float val = default_val;

  float cosmict_10_vtx_z;
  float cosmict_10_flag_shower;
  float cosmict_10_flag_dir_weak;
  float cosmict_10_angle_beam;
  float cosmict_10_length;

  TMVA::Reader reader_cosmict_10;
  reader_cosmict_10.AddVariable("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  reader_cosmict_10.AddVariable("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  reader_cosmict_10.AddVariable("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  reader_cosmict_10.AddVariable("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  reader_cosmict_10.AddVariable("cosmict_10_length",&cosmict_10_length);
      
  reader_cosmict_10.BookMVA( "MyBDT", "input_data_files/weights/cos_tagger_10.weights.xml");
  
  
  if (tagger_info.cosmict_10_length.size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.cosmict_10_length.size();i++){
      cosmict_10_vtx_z = tagger_info.cosmict_10_vtx_z.at(i);
      cosmict_10_flag_shower = tagger_info.cosmict_10_flag_shower.at(i);
      cosmict_10_flag_dir_weak = tagger_info.cosmict_10_flag_dir_weak.at(i);
      cosmict_10_angle_beam = tagger_info.cosmict_10_angle_beam.at(i);
      cosmict_10_length = tagger_info.cosmict_10_length.at(i);

      if (std::isnan(cosmict_10_angle_beam)) cosmict_10_angle_beam = 0;
      
      float tmp_bdt =  reader_cosmict_10.EvaluateMVA("MyBDT");
      if (tmp_bdt < val) val = tmp_bdt;
    }
  }
  
  return val;
}
float WCPPID::NeutrinoID::cal_numu_1_bdt(float default_val){
  float val = default_val;

  float numu_cc_flag_1;
  float numu_cc_1_particle_type;
  float numu_cc_1_length;
  float numu_cc_1_medium_dQ_dx;
  float numu_cc_1_dQ_dx_cut;
  float numu_cc_1_direct_length;
  float numu_cc_1_n_daughter_tracks;
  float numu_cc_1_n_daughter_all;

  TMVA::Reader reader_numu_1;
  reader_numu_1.AddVariable("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  reader_numu_1.AddVariable("numu_cc_1_length",&numu_cc_1_length);
  reader_numu_1.AddVariable("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  reader_numu_1.AddVariable("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  reader_numu_1.AddVariable("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  reader_numu_1.AddVariable("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  reader_numu_1.AddVariable("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
      
  reader_numu_1.BookMVA( "MyBDT", "input_data_files/weights/numu_tagger1.weights.xml");
  
  if (tagger_info.numu_cc_1_particle_type.size()>0){
    val = -1e9;
    for (size_t i=0;i!=tagger_info.numu_cc_1_particle_type.size();i++){
      numu_cc_flag_1 = tagger_info.numu_cc_flag_1.at(i);
      numu_cc_1_particle_type= tagger_info.numu_cc_1_particle_type.at(i);
      numu_cc_1_length= tagger_info.numu_cc_1_length.at(i);
      numu_cc_1_medium_dQ_dx= tagger_info.numu_cc_1_medium_dQ_dx.at(i);
      numu_cc_1_dQ_dx_cut= tagger_info.numu_cc_1_dQ_dx_cut.at(i);
      numu_cc_1_direct_length= tagger_info.numu_cc_1_direct_length.at(i);
      numu_cc_1_n_daughter_tracks= tagger_info.numu_cc_1_n_daughter_tracks.at(i);
      numu_cc_1_n_daughter_all= tagger_info.numu_cc_1_n_daughter_all.at(i);

      if (std::isinf(numu_cc_1_dQ_dx_cut))  numu_cc_1_dQ_dx_cut = 10;
      
      float tmp_bdt =  reader_numu_1.EvaluateMVA("MyBDT");
      if (tmp_bdt > val) val = tmp_bdt;
    }
  }
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_numu_2_bdt(float default_val){
  float numu_cc_2_length;
  float numu_cc_2_total_length;
  float numu_cc_2_n_daughter_tracks;
  float numu_cc_2_n_daughter_all;

  float val = default_val;

  TMVA::Reader reader_numu_2;
  reader_numu_2.AddVariable("numu_cc_2_length",&numu_cc_2_length);
  reader_numu_2.AddVariable("numu_cc_2_total_length",&numu_cc_2_total_length);
  reader_numu_2.AddVariable("numu_cc_2_n_daughter_tracks",&numu_cc_2_n_daughter_tracks);
  reader_numu_2.AddVariable("numu_cc_2_n_daughter_all",&numu_cc_2_n_daughter_all);
  
  reader_numu_2.BookMVA( "MyBDT", "input_data_files/weights/numu_tagger2.weights.xml");
  
  if (tagger_info.numu_cc_2_length.size()>0){
    val = -1e9;
    for (size_t i=0;i!=tagger_info.numu_cc_2_length.size();i++){
      numu_cc_2_length = tagger_info.numu_cc_2_length.at(i);
      numu_cc_2_total_length = tagger_info.numu_cc_2_total_length.at(i);
      numu_cc_2_n_daughter_tracks = tagger_info.numu_cc_2_n_daughter_tracks.at(i);
      numu_cc_2_n_daughter_all = tagger_info.numu_cc_2_n_daughter_all.at(i);
	
      float tmp_bdt =  reader_numu_2.EvaluateMVA("MyBDT");
      if (tmp_bdt > val) val = tmp_bdt;
    }
  }

  return val;
}

float WCPPID::NeutrinoID::cal_numu_3_bdt(float default_val){
  float val = default_val;

  TMVA::Reader reader_numu_3;
  reader_numu_3.AddVariable( "numu_cc_3_particle_type", &tagger_info.numu_cc_3_particle_type);
  reader_numu_3.AddVariable( "numu_cc_3_max_length", &tagger_info.numu_cc_3_max_length);
  reader_numu_3.AddVariable( "numu_cc_3_acc_track_length", &tagger_info.numu_cc_3_acc_track_length);
  reader_numu_3.AddVariable( "numu_cc_3_max_length_all", &tagger_info.numu_cc_3_max_length_all);
  reader_numu_3.AddVariable( "numu_cc_3_max_muon_length", &tagger_info.numu_cc_3_max_muon_length);
  reader_numu_3.AddVariable( "numu_cc_3_n_daughter_tracks", &tagger_info.numu_cc_3_n_daughter_tracks);
  reader_numu_3.AddVariable( "numu_cc_3_n_daughter_all", &tagger_info.numu_cc_3_n_daughter_all);
    
  reader_numu_3.BookMVA( "MyBDT", "input_data_files/weights/numu_tagger3.weights.xml");
  
  val = reader_numu_3.EvaluateMVA("MyBDT");

  return val;
}
