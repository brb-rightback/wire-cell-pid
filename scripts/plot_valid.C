void plot_valid(){
  TChain *T = new TChain("T_tagger","T_tagger");
  T->Add("nue*.root");
  // cosmic tagger ones, one case of cosmics ...
  int cosmic_flag;
  int cosmic_n_solid_tracks;
  double cosmic_energy_main_showers;
  double cosmic_energy_direct_showers;
  double cosmic_energy_indirect_showers;
  int cosmic_n_direct_showers;
  int cosmic_n_indirect_showers;
  int cosmic_n_main_showers;
  int cosmic_filled;
  T->SetBranchAddress("cosmic_flag",&cosmic_flag);
  T->SetBranchAddress("cosmic_energy_main_showers",&cosmic_energy_main_showers);
  T->SetBranchAddress("cosmic_energy_indirect_showers",&cosmic_energy_indirect_showers);
  T->SetBranchAddress("cosmic_energy_direct_showers",&cosmic_energy_direct_showers);
  T->SetBranchAddress("cosmic_n_direct_showers",&cosmic_n_direct_showers);
  T->SetBranchAddress("cosmic_n_indirect_showers",&cosmic_n_indirect_showers);
  T->SetBranchAddress("cosmic_n_solid_tracks",&cosmic_n_solid_tracks);
  T->SetBranchAddress("cosmic_n_main_showers",&cosmic_n_main_showers);
  T->SetBranchAddress("cosmic_filled",&cosmic_filled);

  TH1F *hc1 = new TH1F("hc1","hc1",5,-3,2);
  TH1F *hc2 = new TH1F("hc2","hc2",100,0,2000);
  TH1F *hc3 = new TH1F("hc3","hc3",100,0,2000);
  TH1F *hc4 = new TH1F("hc4","hc4",100,0,2000);
  TH1F *hc5 = new TH1F("hc5","hc5",5,-1,4);
  TH1F *hc6 = new TH1F("hc6","hc6",5,-1,4);
  TH1F *hc7 = new TH1F("hc7","hc7",5,-1,4);
  TH1F *hc8 = new TH1F("hc8","hc8",5,-1,4);
  TH1F *hc9 = new TH1F("hc9","hc9",5,-3,2);
  
  
  // shower gap identification
  int gap_flag;
  int gap_flag_prolong_u;
  int gap_flag_prolong_v;
  int gap_flag_prolong_w;
  int gap_flag_parallel;
  int gap_n_points;
  int gap_n_bad;
  double gap_energy;
  int gap_num_valid_tracks;
  int gap_flag_single_shower;
  int gap_filled;
  T->SetBranchAddress("gap_flag",&gap_flag);
  T->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
  T->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
  T->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
  T->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
  T->SetBranchAddress("gap_n_points",&gap_n_points);
  T->SetBranchAddress("gap_n_bad",&gap_n_bad);
  T->SetBranchAddress("gap_energy",&gap_energy);
  T->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
  T->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
  T->SetBranchAddress("gap_filled",&gap_filled);
  
  TH1F *hg1 = new TH1F("hg1","hg1",5,-3,2);
  TH1F *hg2 = new TH1F("hg2","hg2",5,-3,2);
  TH1F *hg3 = new TH1F("hg3","hg3",5,-3,2);
  TH1F *hg4 = new TH1F("hg4","hg4",5,-3,2);
  TH1F *hg5 = new TH1F("hg5","hg5",5,-3,2);
  TH1F *hg6 = new TH1F("hg6","hg6",100,0,100);
  TH1F *hg7 = new TH1F("hg7","hg7",10,0,10);
  TH1F *hg8 = new TH1F("hg8","hg8",100,0,2000);
  TH1F *hg9 = new TH1F("hg9","hg9",5,0,5);
  TH1F *hg10 = new TH1F("hg10","hg10",5,-3,2);
  TH1F *hg11 = new TH1F("hg11","hg11",5,-3,2);

  int mip_quality_flag;
  double mip_quality_energy;
  int mip_quality_overlap;
  int mip_quality_n_showers;
  int mip_quality_n_tracks;
  int mip_quality_flag_inside_pi0;
  int mip_quality_n_pi0_showers;
  double mip_quality_shortest_length;
  double mip_quality_acc_length;
  double mip_quality_shortest_angle;
  int mip_quality_flag_proton;
  int mip_quality_filled;

  T->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  T->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  T->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  T->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  T->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  T->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  T->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  T->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  T->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  T->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  T->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  T->SetBranchAddress("mip_quality_filled",&mip_quality_filled);

  TH1F *hmq1 = new TH1F("hmq1","hmq1",5,-3,2);
  TH1F *hmq2 = new TH1F("hmq2","hmq2",100,0,2000);
  TH1F *hmq3 = new TH1F("hmq3","hmq3",5,-3,2);
  TH1F *hmq4 = new TH1F("hmq4","hmq4",5,0,5);
  TH1F *hmq5 = new TH1F("hmq5","hmq5",5,0,5);
  TH1F *hmq6 = new TH1F("hmq6","hmq6",5,-3,2);
  TH1F *hmq7 = new TH1F("hmq7","hmq7",5,0,5);
  TH1F *hmq8 = new TH1F("hmq8","hmq8",100,0,300);
  TH1F *hmq9 = new TH1F("hmq9","hmq9",100,0,300);
  TH1F *hmq10 = new TH1F("hmq10","hmq10",100,0,180);
  TH1F *hmq11 = new TH1F("hmq11","hmq11",5,-3,2);
  TH1F *hmq12 = new TH1F("hmq12","hmq12",5,-3,2);

  // mip identification
  int mip_flag;
  double mip_energy;
  int mip_n_end_reduction;    
  int mip_n_first_mip;
  int mip_n_first_non_mip;
  int mip_n_first_non_mip_1;
  int mip_n_first_non_mip_2;
  double mip_vec_dQ_dx_0;
  double mip_vec_dQ_dx_1;
  double mip_max_dQ_dx_sample;
  int mip_n_below_threshold;
  int mip_n_below_zero;
  int mip_n_lowest;
  int mip_n_highest;
  double mip_lowest_dQ_dx;
  double mip_highest_dQ_dx;
  double mip_medium_dQ_dx;
  double mip_stem_length;
  double mip_length_main;
  double mip_length_total;
  double mip_angle_beam;
  double mip_iso_angle;
  int mip_n_vertex;
  int mip_n_good_tracks;
  double mip_E_indirect_max_energy;
  int mip_flag_all_above;
  double mip_min_dQ_dx_5;
  int mip_n_other_vertex; 
  int mip_n_stem_size;
  int mip_flag_stem_trajectory;
  double mip_min_dis;
  int mip_filled;
  
  // extra
  double mip_vec_dQ_dx_2;
  double mip_vec_dQ_dx_3;
  double mip_vec_dQ_dx_4;
  double mip_vec_dQ_dx_5;
  double mip_vec_dQ_dx_6;
  double mip_vec_dQ_dx_7;
  double mip_vec_dQ_dx_8;
  double mip_vec_dQ_dx_9;
  double mip_vec_dQ_dx_10;
  double mip_vec_dQ_dx_11;
  double mip_vec_dQ_dx_12;
  double mip_vec_dQ_dx_13;
  double mip_vec_dQ_dx_14;
  double mip_vec_dQ_dx_15;
  double mip_vec_dQ_dx_16;
  double mip_vec_dQ_dx_17;
  double mip_vec_dQ_dx_18;
  double mip_vec_dQ_dx_19;
  

  // mip
  T->SetBranchAddress("mip_flag",&mip_flag);
  T->SetBranchAddress("mip_energy",&mip_energy);
  T->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
  T->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
  T->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
  T->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
  T->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
  
  T->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
  T->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
  T->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
  T->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
  T->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
  T->SetBranchAddress("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5);
  T->SetBranchAddress("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6);
  T->SetBranchAddress("mip_vec_dQ_dx_7",&mip_vec_dQ_dx_7);
  T->SetBranchAddress("mip_vec_dQ_dx_8",&mip_vec_dQ_dx_8);
  T->SetBranchAddress("mip_vec_dQ_dx_9",&mip_vec_dQ_dx_9);
  T->SetBranchAddress("mip_vec_dQ_dx_10",&mip_vec_dQ_dx_10);
  T->SetBranchAddress("mip_vec_dQ_dx_11",&mip_vec_dQ_dx_11);
  T->SetBranchAddress("mip_vec_dQ_dx_12",&mip_vec_dQ_dx_12);
  T->SetBranchAddress("mip_vec_dQ_dx_13",&mip_vec_dQ_dx_13);
  T->SetBranchAddress("mip_vec_dQ_dx_14",&mip_vec_dQ_dx_14);
  T->SetBranchAddress("mip_vec_dQ_dx_15",&mip_vec_dQ_dx_15);
  T->SetBranchAddress("mip_vec_dQ_dx_16",&mip_vec_dQ_dx_16);
  T->SetBranchAddress("mip_vec_dQ_dx_17",&mip_vec_dQ_dx_17);
  T->SetBranchAddress("mip_vec_dQ_dx_18",&mip_vec_dQ_dx_18);
  T->SetBranchAddress("mip_vec_dQ_dx_19",&mip_vec_dQ_dx_19);
  
  T->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
  T->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
  T->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
  T->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
  T->SetBranchAddress("mip_n_highest",&mip_n_highest);
  
  T->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
  T->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
  T->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
  T->SetBranchAddress("mip_stem_length",&mip_stem_length);
  T->SetBranchAddress("mip_length_main",&mip_length_main);
  T->SetBranchAddress("mip_length_total",&mip_length_total);
  T->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
  T->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
  
  T->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
  T->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
  T->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
  T->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
  T->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
  T->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
  T->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
  T->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
  T->SetBranchAddress("mip_min_dis",&mip_min_dis);
  T->SetBranchAddress("mip_filled",&mip_filled);

  TH1F *hmip1 = new TH1F("hmip1","hmip1",5,-3,2);
  TH1F *hmip2 = new TH1F("hmip2","hmip2",100,0,2000);
  TH1F *hmip3 = new TH1F("hmip3","hmip3",20,0,20);
  TH1F *hmip4 = new TH1F("hmip4","hmip4",20,0,20);
  TH1F *hmip5 = new TH1F("hmip5","hmip5",20,0,20);
  TH1F *hmip6 = new TH1F("hmip6","hmip6",20,0,20);
  TH1F *hmip7 = new TH1F("hmip7","hmip7",20,0,20);

  TH1F *hmip8 = new TH1F("hmip8","hmip8",20,0,6);
  TH1F *hmip9 = new TH1F("hmip9","hmip9",20,0,6);
  TH1F *hmip10 = new TH1F("hmip10","hmip10",20,0,6);
  TH1F *hmip11 = new TH1F("hmip11","hmip11",20,0,6);
  TH1F *hmip12 = new TH1F("hmip12","hmip12",20,0,6);
  TH1F *hmip13 = new TH1F("hmip13","hmip13",20,0,6);
  TH1F *hmip14 = new TH1F("hmip14","hmip14",20,0,6);
  TH1F *hmip15 = new TH1F("hmip15","hmip15",20,0,6);
  TH1F *hmip16 = new TH1F("hmip16","hmip16",20,0,6);
  TH1F *hmip17 = new TH1F("hmip17","hmip17",20,0,6);
  TH1F *hmip18 = new TH1F("hmip18","hmip18",20,0,6);
  TH1F *hmip19 = new TH1F("hmip19","hmip19",20,0,6);
  TH1F *hmip20 = new TH1F("hmip20","hmip20",20,0,6);
  TH1F *hmip21 = new TH1F("hmip21","hmip21",20,0,6);
  TH1F *hmip22 = new TH1F("hmip22","hmip22",20,0,6);
  TH1F *hmip23 = new TH1F("hmip23","hmip23",20,0,6);
  TH1F *hmip24 = new TH1F("hmip24","hmip24",20,0,6);
  TH1F *hmip25 = new TH1F("hmip25","hmip25",20,0,6);
  TH1F *hmip26 = new TH1F("hmip26","hmip26",20,0,6);
  TH1F *hmip27 = new TH1F("hmip27","hmip27",20,0,6);

  TH1F *hmip28 = new TH1F("hmip28","hmip28",20,0,6);
  TH1F *hmip29 = new TH1F("hmip29","hmip29",20,0,20);
  TH1F *hmip30 = new TH1F("hmip30","hmip30",20,0,20);
  TH1F *hmip31 = new TH1F("hmip31","hmip31",20,0,20);
  TH1F *hmip32 = new TH1F("hmip32","hmip32",20,0,20);

  TH1F *hmip33 = new TH1F("hmip33","hmip33",20,0,6);
  TH1F *hmip34 = new TH1F("hmip34","hmip34",20,0,6);
  TH1F *hmip35 = new TH1F("hmip35","hmip35",20,0,6);
  TH1F *hmip36 = new TH1F("hmip36","hmip36",100,0,100);
  TH1F *hmip37 = new TH1F("hmip37","hmip37",100,0,100);
  TH1F *hmip38 = new TH1F("hmip38","hmip38",100,0,100);
  TH1F *hmip39 = new TH1F("hmip39","hmip39",100,0,180);
  TH1F *hmip40 = new TH1F("hmip40","hmip40",100,0,90);

  TH1F *hmip41 = new TH1F("hmip41","hmip41",20,0,20);
  TH1F *hmip42 = new TH1F("hmip42","hmip42",20,0,20);
  TH1F *hmip43 = new TH1F("hmip43","hmip43",100,0,500);
  TH1F *hmip44 = new TH1F("hmip44","hmip44",5,-3,2);
  TH1F *hmip45 = new TH1F("hmip45","hmip45",20,0,6);
  TH1F *hmip46 = new TH1F("hmip46","hmip46",20,0,20);
  TH1F *hmip47 = new TH1F("hmip47","hmip47",20,0,20);
  TH1F *hmip48 = new TH1F("hmip48","hmip48",5, -3, 2);
  TH1F *hmip49 = new TH1F("hmip49","hmip49",20,0,200);
  TH1F *hmip50 = new TH1F("hmip50","hmip50",5,-3,2);

  int pio_flag;
  int pio_mip_id;
  int pio_filled;
  int pio_flag_pio;
  
  int pio_1_flag;
  double pio_1_mass;
  int pio_1_pio_type;
  double pio_1_energy_1;
  double pio_1_energy_2;
  double pio_1_dis_1;
  double pio_1_dis_2;
  
  std::vector<double> *pio_2_v_dis2 = new std::vector<double>;
  std::vector<double> *pio_2_v_angle2 = new std::vector<double>;
  std::vector<double> *pio_2_v_acc_length = new std::vector<double>;
  std::vector<int> *pio_2_v_flag = new std::vector<int>;

  T->SetBranchAddress("pio_flag",&pio_flag);
  T->SetBranchAddress("pio_mip_id",&pio_mip_id);
  T->SetBranchAddress("pio_filled",&pio_filled);
  T->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  T->SetBranchAddress("pio_1_flag",&pio_1_flag);
  T->SetBranchAddress("pio_1_mass",&pio_1_mass);
  T->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
  T->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
  T->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
  T->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
  T->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
  
  T->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  T->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  T->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  T->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);

  TH1F *hpio1 = new TH1F("hpio1","hpio1",5,-3,2);
  TH1F *hpio2 = new TH1F("hpio2","hpio2",5,-3,2);
  TH1F *hpio3 = new TH1F("hpio3","hpio3",5,-3,2);
  TH1F *hpio4 = new TH1F("hpio4","hpio4",5,-3,2);
  TH1F *hpio5 = new TH1F("hpio5","hpio5",100,0,250);
  TH1F *hpio6 = new TH1F("hpio6","hpio6",5,-3,2);
  TH1F *hpio7 = new TH1F("hpio7","hpio7",100,0,1000);
  TH1F *hpio8 = new TH1F("hpio8","hpio8",100,0,1000);
  TH1F *hpio9 = new TH1F("hpio9","hpio9",100,0,100);
  TH1F *hpio10 = new TH1F("hpio10","hpio10",100,0,100);
  TH1F *hpio11 = new TH1F("hpio11","hpio11",100,0,100);
  TH1F *hpio12 = new TH1F("hpio12","hpio12",90,0,90);
  TH1F *hpio13 = new TH1F("hpio13","hpio13",100,0,100);

  TH1F *hpio14 = new TH1F("hpio14","hpio14",5,-3,2);
  TH1F *hpio15 = new TH1F("hpio15","hpio15",5,-3,2);
  
  // stem direction
  int stem_dir_flag;
  int stem_dir_flag_single_shower;
  int stem_dir_filled;
  double stem_dir_angle;
  double stem_dir_energy;
  double stem_dir_angle1;
  double stem_dir_angle2;
  double stem_dir_angle3;
  double stem_dir_ratio;
  
  T->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
  T->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  T->SetBranchAddress("stem_dir_filled",&stem_dir_filled);
  T->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
  T->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
  T->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
  T->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
  T->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
  T->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);
  
  
  TH1F *hsd1 = new TH1F("hsd1","hsd1",5,-3,2);
  TH1F *hsd2 = new TH1F("hsd2","hsd2",5,-3,2);
  TH1F *hsd3 = new TH1F("hsd3","hsd3",5,-3,2);
  TH1F *hsd4 = new TH1F("hsd4","hsd4",90,0,180);
  TH1F *hsd5 = new TH1F("hsd5","hsd5",100,0,2000);
  TH1F *hsd6 = new TH1F("hsd6","hsd6",90,0,180);
  TH1F *hsd7 = new TH1F("hsd7","hsd7",90,0,180);
  TH1F *hsd8 = new TH1F("hsd8","hsd8",90,0,180);
  TH1F *hsd9 = new TH1F("hsd9","hsd9",10,0,1);

  int br_filled;
  int br1_flag;
  
  //bad reconstruction 1_1
  int br1_1_flag;
  int br1_1_shower_type;
  int br1_1_vtx_n_segs;
  double br1_1_energy;
  int br1_1_n_segs;
  int br1_1_flag_sg_topology;
  int br1_1_flag_sg_trajectory;
  double br1_1_sg_length;
  
  // bad reconstruction 1_2
  int br1_2_flag;
  double br1_2_energy;
  int br1_2_n_connected;
  double br1_2_max_length;
  int br1_2_n_connected_1;
  int br1_2_vtx_n_segs;
  int br1_2_n_shower_segs;
  double br1_2_max_length_ratio;
  double br1_2_shower_length;
  
  // bad_reconstruction 1_3
  int br1_3_flag;
  double br1_3_energy;
  int br1_3_n_connected_p;
  double br1_3_max_length_p;
  int br1_3_n_shower_segs;
  int br1_3_flag_sg_topology;
  int br1_3_flag_sg_trajectory;
  int br1_3_n_shower_main_segs;
  double br1_3_sg_length;
  
  T->SetBranchAddress("br_filled",&br_filled);
  T->SetBranchAddress("br1_flag",&br1_flag);
  
  T->SetBranchAddress("br1_1_flag",&br1_1_flag);
  T->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
  T->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
  T->SetBranchAddress("br1_1_energy",&br1_1_energy);
  T->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
  T->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
  T->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
  T->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
  
  T->SetBranchAddress("br1_2_flag",&br1_2_flag);
  T->SetBranchAddress("br1_2_energy",&br1_2_energy);
  T->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
  T->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
  T->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
  T->SetBranchAddress("br1_2_vtx_n_segs",&br1_2_vtx_n_segs);
  T->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
  T->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
  T->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
  
  T->SetBranchAddress("br1_3_flag",&br1_3_flag);
  T->SetBranchAddress("br1_3_energy",&br1_3_energy);
  T->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
  T->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
  T->SetBranchAddress("br1_3_n_shower_segs",&br1_3_n_shower_segs);
  T->SetBranchAddress("br1_3_flag_sg_topology",&br1_3_flag_sg_topology);
  T->SetBranchAddress("br1_3_sg_trajectory",&br1_3_flag_sg_trajectory);
  T->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
  T->SetBranchAddress("br1_3_sg_length",&br1_3_sg_length);

  TH1F *hbr1_1 = new TH1F("hbr1_1","hbr1_1",5,-3,2);
  TH1F *hbr1_2 = new TH1F("hbr1_2","hbr1_2",5,-3,2);
  TH1F *hbr1_3 = new TH1F("hbr1_3","hbr1_3",10,0,20);
  TH1F *hbr1_4 = new TH1F("hbr1_4","hbr1_4",10,0,10);
  TH1F *hbr1_5 = new TH1F("hbr1_5","hbr1_5",100,0,2000);
  TH1F *hbr1_6 = new TH1F("hbr1_6","hbr1_6",10,0,10);
  TH1F *hbr1_7 = new TH1F("hbr1_7","hbr1_7",5,-3,2);
  TH1F *hbr1_8 = new TH1F("hbr1_8","hbr1_8",5,-3,2);
  TH1F *hbr1_9 = new TH1F("hbr1_9","hbr1_9",100,0,300);
  TH1F *hbr1_10 = new TH1F("hbr1_10","hbr1_10",5,-3,2);
  TH1F *hbr1_11 = new TH1F("hbr1_11","hbr1_11",100,0,2000);
  TH1F *hbr1_12 = new TH1F("hbr1_12","hbr1_12",5,0,5);
  TH1F *hbr1_13 = new TH1F("hbr1_13","hbr1_13",100,0,300);
  TH1F *hbr1_14 = new TH1F("hbr1_14","hbr1_14",10,0,30);
  TH1F *hbr1_15 = new TH1F("hbr1_15","hbr1_15",5,0,10);
  TH1F *hbr1_16 = new TH1F("hbr1_16","hbr1_16",20,0,100);
  TH1F *hbr1_17 = new TH1F("hbr1_17","hbr1_17",10,0,1);
  TH1F *hbr1_18 = new TH1F("hbr1_18","hbr1_18",100,0,300);
  TH1F *hbr1_19 = new TH1F("hbr1_19","hbr1_19",5,-3,2);
  TH1F *hbr1_20 = new TH1F("hbr1_20","hbr1_20",100,0,2000);
  TH1F *hbr1_21 = new TH1F("hbr1_21","hbr1_21",10,0,30);
  TH1F *hbr1_22 = new TH1F("hbr1_22","hbr1_22",100,0,300);
  TH1F *hbr1_23 = new TH1F("hbr1_23","hbr1_23",10,0,20);
  TH1F *hbr1_24 = new TH1F("hbr1_24","hbr1_24",5,-3,2);
  TH1F *hbr1_25 = new TH1F("hbr1_25","hbr1_25",5,-3,2);
  TH1F *hbr1_26 = new TH1F("hbr1_26","hbr1_26",5,0,10);
  TH1F *hbr1_27 = new TH1F("hbr1_27","hbr1_27",20,0,100);

  // bad reconstruction 2
  int br2_flag;
  int br2_flag_single_shower;
  int br2_num_valid_tracks;
  double br2_energy;
  double br2_angle1;
  double br2_angle2;
  double br2_angle;
  double br2_angle3;
  int br2_n_shower_main_segs;
  double br2_max_angle;
  double br2_sg_length;
  int br2_flag_sg_trajectory;
  
  T->SetBranchAddress("br2_flag",&br2_flag);
  T->SetBranchAddress("br2_flag_single_shower",&br2_flag_single_shower);
  T->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
  T->SetBranchAddress("br2_energy",&br2_energy);
  T->SetBranchAddress("br2_angle1",&br2_angle1);
  T->SetBranchAddress("br2_angle2",&br2_angle2);
  T->SetBranchAddress("br2_angle",&br2_angle);
  T->SetBranchAddress("br2_angle3",&br2_angle3);
  T->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
  T->SetBranchAddress("br2_max_angle",&br2_max_angle);
  T->SetBranchAddress("br2_sg_length",&br2_sg_length);
  T->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);
  
  TH1F *hbr2_1 = new TH1F("hbr2_1","hbr2_1",5,-3,2);
  TH1F *hbr2_2 = new TH1F("hbr2_2","hbr2_2",5,-3,2);
  TH1F *hbr2_3 = new TH1F("hbr2_3","hbr2_3",5,0,10);
  TH1F *hbr2_4 = new TH1F("hbr2_4","hbr2_4",100,0,2000);
  TH1F *hbr2_5 = new TH1F("hbr2_5","hbr2_5",90,0,180);
  TH1F *hbr2_6 = new TH1F("hbr2_6","hbr2_6",90,0,180);
  TH1F *hbr2_7 = new TH1F("hbr2_7","hbr2_7",90,0,180);
  TH1F *hbr2_8 = new TH1F("hbr2_8","hbr2_8",90,0,180);
  TH1F *hbr2_9 = new TH1F("hbr2_9","hbr2_9",5,0,10);
  TH1F *hbr2_10 = new TH1F("hbr2_10","hbr2_10",90,0,180);
  TH1F *hbr2_11 = new TH1F("hbr2_11","hbr2_11",100,0,300);
  TH1F *hbr2_12 = new TH1F("hbr2_12","hbr2_12",5,-3,2);

  // low energy overlap
  int lol_flag;

  std::vector<double> *lol_1_v_energy = new std::vector<double>;
  std::vector<int> *lol_1_v_vtx_n_segs= new std::vector<int>;
  std::vector<int> *lol_1_v_nseg= new std::vector<int>;
  std::vector<double> *lol_1_v_angle= new std::vector<double>;
  std::vector<int> *lol_1_v_flag= new std::vector<int>;
  
  std::vector<double> *lol_2_v_length= new std::vector<double>;
  std::vector<double> *lol_2_v_angle= new std::vector<double>;
  std::vector<int> *lol_2_v_type= new std::vector<int>;
  std::vector<int> *lol_2_v_vtx_n_segs= new std::vector<int>;
  std::vector<double> *lol_2_v_energy= new std::vector<double>;
  std::vector<double> *lol_2_v_shower_main_length= new std::vector<double>;
  std::vector<int> *lol_2_v_flag_dir_weak= new std::vector<int>;
  std::vector<int> *lol_2_v_flag = new std::vector<int>;
  
  double lol_3_angle_beam;
  int lol_3_n_valid_tracks;
  double lol_3_min_angle;
  int lol_3_vtx_n_segs;
  double lol_3_energy;
  double lol_3_shower_main_length;
  int lol_3_n_out;
  int lol_3_n_sum;    
  int lol_3_flag;

  T->SetBranchAddress("lol_flag",&lol_flag);
  
  T->SetBranchAddress("lol_1_v_energy",&lol_1_v_energy);
  T->SetBranchAddress("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
  T->SetBranchAddress("lol_1_v_nseg",&lol_1_v_nseg);
  T->SetBranchAddress("lol_1_v_angle",&lol_1_v_angle);
  T->SetBranchAddress("lol_1_v_flag",&lol_1_v_flag);
  
  T->SetBranchAddress("lol_2_v_length",&lol_2_v_length);
  T->SetBranchAddress("lol_2_v_angle",&lol_2_v_angle);
  T->SetBranchAddress("lol_2_v_type",&lol_2_v_type);
  T->SetBranchAddress("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
  T->SetBranchAddress("lol_2_v_energy",&lol_2_v_energy);
  T->SetBranchAddress("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
  T->SetBranchAddress("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);
  T->SetBranchAddress("lol_2_v_flag",&lol_2_v_flag);
  
  T->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
  T->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
  T->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
  T->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
  T->SetBranchAddress("lol_3_energy",&lol_3_energy);
  T->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
  T->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
  T->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
  T->SetBranchAddress("lol_3_flag",&lol_3_flag);
  
  TH1F *hlol_1 = new TH1F("hlol_1","hlol_1",5,-3,2);
  TH1F *hlol_2 = new TH1F("hlol_2","hlol_2",100,0,2000);
  TH1F *hlol_3 = new TH1F("hlol_3","hlol_3",5,0,10);
  TH1F *hlol_4 = new TH1F("hlol_4","hlol_4",5,0,10);
  TH1F *hlol_5 = new TH1F("hlol_5","hlol_5",90,0,180);
  TH1F *hlol_6 = new TH1F("hlol_6","hlol_6",5,-3,2);
  
  TH1F *hlol_7 = new TH1F("hlol_7","hlol_7",20,0,100);
  TH1F *hlol_8 = new TH1F("hlol_8","hlol_8",90,0,180);
  TH1F *hlol_9 = new TH1F("hlol_9","hlol_9",10,0,2000);
  TH1F *hlol_10 = new TH1F("hlol_10","hlol_10",5,0,20);
  TH1F *hlol_11 = new TH1F("hlol_11","hlol_11",100,0,2000);
  TH1F *hlol_12 = new TH1F("hlol_12","hlol_12",90,0,200);
  TH1F *hlol_13 = new TH1F("hlol_13","hlol_13",5,-3,2);
  TH1F *hlol_14 = new TH1F("hlol_14","hlol_14",5,-3,2);
  
  TH1F *hlol_15 = new TH1F("hlol_15","hlol_15",90,0,180);
  TH1F *hlol_16 = new TH1F("hlol_16","hlol_16",5,0,20);
  TH1F *hlol_17 = new TH1F("hlol_17","hlol_17",90,0,180);
  TH1F *hlol_18 = new TH1F("hlol_18","hlol_18",90,0,20);
  TH1F *hlol_19 = new TH1F("hlol_19","hlol_19",100,0,2000);
  TH1F *hlol_20 = new TH1F("hlol_20","hlol_20",10,0,200);
  TH1F *hlol_21 = new TH1F("hlol_21","hlol_21",10,0,200);
  TH1F *hlol_22 = new TH1F("hlol_22","hlol_22",10,0,200);
  TH1F *hlol_23 = new TH1F("hlol_23","hlol_23",5,-3,2);

  //bad reconstruction 3
  double br3_1_energy;
  int br3_1_n_shower_segments;
  int br3_1_sg_flag_trajectory;
  double br3_1_sg_direct_length;
  double br3_1_sg_length;
  double br3_1_total_main_length;
  double br3_1_total_length;
  double br3_1_iso_angle;
  int br3_1_sg_flag_topology;
  int br3_1_flag;

  int br3_2_n_ele;
  int br3_2_n_other;
  double br3_2_energy;
  double br3_2_total_main_length;
  double br3_2_total_length;
  int br3_2_other_fid;
  int br3_2_flag;
  
  std::vector<double> *br3_3_v_energy = new std::vector<double>;
  std::vector<double> *br3_3_v_angle = new std::vector<double>;
  std::vector<double> *br3_3_v_dir_length = new std::vector<double>;
  std::vector<double> *br3_3_v_length = new std::vector<double>;
  std::vector<int> *br3_3_v_flag = new std::vector<int>;
  
  double br3_4_acc_length;
  double br3_4_total_length;
  double br3_4_energy;
  int br3_4_flag;
  
  std::vector<double> *br3_5_v_dir_length = new std::vector<double>;
  std::vector<double> *br3_5_v_total_length = new std::vector<double>;
  std::vector<int> *br3_5_v_flag_avoid_muon_check = new std::vector<int>;
  std::vector<int> *br3_5_v_n_seg = new std::vector<int>;
  std::vector<double> *br3_5_v_angle = new std::vector<double>;
  std::vector<double> *br3_5_v_sg_length = new std::vector<double>;
  std::vector<double> *br3_5_v_energy = new std::vector<double>;
  std::vector<int> *br3_5_v_n_main_segs = new std::vector<int>;
  std::vector<int> *br3_5_v_n_segs = new std::vector<int>;
  std::vector<double> *br3_5_v_shower_main_length = new std::vector<double>;
  std::vector<double> *br3_5_v_shower_total_length = new std::vector<double>;
  std::vector<int> *br3_5_v_flag = new std::vector<int>;
  
  std::vector<double> *br3_6_v_angle = new std::vector<double>;
  std::vector<double> *br3_6_v_angle1 = new std::vector<double>;
  std::vector<int> *br3_6_v_flag_shower_trajectory = new std::vector<int>;
  std::vector<double> *br3_6_v_direct_length = new std::vector<double>;
  std::vector<double> *br3_6_v_length = new std::vector<double>;
  std::vector<int> *br3_6_v_n_other_vtx_segs = new std::vector<int>;
  std::vector<double> *br3_6_v_energy = new std::vector<double>;
  std::vector<int> *br3_6_v_flag = new std::vector<int>;
  
  double br3_7_energy;
  double br3_7_min_angle;
  double br3_7_sg_length;
  double br3_7_shower_main_length;
  int br3_7_flag;
  
  double br3_8_max_dQ_dx;
  double br3_8_energy;
  int br3_8_n_main_segs;
  double br3_8_shower_main_length;
  double br3_8_shower_length;
  int br3_8_flag;
  
  int br3_flag;
  
  TH1F *hbr3_1 = new TH1F("hbr3_1","hbr3_1",100,0,2000);
  TH1F *hbr3_2 = new TH1F("hbr3_2","hbr3_2",15,0,30);
  TH1F *hbr3_3 = new TH1F("hbr3_3","hbr3_3",5,-3,2);
  TH1F *hbr3_4 = new TH1F("hbr3_4","hbr3_4",20,0,100);
  TH1F *hbr3_5 = new TH1F("hbr3_5","hbr3_5",20,0,100);
  TH1F *hbr3_6 = new TH1F("hbr3_6","hbr3_6",30,0,120);
  TH1F *hbr3_7 = new TH1F("hbr3_7","hbr3_7",30,0,120);
  TH1F *hbr3_8 = new TH1F("hbr3_8","hbr3_8",90,0,90);
  TH1F *hbr3_9 = new TH1F("hbr3_9","hbr3_9",5,-3,2);
  TH1F *hbr3_10 = new TH1F("hbr3_10","hbr3_10",5,-3,2);

  TH1F *hbr3_11 = new TH1F("hbr3_11","hbr3_11",5,0,10);
  TH1F *hbr3_12 = new TH1F("hbr3_12","hbr3_12",5,0,10);
  TH1F *hbr3_13 = new TH1F("hbr3_13","hbr3_13",5,0,2000);
  TH1F *hbr3_14 = new TH1F("hbr3_14","hbr3_14",20,0,100);
  TH1F *hbr3_15 = new TH1F("hbr3_15","hbr3_15",20,0,100);
  TH1F *hbr3_16 = new TH1F("hbr3_16","hbr3_16",5,-3,2);
  TH1F *hbr3_17 = new TH1F("hbr3_17","hbr3_17",5,-3,2);

  TH1F *hbr3_18 = new TH1F("hbr3_18","hbr3_18",100,0,2000);
  TH1F *hbr3_19 = new TH1F("hbr3_19","hbr3_19",90,0,180);
  TH1F *hbr3_20 = new TH1F("hbr3_20","hbr3_20",20,0,100);
  TH1F *hbr3_21 = new TH1F("hbr3_21","hbr3_21",20,0,100);
  TH1F *hbr3_22 = new TH1F("hbr3_22","hbr3_22",5,-3,2);

  TH1F *hbr3_23 = new TH1F("hbr3_23","hbr3_23",20,0,100);
  TH1F *hbr3_24 = new TH1F("hbr3_24","hbr3_24",20,0,100);
  TH1F *hbr3_25 = new TH1F("hbr3_25","hbr3_25",100,0,2000);
  TH1F *hbr3_26 = new TH1F("hbr3_26","hbr3_26",5,-3,2);

  TH1F *hbr3_27 = new TH1F("hbr3_27","hbr3_27",20,0,100);
  TH1F *hbr3_28 = new TH1F("hbr3_28","hbr3_28",20,0,100);
  TH1F *hbr3_29 = new TH1F("hbr3_29","hbr3_29",5,-3,2);
  TH1F *hbr3_30 = new TH1F("hbr3_30","hbr3_30",5,0,20);
  TH1F *hbr3_31 = new TH1F("hbr3_31","hbr3_31",90,0,180);
  TH1F *hbr3_32 = new TH1F("hbr3_32","hbr3_32",20,0,100);
  TH1F *hbr3_33 = new TH1F("hbr3_33","hbr3_33",100,0,2000);
  TH1F *hbr3_34 = new TH1F("hbr3_34","hbr3_34",20,0,40);
  TH1F *hbr3_35 = new TH1F("hbr3_35","hbr3_35",20,0,40);
  TH1F *hbr3_36 = new TH1F("hbr3_36","hbr3_36",20,0,100);
  TH1F *hbr3_37 = new TH1F("hbr3_37","hbr3_37",20,0,100);
  TH1F *hbr3_38 = new TH1F("hbr3_38","hbr3_38",5,-3,2);

  TH1F *hbr3_39 = new TH1F("hbr3_39","hbr3_39",90,0,180);
  TH1F *hbr3_40 = new TH1F("hbr3_40","hbr3_40",90,0,90);
  TH1F *hbr3_41 = new TH1F("hbr3_41","hbr3_41",5,-3,2);
  TH1F *hbr3_42 = new TH1F("hbr3_42","hbr3_42",20,0,100);
  TH1F *hbr3_43 = new TH1F("hbr3_43","hbr3_43",20,0,100);
  TH1F *hbr3_44 = new TH1F("hbr3_44","hbr3_44",5,0,20);
  TH1F *hbr3_45 = new TH1F("hbr3_45","hbr3_45",100,0,2000);
  TH1F *hbr3_46 = new TH1F("hbr3_46","hbr3_46",5,-3,2);

  TH1F *hbr3_47 = new TH1F("hbr3_47","hbr3_47",100,0,2000);
  TH1F *hbr3_48 = new TH1F("hbr3_48","hbr3_48",90,0,180);
  TH1F *hbr3_49 = new TH1F("hbr3_49","hbr3_49",20,0,100);
  TH1F *hbr3_50 = new TH1F("hbr3_50","hbr3_50",20,0,100);
  TH1F *hbr3_51 = new TH1F("hbr3_51","hbr3_51",5,-3,2);

  TH1F *hbr3_52 = new TH1F("hbr3_52","hbr3_52",40,0,10);
  TH1F *hbr3_53 = new TH1F("hbr3_53","hbr3_53",100,0,2000);
  TH1F *hbr3_54 = new TH1F("hbr3_54","hbr3_53",5,0,20);
  TH1F *hbr3_55 = new TH1F("hbr3_55","hbr3_55",20,0,100);
  TH1F *hbr3_56 = new TH1F("hbr3_56","hbr3_56",20,0,100);
  TH1F *hbr3_57 = new TH1F("hbr3_57","hbr3_57",5,-3,2);

  TH1F *hbr3_58 = new TH1F("hbr3_58","hbr3_1",5,-3,2);

  T->SetBranchAddress("br3_1_energy",&br3_1_energy);
  T->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
  T->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
  T->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
  T->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
  T->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
  T->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
  T->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
  T->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
  T->SetBranchAddress("br3_1_flag",&br3_1_flag);
  
  T->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
  T->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
  T->SetBranchAddress("br3_2_energy",&br3_2_energy);
  T->SetBranchAddress("br3_2_total_main_length",&br3_2_total_main_length);
  T->SetBranchAddress("br3_2_total_length",&br3_2_total_length);
  T->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
  T->SetBranchAddress("br3_2_flag",&br3_2_flag);
  
  T->SetBranchAddress("br3_3_v_energy",&br3_3_v_energy);
  T->SetBranchAddress("br3_3_v_angle",&br3_3_v_angle);
  T->SetBranchAddress("br3_3_v_dir_length",&br3_3_v_dir_length);
  T->SetBranchAddress("br3_3_v_length",&br3_3_v_length);
  T->SetBranchAddress("br3_3_v_flag",&br3_3_v_flag);
  
  T->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
  T->SetBranchAddress("br3_4_total_length", &br3_4_total_length);
  T->SetBranchAddress("br3_4_energy", &br3_4_energy);
  T->SetBranchAddress("br3_4_flag", &br3_4_flag);
  
  T->SetBranchAddress("br3_5_v_dir_length", &br3_5_v_dir_length);
  T->SetBranchAddress("br3_5_v_total_length", &br3_5_v_total_length);
  T->SetBranchAddress("br3_5_v_flag_avoid_muon_check", &br3_5_v_flag_avoid_muon_check);
  T->SetBranchAddress("br3_5_v_n_seg", &br3_5_v_n_seg);
  T->SetBranchAddress("br3_5_v_angle", &br3_5_v_angle);
  T->SetBranchAddress("br3_5_v_sg_length", &br3_5_v_sg_length);
  T->SetBranchAddress("br3_5_v_energy", &br3_5_v_energy);
  T->SetBranchAddress("br3_5_v_n_main_segs", &br3_5_v_n_main_segs);
  T->SetBranchAddress("br3_5_v_n_segs", &br3_5_v_n_segs);
  T->SetBranchAddress("br3_5_v_shower_main_length", &br3_5_v_shower_main_length);
  T->SetBranchAddress("br3_5_v_shower_total_length", &br3_5_v_shower_total_length);
  T->SetBranchAddress("br3_5_v_flag", &br3_5_v_flag);
  
  T->SetBranchAddress("br3_6_v_angle",&br3_6_v_angle);
  T->SetBranchAddress("br3_6_v_angle1",&br3_6_v_angle1);
  T->SetBranchAddress("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
  T->SetBranchAddress("br3_6_v_direct_length",&br3_6_v_direct_length);
  T->SetBranchAddress("br3_6_v_length",&br3_6_v_length);
  T->SetBranchAddress("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
  T->SetBranchAddress("br3_6_v_energy",&br3_6_v_energy);
  T->SetBranchAddress("br3_6_v_flag",&br3_6_v_flag);
  
  T->SetBranchAddress("br3_7_energy",&br3_7_energy);
  T->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
  T->SetBranchAddress("br3_7_sg_length",&br3_7_sg_length);
  T->SetBranchAddress("br3_7_main_length",&br3_7_shower_main_length);
  T->SetBranchAddress("br3_7_flag",&br3_7_flag);
  
  T->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
  T->SetBranchAddress("br3_8_energy",&br3_8_energy);
  T->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
  T->SetBranchAddress("br3_8_shower_main_length",&br3_8_shower_main_length);
  T->SetBranchAddress("br3_8_shower_length",&br3_8_shower_length);
  T->SetBranchAddress("br3_8_flag",&br3_8_flag);
  
  T->SetBranchAddress("br3_flag",&br3_flag);


  // BR 4
  double br4_1_shower_main_length;
  double br4_1_shower_total_length;
  double br4_1_min_dis;
  double br4_1_energy;
  int br4_1_flag_avoid_muon_check;
  int br4_1_n_vtx_segs;
  int br4_1_n_main_segs;
  int br4_1_flag;
  
  double br4_2_ratio_45;
  double br4_2_ratio_35;
  double br4_2_ratio_25;
  double br4_2_ratio_15;
  double br4_2_energy;
  double br4_2_ratio1_45;
  double br4_2_ratio1_35;
  double br4_2_ratio1_25;
  double br4_2_ratio1_15;
  double br4_2_iso_angle;
  double br4_2_iso_angle1;
  double br4_2_angle;
  int br4_2_flag;
  
  int br4_flag;
  
  T->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
  T->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
  T->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
  T->SetBranchAddress("br4_1_energy", &br4_1_energy);
  T->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
  T->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
  T->SetBranchAddress("br4_1_n_main_segs", &br4_1_n_main_segs);
  T->SetBranchAddress("br4_1_flag", &br4_1_flag);
  
  T->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
  T->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
  T->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
  T->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
  T->SetBranchAddress("br4_2_energy",   &br4_2_energy);
  T->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
  T->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
  T->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
  T->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
  T->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
  T->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
  T->SetBranchAddress("br4_2_angle", &br4_2_angle);
  T->SetBranchAddress("br4_2_flag", &br4_2_flag);
  
  T->SetBranchAddress("br4_flag", &br4_flag);
  
  TH1F *hbr4_1 = new TH1F("hbr4_1","hbr4_1",20,0,100);
  TH1F *hbr4_2 = new TH1F("hbr4_2","hbr4_2",20,0,100);
  TH1F *hbr4_3 = new TH1F("hbr4_3","hbr4_3",20,0,100);
  TH1F *hbr4_4 = new TH1F("hbr4_4","hbr4_4",100,0,2000);
  TH1F *hbr4_5 = new TH1F("hbr4_5","hbr4_5",5,-3,2);
  TH1F *hbr4_6 = new TH1F("hbr4_6","hbr4_6",20,0,20);
  TH1F *hbr4_7 = new TH1F("hbr4_7","hbr4_7",20,0,20);
  TH1F *hbr4_8 = new TH1F("hbr4_8","hbr4_8",5,-3,2);

  TH1F *hbr4_9 = new TH1F("hbr4_9","hbr4_9",20,0,1);
  TH1F *hbr4_10 = new TH1F("hbr4_10","hbr4_10",20,0,1);
  TH1F *hbr4_11 = new TH1F("hbr4_11","hbr4_11",20,0,1);
  TH1F *hbr4_12 = new TH1F("hbr4_12","hbr4_12",20,0,1);
  TH1F *hbr4_13 = new TH1F("hbr4_13","hbr4_13",100,0,2000);
  TH1F *hbr4_14 = new TH1F("hbr4_14","hbr4_14",20,0,1);
  TH1F *hbr4_15 = new TH1F("hbr4_15","hbr4_15",20,0,1);
  TH1F *hbr4_16 = new TH1F("hbr4_16","hbr4_16",20,0,1);
  TH1F *hbr4_17 = new TH1F("hbr4_17","hbr4_17",20,0,1);
  TH1F *hbr4_18 = new TH1F("hbr4_18","hbr4_18",45,0,90);
  TH1F *hbr4_19 = new TH1F("hbr4_19","hbr4_19",45,0,90);
  TH1F *hbr4_20 = new TH1F("hbr4_20","hbr4_20",90,0,180);
  TH1F *hbr4_21 = new TH1F("hbr4_21","hbr4_21",5,-3,2);

  TH1F *hbr4_22 = new TH1F("hbr4_22","hbr4_22",5,-3,2);
    
  int hol_1_n_valid_tracks;
  double hol_1_min_angle;
  double hol_1_energy;
  int hol_1_flag_all_shower;
  double hol_1_min_length;
  int hol_1_flag;
  
  double hol_2_min_angle;
  double hol_2_medium_dQ_dx;
  int hol_2_ncount;
  double hol_2_energy;
  int hol_2_flag;
  
  int hol_flag;

  T->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
  T->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
  T->SetBranchAddress("hol_1_energy", &hol_1_energy);
  T->SetBranchAddress("hol_1_all_shower", &hol_1_flag_all_shower);
  T->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
  T->SetBranchAddress("hol_1_flag", &hol_1_flag);
  
  T->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
  T->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
  T->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
  T->SetBranchAddress("hol_2_energy", &hol_2_energy);
  T->SetBranchAddress("hol_2_flag", &hol_2_flag);
  
  T->SetBranchAddress("hol_flag", &hol_flag);

  TH1F *hol_1 = new TH1F("hol_1","hol_1",5,0,10);
  TH1F *hol_2 = new TH1F("hol_2","hol_2",90,0,200);
  TH1F *hol_3 = new TH1F("hol_3","hol_3",100,0,2000);
  TH1F *hol_4 = new TH1F("hol_4","hol_4",5,-3,2);
  TH1F *hol_5 = new TH1F("hol_5","hol_5",20,0,100);
  TH1F *hol_6 = new TH1F("hol_6","hol_6",5,-3,2);

  TH1F *hol_7 = new TH1F("hol_7","hol_7",90,0,90);
  TH1F *hol_8 = new TH1F("hol_8","hol_8",20,0,5);
  TH1F *hol_9 = new TH1F("hol_9","hol_9",5,0,10);
  TH1F *hol_10 = new TH1F("hol_10","hol_10",100,0,2000);
  TH1F *hol_11 = new TH1F("hol_11","hol_11",5,-3,2);
  
  TH1F *hol_12 = new TH1F("hol_12","hol_12",5,-3,2);
  
  // vertex inside shower
  int vis_1_filled;
  int vis_1_n_vtx_segs;
  double vis_1_energy;
  int vis_1_num_good_tracks;
  double vis_1_max_angle;
  double vis_1_max_shower_angle;
  double vis_1_tmp_length1;
  double vis_1_tmp_length2;
  double vis_1_particle_type;
  int vis_1_flag;
  
  int vis_2_filled;
  int vis_2_n_vtx_segs;
  double vis_2_min_angle;
  int vis_2_min_weak_track;
  double vis_2_angle_beam;
  double vis_2_min_angle1;
  double vis_2_iso_angle1;
  double vis_2_min_medium_dQ_dx;
  double vis_2_min_length;
  double vis_2_sg_length;
  double vis_2_max_angle;
  int vis_2_max_weak_track;
  int vis_2_flag;
  
  int vis_flag;

  T->SetBranchAddress("vis_1_filled",&vis_1_filled);
  T->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
  T->SetBranchAddress("vis_1_energy",&vis_1_energy);
  T->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
  T->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
  T->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
  T->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
  T->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
  T->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);
  T->SetBranchAddress("vis_1_flag",&vis_1_flag);
  
  T->SetBranchAddress("vis_2_filled",&vis_2_filled);
  T->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
  T->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
  T->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
  T->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
  T->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
  T->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
  T->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
  T->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
  T->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
  T->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
  T->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
  T->SetBranchAddress("vis_2_flab",&vis_2_flag);
  
  T->SetBranchAddress("vis_flag",&vis_flag);

  TH1F *hvis_1 = new TH1F("hvis_1","hvis_1",5,-3,2);
  TH1F *hvis_2 = new TH1F("hvis_2","hvis_2",5,0,10);
  TH1F *hvis_3 = new TH1F("hvis_3","hvis_3",100,0,2000);
  TH1F *hvis_4 = new TH1F("hvis_4","hvis_4",5,0,10);
  TH1F *hvis_5 = new TH1F("hvis_5","hvis_5",90,0,180);
  TH1F *hvis_6 = new TH1F("hvis_6","hvis_6",20,0,100);
  TH1F *hvis_7 = new TH1F("hvis_7","hvis_7",20,0,100);
  TH1F *hvis_8 = new TH1F("hvis_8","hvis_8",20,0,100);
  TH1F *hvis_9 = new TH1F("hvis_9","hvis_9",5,-3,2);

  TH1F *hvis_10 = new TH1F("hvis_10","hvis_10",5,-3,2);
  TH1F *hvis_11 = new TH1F("hvis_11","hvis_11",5,0,20);
  TH1F *hvis_12 = new TH1F("hvis_12","hvis_12",90,0,180);
  TH1F *hvis_13 = new TH1F("hvis_13","hvis_13",5,-3,2);
  TH1F *hvis_14 = new TH1F("hvis_14","hvis_14",90,0,180);
  TH1F *hvis_15 = new TH1F("hvis_15","hvis_15",90,0,180);
  TH1F *hvis_16 = new TH1F("hvis_16","hvis_16",90,0,180);
  TH1F *hvis_17 = new TH1F("hvis_17","hvis_17",10,0,5);
  TH1F *hvis_18 = new TH1F("hvis_18","hvis_18",20,0,100);
  TH1F *hvis_19 = new TH1F("hvis_19","hvis_19",20,0,100);
  TH1F *hvis_20 = new TH1F("hvis_20","hvis_20",90,0,180);
  TH1F *hvis_21 = new TH1F("hvis_21","hvis_21",5,-3,2);
  TH1F *hvis_22 = new TH1F("hvis_22","hvis_22",5,-3,2);

  TH1F *hvis_23 = new TH1F("hvis_23","hvis_23",5,-3,2);
  
  double stem_len_energy;
  double stem_len_length;
  int stem_len_flag_avoid_muon_check;
  int stem_len_num_daughters;
  double stem_len_daughter_length;
  int stem_len_flag;

  T->SetBranchAddress("stem_len_energy", &stem_len_energy);
  T->SetBranchAddress("stem_len_length", &stem_len_length);
  T->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
  T->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
  T->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);
  T->SetBranchAddress("stem_len_flag", &stem_len_flag);

  TH1F *hstl_1 = new TH1F("hstl_1","hstl_1", 100, 0, 2000);
  TH1F *hstl_2 = new TH1F("hstl_2","hstl_2", 20, 0, 100);
  TH1F *hstl_3 = new TH1F("hstl_3","hstl_3", 5, -3, 2);
  TH1F *hstl_4 = new TH1F("hstl_4","hstl_4", 20, 0, 20);
  TH1F *hstl_5 = new TH1F("hstl_5","hstl_5", 100, 0, 200);
  TH1F *hstl_6 = new TH1F("hstl_6","hstl_6", 5, -3, 2);

  int brm_n_mu_segs;
  double brm_Ep;
  double brm_energy;
  double brm_acc_length;
  double brm_shower_total_length;
  double brm_connected_length;
  int brm_n_size;
  double brm_acc_direct_length;
  int brm_n_shower_main_segs;
  int brm_n_mu_main;
  int brm_flag;

  T->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
  T->SetBranchAddress("brm_Ep",&brm_Ep);
  T->SetBranchAddress("brm_energy",&brm_energy);
  T->SetBranchAddress("brm_acc_length",&brm_acc_length);
  T->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
  T->SetBranchAddress("brm_connected_length",&brm_connected_length);
  T->SetBranchAddress("brm_n_size",&brm_n_size);
  T->SetBranchAddress("brm_nacc_direct_length",&brm_acc_direct_length);
  T->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
  T->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
  T->SetBranchAddress("brm_n_flag",&brm_flag);

  TH1F *hbrm_1 = new TH1F("brm_1","brm_1",20,0,20);
  TH1F *hbrm_2 = new TH1F("brm_2","brm_2",100,0,2000);
  TH1F *hbrm_3 = new TH1F("brm_3","brm_3",100,0,2000);
  TH1F *hbrm_4 = new TH1F("brm_4","brm_4",20,0,100);
  TH1F *hbrm_5 = new TH1F("brm_5","brm_5",20,0,100);
  TH1F *hbrm_6 = new TH1F("brm_6","brm_6",20,0,100);
  TH1F *hbrm_7 = new TH1F("brm_7","brm_7",20,0,20);
  TH1F *hbrm_8 = new TH1F("brm_8","brm_8",20,0,100);
  TH1F *hbrm_9 = new TH1F("brm_9","brm_9",20,0,20);
  TH1F *hbrm_10 = new TH1F("brm_10","brm_10",20,0,20);
  TH1F *hbrm_11 = new TH1F("brm_11","brm_11",5,-3,2);

  // compare with muon
  double cme_mu_energy;
  double cme_energy;
  double cme_mu_length;
  double cme_length;
  double cme_angle_beam;
  int cme_flag;

  T->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
  T->SetBranchAddress("cme_energy",&cme_energy);
  T->SetBranchAddress("cme_mu_length",&cme_mu_length);
  T->SetBranchAddress("cme_length",&cme_length);
  T->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
  T->SetBranchAddress("cme_flag",&cme_flag);

  TH1F *hcme_1 = new TH1F("hcme_1","hcme_1",100,0,2000);
  TH1F *hcme_2 = new TH1F("hcme_2","hcme_2",100,0,2000);
  TH1F *hcme_3 = new TH1F("hcme_3","hcme_3",100,0,200);
  TH1F *hcme_4 = new TH1F("hcme_4","hcme_4",100,0,200);
  TH1F *hcme_5 = new TH1F("hcme_5","hcme_5",100,0,180);
  TH1F *hcme_6 = new TH1F("hcme_6","hcme_6",5,-3,2);
  
  // angular cut
  double anc_energy;
  double anc_angle;
  double anc_max_angle;
  double anc_max_length;
  double anc_acc_forward_length;
  double anc_acc_backward_length;
  double anc_acc_forward_length1;
  double anc_shower_main_length;
  double anc_shower_total_length;
  int anc_flag_main_outside;
  int anc_flag;

  T->SetBranchAddress("anc_energy",&anc_energy);
  T->SetBranchAddress("anc_angle",&anc_angle);
  T->SetBranchAddress("anc_max_angle",&anc_max_angle);
  T->SetBranchAddress("anc_max_length",&anc_max_length);
  T->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
  T->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
  T->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
  T->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
  T->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
  T->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
  T->SetBranchAddress("anc_flag",&anc_flag);

  TH1F *hanc_1 = new TH1F("hanc_1","hanc_1",100,0,2000);
  TH1F *hanc_2 = new TH1F("hanc_2","hanc_2",100,0,180);
  TH1F *hanc_3 = new TH1F("hanc_3","hanc_3",100,0,180);
  TH1F *hanc_4 = new TH1F("hanc_4","hanc_4",20,0,200);
  TH1F *hanc_5 = new TH1F("hanc_5","hanc_5",20,0,200);
  TH1F *hanc_6 = new TH1F("hanc_6","hanc_6",20,0,200);
  TH1F *hanc_7 = new TH1F("hanc_7","hanc_7",20,0,200);
  TH1F *hanc_8 = new TH1F("hanc_8","hanc_8",20,0,200);
  TH1F *hanc_9 = new TH1F("hanc_9","hanc_9",20,0,200);
  TH1F *hanc_10 = new TH1F("hanc_10","hanc_10",5,-3,2);
  TH1F *hanc_11 = new TH1F("hanc_11","hanc_11",5,-3,2);
  
  // low-energy michel
  double lem_shower_total_length;
  double lem_shower_main_length;
  int lem_n_3seg;
  double lem_e_charge;
  double lem_e_dQdx;
  int lem_shower_num_segs;
  int lem_shower_num_main_segs;
  int lem_flag;
  
  T->SetBranchAddress("lem_shower_total_length",&lem_shower_total_length);
  T->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
  T->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
  T->SetBranchAddress("lem_e_charge",&lem_e_charge);
  T->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
  T->SetBranchAddress("lem_shower_num_segs",&lem_shower_num_segs);
  T->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  T->SetBranchAddress("lem_flag",&lem_flag);

  TH1F *hlem_1 = new TH1F("hlem_1","hlem_1",20,0,200);
  TH1F *hlem_2 = new TH1F("hlem_2","hlem_2",20,0,200);
  TH1F *hlem_3 = new TH1F("hlem_3","hlem_3",20,0,20);
  TH1F *hlem_4 = new TH1F("hlem_4","hlem_4",20,0,400);
  TH1F *hlem_5 = new TH1F("hlem_5","hlem_5",20,0,400);
  TH1F *hlem_6 = new TH1F("hlem_6","hlem_6",20,0,50);
  TH1F *hlem_7 = new TH1F("hlem_7","hlem_7",20,0,50);
  TH1F *hlem_8 = new TH1F("hlem_8","hlem_8",5,-3,2);

  // shower to wall
  double stw_1_energy;
  double stw_1_dis;
  double stw_1_dQ_dx;
  int stw_1_flag_single_shower;
  int stw_1_n_pi0;
  int stw_1_num_valid_tracks;
  int stw_1_flag;
  
  std::vector<double> *stw_2_v_medium_dQ_dx = new std::vector<double>;
  std::vector<double> *stw_2_v_energy = new std::vector<double>;
  std::vector<double> *stw_2_v_angle = new std::vector<double>;
  std::vector<double> *stw_2_v_dir_length = new std::vector<double>;
  std::vector<double> *stw_2_v_max_dQ_dx = new std::vector<double>;
  std::vector<int> *stw_2_v_flag = new std::vector<int>;
  
  std::vector<double> *stw_3_v_angle = new std::vector<double>;
  std::vector<double> *stw_3_v_dir_length = new std::vector<double>;
  std::vector<double> *stw_3_v_energy = new std::vector<double>;
  std::vector<double> *stw_3_v_medium_dQ_dx = new std::vector<double>;
  std::vector<int> *stw_3_v_flag = new std::vector<int>;
  
  std::vector<double> *stw_4_v_angle = new std::vector<double>;
  std::vector<double> *stw_4_v_dis = new std::vector<double>;
  std::vector<double> *stw_4_v_energy = new std::vector<double>;
  std::vector<int> *stw_4_v_flag = new std::vector<int>;
  
  int stw_flag;
  
  T->SetBranchAddress("stw_1_energy",&stw_1_energy);
  T->SetBranchAddress("stw_1_dis",&stw_1_dis);
  T->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
  T->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
  T->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
  T->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
  T->SetBranchAddress("stw_1_flag",&stw_1_flag);
  
  T->SetBranchAddress("stw_2_v_medium_dQ_dx", &stw_2_v_medium_dQ_dx);
  T->SetBranchAddress("stw_2_v_energy", &stw_2_v_energy);
  T->SetBranchAddress("stw_2_v_angle", &stw_2_v_angle);
  T->SetBranchAddress("stw_2_v_dir_length", &stw_2_v_dir_length);
  T->SetBranchAddress("stw_2_v_max_dQ_dx", &stw_2_v_max_dQ_dx);
  T->SetBranchAddress("stw_2_v_flag", &stw_2_v_flag);
  
  T->SetBranchAddress("stw_3_v_angle",&stw_3_v_angle);
  T->SetBranchAddress("stw_3_v_dir_length",&stw_3_v_dir_length);
  T->SetBranchAddress("stw_3_v_energy",&stw_3_v_energy);
  T->SetBranchAddress("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);
  T->SetBranchAddress("stw_3_v_flag",&stw_3_v_flag);
  
  T->SetBranchAddress("stw_4_v_angle",&stw_4_v_angle);
  T->SetBranchAddress("stw_4_v_dis",&stw_4_v_dis);
  T->SetBranchAddress("stw_4_v_energy",&stw_4_v_energy);
  T->SetBranchAddress("stw_4_v_flag",&stw_4_v_flag);
  
  T->SetBranchAddress("stw_flag", &stw_flag);

  TH1F *hstw_1 = new TH1F("hstw_1","hstw_1",100,0,2000);
  TH1F *hstw_2 = new TH1F("hstw_2","hstw_2",20,0,300);
  TH1F *hstw_3 = new TH1F("hstw_3","hstw_3",20,0,5);
  TH1F *hstw_4 = new TH1F("hstw_4","hstw_4",5,-3,2);
  TH1F *hstw_5 = new TH1F("hstw_5","hstw_5",20,0,20);
  TH1F *hstw_6 = new TH1F("hstw_6","hstw_6",20,0,20);
  TH1F *hstw_7 = new TH1F("hstw_7","hstw_7",5,-3,2);

  TH1F *hstw_8 = new TH1F("hstw_8","hstw_8",20,0,5);
  TH1F *hstw_9 = new TH1F("hstw_9","hstw_9",100,0,2000);
  TH1F *hstw_10 = new TH1F("hstw_10","hstw_10",100,0,180);
  TH1F *hstw_11 = new TH1F("hstw_11","hstw_11",20,0,100);
  TH1F *hstw_12 = new TH1F("hstw_12","hstw_12",20, 0, 10);
  TH1F *hstw_13 = new TH1F("hstw_13","hstw_13",5, -3, 2);
  
  TH1F *hstw_14 = new TH1F("hstw_14","hstw_14",100,0,180);
  TH1F *hstw_15 = new TH1F("hstw_15","hstw_15",20,0,100);
  TH1F *hstw_16 = new TH1F("hstw_16","hstw_16",100,0,2000);
  TH1F *hstw_17 = new TH1F("hstw_17","hstw_17",100,0,5);
  TH1F *hstw_18 = new TH1F("hstw_18","hstw_18",5, -3, 2);

  TH1F *hstw_19 = new TH1F("hstw_19","hstw_19",100,0,180);
  TH1F *hstw_20 = new TH1F("hstw_20","hstw_20",20,0,100);
  TH1F *hstw_21 = new TH1F("hstw_21","hstw_21",100,0,2000);
  TH1F *hstw_22 = new TH1F("hstw_22","hstw_22",5, -3, 2);

  TH1F *hstw_23 = new TH1F("hstw_23","hstw_23",5, -3,2);
  
  // single photon  cases ...
  int spt_flag_single_shower;
  double spt_energy;
  double spt_shower_main_length;
  double spt_shower_total_length;
  double spt_angle_beam;
  double spt_angle_vertical;
  double spt_max_dQ_dx;
  double spt_angle_beam_1;
  double spt_angle_drift;
  double spt_angle_drift_1;
  int spt_num_valid_tracks;
  double spt_n_vtx_segs;
  double spt_max_length;
  int spt_flag;

  T->SetBranchAddress("spt_flag_single_shower", &spt_flag_single_shower);
  T->SetBranchAddress("spt_energy", &spt_energy);
  T->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
  T->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
  T->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
  T->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
  T->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
  T->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
  T->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
  T->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
  T->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
  T->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
  T->SetBranchAddress("spt_max_length", &spt_max_length);
  T->SetBranchAddress("spt_flag", &spt_flag);

  TH1F *hspt_1 = new TH1F("hspt_1","hspt_1",5,-3,2);
  TH1F *hspt_2 = new TH1F("hspt_2","hspt_2",100,0,2000);
  TH1F *hspt_3 = new TH1F("hspt_3","hspt_3",20,0,200);
  TH1F *hspt_4 = new TH1F("hspt_4","hspt_4",20,0,200);
  TH1F *hspt_5 = new TH1F("hspt_5","hspt_5",60,0,180);
  TH1F *hspt_6 = new TH1F("hspt_6","hspt_6",60,0,180);
  TH1F *hspt_7 = new TH1F("hspt_7","hspt_7",5,0,10);
  TH1F *hspt_8 = new TH1F("hspt_8","hspt_8",60,0,180);
  TH1F *hspt_9 = new TH1F("hspt_9","hspt_9",60,0,180);
  TH1F *hspt_10 = new TH1F("hspt_10","hspt_10",60,0,180);
  TH1F *hspt_11 = new TH1F("hspt_11","hspt_11",50,0,20);
  TH1F *hspt_12 = new TH1F("hspt_12","hspt_12",5,0,20);
  TH1F *hspt_13 = new TH1F("hspt_13","hspt_13",20,0,200);
  TH1F *hspt_14 = new TH1F("hspt_14","hspt_14",5,-3,2);

  // multiple gamma I
  double mgo_energy;
  double mgo_max_energy;
  double mgo_total_energy;
  int mgo_n_showers;
  double mgo_max_energy_1;
  double mgo_max_energy_2;
  double mgo_total_other_energy;
  int mgo_n_total_showers;
  double mgo_total_other_energy_1;
  int mgo_flag;
  
  T->SetBranchAddress("mgo_energy",&mgo_energy);
  T->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
  T->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
  T->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
  T->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
  T->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
  T->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
  T->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
  T->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
  T->SetBranchAddress("mgo_flag",&mgo_flag);

  TH1F *hmgo_1 = new TH1F("hmgo_1","hmgo_1",100,0,2000);
  TH1F *hmgo_2 = new TH1F("hmgo_2","hmgo_2",100,0,2000);
  TH1F *hmgo_3 = new TH1F("hmgo_3","hmgo_3",100,0,2000);
  TH1F *hmgo_4 = new TH1F("hmgo_4","hmgo_4",20,0,20);
  TH1F *hmgo_5 = new TH1F("hmgo_5","hmgo_5",100,0,2000);
  TH1F *hmgo_6 = new TH1F("hmgo_6","hmgo_6",100,0,2000);
  TH1F *hmgo_7 = new TH1F("hmgo_7","hmgo_7",100,0,2000);
  TH1F *hmgo_8 = new TH1F("hmgo_8","hmgo_8",20,0,20);
  TH1F *hmgo_9 = new TH1F("hmgo_9","hmgo_9",100,0,2000);
  TH1F *hmgo_10 = new TH1F("hmgo_10","hmgo_10",5, -3, 2);
  
  // multiple gamma II
  int mgt_flag_single_shower;
  double mgt_max_energy;
  double mgt_energy;
  double mgt_total_other_energy;
  double mgt_max_energy_1;
  double mgt_e_indirect_max_energy;
  double mgt_e_direct_max_energy;
  int mgt_n_direct_showers;
  double mgt_e_direct_total_energy;
  int mgt_flag_indirect_max_pio;
  double mgt_e_indirect_total_energy;
  int mgt_flag;

  T->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
  T->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
  T->SetBranchAddress("mgt_energy",&mgt_energy);
  T->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
  T->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
  T->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
  T->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
  T->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
  T->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
  T->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
  T->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
  T->SetBranchAddress("mgt_flag",&mgt_flag);

  TH1F *hmgt_1 = new TH1F("hmgt_1","hmgt_1",5,-3,2);
  TH1F *hmgt_2 = new TH1F("hmgt_2","hmgt_2",100,0,2000);
  TH1F *hmgt_3 = new TH1F("hmgt_3","hmgt_3",100,0,2000);
  TH1F *hmgt_4 = new TH1F("hmgt_4","hmgt_4",100,0,2000);
  TH1F *hmgt_5 = new TH1F("hmgt_5","hmgt_5",100,0,2000);
  TH1F *hmgt_6 = new TH1F("hmgt_6","hmgt_6",100,0,2000);
  TH1F *hmgt_7 = new TH1F("hmgt_7","hmgt_7",100,0,2000);
  TH1F *hmgt_8 = new TH1F("hmgt_8","hmgt_8",20,0,20);
  TH1F *hmgt_9 = new TH1F("hmgt_9","hmgt_9",100,0,2000);
  TH1F *hmgt_10 = new TH1F("hmgt_10","hmgt_10",5,-3,2);
  TH1F *hmgt_11 = new TH1F("hmgt_11","hmgt_11",100,0,2000);
  TH1F *hmgt_12 = new TH1F("hmgt_12","hmgt_12",5,-3,2);

  
  //single shower
  std::vector<double> *sig_1_v_angle = new std::vector<double>;
  std::vector<int> *sig_1_v_flag_single_shower= new std::vector<int>;
  std::vector<double> *sig_1_v_energy= new std::vector<double>;
  std::vector<double> *sig_1_v_energy_1= new std::vector<double>;
  std::vector<int> *sig_1_v_flag= new std::vector<int>;
  
  std::vector<double> *sig_2_v_energy= new std::vector<double>;
  std::vector<double> *sig_2_v_shower_angle= new std::vector<double>;
  std::vector<int> *sig_2_v_flag_single_shower= new std::vector<int>;
  std::vector<double> *sig_2_v_medium_dQ_dx= new std::vector<double>;
  std::vector<double> *sig_2_v_start_dQ_dx= new std::vector<double>;
  std::vector<int> *sig_2_v_flag= new std::vector<int>;
  
  int sig_flag;

  T->SetBranchAddress("sig_1_v_angle",&sig_1_v_angle);
  T->SetBranchAddress("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
  T->SetBranchAddress("sig_1_v_energy",&sig_1_v_energy);
  T->SetBranchAddress("sig_1_v_energy_1",&sig_1_v_energy_1);
  T->SetBranchAddress("sig_1_v_flag",&sig_1_v_flag);
  
  T->SetBranchAddress("sig_2_v_energy",&sig_2_v_energy);
  T->SetBranchAddress("sig_2_v_shower_angle",&sig_2_v_shower_angle);
  T->SetBranchAddress("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
  T->SetBranchAddress("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);
  T->SetBranchAddress("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);
  T->SetBranchAddress("sig_2_v_flag",&sig_2_v_flag);
  
  T->SetBranchAddress("sig_flag",&sig_flag);

  TH1F *hsig_1 = new TH1F("hsig_1", "hsig_1",90,0,180);
  TH1F *hsig_2 = new TH1F("hsig_2", "hsig_2",5,-3,2);
  TH1F *hsig_3 = new TH1F("hsig_3", "hsig_3",90,0,2000);
  TH1F *hsig_4 = new TH1F("hsig_4", "hsig_4",90,0,2000);
  TH1F *hsig_5 = new TH1F("hsig_5", "hsig_5",5,-3,2);
  TH1F *hsig_6 = new TH1F("hsig_6", "hsig_6",90,0,2000);
  TH1F *hsig_7 = new TH1F("hsig_7", "hsig_7",90,0,180);
  TH1F *hsig_8 = new TH1F("hsig_8", "hsig_8",5,-3,2);
  TH1F *hsig_9 = new TH1F("hsig_9", "hsig_9",90,0,10);
  TH1F *hsig_10 = new TH1F("hsig_10", "hsig_10",90,0,10);
  TH1F *hsig_11 = new TH1F("hsig_11", "hsig_11",5,-3,2);
  TH1F *hsig_12 = new TH1F("hsig_12", "hsig_12",5,-3,2);

  // track overclustering ..
  std::vector<int> *tro_1_v_particle_type= new std::vector<int>;
  std::vector<int> *tro_1_v_flag_dir_weak= new std::vector<int>;
  std::vector<double> *tro_1_v_min_dis = new std::vector<double>;
  std::vector<double> *tro_1_v_sg1_length = new std::vector<double>;
  std::vector<double> *tro_1_v_shower_main_length = new std::vector<double>;
  std::vector<int> *tro_1_v_max_n_vtx_segs= new std::vector<int>;
  std::vector<double> *tro_1_v_tmp_length = new std::vector<double>;
  std::vector<double> *tro_1_v_medium_dQ_dx = new std::vector<double>;
  std::vector<double> *tro_1_v_dQ_dx_cut = new std::vector<double>;
  std::vector<int> *tro_1_v_flag_shower_topology= new std::vector<int>;
  std::vector<int> *tro_1_v_flag= new std::vector<int>;
  
  std::vector<double> *tro_2_v_energy = new std::vector<double>;
  std::vector<double> *tro_2_v_stem_length = new std::vector<double>;
  std::vector<double> *tro_2_v_iso_angle = new std::vector<double>;
  std::vector<double> *tro_2_v_max_length = new std::vector<double>;
  std::vector<double> *tro_2_v_angle = new std::vector<double>;
  std::vector<int> *tro_2_v_flag= new std::vector<int>;
  
  double tro_3_stem_length;
  int tro_3_n_muon_segs;
  double tro_3_energy;
  int tro_3_flag;
  
  std::vector<double> *tro_4_v_dir2_mag = new std::vector<double>;
  std::vector<double> *tro_4_v_angle = new std::vector<double>;
  std::vector<double> *tro_4_v_angle1 = new std::vector<double>;
  std::vector<double> *tro_4_v_angle2 = new std::vector<double>;
  std::vector<double> *tro_4_v_length = new std::vector<double>;
  std::vector<double> *tro_4_v_length1 = new std::vector<double>;
  std::vector<double> *tro_4_v_medium_dQ_dx = new std::vector<double>;
  std::vector<double> *tro_4_v_end_dQ_dx = new std::vector<double>;
  std::vector<double> *tro_4_v_energy = new std::vector<double>;
  std::vector<double> *tro_4_v_shower_main_length = new std::vector<double>;
  std::vector<int> *tro_4_v_flag_shower_trajectory= new std::vector<int>;
  std::vector<int> *tro_4_v_flag= new std::vector<int>;
  
  std::vector<double> *tro_5_v_max_angle = new std::vector<double>;
  std::vector<double> *tro_5_v_min_angle = new std::vector<double>;
  std::vector<double> *tro_5_v_max_length = new std::vector<double>;
  std::vector<double> *tro_5_v_iso_angle = new std::vector<double>;
  std::vector<int> *tro_5_v_n_vtx_segs= new std::vector<int>;
  std::vector<int> *tro_5_v_min_count= new std::vector<int>;
  std::vector<int> *tro_5_v_max_count= new std::vector<int>;
  std::vector<double> *tro_5_v_energy = new std::vector<double>;
  std::vector<int> *tro_5_v_flag = new std::vector<int>;
  
  int tro_flag;
  
  T->SetBranchAddress("tro_1_v_particle_type",&tro_1_v_particle_type);
  T->SetBranchAddress("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
  T->SetBranchAddress("tro_1_v_min_dis",&tro_1_v_min_dis);
  T->SetBranchAddress("tro_1_v_sg1_length",&tro_1_v_sg1_length);
  T->SetBranchAddress("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
  T->SetBranchAddress("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
  T->SetBranchAddress("tro_1_v_tmp_length",&tro_1_v_tmp_length);
  T->SetBranchAddress("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
  T->SetBranchAddress("tro_1_v_dQ_dx_cut",&tro_1_v_dQ_dx_cut);
  T->SetBranchAddress("tro_1_v_flag_shower_topology",&tro_1_v_flag_shower_topology);
  T->SetBranchAddress("tro_1_v_flag",&tro_1_v_flag);
  
  T->SetBranchAddress("tro_2_v_energy",&tro_2_v_energy);
  T->SetBranchAddress("tro_2_v_stem_length",&tro_2_v_stem_length);
  T->SetBranchAddress("tro_2_v_iso_angle",&tro_2_v_iso_angle);
  T->SetBranchAddress("tro_2_v_max_length",&tro_2_v_max_length);
  T->SetBranchAddress("tro_2_v_angle",&tro_2_v_angle);
  T->SetBranchAddress("tro_2_v_flag",&tro_2_v_flag);
  
  T->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
  T->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
  T->SetBranchAddress("tro_3_energy",&tro_3_energy);
  T->SetBranchAddress("tro_3_flag",&tro_3_flag);
  
  T->SetBranchAddress("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
  T->SetBranchAddress("tro_4_v_angle",&tro_4_v_angle);
  T->SetBranchAddress("tro_4_v_angle1",&tro_4_v_angle1);
  T->SetBranchAddress("tro_4_v_angle2",&tro_4_v_angle2);
  T->SetBranchAddress("tro_4_v_length",&tro_4_v_length);
  T->SetBranchAddress("tro_4_v_length1",&tro_4_v_length1);
  T->SetBranchAddress("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
  T->SetBranchAddress("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
  T->SetBranchAddress("tro_4_v_energy",&tro_4_v_energy);
  T->SetBranchAddress("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
  T->SetBranchAddress("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);
  T->SetBranchAddress("tro_4_v_flag",&tro_4_v_flag);
  
  T->SetBranchAddress("tro_5_v_max_angle",&tro_5_v_max_angle);
  T->SetBranchAddress("tro_5_v_min_angle",&tro_5_v_min_angle);
  T->SetBranchAddress("tro_5_v_max_length",&tro_5_v_max_length);
  T->SetBranchAddress("tro_5_v_iso_angle",&tro_5_v_iso_angle);
  T->SetBranchAddress("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
  T->SetBranchAddress("tro_5_v_min_count",&tro_5_v_min_count);
  T->SetBranchAddress("tro_5_v_max_count",&tro_5_v_max_count);
  T->SetBranchAddress("tro_5_v_energy",&tro_5_v_energy);
  T->SetBranchAddress("tro_5_v_flag",&tro_5_v_flag);
  
  T->SetBranchAddress("tro_flag",&tro_flag);
  
  TH1F *htro_1 = new TH1F("htro_1","htro_1",20,0,100);
  TH1F *htro_2 = new TH1F("htro_2","htro_2",5,-3,2);
  TH1F *htro_3 = new TH1F("htro_3","htro_3",20,0,1000);
  TH1F *htro_4 = new TH1F("htro_4","htro_4",20,0,100);
  TH1F *htro_5 = new TH1F("htro_5","htro_5",20,0,100);
  TH1F *htro_6 = new TH1F("htro_6","htro_6",20,0,10);
  TH1F *htro_7 = new TH1F("htro_7","htro_7",20,0,100);
  TH1F *htro_8 = new TH1F("htro_8","htro_8",20,0,10);
  TH1F *htro_9 = new TH1F("htro_9","htro_9",20,0,10);
  TH1F *htro_10 = new TH1F("htro_10","htro_10",5,-3,2);
  TH1F *htro_11 = new TH1F("htro_11","htro_11",5,-3,2);

  TH1F *htro_12 = new TH1F("htro_12","htro_12",100,0,2000);
  TH1F *htro_13 = new TH1F("htro_13","htro_13",20,0,100);
  TH1F *htro_14 = new TH1F("htro_14","htro_14",45,0,90);
  TH1F *htro_15 = new TH1F("htro_15","htro_15",20,0,100);
  TH1F *htro_16 = new TH1F("htro_16","htro_16",45,0,180);
  TH1F *htro_17 = new TH1F("htro_17","htro_17",5, -3, 2);

  TH1F *htro_18 = new TH1F("htro_18","htro_18",20,0,100);
  TH1F *htro_19 = new TH1F("htro_19","htro_19",20,0,100);
  TH1F *htro_20 = new TH1F("htro_20","htro_20",100,0,2000);
  TH1F *htro_21 = new TH1F("htro_21","htro_21",5, -3, 2);

  TH1F *htro_22 = new TH1F("htro_22","htro_22",20,0,1000);
  TH1F *htro_23 = new TH1F("htro_23","htro_23",100,0,180);
  TH1F *htro_24 = new TH1F("htro_24","htro_24",100,0,180);
  TH1F *htro_25 = new TH1F("htro_25","htro_25",100,0,180);
  TH1F *htro_26 = new TH1F("htro_26","htro_26",20,0,100);
  TH1F *htro_27 = new TH1F("htro_27","htro_27",20,0,100);
  TH1F *htro_28 = new TH1F("htro_28","htro_28",100,0,10);
  TH1F *htro_29 = new TH1F("htro_29","htro_29",100,0,10);
  TH1F *htro_30 = new TH1F("htro_30","htro_30",100,0,2000);
  TH1F *htro_31 = new TH1F("htro_31","htro_31",100,0,100);
  TH1F *htro_32 = new TH1F("htro_32","htro_32",5, -3, 2);
  TH1F *htro_33 = new TH1F("htro_33","htro_33",5, -3, 2);

  TH1F *htro_34 = new TH1F("htro_34","htro_34",100,0,180);
  TH1F *htro_35 = new TH1F("htro_35","htro_35",100,0,180);
  TH1F *htro_36 = new TH1F("htro_36","htro_36",100,0,100);
  TH1F *htro_37 = new TH1F("htro_37","htro_37",100,0,180);
  TH1F *htro_38 = new TH1F("htro_38","htro_38",20,0,20);
  TH1F *htro_39 = new TH1F("htro_39","htro_39",20,0,20);
  TH1F *htro_40 = new TH1F("htro_40","htro_40",20,0,20);
  TH1F *htro_41 = new TH1F("htro_41","htro_41",100,0,2000);
  TH1F *htro_42 = new TH1F("htro_42","htro_42",5,-3,2);

  TH1F *htro_43 = new TH1F("htro_43","htro_43",5, -3, 2);
  
  

  
  for (Int_t j=0;j!=T->GetEntries(); j++){
    T->GetEntry(j);
    hc1->Fill(cosmic_flag);
    hc2->Fill(cosmic_energy_direct_showers);
    hc3->Fill(cosmic_energy_main_showers);
    hc4->Fill(cosmic_energy_indirect_showers);
    hc5->Fill(cosmic_n_direct_showers);
    hc6->Fill(cosmic_n_indirect_showers);
    hc7->Fill(cosmic_n_main_showers);
    hc8->Fill(cosmic_n_solid_tracks);
    hc9->Fill(cosmic_filled);

    hg1->Fill(gap_flag);
    hg2->Fill(gap_flag_prolong_u);
    hg3->Fill(gap_flag_prolong_v);
    hg4->Fill(gap_flag_prolong_w);
    hg5->Fill(gap_flag_parallel);
    hg6->Fill(gap_n_points);
    hg7->Fill(gap_n_bad);
    hg8->Fill(gap_energy);
    hg9->Fill(gap_num_valid_tracks);
    hg10->Fill(gap_flag_single_shower);
    hg11->Fill(gap_filled);

    hmq1->Fill(mip_quality_flag);
    hmq2->Fill(mip_quality_energy);
    hmq3->Fill(mip_quality_overlap);
    hmq4->Fill(mip_quality_n_showers);
    hmq5->Fill(mip_quality_n_tracks);
    hmq6->Fill(mip_quality_flag_inside_pi0);
    hmq7->Fill(mip_quality_n_pi0_showers);
    hmq8->Fill(mip_quality_shortest_length);
    hmq9->Fill(mip_quality_acc_length);
    hmq10->Fill(mip_quality_shortest_angle);
    hmq11->Fill(mip_quality_flag_proton);
    hmq12->Fill(mip_quality_filled);

    hmip1->Fill(mip_flag);
    hmip2->Fill(mip_energy);
    hmip3->Fill(mip_n_first_mip);
    hmip4->Fill(mip_n_first_non_mip);
    hmip5->Fill(mip_n_first_non_mip_1);
    hmip6->Fill(mip_n_first_non_mip_2);
    hmip7->Fill(mip_n_end_reduction);

    hmip8->Fill(mip_vec_dQ_dx_0);
    hmip9->Fill(mip_vec_dQ_dx_1);
    hmip10->Fill(mip_vec_dQ_dx_2);
    hmip11->Fill(mip_vec_dQ_dx_3);
    hmip12->Fill(mip_vec_dQ_dx_4);
    hmip13->Fill(mip_vec_dQ_dx_5);
    hmip14->Fill(mip_vec_dQ_dx_6);
    hmip15->Fill(mip_vec_dQ_dx_7);
    hmip16->Fill(mip_vec_dQ_dx_8);
    hmip17->Fill(mip_vec_dQ_dx_9);
    hmip18->Fill(mip_vec_dQ_dx_10);
    hmip19->Fill(mip_vec_dQ_dx_11);
    hmip20->Fill(mip_vec_dQ_dx_12);
    hmip21->Fill(mip_vec_dQ_dx_13);
    hmip22->Fill(mip_vec_dQ_dx_14);
    hmip23->Fill(mip_vec_dQ_dx_15);
    hmip24->Fill(mip_vec_dQ_dx_16);
    hmip25->Fill(mip_vec_dQ_dx_17);
    hmip26->Fill(mip_vec_dQ_dx_18);
    hmip27->Fill(mip_vec_dQ_dx_19);

    hmip28->Fill(mip_max_dQ_dx_sample);
    hmip29->Fill(mip_n_below_threshold);
    hmip30->Fill(mip_n_below_zero);
    hmip31->Fill(mip_n_lowest);
    hmip32->Fill(mip_n_highest);

    hmip33->Fill(mip_lowest_dQ_dx);
    hmip34->Fill(mip_highest_dQ_dx);
    hmip35->Fill(mip_medium_dQ_dx);
    hmip36->Fill(mip_stem_length);
    hmip37->Fill(mip_length_main);
    hmip38->Fill(mip_length_total);
    hmip39->Fill(mip_angle_beam);
    hmip40->Fill(mip_iso_angle);

    hmip41->Fill(mip_n_vertex);
    hmip42->Fill(mip_n_good_tracks);
    hmip43->Fill(mip_E_indirect_max_energy);
    hmip44->Fill(mip_flag_all_above);
    hmip45->Fill(mip_min_dQ_dx_5);
    hmip46->Fill(mip_n_other_vertex);
    hmip47->Fill(mip_n_stem_size);
    hmip48->Fill(mip_flag_stem_trajectory);
    hmip49->Fill(mip_min_dis);
    hmip50->Fill(mip_filled);

    hpio1->Fill(pio_flag);
    hpio2->Fill(pio_mip_id);
    hpio3->Fill(pio_filled);
    hpio4->Fill(pio_flag_pio);
    
    hpio5->Fill(pio_1_mass);
    hpio6->Fill(pio_1_pio_type);
    hpio7->Fill(pio_1_energy_1);
    hpio8->Fill(pio_1_energy_2);
    hpio9->Fill(pio_1_dis_1);
    hpio10->Fill(pio_1_dis_2);
    hpio14->Fill(pio_1_flag);
		 
    for (size_t i=0;i!=pio_2_v_dis2->size();i++){
      hpio11->Fill(pio_2_v_dis2->at(i));
      hpio12->Fill(pio_2_v_angle2->at(i));
      hpio13->Fill(pio_2_v_acc_length->at(i));
      hpio15->Fill(pio_2_v_flag->at(i));
    }

    hsd1->Fill(stem_dir_flag);
    hsd2->Fill(stem_dir_flag_single_shower);
    hsd3->Fill(stem_dir_filled);
    hsd4->Fill(stem_dir_angle);
    hsd5->Fill(stem_dir_energy);
    hsd6->Fill(stem_dir_angle1);
    hsd7->Fill(stem_dir_angle2);
    hsd8->Fill(stem_dir_angle3);
    hsd9->Fill(stem_dir_ratio);
    
    hbr1_1->Fill(br_filled);
    hbr1_2->Fill(br1_1_flag);
    hbr1_3->Fill(br1_1_shower_type);
    hbr1_4->Fill(br1_1_vtx_n_segs);
    hbr1_5->Fill(br1_1_energy);
    hbr1_6->Fill(br1_1_n_segs);
    hbr1_7->Fill(br1_1_flag_sg_topology);
    hbr1_8->Fill(br1_1_flag_sg_trajectory);
    hbr1_9->Fill(br1_1_sg_length);
    hbr1_10->Fill(br1_2_flag);
    hbr1_11->Fill(br1_2_energy);
    hbr1_12->Fill(br1_2_n_connected);
    hbr1_13->Fill(br1_2_max_length);
    hbr1_14->Fill(br1_2_n_connected_1);
    hbr1_15->Fill(br1_2_vtx_n_segs);
    hbr1_16->Fill(br1_2_n_shower_segs);
    hbr1_17->Fill(br1_2_max_length_ratio);
    hbr1_18->Fill(br1_2_shower_length);
    hbr1_19->Fill(br1_3_flag);
    hbr1_20->Fill(br1_3_energy);
    hbr1_21->Fill(br1_3_n_connected_p);
    hbr1_22->Fill(br1_3_max_length_p);
    hbr1_23->Fill(br1_3_n_shower_segs);
    hbr1_24->Fill(br1_3_flag_sg_topology);
    hbr1_25->Fill(br1_3_flag_sg_trajectory);
    hbr1_26->Fill(br1_3_n_shower_main_segs);
    hbr1_27->Fill(br1_3_sg_length);

    hbr2_1->Fill(br2_flag);
    hbr2_2->Fill(br2_flag_single_shower);
    hbr2_3->Fill(br2_num_valid_tracks);
    hbr2_4->Fill(br2_energy);
    hbr2_5->Fill(br2_angle1);
    hbr2_6->Fill(br2_angle2);
    hbr2_7->Fill(br2_angle);
    hbr2_8->Fill(br2_angle3);
    hbr2_9->Fill(br2_n_shower_main_segs);
    hbr2_10->Fill(br2_max_angle);
    hbr2_11->Fill(br2_sg_length);
    hbr2_12->Fill(br2_flag_sg_trajectory);

    hlol_1->Fill(lol_flag);
    for (size_t i=0; i!= lol_1_v_energy->size(); i++){
      hlol_2->Fill(lol_1_v_energy->at(i));
      hlol_3->Fill(lol_1_v_vtx_n_segs->at(i));
      hlol_4->Fill(lol_1_v_nseg->at(i));
      hlol_5->Fill(lol_1_v_angle->at(i));
      hlol_6->Fill(lol_1_v_flag->at(i));
    }
    
    for (size_t i=0; i!= lol_2_v_length->size();i++){
      hlol_7->Fill(lol_2_v_length->at(i));
      hlol_8->Fill(lol_2_v_angle->at(i));
      hlol_9->Fill(lol_2_v_type->at(i));
      hlol_10->Fill(lol_2_v_vtx_n_segs->at(i));
      hlol_11->Fill(lol_2_v_energy->at(i));
      hlol_12->Fill(lol_2_v_shower_main_length->at(i));
      hlol_13->Fill(lol_2_v_flag_dir_weak->at(i));
      hlol_14->Fill(lol_2_v_flag->at(i));
    }
    hlol_15->Fill(lol_3_angle_beam);
    hlol_16->Fill(lol_3_n_valid_tracks);
    hlol_17->Fill(lol_3_min_angle);
    hlol_18->Fill(lol_3_vtx_n_segs);
    hlol_19->Fill(lol_3_energy);
    hlol_20->Fill(lol_3_shower_main_length);
    hlol_21->Fill(lol_3_n_out);
    hlol_22->Fill(lol_3_n_sum);
    hlol_23->Fill(lol_3_flag);
    
    

    hbr3_1->Fill(br3_1_energy);
    hbr3_2->Fill(br3_1_n_shower_segments);
    hbr3_3->Fill(br3_1_sg_flag_trajectory);
    hbr3_4->Fill(br3_1_sg_direct_length);
    hbr3_5->Fill(br3_1_sg_length);
    hbr3_6->Fill(br3_1_total_main_length);
    hbr3_7->Fill(br3_1_total_length);
    hbr3_8->Fill(br3_1_iso_angle);
    hbr3_9->Fill(br3_1_sg_flag_topology);
    hbr3_10->Fill(br3_1_flag);

    hbr3_11->Fill(br3_2_n_ele);
    hbr3_12->Fill(br3_2_n_other);
    hbr3_13->Fill(br3_2_energy);
    hbr3_14->Fill(br3_2_total_main_length);
    hbr3_15->Fill(br3_2_total_length);
    hbr3_16->Fill(br3_2_other_fid);
    hbr3_17->Fill(br3_2_flag);

    for (size_t i=0;i!=br3_3_v_energy->size();i++){
      hbr3_18->Fill(br3_3_v_energy->at(i));
      hbr3_19->Fill(br3_3_v_angle->at(i));
      hbr3_20->Fill(br3_3_v_dir_length->at(i));
      hbr3_21->Fill(br3_3_v_length->at(i));
      hbr3_22->Fill(br3_3_v_flag->at(i));
    }

    hbr3_23->Fill(br3_4_acc_length);
    hbr3_24->Fill(br3_4_total_length);
    hbr3_25->Fill(br3_4_energy);
    hbr3_26->Fill(br3_4_flag);

    for (size_t i=0;i!=br3_5_v_dir_length->size();i++){
      hbr3_27->Fill(br3_5_v_dir_length->at(i));
      hbr3_28->Fill(br3_5_v_total_length->at(i));
      hbr3_29->Fill(br3_5_v_flag_avoid_muon_check->at(i));
      hbr3_30->Fill(br3_5_v_n_seg->at(i));
      hbr3_31->Fill(br3_5_v_angle->at(i));
      hbr3_32->Fill(br3_5_v_sg_length->at(i));
      hbr3_33->Fill(br3_5_v_energy->at(i));
      hbr3_34->Fill(br3_5_v_n_main_segs->at(i));
      hbr3_35->Fill(br3_5_v_n_segs->at(i));
      hbr3_36->Fill(br3_5_v_shower_main_length->at(i));
      hbr3_37->Fill(br3_5_v_shower_total_length->at(i));
      hbr3_38->Fill(br3_5_v_flag->at(i));
    }

    for (size_t i=0;i!=br3_6_v_angle->size();i++){
      hbr3_39->Fill(br3_6_v_angle->at(i));
      hbr3_40->Fill(br3_6_v_angle1->at(i));
      hbr3_41->Fill(br3_6_v_flag_shower_trajectory->at(i));
      hbr3_42->Fill(br3_6_v_direct_length->at(i));
      hbr3_43->Fill(br3_6_v_length->at(i));
      hbr3_44->Fill(br3_6_v_n_other_vtx_segs->at(i));
      hbr3_45->Fill(br3_6_v_energy->at(i));
      hbr3_46->Fill(br3_6_v_flag->at(i));
    }

    hbr3_47->Fill(br3_7_energy);
    hbr3_48->Fill(br3_7_min_angle);
    hbr3_49->Fill(br3_7_sg_length);
    hbr3_50->Fill(br3_7_shower_main_length);
    hbr3_51->Fill(br3_7_flag);

    hbr3_52->Fill(br3_8_max_dQ_dx);
    hbr3_53->Fill(br3_8_energy);
    hbr3_54->Fill(br3_8_n_main_segs);
    hbr3_55->Fill(br3_8_shower_main_length);
    hbr3_56->Fill(br3_8_shower_length);
    hbr3_57->Fill(br3_8_flag);

    hbr3_58->Fill(br3_flag);


    hbr4_1->Fill(br4_1_shower_main_length);
    hbr4_2->Fill(br4_1_shower_total_length);
    hbr4_3->Fill(br4_1_min_dis);
    hbr4_4->Fill(br4_1_energy);
    hbr4_5->Fill(br4_1_flag_avoid_muon_check);
    hbr4_6->Fill(br4_1_n_vtx_segs);
    hbr4_7->Fill(br4_1_n_main_segs);
    hbr4_8->Fill(br4_1_flag);

    hbr4_9->Fill(br4_2_ratio_45);
    hbr4_10->Fill(br4_2_ratio_35);
    hbr4_11->Fill(br4_2_ratio_25);
    hbr4_12->Fill(br4_2_ratio_15);
    hbr4_13->Fill(br4_2_energy);
    hbr4_14->Fill(br4_2_ratio1_45);
    hbr4_15->Fill(br4_2_ratio1_35);
    hbr4_16->Fill(br4_2_ratio1_25);
    hbr4_17->Fill(br4_2_ratio1_15);
    hbr4_18->Fill(br4_2_iso_angle);
    hbr4_19->Fill(br4_2_iso_angle1);
    hbr4_20->Fill(br4_2_angle);
    hbr4_21->Fill(br4_2_flag);

    hbr4_22->Fill(br4_flag);

    hol_1->Fill(hol_1_n_valid_tracks);
    hol_2->Fill(hol_1_min_angle);
    hol_3->Fill(hol_1_energy);
    hol_4->Fill(hol_1_flag_all_shower);
    hol_5->Fill(hol_1_min_length);
    hol_6->Fill(hol_1_flag);

    hol_7->Fill(hol_2_min_angle);
    hol_8->Fill(hol_2_medium_dQ_dx);
    hol_9->Fill(hol_2_ncount);
    hol_10->Fill(hol_2_energy);
    hol_11->Fill(hol_2_flag);

    hol_12->Fill(hol_flag);

    hvis_1->Fill(vis_1_filled);
    hvis_2->Fill(vis_1_n_vtx_segs);
    hvis_3->Fill(vis_1_energy);
    hvis_4->Fill(vis_1_num_good_tracks);
    hvis_5->Fill(vis_1_max_angle);
    hvis_6->Fill(vis_1_tmp_length1);
    hvis_7->Fill(vis_1_tmp_length2);
    hvis_8->Fill(vis_1_particle_type);
    hvis_9->Fill(vis_1_flag);

    hvis_10->Fill(vis_2_filled);
    hvis_11->Fill(vis_2_n_vtx_segs);
    hvis_12->Fill(vis_2_min_angle);
    hvis_13->Fill(vis_2_min_weak_track);
    hvis_14->Fill(vis_2_angle_beam);
    hvis_15->Fill(vis_2_min_angle1);
    hvis_16->Fill(vis_2_iso_angle1);
    hvis_17->Fill(vis_2_min_medium_dQ_dx);
    hvis_18->Fill(vis_2_min_length);
    hvis_19->Fill(vis_2_sg_length);
    hvis_20->Fill(vis_2_max_angle);
    hvis_21->Fill(vis_2_max_weak_track);
    hvis_22->Fill(vis_2_flag);

    hvis_23->Fill(vis_flag);


    hstl_1->Fill(stem_len_energy);
    hstl_2->Fill(stem_len_length);
    hstl_3->Fill(stem_len_flag_avoid_muon_check);
    hstl_4->Fill(stem_len_num_daughters);
    hstl_5->Fill(stem_len_daughter_length);
    hstl_6->Fill(stem_len_flag);


    hbrm_1->Fill(brm_n_mu_segs);
    hbrm_2->Fill(brm_Ep);
    hbrm_3->Fill(brm_energy);
    hbrm_4->Fill(brm_acc_length);
    hbrm_5->Fill(brm_shower_total_length);
    hbrm_6->Fill(brm_connected_length);
    hbrm_7->Fill(brm_n_size);
    hbrm_8->Fill(brm_acc_direct_length);
    hbrm_9->Fill(brm_n_shower_main_segs);
    hbrm_10->Fill(brm_n_mu_main);
    hbrm_11->Fill(brm_flag);

    hcme_1->Fill(cme_mu_energy);
    hcme_2->Fill(cme_energy);
    hcme_3->Fill(cme_mu_length);
    hcme_4->Fill(cme_length);
    hcme_5->Fill(cme_angle_beam);
    hcme_6->Fill(cme_flag);


    hanc_1->Fill(anc_energy);
    hanc_2->Fill(anc_angle);
    hanc_3->Fill(anc_max_angle);
    hanc_4->Fill(anc_max_length);
    hanc_5->Fill(anc_acc_forward_length);
    hanc_6->Fill(anc_acc_backward_length);
    hanc_7->Fill(anc_acc_forward_length1);
    hanc_8->Fill(anc_shower_main_length);
    hanc_9->Fill(anc_shower_total_length);
    hanc_10->Fill(anc_flag_main_outside);
    hanc_11->Fill(anc_flag);

    hlem_1->Fill(lem_shower_total_length);
    hlem_2->Fill(lem_shower_main_length);
    hlem_3->Fill(lem_n_3seg);
    hlem_4->Fill(lem_e_charge);
    hlem_5->Fill(lem_e_dQdx);
    hlem_6->Fill(lem_shower_num_segs);
    hlem_7->Fill(lem_shower_num_main_segs);
    hlem_8->Fill(lem_flag);


    hstw_1->Fill(stw_1_energy);
    hstw_2->Fill(stw_1_dis);
    hstw_3->Fill(stw_1_dQ_dx);
    hstw_4->Fill(stw_1_flag_single_shower);
    hstw_5->Fill(stw_1_n_pi0);
    hstw_6->Fill(stw_1_num_valid_tracks);
    hstw_7->Fill(stw_1_flag);

    for (size_t i=0;i!=stw_2_v_medium_dQ_dx->size();i++){
      hstw_8->Fill(stw_2_v_medium_dQ_dx->at(i));
      hstw_9->Fill(stw_2_v_energy->at(i));
      hstw_10->Fill(stw_2_v_angle->at(i));
      hstw_11->Fill(stw_2_v_dir_length->at(i));
      hstw_12->Fill(stw_2_v_max_dQ_dx->at(i));
      hstw_13->Fill(stw_2_v_flag->at(i));
    }

    for (size_t i=0;i!=stw_3_v_angle->size(); i++){
      hstw_14->Fill(stw_3_v_angle->at(i));
      hstw_15->Fill(stw_3_v_dir_length->at(i));
      hstw_16->Fill(stw_3_v_energy->at(i));
      hstw_17->Fill(stw_3_v_medium_dQ_dx->at(i));
      hstw_18->Fill(stw_3_v_flag->at(i));
    }

    for (size_t i=0;i!= stw_4_v_angle->size();i++){
      hstw_19->Fill(stw_4_v_angle->at(i));
      hstw_20->Fill(stw_4_v_dis->at(i));
      hstw_21->Fill(stw_4_v_energy->at(i));
      hstw_22->Fill(stw_4_v_flag->at(i));
    }

    hstw_23->Fill(stw_flag);
    
    hspt_1->Fill(spt_flag_single_shower);
    hspt_2->Fill(spt_energy);
    hspt_3->Fill(spt_shower_main_length);
    hspt_4->Fill(spt_shower_total_length);
    hspt_5->Fill(spt_angle_beam);
    hspt_6->Fill(spt_angle_vertical);
    hspt_7->Fill(spt_max_dQ_dx);
    hspt_8->Fill(spt_angle_beam_1);
    hspt_9->Fill(spt_angle_drift);
    hspt_10->Fill(spt_angle_drift_1);
    hspt_11->Fill(spt_num_valid_tracks);
    hspt_12->Fill(spt_n_vtx_segs);
    hspt_13->Fill(spt_max_length);
    hspt_14->Fill(spt_flag);

    hmgo_1->Fill(mgo_energy);
    hmgo_2->Fill(mgo_max_energy);
    hmgo_3->Fill(mgo_total_energy);
    hmgo_4->Fill(mgo_n_showers);
    hmgo_5->Fill(mgo_max_energy_1);
    hmgo_6->Fill(mgo_max_energy_2);
    hmgo_7->Fill(mgo_total_other_energy);
    hmgo_8->Fill(mgo_n_total_showers);
    hmgo_9->Fill(mgo_total_other_energy_1);
    hmgo_10->Fill(mgo_flag);

    hmgt_1->Fill(mgt_flag_single_shower);
    hmgt_2->Fill(mgt_max_energy);
    hmgt_3->Fill(mgt_energy);
    hmgt_4->Fill(mgt_total_other_energy);
    hmgt_5->Fill(mgt_max_energy_1);
    hmgt_6->Fill(mgt_e_indirect_max_energy);
    hmgt_7->Fill(mgt_e_direct_max_energy);
    hmgt_8->Fill(mgt_n_direct_showers);
    hmgt_9->Fill(mgt_e_direct_total_energy);
    hmgt_10->Fill(mgt_flag_indirect_max_pio);
    hmgt_11->Fill(mgt_e_indirect_total_energy);
    hmgt_12->Fill(mgt_flag);

    for (size_t i=0;i!=sig_1_v_angle->size();i++){
      hsig_1->Fill(sig_1_v_angle->at(i));
      hsig_2->Fill(sig_1_v_flag_single_shower->at(i));
      hsig_3->Fill(sig_1_v_energy->at(i));
      hsig_4->Fill(sig_1_v_energy_1->at(i));
      hsig_5->Fill(sig_1_v_flag->at(i));
    }

    for (size_t i=0;i!=sig_2_v_energy->size();i++){
      hsig_6->Fill(sig_2_v_energy->at(i));
      hsig_7->Fill(sig_2_v_shower_angle->at(i));
      hsig_8->Fill(sig_2_v_flag_single_shower->at(i));
      hsig_9->Fill(sig_2_v_medium_dQ_dx->at(i));
      hsig_10->Fill(sig_2_v_start_dQ_dx->at(i));
      hsig_11->Fill(sig_2_v_flag->at(i));
    }

    hsig_12->Fill(sig_flag);

    for (size_t i=0;i!=tro_1_v_particle_type->size();i++){
      htro_1->Fill(tro_1_v_particle_type->at(i));
      htro_2->Fill(tro_1_v_flag_dir_weak->at(i));
      htro_3->Fill(tro_1_v_min_dis->at(i));
      htro_4->Fill(tro_1_v_sg1_length->at(i));
      htro_5->Fill(tro_1_v_shower_main_length->at(i));
      htro_6->Fill(tro_1_v_max_n_vtx_segs->at(i));
      htro_7->Fill(tro_1_v_tmp_length->at(i));
      htro_8->Fill(tro_1_v_medium_dQ_dx->at(i));
      htro_9->Fill(tro_1_v_dQ_dx_cut->at(i));
      htro_10->Fill(tro_1_v_flag_shower_topology->at(i));
      htro_11->Fill(tro_1_v_flag->at(i));
    }

    for (size_t i=0;i!=tro_2_v_energy->size();i++){
      htro_12->Fill(tro_2_v_energy->at(i));
      htro_13->Fill(tro_2_v_stem_length->at(i));
      htro_14->Fill(tro_2_v_iso_angle->at(i));
      htro_15->Fill(tro_2_v_max_length->at(i));
      htro_16->Fill(tro_2_v_angle->at(i));
      htro_17->Fill(tro_2_v_flag->at(i));
    }

    htro_18->Fill(tro_3_stem_length);
    htro_19->Fill(tro_3_n_muon_segs);
    htro_20->Fill(tro_3_energy);
    htro_21->Fill(tro_3_flag);

    for (size_t i=0;i!=tro_4_v_dir2_mag->size();i++){
      htro_22->Fill(tro_4_v_dir2_mag->at(i));
      htro_23->Fill(tro_4_v_angle->at(i));
      htro_24->Fill(tro_4_v_angle1->at(i));
      htro_25->Fill(tro_4_v_angle2->at(i));
      htro_26->Fill(tro_4_v_length->at(i));
      htro_27->Fill(tro_4_v_length1->at(i));
      htro_28->Fill(tro_4_v_medium_dQ_dx->at(i));
      htro_29->Fill(tro_4_v_end_dQ_dx->at(i));
      htro_30->Fill(tro_4_v_energy->at(i));
      htro_31->Fill(tro_4_v_shower_main_length->at(i));
      htro_32->Fill(tro_4_v_flag_shower_trajectory->at(i));
      htro_33->Fill(tro_4_v_flag->at(i));
    }

    for (size_t i=0; i!= tro_5_v_max_angle->size(); i++){
      htro_34->Fill(tro_5_v_max_angle->at(i));
      htro_35->Fill(tro_5_v_min_angle->at(i));
      htro_36->Fill(tro_5_v_max_length->at(i));
      htro_37->Fill(tro_5_v_iso_angle->at(i));
      htro_38->Fill(tro_5_v_n_vtx_segs->at(i));
      htro_39->Fill(tro_5_v_min_count->at(i));
      htro_40->Fill(tro_5_v_max_count->at(i));
      htro_41->Fill(tro_5_v_energy->at(i));
      htro_42->Fill(tro_5_v_flag->at(i));
    }

    htro_43->Fill(tro_flag);
    
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(5,5);
  c1->cd(1); hc1->Draw();
  c1->cd(2); hc2->Draw();
  c1->cd(3); hc3->Draw();
  c1->cd(4); hc4->Draw();
  c1->cd(5); hc5->Draw();
  c1->cd(6); hc6->Draw();
  c1->cd(7); hc7->Draw();
  c1->cd(8); hc8->Draw();
  c1->cd(9); hc9->Draw();

  c1->cd(10); hg1->Draw();
  c1->cd(11); hg2->Draw();
  c1->cd(12); hg3->Draw();
  c1->cd(13); hg4->Draw();
  c1->cd(14); hg5->Draw();
  c1->cd(15); hg6->Draw();
  c1->cd(16); hg7->Draw();
  c1->cd(17); hg8->Draw();
  c1->cd(18); hg9->Draw();
  c1->cd(19); hg10->Draw();
  c1->cd(20); hg11->Draw();
    
  c1->cd(21); hmq1->Draw();
  c1->cd(22); hmq2->Draw();
  c1->cd(23); hmq3->Draw();
  c1->cd(24); hmq4->Draw();
  c1->cd(25); hmq5->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  c2->Divide(5,5);
  c2->cd(1); hmq6->Draw();
  c2->cd(2); hmq7->Draw();
  c2->cd(3); hmq8->Draw();
  c2->cd(4); hmq9->Draw();
  c2->cd(5); hmq10->Draw();
  c2->cd(6); hmq11->Draw();
  c2->cd(7); hmq12->Draw();

  c2->cd(8); hmip1->Draw();
  c2->cd(9); hmip2->Draw();
  c2->cd(10); hmip3->Draw();
  c2->cd(11); hmip4->Draw();
  c2->cd(12); hmip5->Draw();
  c2->cd(13); hmip6->Draw();
  c2->cd(14); hmip7->Draw();
  c2->cd(15); hmip8->Draw();
  c2->cd(16); hmip9->Draw();
  c2->cd(17); hmip10->Draw();
  c2->cd(18); hmip11->Draw();
  c2->cd(19); hmip12->Draw();
  c2->cd(20); hmip13->Draw();
  c2->cd(21); hmip14->Draw();
  c2->cd(22); hmip15->Draw();
  c2->cd(23); hmip16->Draw();
  c2->cd(24); hmip17->Draw();
  c2->cd(25); hmip18->Draw();


  TCanvas *c3 = new TCanvas("c3","c3",1200,900);
  c3->Divide(5,5);
  c3->cd(1); hmip19->Draw();
  c3->cd(2); hmip20->Draw();
  c3->cd(3); hmip21->Draw();
  c3->cd(4); hmip22->Draw();
  c3->cd(5); hmip23->Draw();
  c3->cd(6); hmip24->Draw();
  c3->cd(7); hmip25->Draw();

  c3->cd(8); hmip26->Draw();
  c3->cd(9); hmip27->Draw();
  c3->cd(10); hmip28->Draw();
  c3->cd(11); hmip29->Draw();
  c3->cd(12); hmip30->Draw();
  c3->cd(13); hmip31->Draw();
  c3->cd(14); hmip32->Draw();
  c3->cd(15); hmip32->Draw();
  c3->cd(16); hmip33->Draw();
  c3->cd(17); hmip34->Draw();
  c3->cd(18); hmip35->Draw();
  c3->cd(19); hmip36->Draw();
  c3->cd(20); hmip37->Draw();
  c3->cd(21); hmip38->Draw();
  c3->cd(22); hmip39->Draw();
  c3->cd(23); hmip40->Draw();
  c3->cd(24); hmip41->Draw();
  c3->cd(25); hmip42->Draw();

  TCanvas *c4 = new TCanvas("c4","c4",1200,900);
  c4->Divide(5,5);
  c4->cd(1); hmip43->Draw();
  c4->cd(2); hmip44->Draw();
  c4->cd(3); hmip45->Draw();
  c4->cd(4); hmip46->Draw();
  c4->cd(5); hmip47->Draw();
  c4->cd(6); hmip48->Draw();
  c4->cd(7); hmip49->Draw();
  c4->cd(8); hmip50->Draw();

  c4->cd(9); hpio1->Draw();
  c4->cd(10); hpio2->Draw();
  c4->cd(11); hpio3->Draw();
  c4->cd(12); hpio4->Draw();
  c4->cd(13); hpio5->Draw();
  c4->cd(14); hpio6->Draw();
  c4->cd(15); hpio7->Draw();
  c4->cd(16); hpio8->Draw();
  c4->cd(17); hpio9->Draw();
  c4->cd(18); hpio10->Draw();
  c4->cd(19); hpio11->Draw();
  c4->cd(20); hpio12->Draw();
  c4->cd(21); hpio13->Draw();

  c4->cd(22); hsd1->Draw();
  c4->cd(23); hsd2->Draw();
  c4->cd(24); hsd3->Draw();
  c4->cd(25); hsd4->Draw();

  TCanvas *c5 = new TCanvas("c5","c5",1200,900);
  c5->Divide(5,5);
  
  c5->cd(1); hsd5->Draw();
  c5->cd(2); hsd6->Draw();
  c5->cd(3); hsd7->Draw();
  c5->cd(4); hsd8->Draw();
  c5->cd(5); hsd9->Draw();
  
  c5->cd(6); hbr1_1->Draw();
  c5->cd(7); hbr1_2->Draw();
  c5->cd(8); hbr1_3->Draw();
  c5->cd(9); hbr1_4->Draw();
  c5->cd(10); hbr1_5->Draw();
  c5->cd(11); hbr1_6->Draw();
  c5->cd(12); hbr1_7->Draw();
  c5->cd(13); hbr1_8->Draw();
  c5->cd(14); hbr1_9->Draw();
  c5->cd(15); hbr1_10->Draw();
  c5->cd(16); hbr1_11->Draw();
  c5->cd(17); hbr1_12->Draw();
  c5->cd(18); hbr1_13->Draw();
  c5->cd(19); hbr1_14->Draw();
  c5->cd(20); hbr1_15->Draw();
  c5->cd(21); hbr1_16->Draw();
  c5->cd(22); hbr1_17->Draw();
  c5->cd(23); hbr1_18->Draw();
  c5->cd(24); hbr1_19->Draw();
  c5->cd(25); hbr1_20->Draw();

  TCanvas *c6 = new TCanvas("c6","c6",1200,900);
  c6->Divide(5,5);
  c6->cd(1); hbr1_21->Draw();
  c6->cd(2); hbr1_22->Draw();
  c6->cd(3); hbr1_23->Draw();
  c6->cd(4); hbr1_24->Draw();
  c6->cd(5); hbr1_25->Draw();
  c6->cd(6); hbr1_26->Draw();
  c6->cd(7); hbr1_27->Draw();

  c6->cd(8); hbr2_1->Draw();
  c6->cd(9); hbr2_2->Draw();
  c6->cd(10); hbr2_3->Draw();
  c6->cd(11); hbr2_4->Draw();
  c6->cd(12); hbr2_5->Draw();
  c6->cd(13); hbr2_6->Draw();
  c6->cd(14); hbr2_7->Draw();
  c6->cd(15); hbr2_8->Draw();
  c6->cd(16); hbr2_9->Draw();
  c6->cd(17); hbr2_10->Draw();
  c6->cd(18); hbr2_11->Draw();
  c6->cd(19); hbr2_12->Draw();

  c6->cd(20); hlol_1->Draw();
  c6->cd(21); hlol_2->Draw();
  c6->cd(22); hlol_3->Draw();
  c6->cd(23); hlol_4->Draw();
  c6->cd(24); hlol_5->Draw();
  c6->cd(25); hlol_6->Draw();

  TCanvas *c7 = new TCanvas("c7","c7",1200,900);
  c7->Divide(5,5);
  c7->cd(1); hlol_7->Draw();
  c7->cd(2); hlol_8->Draw();
  c7->cd(3); hlol_9->Draw();
  c7->cd(4); hlol_10->Draw();
  c7->cd(5); hlol_11->Draw();
  c7->cd(6); hlol_12->Draw();
  c7->cd(7); hlol_13->Draw();
  c7->cd(8); hlol_14->Draw();
  c7->cd(9); hlol_15->Draw();
  c7->cd(10); hlol_16->Draw();
  c7->cd(11); hlol_17->Draw();
  c7->cd(12); hlol_18->Draw();
  c7->cd(13); hlol_19->Draw();
  c7->cd(14); hlol_20->Draw();
  c7->cd(15); hlol_21->Draw();
  c7->cd(16); hlol_22->Draw();

  c7->cd(17); hbr3_1->Draw();
  c7->cd(18); hbr3_2->Draw();
  c7->cd(19); hbr3_3->Draw();
  c7->cd(20); hbr3_4->Draw();
  c7->cd(21); hbr3_5->Draw();
  c7->cd(22); hbr3_6->Draw();
  c7->cd(23); hbr3_7->Draw();
  c7->cd(24); hbr3_8->Draw();
  c7->cd(25); hbr3_9->Draw();

  TCanvas *c8 = new TCanvas("c8","c8",1200,900);
  c8->Divide(5,5);
  c8->cd(1); hbr3_10->Draw();
  c8->cd(2); hbr3_11->Draw();
  c8->cd(3); hbr3_12->Draw();
  c8->cd(4); hbr3_13->Draw();
  c8->cd(5); hbr3_14->Draw();
  c8->cd(6); hbr3_15->Draw();
  c8->cd(7); hbr3_16->Draw();
  c8->cd(8); hbr3_17->Draw();
  c8->cd(9); hbr3_18->Draw();
  c8->cd(10); hbr3_19->Draw();
  c8->cd(11); hbr3_20->Draw();
  c8->cd(12); hbr3_21->Draw();
  c8->cd(13); hbr3_22->Draw();
  c8->cd(14); hbr3_23->Draw();
  c8->cd(15); hbr3_24->Draw();
  c8->cd(16); hbr3_25->Draw();
  c8->cd(17); hbr3_26->Draw();
  c8->cd(18); hbr3_27->Draw();
  c8->cd(19); hbr3_28->Draw();
  c8->cd(20); hbr3_29->Draw();
  c8->cd(21); hbr3_30->Draw();
  c8->cd(22); hbr3_31->Draw();
  c8->cd(23); hbr3_32->Draw();
  c8->cd(24); hbr3_33->Draw();
  c8->cd(25); hbr3_34->Draw();
  
  TCanvas *c9 = new TCanvas("c9","c9",1200,900);
  c9->Divide(5,5);
  c9->cd(1); hbr3_35->Draw();
  c9->cd(2); hbr3_36->Draw();
  c9->cd(3); hbr3_37->Draw();
  c9->cd(4); hbr3_38->Draw();
  c9->cd(5); hbr3_39->Draw();
  c9->cd(6); hbr3_40->Draw();
  c9->cd(7); hbr3_41->Draw();
  c9->cd(8); hbr3_42->Draw();
  c9->cd(9); hbr3_43->Draw();
  c9->cd(10); hbr3_44->Draw();
  c9->cd(11); hbr3_45->Draw();
  c9->cd(12); hbr3_46->Draw();
  c9->cd(13); hbr3_47->Draw();
  c9->cd(14); hbr3_48->Draw();
  c9->cd(15); hbr3_49->Draw();
  c9->cd(16); hbr3_50->Draw();
  c9->cd(17); hbr3_51->Draw();
  c9->cd(18); hbr3_52->Draw();
  c9->cd(19); hbr3_53->Draw();
  c9->cd(20); hbr3_54->Draw();
  c9->cd(21); hbr3_55->Draw();
  c9->cd(22); hbr3_56->Draw();
  c9->cd(23); hbr3_57->Draw();
  c9->cd(24); hbr3_58->Draw();

  c9->cd(25); hbr4_1->Draw();

  TCanvas *c10 = new TCanvas("c10","c10",1200,900);
  c10->Divide(5,5);
  c10->cd(1); hbr4_2->Draw();
  c10->cd(2); hbr4_3->Draw();
  c10->cd(3); hbr4_4->Draw();
  c10->cd(4); hbr4_5->Draw();
  c10->cd(5); hbr4_6->Draw();
  c10->cd(6); hbr4_7->Draw();
  c10->cd(7); hbr4_8->Draw();
  c10->cd(8); hbr4_9->Draw();
  c10->cd(9); hbr4_10->Draw();
  c10->cd(10); hbr4_11->Draw();

  c10->cd(11); hbr4_12->Draw();
  c10->cd(12); hbr4_13->Draw();
  c10->cd(13); hbr4_14->Draw();
  c10->cd(14); hbr4_15->Draw();
  c10->cd(15); hbr4_16->Draw();
  c10->cd(16); hbr4_17->Draw();
  c10->cd(17); hbr4_18->Draw();
  c10->cd(18); hbr4_19->Draw();
  c10->cd(19); hbr4_20->Draw();
  c10->cd(20); hbr4_21->Draw();
  c10->cd(21); hbr4_22->Draw();

  c10->cd(22); hol_1->Draw();
  c10->cd(23); hol_2->Draw();
  c10->cd(24); hol_3->Draw();
  c10->cd(25); hol_4->Draw();

  TCanvas *c11 = new TCanvas("c11","c11",1200,900);
  c11->Divide(5,5);
  c11->cd(1); hol_5->Draw();
  c11->cd(2); hol_6->Draw();
  c11->cd(3); hol_7->Draw();
  c11->cd(4); hol_8->Draw();
  c11->cd(5); hol_9->Draw();
  c11->cd(6); hol_10->Draw();
  c11->cd(7); hol_11->Draw();
  c11->cd(8); hol_12->Draw();

  c11->cd(9); hvis_1->Draw();
  c11->cd(10); hvis_2->Draw();
  c11->cd(11); hvis_3->Draw();
  c11->cd(12); hvis_4->Draw();
  c11->cd(13); hvis_5->Draw();
  c11->cd(14); hvis_6->Draw();
  c11->cd(15); hvis_7->Draw();
  c11->cd(16); hvis_8->Draw();
  c11->cd(17); hvis_9->Draw();
  c11->cd(18); hvis_10->Draw();
  c11->cd(19); hvis_11->Draw();
  c11->cd(20); hvis_12->Draw();
  c11->cd(21); hvis_13->Draw();
  c11->cd(22); hvis_14->Draw();
  c11->cd(23); hvis_15->Draw();
  c11->cd(24); hvis_16->Draw();
  c11->cd(25); hvis_17->Draw();

  TCanvas *c12 = new TCanvas("c12","c12",1200,900);
  c12->Divide(5,5);
  c12->cd(1); hvis_18->Draw();
  c12->cd(2); hvis_19->Draw();
  c12->cd(3); hvis_20->Draw();
  c12->cd(4); hvis_21->Draw();
  c12->cd(5); hvis_22->Draw();
  c12->cd(6); hvis_23->Draw();

  c12->cd(7); hstl_1->Draw();
  c12->cd(8); hstl_2->Draw();
  c12->cd(9); hstl_3->Draw();
  c12->cd(10); hstl_4->Draw();
  c12->cd(11); hstl_5->Draw();
  c12->cd(12); hstl_6->Draw();

  c12->cd(13); hbrm_1->Draw();
  c12->cd(14); hbrm_2->Draw();
  c12->cd(15); hbrm_3->Draw();
  c12->cd(16); hbrm_4->Draw();
  c12->cd(17); hbrm_5->Draw();
  c12->cd(18); hbrm_6->Draw();
  c12->cd(19); hbrm_7->Draw();
  c12->cd(20); hbrm_8->Draw();
  c12->cd(21); hbrm_9->Draw();
  c12->cd(22); hbrm_10->Draw();
  c12->cd(23); hbrm_11->Draw();

  c12->cd(24); hcme_1->Draw();
  c12->cd(25); hcme_2->Draw();

  TCanvas *c13 = new TCanvas("c13","c13",1200,900);
  c13->Divide(5,5);
  c13->cd(1); hcme_3->Draw();
  c13->cd(2); hcme_4->Draw();
  c13->cd(3); hcme_5->Draw();
  c13->cd(4); hcme_6->Draw();

  c13->cd(5); hanc_1->Draw();
  c13->cd(6); hanc_2->Draw();
  c13->cd(7); hanc_3->Draw();
  c13->cd(8); hanc_4->Draw();
  c13->cd(9); hanc_5->Draw();
  c13->cd(10); hanc_6->Draw();
  c13->cd(11); hanc_7->Draw();
  c13->cd(12); hanc_8->Draw();
  c13->cd(13); hanc_9->Draw();
  c13->cd(14); hanc_10->Draw();
  c13->cd(15); hanc_11->Draw();

  c13->cd(16); hlem_1->Draw();
  c13->cd(17); hlem_2->Draw();
  c13->cd(18); hlem_3->Draw();
  c13->cd(19); hlem_4->Draw();
  c13->cd(20); hlem_5->Draw();
  c13->cd(21); hlem_6->Draw();
  c13->cd(22); hlem_7->Draw();
  c13->cd(23); hlem_8->Draw();

  c13->cd(24); hstw_1->Draw();
  c13->cd(25); hstw_2->Draw();

  TCanvas *c14 = new TCanvas("c14","c14",1200,900);
  c14->Divide(5,5);
  c14->cd(1); hstw_3->Draw();
  c14->cd(2); hstw_4->Draw();
  c14->cd(3); hstw_5->Draw();
  c14->cd(4); hstw_6->Draw();
  c14->cd(5); hstw_7->Draw();
  c14->cd(6); hstw_8->Draw();
  c14->cd(7); hstw_9->Draw();
  c14->cd(8); hstw_10->Draw();
  c14->cd(9); hstw_11->Draw();
  c14->cd(10); hstw_12->Draw();

  c14->cd(11); hstw_13->Draw();
  c14->cd(12); hstw_14->Draw();
  c14->cd(13); hstw_15->Draw();
  c14->cd(14); hstw_16->Draw();
  c14->cd(15); hstw_17->Draw();
  c14->cd(16); hstw_18->Draw();
  c14->cd(17); hstw_19->Draw();
  c14->cd(18); hstw_20->Draw();
  c14->cd(19); hstw_21->Draw();
  c14->cd(20); hstw_22->Draw();
  c14->cd(21); hstw_23->Draw();

  c14->cd(22); hspt_1->Draw();
  c14->cd(23); hspt_2->Draw();
  c14->cd(24); hspt_3->Draw();
  c14->cd(25); hspt_4->Draw();

  TCanvas *c15 = new TCanvas("c15","c15",1200,900);
  c15->Divide(5,5);
  c15->cd(1); hspt_5->Draw();
  c15->cd(2); hspt_6->Draw();
  c15->cd(3); hspt_7->Draw();
  c15->cd(4); hspt_8->Draw();
  c15->cd(5); hspt_9->Draw();
  c15->cd(6); hspt_10->Draw();
  c15->cd(7); hspt_11->Draw();
  c15->cd(8); hspt_12->Draw();
  c15->cd(9); hspt_13->Draw();
  c15->cd(10); hspt_14->Draw();

  c15->cd(11); hmgo_1->Draw();
  c15->cd(12); hmgo_2->Draw();
  c15->cd(13); hmgo_3->Draw();
  c15->cd(14); hmgo_4->Draw();
  c15->cd(15); hmgo_5->Draw();
  c15->cd(16); hmgo_6->Draw();
  c15->cd(17); hmgo_7->Draw();
  c15->cd(18); hmgo_8->Draw();
  c15->cd(19); hmgo_9->Draw();
  c15->cd(20); hmgo_10->Draw();

  c15->cd(21); hmgt_1->Draw();
  c15->cd(22); hmgt_2->Draw();
  c15->cd(23); hmgt_3->Draw();
  c15->cd(24); hmgt_4->Draw();
  c15->cd(25); hmgt_5->Draw();

  TCanvas *c16 = new TCanvas("c16","c16",1200,900);
  c16->Divide(5,5);
  c16->cd(1); hmgt_6->Draw();
  c16->cd(2); hmgt_7->Draw();
  c16->cd(3); hmgt_8->Draw();
  c16->cd(4); hmgt_9->Draw();
  c16->cd(5); hmgt_10->Draw();
  c16->cd(6); hmgt_11->Draw();
  c16->cd(7); hmgt_12->Draw();

  c16->cd(8); hsig_1->Draw();
  c16->cd(9); hsig_2->Draw();
  c16->cd(10); hsig_3->Draw();
  c16->cd(11); hsig_4->Draw();
  c16->cd(12); hsig_5->Draw();
  c16->cd(13); hsig_6->Draw();
  c16->cd(14); hsig_7->Draw();
  c16->cd(15); hsig_8->Draw();
  c16->cd(16); hsig_9->Draw();
  c16->cd(17); hsig_10->Draw();
  c16->cd(18); hsig_11->Draw();
  c16->cd(19); hsig_12->Draw();

  // additional ones ...
  c16->cd(20); hpio14->Draw();
  c16->cd(21); hpio15->Draw();

  c16->cd(22); htro_1->Draw();
  c16->cd(23); htro_2->Draw();
  c16->cd(24); htro_3->Draw();
  c16->cd(25); htro_4->Draw();
  
  TCanvas *c17 = new TCanvas("c17","c17",1200,900);
  c17->Divide(5,5);
  c17->cd(1); htro_5->Draw();
  c17->cd(2); htro_6->Draw();
  c17->cd(3); htro_7->Draw();
  c17->cd(4); htro_8->Draw();
  c17->cd(5); htro_9->Draw();
  c17->cd(6); htro_10->Draw();
  c17->cd(7); htro_11->Draw();
  c17->cd(8); htro_12->Draw();
  c17->cd(9); htro_13->Draw();
  c17->cd(10); htro_14->Draw();

  c17->cd(11); htro_15->Draw();
  c17->cd(12); htro_16->Draw();
  c17->cd(13); htro_17->Draw();
  c17->cd(14); htro_18->Draw();
  c17->cd(15); htro_19->Draw();
  c17->cd(16); htro_20->Draw();
  c17->cd(17); htro_21->Draw();
  c17->cd(18); htro_22->Draw();
  c17->cd(19); htro_23->Draw();
  c17->cd(20); htro_24->Draw();

  c17->cd(21); htro_25->Draw();
  c17->cd(22); htro_26->Draw();
  c17->cd(23); htro_27->Draw();
  c17->cd(24); htro_28->Draw();
  c17->cd(25); htro_29->Draw();

  TCanvas *c18 = new TCanvas("c18","c18",1200,900);
  c18->Divide(5,5);
  c18->cd(1); htro_30->Draw();
  c18->cd(2); htro_31->Draw();
  c18->cd(3); htro_32->Draw();
  c18->cd(4); htro_33->Draw();
  c18->cd(5); htro_34->Draw();
  c18->cd(6); htro_35->Draw();
  c18->cd(7); htro_36->Draw();
  c18->cd(8); htro_37->Draw();
  c18->cd(9); htro_38->Draw();
  c18->cd(10); htro_39->Draw();

  c18->cd(11); htro_40->Draw();
  c18->cd(12); htro_41->Draw();
  c18->cd(13); htro_42->Draw();
  c18->cd(14); htro_43->Draw();
}
