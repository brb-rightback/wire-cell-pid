#include <vector>

void print_valid(TString run_info = "7019_5_264"){
  TChain *T = new TChain("T_tagger","T_tagger");
  TString filename = "nue_" + run_info + ".root";
  T->AddFile(filename);

  std::vector<int> flag_info_vec;
  
  // cosmic information
  float cosmic_flag;
  float cosmic_n_solid_tracks;
  float cosmic_energy_main_showers;
  float cosmic_energy_direct_showers;
  float cosmic_energy_indirect_showers;
  float cosmic_n_direct_showers;
  float cosmic_n_indirect_showers;
  float cosmic_n_main_showers;
  float cosmic_filled;
  T->SetBranchAddress("cosmic_flag",&cosmic_flag);
  T->SetBranchAddress("cosmic_energy_main_showers",&cosmic_energy_main_showers);
  T->SetBranchAddress("cosmic_energy_indirect_showers",&cosmic_energy_indirect_showers);
  T->SetBranchAddress("cosmic_energy_direct_showers",&cosmic_energy_direct_showers);
  T->SetBranchAddress("cosmic_n_direct_showers",&cosmic_n_direct_showers);
  T->SetBranchAddress("cosmic_n_indirect_showers",&cosmic_n_indirect_showers);
  T->SetBranchAddress("cosmic_n_solid_tracks",&cosmic_n_solid_tracks);
  T->SetBranchAddress("cosmic_n_main_showers",&cosmic_n_main_showers);
  T->SetBranchAddress("cosmic_filled",&cosmic_filled);

   // shower gap identification
  float gap_flag;
  float gap_flag_prolong_u;
  float gap_flag_prolong_v;
  float gap_flag_prolong_w;
  float gap_flag_parallel;
  float gap_n_points;
  float gap_n_bad;
  float gap_energy;
  float gap_num_valid_tracks;
  float gap_flag_single_shower;
  float gap_filled;
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

  float mip_quality_flag;
  float mip_quality_energy;
  float mip_quality_overlap;
  float mip_quality_n_showers;
  float mip_quality_n_tracks;
  float mip_quality_flag_inside_pi0;
  float mip_quality_n_pi0_showers;
  float mip_quality_shortest_length;
  float mip_quality_acc_length;
  float mip_quality_shortest_angle;
  float mip_quality_flag_proton;
  float mip_quality_filled;

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

   // mip identification
  float mip_flag;
  float mip_energy;
  float mip_n_end_reduction;    
  float mip_n_first_mip;
  float mip_n_first_non_mip;
  float mip_n_first_non_mip_1;
  float mip_n_first_non_mip_2;
  float mip_vec_dQ_dx_0;
  float mip_vec_dQ_dx_1;
  float mip_max_dQ_dx_sample;
  float mip_n_below_threshold;
  float mip_n_below_zero;
  float mip_n_lowest;
  float mip_n_highest;
  float mip_lowest_dQ_dx;
  float mip_highest_dQ_dx;
  float mip_medium_dQ_dx;
  float mip_stem_length;
  float mip_length_main;
  float mip_length_total;
  float mip_angle_beam;
  float mip_iso_angle;
  float mip_n_vertex;
  float mip_n_good_tracks;
  float mip_E_indirect_max_energy;
  float mip_flag_all_above;
  float mip_min_dQ_dx_5;
  float mip_n_other_vertex; 
  float mip_n_stem_size;
  float mip_flag_stem_trajectory;
  float mip_min_dis;
  float mip_filled;
  
  // extra
  float mip_vec_dQ_dx_2;
  float mip_vec_dQ_dx_3;
  float mip_vec_dQ_dx_4;
  float mip_vec_dQ_dx_5;
  float mip_vec_dQ_dx_6;
  float mip_vec_dQ_dx_7;
  float mip_vec_dQ_dx_8;
  float mip_vec_dQ_dx_9;
  float mip_vec_dQ_dx_10;
  float mip_vec_dQ_dx_11;
  float mip_vec_dQ_dx_12;
  float mip_vec_dQ_dx_13;
  float mip_vec_dQ_dx_14;
  float mip_vec_dQ_dx_15;
  float mip_vec_dQ_dx_16;
  float mip_vec_dQ_dx_17;
  float mip_vec_dQ_dx_18;
  float mip_vec_dQ_dx_19;
  

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

  
  float pio_flag;
  float pio_mip_id;
  float pio_filled;
  float pio_flag_pio;
  
  float pio_1_flag;
  float pio_1_mass;
  float pio_1_pio_type;
  float pio_1_energy_1;
  float pio_1_energy_2;
  float pio_1_dis_1;
  float pio_1_dis_2;
  
  std::vector<float> *pio_2_v_dis2 = new std::vector<float>;
  std::vector<float> *pio_2_v_angle2 = new std::vector<float>;
  std::vector<float> *pio_2_v_acc_length = new std::vector<float>;
  std::vector<float> *pio_2_v_flag = new std::vector<float>;

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

   // stem direction
  float stem_dir_flag;
  float stem_dir_flag_single_shower;
  float stem_dir_filled;
  float stem_dir_angle;
  float stem_dir_energy;
  float stem_dir_angle1;
  float stem_dir_angle2;
  float stem_dir_angle3;
  float stem_dir_ratio;
  
  T->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
  T->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  T->SetBranchAddress("stem_dir_filled",&stem_dir_filled);
  T->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
  T->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
  T->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
  T->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
  T->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
  T->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);

  float br_filled;
  float br1_flag;
  
  //bad reconstruction 1_1
  float br1_1_flag;
  float br1_1_shower_type;
  float br1_1_vtx_n_segs;
  float br1_1_energy;
  float br1_1_n_segs;
  float br1_1_flag_sg_topology;
  float br1_1_flag_sg_trajectory;
  float br1_1_sg_length;
  
  // bad reconstruction 1_2
  float br1_2_flag;
  float br1_2_energy;
  float br1_2_n_connected;
  float br1_2_max_length;
  float br1_2_n_connected_1;
  float br1_2_vtx_n_segs;
  float br1_2_n_shower_segs;
  float br1_2_max_length_ratio;
  float br1_2_shower_length;
  
  // bad_reconstruction 1_3
  float br1_3_flag;
  float br1_3_energy;
  float br1_3_n_connected_p;
  float br1_3_max_length_p;
  float br1_3_n_shower_segs;
  float br1_3_flag_sg_topology;
  float br1_3_flag_sg_trajectory;
  float br1_3_n_shower_main_segs;
  float br1_3_sg_length;
  
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
  T->SetBranchAddress("br1_3_flag_sg_trajectory",&br1_3_flag_sg_trajectory);
  T->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
  T->SetBranchAddress("br1_3_sg_length",&br1_3_sg_length);

   // bad reconstruction 2
  float br2_flag;
  float br2_flag_single_shower;
  float br2_num_valid_tracks;
  float br2_energy;
  float br2_angle1;
  float br2_angle2;
  float br2_angle;
  float br2_angle3;
  float br2_n_shower_main_segs;
  float br2_max_angle;
  float br2_sg_length;
  float br2_flag_sg_trajectory;
  
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

  // low energy overlap
  float lol_flag;

  std::vector<float> *lol_1_v_energy = new std::vector<float>;
  std::vector<float> *lol_1_v_vtx_n_segs= new std::vector<float>;
  std::vector<float> *lol_1_v_nseg= new std::vector<float>;
  std::vector<float> *lol_1_v_angle= new std::vector<float>;
  std::vector<float> *lol_1_v_flag= new std::vector<float>;
  
  std::vector<float> *lol_2_v_length= new std::vector<float>;
  std::vector<float> *lol_2_v_angle= new std::vector<float>;
  std::vector<float> *lol_2_v_type= new std::vector<float>;
  std::vector<float> *lol_2_v_vtx_n_segs= new std::vector<float>;
  std::vector<float> *lol_2_v_energy= new std::vector<float>;
  std::vector<float> *lol_2_v_shower_main_length= new std::vector<float>;
  std::vector<float> *lol_2_v_flag_dir_weak= new std::vector<float>;
  std::vector<float> *lol_2_v_flag = new std::vector<float>;
  
  float lol_3_angle_beam;
  float lol_3_n_valid_tracks;
  float lol_3_min_angle;
  float lol_3_vtx_n_segs;
  float lol_3_energy;
  float lol_3_shower_main_length;
  float lol_3_n_out;
  float lol_3_n_sum;    
  float lol_3_flag;

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
  
   //bad reconstruction 3
  float br3_1_energy;
  float br3_1_n_shower_segments;
  float br3_1_sg_flag_trajectory;
  float br3_1_sg_direct_length;
  float br3_1_sg_length;
  float br3_1_total_main_length;
  float br3_1_total_length;
  float br3_1_iso_angle;
  float br3_1_sg_flag_topology;
  float br3_1_flag;

  float br3_2_n_ele;
  float br3_2_n_other;
  float br3_2_energy;
  float br3_2_total_main_length;
  float br3_2_total_length;
  float br3_2_other_fid;
  float br3_2_flag;
  
  std::vector<float> *br3_3_v_energy = new std::vector<float>;
  std::vector<float> *br3_3_v_angle = new std::vector<float>;
  std::vector<float> *br3_3_v_dir_length = new std::vector<float>;
  std::vector<float> *br3_3_v_length = new std::vector<float>;
  std::vector<float> *br3_3_v_flag = new std::vector<float>;
  
  float br3_4_acc_length;
  float br3_4_total_length;
  float br3_4_energy;
  float br3_4_flag;
  
  std::vector<float> *br3_5_v_dir_length = new std::vector<float>;
  std::vector<float> *br3_5_v_total_length = new std::vector<float>;
  std::vector<float> *br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  std::vector<float> *br3_5_v_n_seg = new std::vector<float>;
  std::vector<float> *br3_5_v_angle = new std::vector<float>;
  std::vector<float> *br3_5_v_sg_length = new std::vector<float>;
  std::vector<float> *br3_5_v_energy = new std::vector<float>;
  std::vector<float> *br3_5_v_n_main_segs = new std::vector<float>;
  std::vector<float> *br3_5_v_n_segs = new std::vector<float>;
  std::vector<float> *br3_5_v_shower_main_length = new std::vector<float>;
  std::vector<float> *br3_5_v_shower_total_length = new std::vector<float>;
  std::vector<float> *br3_5_v_flag = new std::vector<float>;
  
  std::vector<float> *br3_6_v_angle = new std::vector<float>;
  std::vector<float> *br3_6_v_angle1 = new std::vector<float>;
  std::vector<float> *br3_6_v_flag_shower_trajectory = new std::vector<float>;
  std::vector<float> *br3_6_v_direct_length = new std::vector<float>;
  std::vector<float> *br3_6_v_length = new std::vector<float>;
  std::vector<float> *br3_6_v_n_other_vtx_segs = new std::vector<float>;
  std::vector<float> *br3_6_v_energy = new std::vector<float>;
  std::vector<float> *br3_6_v_flag = new std::vector<float>;
  
  float br3_7_energy;
  float br3_7_min_angle;
  float br3_7_sg_length;
  float br3_7_shower_main_length;
  float br3_7_flag;
  
  float br3_8_max_dQ_dx;
  float br3_8_energy;
  float br3_8_n_main_segs;
  float br3_8_shower_main_length;
  float br3_8_shower_length;
  float br3_8_flag;
  
  float br3_flag;

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
  float br4_1_shower_main_length;
  float br4_1_shower_total_length;
  float br4_1_min_dis;
  float br4_1_energy;
  float br4_1_flag_avoid_muon_check;
  float br4_1_n_vtx_segs;
  float br4_1_n_main_segs;
  float br4_1_flag;
  
  float br4_2_ratio_45;
  float br4_2_ratio_35;
  float br4_2_ratio_25;
  float br4_2_ratio_15;
  float br4_2_energy;
  float br4_2_ratio1_45;
  float br4_2_ratio1_35;
  float br4_2_ratio1_25;
  float br4_2_ratio1_15;
  float br4_2_iso_angle;
  float br4_2_iso_angle1;
  float br4_2_angle;
  float br4_2_flag;
  
  float br4_flag;
  
  T->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
  T->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
  T->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
  T->SetBranchAddress("br4_1_energy", &br4_1_energy);
  T->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
  T->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
  T->SetBranchAddress("br4_1_br4_1_n_main_segs", &br4_1_n_main_segs);
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
  
  float hol_1_n_valid_tracks;
  float hol_1_min_angle;
  float hol_1_energy;
  float hol_1_flag_all_shower;
  float hol_1_min_length;
  float hol_1_flag;
  
  float hol_2_min_angle;
  float hol_2_medium_dQ_dx;
  float hol_2_ncount;
  float hol_2_energy;
  float hol_2_flag;
  
  float hol_flag;

  T->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
  T->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
  T->SetBranchAddress("hol_1_energy", &hol_1_energy);
  T->SetBranchAddress("hol_1_flag_all_shower", &hol_1_flag_all_shower);
  T->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
  T->SetBranchAddress("hol_1_flag", &hol_1_flag);
  
  T->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
  T->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
  T->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
  T->SetBranchAddress("hol_2_energy", &hol_2_energy);
  T->SetBranchAddress("hol_2_flag", &hol_2_flag);
  
  T->SetBranchAddress("hol_flag", &hol_flag);

   // vertex inside shower
  float vis_1_filled;
  float vis_1_n_vtx_segs;
  float vis_1_energy;
  float vis_1_num_good_tracks;
  float vis_1_max_angle;
  float vis_1_max_shower_angle;
  float vis_1_tmp_length1;
  float vis_1_tmp_length2;
  float vis_1_particle_type;
  float vis_1_flag;
  
  float vis_2_filled;
  float vis_2_n_vtx_segs;
  float vis_2_min_angle;
  float vis_2_min_weak_track;
  float vis_2_angle_beam;
  float vis_2_min_angle1;
  float vis_2_iso_angle1;
  float vis_2_min_medium_dQ_dx;
  float vis_2_min_length;
  float vis_2_sg_length;
  float vis_2_max_angle;
  float vis_2_max_weak_track;
  float vis_2_flag;
  
  float vis_flag;

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
  T->SetBranchAddress("vis_2_flag",&vis_2_flag);
  
  T->SetBranchAddress("vis_flag",&vis_flag);
  
  float stem_len_energy;
  float stem_len_length;
  float stem_len_flag_avoid_muon_check;
  float stem_len_num_daughters;
  float stem_len_daughter_length;
  float stem_len_flag;

  T->SetBranchAddress("stem_len_energy", &stem_len_energy);
  T->SetBranchAddress("stem_len_length", &stem_len_length);
  T->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
  T->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
  T->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);
  T->SetBranchAddress("stem_len_flag", &stem_len_flag);

  float brm_n_mu_segs;
  float brm_Ep;
  float brm_energy;
  float brm_acc_length;
  float brm_shower_total_length;
  float brm_connected_length;
  float brm_n_size;
  float brm_acc_direct_length;
  float brm_n_shower_main_segs;
  float brm_n_mu_main;
  float brm_flag;

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
  T->SetBranchAddress("brm_flag",&brm_flag);

  // compare with muon
  float cme_mu_energy;
  float cme_energy;
  float cme_mu_length;
  float cme_length;
  float cme_angle_beam;
  float cme_flag;

  T->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
  T->SetBranchAddress("cme_energy",&cme_energy);
  T->SetBranchAddress("cme_mu_length",&cme_mu_length);
  T->SetBranchAddress("cme_length",&cme_length);
  T->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
  T->SetBranchAddress("cme_flag",&cme_flag);

 // angular cut
  float anc_energy;
  float anc_angle;
  float anc_max_angle;
  float anc_max_length;
  float anc_acc_forward_length;
  float anc_acc_backward_length;
  float anc_acc_forward_length1;
  float anc_shower_main_length;
  float anc_shower_total_length;
  float anc_flag_main_outside;
  float anc_flag;

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
  
  // low-energy michel
  float lem_shower_total_length;
  float lem_shower_main_length;
  float lem_n_3seg;
  float lem_e_charge;
  float lem_e_dQdx;
  float lem_shower_num_segs;
  float lem_shower_num_main_segs;
  float lem_flag;
  
  T->SetBranchAddress("lem_shower_total_length",&lem_shower_total_length);
  T->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
  T->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
  T->SetBranchAddress("lem_e_charge",&lem_e_charge);
  T->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
  T->SetBranchAddress("lem_shower_num_segs",&lem_shower_num_segs);
  T->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  T->SetBranchAddress("lem_flag",&lem_flag);

   // shower to wall
  float stw_1_energy;
  float stw_1_dis;
  float stw_1_dQ_dx;
  float stw_1_flag_single_shower;
  float stw_1_n_pi0;
  float stw_1_num_valid_tracks;
  float stw_1_flag;
  
  std::vector<float> *stw_2_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *stw_2_v_energy = new std::vector<float>;
  std::vector<float> *stw_2_v_angle = new std::vector<float>;
  std::vector<float> *stw_2_v_dir_length = new std::vector<float>;
  std::vector<float> *stw_2_v_max_dQ_dx = new std::vector<float>;
  std::vector<float> *stw_2_v_flag = new std::vector<float>;
  
  std::vector<float> *stw_3_v_angle = new std::vector<float>;
  std::vector<float> *stw_3_v_dir_length = new std::vector<float>;
  std::vector<float> *stw_3_v_energy = new std::vector<float>;
  std::vector<float> *stw_3_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *stw_3_v_flag = new std::vector<float>;
  
  std::vector<float> *stw_4_v_angle = new std::vector<float>;
  std::vector<float> *stw_4_v_dis = new std::vector<float>;
  std::vector<float> *stw_4_v_energy = new std::vector<float>;
  std::vector<float> *stw_4_v_flag = new std::vector<float>;
  
  float stw_flag;
  
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

   // single photon  cases ...
  float spt_flag_single_shower;
  float spt_energy;
  float spt_shower_main_length;
  float spt_shower_total_length;
  float spt_angle_beam;
  float spt_angle_vertical;
  float spt_max_dQ_dx;
  float spt_angle_beam_1;
  float spt_angle_drift;
  float spt_angle_drift_1;
  float spt_num_valid_tracks;
  float spt_n_vtx_segs;
  float spt_max_length;
  float spt_flag;

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

  // multiple gamma I
  float mgo_energy;
  float mgo_max_energy;
  float mgo_total_energy;
  float mgo_n_showers;
  float mgo_max_energy_1;
  float mgo_max_energy_2;
  float mgo_total_other_energy;
  float mgo_n_total_showers;
  float mgo_total_other_energy_1;
  float mgo_flag;
  
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
  
  // multiple gamma II
  float mgt_flag_single_shower;
  float mgt_max_energy;
  float mgt_energy;
  float mgt_total_other_energy;
  float mgt_max_energy_1;
  float mgt_e_indirect_max_energy;
  float mgt_e_direct_max_energy;
  float mgt_n_direct_showers;
  float mgt_e_direct_total_energy;
  float mgt_flag_indirect_max_pio;
  float mgt_e_indirect_total_energy;
  float mgt_flag;

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

   //single shower
  std::vector<float> *sig_1_v_angle = new std::vector<float>;
  std::vector<float> *sig_1_v_flag_single_shower= new std::vector<float>;
  std::vector<float> *sig_1_v_energy= new std::vector<float>;
  std::vector<float> *sig_1_v_energy_1= new std::vector<float>;
  std::vector<float> *sig_1_v_flag= new std::vector<float>;
  
  std::vector<float> *sig_2_v_energy= new std::vector<float>;
  std::vector<float> *sig_2_v_shower_angle= new std::vector<float>;
  std::vector<float> *sig_2_v_flag_single_shower= new std::vector<float>;
  std::vector<float> *sig_2_v_medium_dQ_dx= new std::vector<float>;
  std::vector<float> *sig_2_v_start_dQ_dx= new std::vector<float>;
  std::vector<float> *sig_2_v_flag= new std::vector<float>;
  
  float sig_flag;

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

  // track overclustering ..
  std::vector<float> *tro_1_v_particle_type= new std::vector<float>;
  std::vector<float> *tro_1_v_flag_dir_weak= new std::vector<float>;
  std::vector<float> *tro_1_v_min_dis = new std::vector<float>;
  std::vector<float> *tro_1_v_sg1_length = new std::vector<float>;
  std::vector<float> *tro_1_v_shower_main_length = new std::vector<float>;
  std::vector<float> *tro_1_v_max_n_vtx_segs= new std::vector<float>;
  std::vector<float> *tro_1_v_tmp_length = new std::vector<float>;
  std::vector<float> *tro_1_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *tro_1_v_dQ_dx_cut = new std::vector<float>;
  std::vector<float> *tro_1_v_flag_shower_topology= new std::vector<float>;
  std::vector<float> *tro_1_v_flag= new std::vector<float>;
  
  std::vector<float> *tro_2_v_energy = new std::vector<float>;
  std::vector<float> *tro_2_v_stem_length = new std::vector<float>;
  std::vector<float> *tro_2_v_iso_angle = new std::vector<float>;
  std::vector<float> *tro_2_v_max_length = new std::vector<float>;
  std::vector<float> *tro_2_v_angle = new std::vector<float>;
  std::vector<float> *tro_2_v_flag= new std::vector<float>;
  
  float tro_3_stem_length;
  float tro_3_n_muon_segs;
  float tro_3_energy;
  float tro_3_flag;
  
  std::vector<float> *tro_4_v_dir2_mag = new std::vector<float>;
  std::vector<float> *tro_4_v_angle = new std::vector<float>;
  std::vector<float> *tro_4_v_angle1 = new std::vector<float>;
  std::vector<float> *tro_4_v_angle2 = new std::vector<float>;
  std::vector<float> *tro_4_v_length = new std::vector<float>;
  std::vector<float> *tro_4_v_length1 = new std::vector<float>;
  std::vector<float> *tro_4_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *tro_4_v_end_dQ_dx = new std::vector<float>;
  std::vector<float> *tro_4_v_energy = new std::vector<float>;
  std::vector<float> *tro_4_v_shower_main_length = new std::vector<float>;
  std::vector<float> *tro_4_v_flag_shower_trajectory= new std::vector<float>;
  std::vector<float> *tro_4_v_flag= new std::vector<float>;
  
  std::vector<float> *tro_5_v_max_angle = new std::vector<float>;
  std::vector<float> *tro_5_v_min_angle = new std::vector<float>;
  std::vector<float> *tro_5_v_max_length = new std::vector<float>;
  std::vector<float> *tro_5_v_iso_angle = new std::vector<float>;
  std::vector<float> *tro_5_v_n_vtx_segs= new std::vector<float>;
  std::vector<float> *tro_5_v_min_count= new std::vector<float>;
  std::vector<float> *tro_5_v_max_count= new std::vector<float>;
  std::vector<float> *tro_5_v_energy = new std::vector<float>;
  std::vector<float> *tro_5_v_flag = new std::vector<float>;
  
  float tro_flag;
  
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

  // cosmic tagger
  float cosmict_flag_1; // fiducial volume vertex
  float cosmict_flag_2;  // single muon
  float cosmict_flag_3;  // single muon (long)
  float cosmict_flag_4;  // kinematics muon
  float cosmict_flag_5; // kinematics muon (long)
  float cosmict_flag_6; // special ...
  float cosmict_flag_7;  // muon+ michel
  float cosmict_flag_8;  // muon + michel + special
  float cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
  std::vector<float> *cosmict_flag_10 = new std::vector<float>;  // front upstream (dirt)
  float cosmict_flag;
  
  // single muon
  float cosmict_2_filled;
  float cosmict_2_particle_type;
  float cosmict_2_n_muon_tracks;
  float cosmict_2_total_shower_length;
  float cosmict_2_flag_inside;
  float cosmict_2_angle_beam;
  float cosmict_2_flag_dir_weak;
  float cosmict_2_dQ_dx_end;
  float cosmict_2_dQ_dx_front;
  float cosmict_2_theta;
  float cosmict_2_phi;
  float cosmict_2_valid_tracks;
  
  // signel muon (long)
  float cosmict_3_filled;
  float cosmict_3_flag_inside;
  float cosmict_3_angle_beam;
  float cosmict_3_flag_dir_weak;
  float cosmict_3_dQ_dx_end;
  float cosmict_3_dQ_dx_front;
  float cosmict_3_theta;
  float cosmict_3_phi;
  float cosmict_3_valid_tracks;
  
  // kinematics muon
  float cosmict_4_filled;
  float cosmict_4_flag_inside;
  float cosmict_4_angle_beam;
  float cosmict_4_connected_showers;  // need to be careful about the nueCC ...
  
  // kinematics muon (long)
  float cosmict_5_filled;
  float cosmict_5_flag_inside;
  float cosmict_5_angle_beam;
  float cosmict_5_connected_showers;
  
  // special
  float cosmict_6_filled;
  float cosmict_6_flag_dir_weak;
  float cosmict_6_flag_inside;
  float cosmict_6_angle;
  
  // muon + michel
  float cosmict_7_filled;
  float cosmict_7_flag_sec;
  float cosmict_7_n_muon_tracks;
  float cosmict_7_total_shower_length;
  float cosmict_7_flag_inside;
  float cosmict_7_angle_beam;
  float cosmict_7_flag_dir_weak;
  float cosmict_7_dQ_dx_end;
  float cosmict_7_dQ_dx_front;
  float cosmict_7_theta;
  float cosmict_7_phi;
  
  // muon + michel + special
  float cosmict_8_filled;
  float cosmict_8_flag_out;
  float cosmict_8_muon_length;
  float cosmict_8_acc_length;
  
  // front upstream (dirt)
  std::vector<float> *cosmict_10_flag_inside= new std::vector<float>;
  std::vector<float> *cosmict_10_vtx_z= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_shower= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_dir_weak= new std::vector<float>;
  std::vector<float> *cosmict_10_angle_beam= new std::vector<float>;
  std::vector<float> *cosmict_10_length = new std::vector<float>;

  T->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  T->SetBranchAddress("cosmict_flag_2",&cosmict_flag_2);
  T->SetBranchAddress("cosmict_flag_3",&cosmict_flag_3);
  T->SetBranchAddress("cosmict_flag_4",&cosmict_flag_4);
  T->SetBranchAddress("cosmict_flag_5",&cosmict_flag_5);
  T->SetBranchAddress("cosmict_flag_6",&cosmict_flag_6);
  T->SetBranchAddress("cosmict_flag_7",&cosmict_flag_7);
  T->SetBranchAddress("cosmict_flag_8",&cosmict_flag_8);
  T->SetBranchAddress("cosmict_flag_9",&cosmict_flag_9);
  T->SetBranchAddress("cosmict_flag_10",&cosmict_flag_10);
  T->SetBranchAddress("cosmict_flag",&cosmict_flag);
  
  T->SetBranchAddress("cosmict_2_filled",&cosmict_2_filled);
  T->SetBranchAddress("cosmict_2_particle_type",&cosmict_2_particle_type);
  T->SetBranchAddress("cosmict_2_n_muon_tracks",&cosmict_2_n_muon_tracks);
  T->SetBranchAddress("cosmict_2_total_shower_length",&cosmict_2_total_shower_length);
  T->SetBranchAddress("cosmict_2_flag_inside",&cosmict_2_flag_inside);
  T->SetBranchAddress("cosmict_2_angle_beam",&cosmict_2_angle_beam);
  T->SetBranchAddress("cosmict_2_flag_dir_weak",&cosmict_2_flag_dir_weak);
  T->SetBranchAddress("cosmict_2_dQ_dx_end",&cosmict_2_dQ_dx_end);
  T->SetBranchAddress("cosmict_2_dQ_dx_front",&cosmict_2_dQ_dx_front);
  T->SetBranchAddress("cosmict_2_theta",&cosmict_2_theta);
  T->SetBranchAddress("cosmict_2_phi",&cosmict_2_phi);
  T->SetBranchAddress("cosmict_2_valid_tracks",&cosmict_2_valid_tracks);
  
  T->SetBranchAddress("cosmict_3_filled",&cosmict_3_filled);
  T->SetBranchAddress("cosmict_3_flag_inside",&cosmict_3_flag_inside);
  T->SetBranchAddress("cosmict_3_angle_beam",&cosmict_3_angle_beam);
  T->SetBranchAddress("cosmict_3_flag_dir_weak",&cosmict_3_flag_dir_weak);
  T->SetBranchAddress("cosmict_3_dQ_dx_end",&cosmict_3_dQ_dx_end);
  T->SetBranchAddress("cosmict_3_dQ_dx_front",&cosmict_3_dQ_dx_front);
  T->SetBranchAddress("cosmict_3_theta",&cosmict_3_theta);
  T->SetBranchAddress("cosmict_3_phi",&cosmict_3_phi);
  T->SetBranchAddress("cosmict_3_valid_tracks",&cosmict_3_valid_tracks);
  
  T->SetBranchAddress("cosmict_4_filled",&cosmict_4_filled);
  T->SetBranchAddress("cosmict_4_flag_inside",&cosmict_4_flag_inside);
  T->SetBranchAddress("cosmict_4_angle_beam",&cosmict_4_angle_beam);
  T->SetBranchAddress("cosmict_4_connected_showers",&cosmict_4_connected_showers);
  
  T->SetBranchAddress("cosmict_5_filled",&cosmict_5_filled);
  T->SetBranchAddress("cosmict_5_flag_inside",&cosmict_5_flag_inside);
  T->SetBranchAddress("cosmict_5_angle_beam",&cosmict_5_angle_beam);
  T->SetBranchAddress("cosmict_5_connected_showers",&cosmict_5_connected_showers);
  
  T->SetBranchAddress("cosmict_6_filled",&cosmict_6_filled);
  T->SetBranchAddress("cosmict_6_flag_dir_weak",&cosmict_6_flag_dir_weak);
  T->SetBranchAddress("cosmict_6_flag_inside",&cosmict_6_flag_inside);
  T->SetBranchAddress("cosmict_6_angle",&cosmict_6_angle);
  
  
  T->SetBranchAddress("cosmict_7_filled",&cosmict_7_filled);
  T->SetBranchAddress("cosmict_7_flag_sec",&cosmict_7_flag_sec);
  T->SetBranchAddress("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks);
  T->SetBranchAddress("cosmict_7_total_shower_length",&cosmict_7_total_shower_length);
  T->SetBranchAddress("cosmict_7_flag_inside",&cosmict_7_flag_inside);
  T->SetBranchAddress("cosmict_7_angle_beam",&cosmict_7_angle_beam);
  T->SetBranchAddress("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak);
  T->SetBranchAddress("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end);
  T->SetBranchAddress("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front);
  T->SetBranchAddress("cosmict_7_theta",&cosmict_7_theta);
  T->SetBranchAddress("cosmict_7_phi",&cosmict_7_phi);
  
  T->SetBranchAddress("cosmict_8_filled",&cosmict_8_filled);
  T->SetBranchAddress("cosmict_8_flag_out",&cosmict_8_flag_out);
  T->SetBranchAddress("cosmict_8_muon_length",&cosmict_8_muon_length);
  T->SetBranchAddress("cosmict_8_acc_length",&cosmict_8_acc_length);
  
  T->SetBranchAddress("cosmict_10_flag_inside",&cosmict_10_flag_inside);
  T->SetBranchAddress("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  T->SetBranchAddress("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  T->SetBranchAddress("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  T->SetBranchAddress("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  T->SetBranchAddress("cosmict_10_length",&cosmict_10_length);
  
  
  
  // numu tagger
  float numu_cc_flag;
  
  // single muon connected to neutrino vertex
  std::vector<float> *numu_cc_flag_1= new std::vector<float>;
  std::vector<float> *numu_cc_1_particle_type= new std::vector<float>;
  std::vector<float> *numu_cc_1_length= new std::vector<float>;
  std::vector<float> *numu_cc_1_medium_dQ_dx= new std::vector<float>;
  std::vector<float> *numu_cc_1_dQ_dx_cut= new std::vector<float>;
  std::vector<float> *numu_cc_1_direct_length= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_tracks= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_all= new std::vector<float>;
  
  // long muon connected to neutrino vertex
  std::vector<float> *numu_cc_flag_2= new std::vector<float>;
  std::vector<float> *numu_cc_2_length= new std::vector<float>;
  std::vector<float> *numu_cc_2_total_length= new std::vector<float>;
  std::vector<float> *numu_cc_2_n_daughter_tracks= new std::vector<float>;
  std::vector<float> *numu_cc_2_n_daughter_all = new std::vector<float>;
  
  // any muon ...
  float numu_cc_flag_3;
  float numu_cc_3_particle_type;
  float numu_cc_3_max_length;
  float numu_cc_3_acc_track_length;
  float numu_cc_3_max_length_all;
  float numu_cc_3_max_muon_length;
  float numu_cc_3_n_daughter_tracks;
  float numu_cc_3_n_daughter_all;

  T->SetBranchAddress("numu_cc_flag",&numu_cc_flag);
  
  T->SetBranchAddress("numu_cc_flag_1",&numu_cc_flag_1);
  T->SetBranchAddress("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  T->SetBranchAddress("numu_cc_1_length",&numu_cc_1_length);
  T->SetBranchAddress("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  T->SetBranchAddress("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  T->SetBranchAddress("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  T->SetBranchAddress("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  T->SetBranchAddress("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
  
  T->SetBranchAddress("numu_cc_flag_2",&numu_cc_flag_2);
  T->SetBranchAddress("numu_cc_2_length",&numu_cc_2_length);
  T->SetBranchAddress("numu_cc_2_total_length",&numu_cc_2_total_length);
  T->SetBranchAddress("numu_cc_2_n_daughter_tracks",&numu_cc_2_n_daughter_tracks);
  T->SetBranchAddress("numu_cc_2_n_daughter_all",&numu_cc_2_n_daughter_all);
  
  T->SetBranchAddress("numu_cc_flag_3",&numu_cc_flag_3);
  T->SetBranchAddress("numu_cc_3_particle_type",&numu_cc_3_particle_type);
  T->SetBranchAddress("numu_cc_3_max_length",&numu_cc_3_max_length);
  T->SetBranchAddress("numu_cc_3_track_length",&numu_cc_3_acc_track_length);
  T->SetBranchAddress("numu_cc_3_max_length_all",&numu_cc_3_max_length_all);
  T->SetBranchAddress("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length);
  T->SetBranchAddress("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks);
  T->SetBranchAddress("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all);
  
  
  T->GetEntry(0);  // total 47 of them ...
  if (!cosmic_flag) flag_info_vec.push_back(0);       // 1  checked
  if (!gap_flag) flag_info_vec.push_back(1);          // 1  checked
  if (!mip_quality_flag) flag_info_vec.push_back(2);  // 1  checked
  if (!mip_flag) flag_info_vec.push_back(3);          // 1  checked
  if (!pio_flag) flag_info_vec.push_back(4);          // 2  checked
  if (!stem_dir_flag) flag_info_vec.push_back(5);     // 1  checked
  if (!br1_flag) flag_info_vec.push_back(6);          // 3  checked
  if (!br2_flag) flag_info_vec.push_back(7);          // 1  checked
  if (!lol_flag) flag_info_vec.push_back(8);          // 3  checked 
  if (!br3_flag) flag_info_vec.push_back(9);          // 8  checked
  if (!br4_flag) flag_info_vec.push_back(10);         // 2  checked
  if (!hol_flag) flag_info_vec.push_back(11);         // 2  checked
  if (!vis_flag) flag_info_vec.push_back(12);         // 2  checked
  if (!stem_len_flag) flag_info_vec.push_back(13);    // 1  checked
  if (!brm_flag) flag_info_vec.push_back(14);         // 1  checked
  if (!cme_flag) flag_info_vec.push_back(15);         // 1  checked
  if (!anc_flag) flag_info_vec.push_back(16);         // 1  checked
  if (!lem_flag) flag_info_vec.push_back(17);         // 1  checked
  if (!stw_flag) flag_info_vec.push_back(18);         // 4  checked
  if (!spt_flag) flag_info_vec.push_back(19);         // 1  checked
  if (!mgo_flag) flag_info_vec.push_back(20);         // 1  checked
  if (!mgt_flag) flag_info_vec.push_back(21);         // 1  checked
  if (!sig_flag) flag_info_vec.push_back(22);         // 2  checked
  if (!tro_flag) flag_info_vec.push_back(23);         // 5  checked

  std::cout << flag_info_vec.size() << std::endl;

  bool cosmic_print = false; 
  bool gap_print = false;
  bool mip_quality_print = false;
  bool mip_print = false;
  bool pio_print = false;
  bool stem_dir_print = false;
  bool br1_print = false;
  bool br2_print = false;
  bool lol_print = false;
  bool br3_print = false;
  bool br4_print = false;
  bool hol_print = false;
  bool vis_print = false;
  bool stem_len_print = false;
  bool brm_print = false;
  bool cme_print = false;
  bool anc_print = false;
  bool lem_print = false;
  bool stw_print = false;
  bool spt_print = false;
  bool mgo_print = false;
  bool mgt_print = false;
  bool sig_print = false;
  bool tro_print = false;
  bool cosmict_print = false;
  bool numu_cc_print = true;
  
  // Now print out the information ...
  if (numu_cc_print){
    std::cout << "numu_cc_general: " << numu_cc_flag << std::endl;
    for (size_t i=0;i!=numu_cc_flag_1->size();i++){
      std::cout << "numu_cc_1: " << numu_cc_flag_1->at(i) << " " << numu_cc_1_particle_type->at(i) << " " << numu_cc_1_length->at(i) << " " << numu_cc_1_medium_dQ_dx->at(i) << " " << numu_cc_1_dQ_dx_cut->at(i) << " " << numu_cc_1_direct_length->at(i) << " " << numu_cc_1_n_daughter_tracks->at(i) << " " << numu_cc_1_n_daughter_all->at(i) << std::endl;
    }
    for (size_t i=0; i!= numu_cc_flag_2->size(); i++){
      std::cout << "numu_cc_2: " << numu_cc_flag_2->at(i) << " " << numu_cc_2_length->at(i) << " " << numu_cc_2_total_length->at(i) << " " << numu_cc_2_n_daughter_tracks->at(i) << " " << numu_cc_2_n_daughter_all->at(i) << std::endl;
    }

    std::cout << "numu_cc_3: " << numu_cc_flag_3 << " " << numu_cc_3_particle_type << " " << numu_cc_3_max_length << " " << numu_cc_3_acc_track_length << " " << numu_cc_3_max_length_all << " "<< numu_cc_3_max_muon_length << " " << numu_cc_3_n_daughter_tracks << " " << numu_cc_3_n_daughter_all << std::endl;
    
  }
  
  if (cosmict_print){
    std::cout << "cosmict_general: " << cosmict_flag << " " << cosmict_flag_1 << " " << cosmict_flag_9 << std::endl;
    std::cout << "cosmcit_2: " << cosmict_flag_2 << " " << cosmict_2_filled << " " << cosmict_2_particle_type << " " << cosmict_2_n_muon_tracks << " " << cosmict_2_total_shower_length << " " << cosmict_2_flag_inside << " " << cosmict_2_angle_beam << " " << cosmict_2_flag_dir_weak <<  " " << cosmict_2_dQ_dx_end << " " << cosmict_2_dQ_dx_front << " " << cosmict_2_theta << " " << cosmict_2_phi << " " << cosmict_2_valid_tracks << std::endl;
    std::cout << "cosmict_3: " << cosmict_flag_3 << " " << cosmict_3_filled << " " << cosmict_3_flag_inside << " " << cosmict_3_angle_beam << " " << cosmict_3_flag_dir_weak << " " << cosmict_3_dQ_dx_end << " " << cosmict_3_dQ_dx_front << " " << cosmict_3_theta << " " << cosmict_3_phi << " " << cosmict_3_valid_tracks << std::endl;
    std::cout << "cosmict_4: " << cosmict_flag_4 << " " << cosmict_4_filled << " " << cosmict_4_flag_inside << " " << cosmict_4_angle_beam << " " << cosmict_4_connected_showers << std::endl;
    std::cout << "cosmict_5: " << cosmict_flag_5 << " " << cosmict_5_filled << " " << cosmict_5_flag_inside << " " << cosmict_5_angle_beam << " " << cosmict_5_connected_showers << std::endl;
    std::cout << "cosmict_6: " << cosmict_flag_6 << " " << cosmict_6_filled << " " << cosmict_6_flag_dir_weak << " " << cosmict_6_flag_inside << " " << cosmict_6_angle << std::endl;
    std::cout << "cosmict_7: " << cosmict_flag_7 << " " << cosmict_7_filled << " " << cosmict_7_flag_sec << " " << cosmict_7_n_muon_tracks << " " << cosmict_7_total_shower_length << " " << cosmict_7_flag_inside << " " << cosmict_7_angle_beam << " " << cosmict_7_flag_dir_weak << " " << cosmict_7_dQ_dx_end << " " << cosmict_7_dQ_dx_front << " " << cosmict_7_theta << " " << cosmict_7_phi << std::endl;
    std::cout << "cosmict_8: " << cosmict_flag_8 << " " << cosmict_8_filled << " " << cosmict_8_flag_out << " " << cosmict_8_muon_length << " " << cosmict_8_acc_length << std::endl;
    for (size_t i=0;i!=cosmict_flag_10->size();i++){
      std::cout << "cosmict_10: " << cosmict_flag_10->at(i) << " " << cosmict_10_flag_inside->at(i) << " " << cosmict_10_vtx_z->at(i) << " " << cosmict_10_flag_shower->at(i) << " " << cosmict_10_angle_beam->at(i) << " " << cosmict_10_length->at(i) << std::endl;
    }
  }

  if (tro_print){
    for (size_t i=0;i!=tro_1_v_particle_type->size();i++){
      std::cout << "tro_1: " << tro_1_v_particle_type->at(i) << " " << tro_1_v_flag_dir_weak->at(i) << " " << tro_1_v_min_dis->at(i) << " " << tro_1_v_sg1_length->at(i) << " " << tro_1_v_shower_main_length->at(i) << " " << tro_1_v_max_n_vtx_segs->at(i) << " " << tro_1_v_tmp_length->at(i) << " " << tro_1_v_medium_dQ_dx->at(i) << " " << tro_1_v_dQ_dx_cut->at(i) << " " << tro_1_v_flag_shower_topology->at(i) << " " << tro_1_v_flag->at(i) << std::endl;
    }
    for (size_t i=0;i!=tro_2_v_energy->size();i++){
      std::cout << "tro_2: " << tro_2_v_energy->at(i) << " " << tro_2_v_stem_length->at(i) << " " << tro_2_v_iso_angle->at(i) << " " << tro_2_v_max_length->at(i) << " " << tro_2_v_angle->at(i) << " " << tro_2_v_flag->at(i) << std::endl;
    }
    std::cout << "tro_3: " << tro_3_stem_length << " " << tro_3_n_muon_segs << " " << tro_3_energy << " "<< tro_3_flag << std::endl;
    for (size_t i=0;i!=tro_4_v_dir2_mag->size();i++){
      std::cout << "tro_4: " << tro_4_v_dir2_mag->at(i) << " " << tro_4_v_angle->at(i) << " " << tro_4_v_angle1->at(i) << " " << tro_4_v_angle2->at(i) << " " << tro_4_v_length->at(i) << " " << tro_4_v_length1->at(i) << " " << tro_4_v_medium_dQ_dx->at(i) << " " << tro_4_v_end_dQ_dx->at(i) << " " << tro_4_v_energy->at(i) << " " << tro_4_v_shower_main_length->at(i) << " " << tro_4_v_flag_shower_trajectory->at(i) << " " << tro_4_v_flag->at(i) << std::endl;
    }
    for (size_t i=0;i!=tro_5_v_max_angle->size();i++){
      std::cout << "tro_5: " << tro_5_v_max_angle->at(i) << " " << tro_5_v_min_angle->at(i) << " " << tro_5_v_max_length->at(i) << " " << tro_5_v_iso_angle->at(i) << " " << tro_5_v_n_vtx_segs->at(i) << " " << tro_5_v_min_count->at(i) << " " << tro_5_v_max_count->at(i) << " " << tro_5_v_energy->at(i) << " " << tro_5_v_flag->at(i) << std::endl;
    }
    
    
  }

  
  if (sig_print){
    for (size_t i=0;i!=sig_1_v_angle->size();i++){
      std::cout << "sig_1: " << sig_1_v_angle->at(i) << " " << sig_1_v_flag_single_shower->at(i) << " " << sig_1_v_energy->at(i) << " " << sig_1_v_energy_1->at(i) << " " << sig_1_v_flag->at(i) << std::endl;
    }
    for (size_t i=0;i!=sig_2_v_energy->size();i++){
      std::cout << "sig_2: " << sig_2_v_energy->at(i) << " " << sig_2_v_shower_angle->at(i) << " " << sig_2_v_flag_single_shower->at(i) << " " << sig_2_v_medium_dQ_dx->at(i) << " " << sig_2_v_start_dQ_dx->at(i) << " " << sig_2_v_flag->at(i) << std::endl;
    }
  }
  
  if (mgt_print){
    std::cout << "mgt: " << mgt_flag_single_shower << " " << mgt_max_energy << " " << mgt_energy << " " << mgt_total_other_energy << " " << mgt_max_energy_1 << " " << mgt_e_indirect_max_energy << " " << mgt_e_direct_max_energy << " " << mgt_n_direct_showers << " " << mgt_e_direct_total_energy << " " << mgt_flag_indirect_max_pio << " " << mgt_e_indirect_total_energy << " " << mgt_flag << std::endl;
  }

  if (mgo_print){
    std::cout << "mgo: " << mgo_energy << " " << mgo_max_energy << " " << mgo_total_energy << " " << mgo_n_showers << " " << mgo_max_energy_1 << " " << mgo_max_energy_2 << " " << mgo_total_other_energy << " " << mgo_n_total_showers << " " << mgo_total_other_energy_1 << " " << mgo_flag << std::endl;
  }



  if (spt_print){
    std::cout << "spt: " << spt_flag_single_shower << " " << spt_energy << " " << spt_shower_main_length << " " << spt_shower_total_length << " " << spt_angle_beam << " " << spt_angle_vertical << " " << spt_max_dQ_dx << " " << spt_angle_beam_1 << " " << spt_angle_drift << " " << spt_angle_drift_1 << " " << spt_num_valid_tracks << " " << spt_n_vtx_segs << " " << spt_max_length << " " << spt_flag << std::endl;
  }


  if (stw_print){
    std::cout << "stw_1: " << stw_1_energy << " " << stw_1_dis << " " << stw_1_dQ_dx << " " << stw_1_flag_single_shower << " " << stw_1_n_pi0 << " " << stw_1_num_valid_tracks << " " << stw_1_flag << std::endl;
    for (size_t i=0;i!=stw_2_v_medium_dQ_dx->size();i++){
      std::cout << "stw_2: " << stw_2_v_medium_dQ_dx->at(i) << " " << stw_2_v_energy->at(i) << " " << stw_2_v_angle->at(i) << " " << stw_2_v_dir_length->at(i) << " " << stw_2_v_max_dQ_dx->at(i) << " " << stw_2_v_flag->at(i) << std::endl;
    }
    for (size_t i=0;i!=stw_3_v_angle->size();i++){
      std::cout << "stw_3: " << stw_3_v_angle->at(i) << " " << stw_3_v_dir_length->at(i) << " " << stw_3_v_energy->at(i) << " " << stw_3_v_medium_dQ_dx->at(i) << " " << stw_3_v_flag->at(i) << std::endl;
    }
    for (size_t i=0;i!=stw_4_v_angle->size();i++){
      std::cout << "stw_4: " << stw_4_v_angle->at(i) << " " << stw_4_v_dis->at(i) << " " <<stw_4_v_energy->at(i) << " " << stw_4_v_flag->at(i) << std::endl;
    }
    
    
  }

  if (lem_print){
    std::cout << "lem: " << lem_shower_total_length << " " << lem_shower_main_length << " " << lem_n_3seg << " "<< lem_e_charge << " " << lem_e_dQdx << " " << lem_shower_num_segs << " " << lem_shower_num_main_segs << " " << lem_flag << std::endl;
  }

  
  if (anc_print){
    std::cout << "anc: " << anc_energy << " " << anc_angle << " " << anc_max_angle << " " << anc_max_length << " " << anc_acc_forward_length << " " << anc_acc_backward_length << " " << anc_acc_forward_length1 << " " << anc_shower_main_length << " " << anc_shower_total_length << " " << anc_flag_main_outside << " " << anc_flag << std::endl;
  }

  if (cme_print){
    std::cout << "cme: " << cme_mu_energy << " " << cme_energy << " " << cme_mu_length << " " << cme_length << " " << cme_angle_beam << " " << cme_flag << std::endl;
  }

  if (brm_print){
    std::cout << "brm: " << brm_n_mu_segs << " " << brm_Ep << " " << brm_energy << " " << brm_acc_length << " " << brm_shower_total_length << " " << brm_connected_length << " " << brm_n_size << " " << brm_acc_direct_length << " " << brm_n_shower_main_segs << " " << brm_n_mu_main << " " << brm_flag << std::endl;
  }

  if (cosmic_print){
    std::cout << "cosmic: " << cosmic_flag << " " << cosmic_n_solid_tracks << " " << cosmic_energy_main_showers << " " << cosmic_energy_direct_showers << " " << cosmic_energy_indirect_showers << " " << cosmic_n_direct_showers << " " << cosmic_n_indirect_showers << " " << cosmic_n_main_showers << " " << cosmic_filled << std::endl;
  }

  if (gap_print){
    std::cout << "gap: " << gap_flag << " " << gap_flag_prolong_u << " " << gap_flag_prolong_v << " " << gap_flag_prolong_w << " " << gap_flag_parallel << " " << gap_n_points << " " << gap_n_bad << " " << gap_energy << " " << gap_num_valid_tracks << " " << gap_flag_single_shower << " " << gap_filled << std::endl;
  }

  if (mip_quality_print){
    std::cout << "mip_quality: " << mip_quality_flag << " " << mip_quality_energy << " " << mip_quality_overlap << " " << mip_quality_n_showers << " " << mip_quality_n_tracks << " " << mip_quality_flag_inside_pi0 << " " << mip_quality_n_pi0_showers << " " << mip_quality_shortest_length << " " << mip_quality_acc_length << " " << mip_quality_shortest_angle << " " << mip_quality_flag_proton << " " << mip_quality_filled << std::endl;
  }

  if (mip_print){
    std::cout << "mip_1: " << mip_flag << " " << mip_energy << " " << mip_n_end_reduction << " " << mip_n_first_mip << " " << mip_n_first_non_mip << " " << mip_n_first_non_mip_1 << " " << mip_n_first_non_mip_2 << " " << mip_vec_dQ_dx_0 << " " << mip_vec_dQ_dx_1 << " " << mip_max_dQ_dx_sample << std::endl;
    std::cout << "mip_2: " << mip_n_below_threshold << " " << mip_n_below_zero << " " << mip_n_lowest << " " << mip_n_highest << " " << mip_lowest_dQ_dx << " " << mip_highest_dQ_dx << " " << mip_medium_dQ_dx << " " << mip_stem_length << " " << mip_length_main << " " << mip_length_total << " " << std::endl;
    std::cout << "mip_3: " << mip_angle_beam << " " << mip_iso_angle << " " << mip_n_vertex << " " << mip_n_good_tracks << " " << mip_E_indirect_max_energy << " " << mip_flag_all_above << " " << mip_min_dQ_dx_5 << " " << mip_n_other_vertex << " " << mip_n_stem_size << " " << mip_flag_stem_trajectory << " " << mip_min_dis << " " << mip_filled << std::endl;
    std::cout << "mip_4: " << mip_vec_dQ_dx_2 << " " << mip_vec_dQ_dx_3 << " " << mip_vec_dQ_dx_4 << " " << mip_vec_dQ_dx_5 << " " << mip_vec_dQ_dx_6 << " " << mip_vec_dQ_dx_7 << " " << mip_vec_dQ_dx_8 << " " << mip_vec_dQ_dx_9 << " " << mip_vec_dQ_dx_10 << std::endl;
    std::cout << "mip_5: " <<  mip_vec_dQ_dx_11 << " " << mip_vec_dQ_dx_12 << " " << mip_vec_dQ_dx_13 << " " << mip_vec_dQ_dx_14 << " " << mip_vec_dQ_dx_15 << " " << mip_vec_dQ_dx_16 << " " << mip_vec_dQ_dx_17 << " " << mip_vec_dQ_dx_18 << " " << mip_vec_dQ_dx_19 << std::endl;
  }

  if (pio_print){
    std::cout << "pio: " << pio_flag << " " << pio_mip_id << " " << pio_filled << " " << pio_flag_pio << std::endl;
    std::cout << "pio_1: " << pio_1_flag << " " << pio_1_mass << " " << pio_1_pio_type << " " << pio_1_energy_1 << " " << pio_1_energy_2 << " " << pio_1_dis_1 << " " << pio_1_dis_2 << std::endl;
    for (size_t i=0;i!=pio_2_v_dis2->size();i++){
      std::cout << "pio_2: " << pio_2_v_dis2->at(i) << " " << pio_2_v_angle2->at(i) << " " << pio_2_v_acc_length->at(i) << " " << pio_2_v_flag->at(i) << std::endl;
    }
  }

  if (stem_dir_print){
    std::cout << "stem_dir: " << stem_dir_flag << " " << stem_dir_flag_single_shower << " " << stem_dir_filled << " " << stem_dir_angle << " " << stem_dir_energy << " " << stem_dir_angle1 << " " << stem_dir_angle2 << " " << stem_dir_angle3 << " " << stem_dir_ratio << std::endl;
  }
  
  if (br1_print){
    std::cout << "br1_1: " << br1_1_flag << " " << br1_1_shower_type << " " << br1_1_vtx_n_segs << " " << br1_1_energy << " " << br1_1_n_segs << " " << br1_1_flag_sg_topology << " " << br1_1_flag_sg_trajectory << " "<< br1_1_sg_length << std::endl;
    std::cout << "br1_2: " << br1_2_flag << " " << br1_2_energy << " " << br1_2_n_connected << " " << br1_2_max_length << " " << br1_2_n_connected_1 << " " << br1_2_vtx_n_segs << " " << br1_2_n_shower_segs << " " << br1_2_max_length_ratio << " " << br1_2_shower_length << std::endl;
    std::cout << "br1_3: " << br1_3_flag << " " << br1_3_energy << " " << br1_3_n_connected_p << " " << br1_3_max_length_p << " " << br1_3_n_shower_segs << "  " << br1_3_flag_sg_topology << " " << br1_3_flag_sg_trajectory << " " << br1_3_n_shower_main_segs << " " << br1_3_sg_length << std::endl;
  }

  if (br2_print){
    std::cout << "br2: " << br2_flag << " " << br2_flag_single_shower << " " << br2_num_valid_tracks << " " << br2_energy << " " << br2_angle1 << " " << br2_angle2 << " " << br2_angle << " " << br2_angle3 << " " << br2_n_shower_main_segs << " " << br2_max_angle << " " << br2_sg_length << " " << br2_flag_sg_trajectory << std::endl;
  }

  if (lol_print){
    for (size_t i=0; i!= lol_1_v_energy->size(); i++){
      std::cout << "lol_1: " << lol_1_v_energy->at(i) << " " << lol_1_v_vtx_n_segs->at(i) << " " << lol_1_v_nseg->at(i) << " " << lol_1_v_angle->at(i) << " " << lol_1_v_flag->at(i) << std::endl;
    }
    for (size_t i=0; i!= lol_2_v_length->size(); i++){
      std::cout << "lol_2: " << lol_2_v_length->at(i) << " " << lol_2_v_angle->at(i) << " " << lol_2_v_type->at(i) << " " << lol_2_v_vtx_n_segs->at(i) << " " << lol_2_v_energy->at(i) << " " << lol_2_v_shower_main_length->at(i) << " " << lol_2_v_flag_dir_weak->at(i) << " " << lol_2_v_flag->at(i) <<std::endl;
    }
    std::cout << "lol_3: " << lol_3_angle_beam << " " << lol_3_n_valid_tracks << " " << lol_3_min_angle << " " << lol_3_vtx_n_segs << " " << lol_3_energy << " " << lol_3_shower_main_length << " " << lol_3_n_out << " " << lol_3_n_sum << " " << lol_3_flag << std::endl;
  }

  if (br3_print){
    std::cout << "br3_1: " << br3_1_energy << " " << br3_1_n_shower_segments << " " << br3_1_sg_flag_trajectory << " " << br3_1_sg_direct_length << " " << br3_1_sg_length << " " << br3_1_total_main_length << " " << br3_1_total_length << " " << br3_1_iso_angle << " " << br3_1_sg_flag_topology << " " << br3_1_flag << std::endl;
    std::cout << "br3_2: " << br3_2_n_ele << " " << br3_2_n_other << " " << br3_2_energy << " " << br3_2_total_main_length << " " << br3_2_total_length << " " << br3_2_other_fid << " " << br3_2_flag << std::endl;
    for (size_t i=0;i!=br3_3_v_energy->size();i++){
      std::cout << "br3_3: " << br3_3_v_energy->at(i) << " " << br3_3_v_angle->at(i) << " " << br3_3_v_dir_length->at(i) << " " << br3_3_v_length->at(i) << " " << br3_3_v_flag->at(i) << std::endl;
    }
    std::cout << "br3_4: " << br3_4_acc_length << " " << br3_4_total_length << " " << br3_4_energy << " " << br3_4_flag << std::endl;
    for (size_t i=0;i!=br3_5_v_dir_length->size();i++){
      std::cout << "br3_5: " << br3_5_v_dir_length->at(i) << " " << br3_5_v_total_length->at(i) << " " << br3_5_v_flag_avoid_muon_check->at(i) << " " << br3_5_v_n_seg->at(i) << " " << br3_5_v_angle->at(i) << " " << br3_5_v_sg_length->at(i) << " " << br3_5_v_energy->at(i) << " " << br3_5_v_n_main_segs->at(i) << " " << br3_5_v_n_segs->at(i) << " " << br3_5_v_shower_main_length->at(i) << " " << br3_5_v_shower_total_length->at(i) << " " << br3_5_v_flag->at(i) << std::endl;
    }
    for (size_t i=0;i!=br3_6_v_angle->size();i++){
      std::cout << "br3_6: " << br3_6_v_angle->at(i) << " " << br3_6_v_angle1->at(i) << " " << br3_6_v_flag_shower_trajectory->at(i) << " " << br3_6_v_direct_length->at(i) << " " << br3_6_v_length->at(i) << " " << br3_6_v_n_other_vtx_segs->at(i) << " " << br3_6_v_energy->at(i) << " " << br3_6_v_flag->at(i) << std::endl;
    }
    std::cout << "br3_7: " << br3_7_energy << " " << br3_7_min_angle << " " << br3_7_sg_length << " " << br3_7_shower_main_length << " " << br3_7_flag << std::endl;
    std::cout << "br3_8: " << br3_8_max_dQ_dx << " " << br3_8_energy << " " << br3_8_n_main_segs << " " << br3_8_shower_main_length << " " << br3_8_shower_length << " " << br3_8_flag << std::endl;
  }

  if (br4_print){
    std::cout << "br4_1: " << br4_1_shower_main_length << " " << br4_1_shower_total_length << " " << br4_1_min_dis << " " << br4_1_energy << " " << br4_1_flag_avoid_muon_check << " " << br4_1_n_vtx_segs << " " << br4_1_n_main_segs << " " << br4_1_flag << std::endl;
    std::cout << "br4_2: " << br4_2_ratio_45 << " " << br4_2_ratio_35 << " " << br4_2_ratio_25 << " " << br4_2_ratio_15 << " " << br4_2_energy << " " << br4_2_ratio1_45 << " " << br4_2_ratio1_35 << " "<< br4_2_ratio1_25 << " " << br4_2_ratio1_15 << " " << br4_2_iso_angle << " " << br4_2_iso_angle1 << " " << br4_2_angle << " " << br4_2_flag << std::endl;
  }

  if (hol_print){
    std::cout << "hol_1: " << hol_1_n_valid_tracks << " " << hol_1_min_angle << " " << hol_1_energy << " " << hol_1_flag_all_shower << " " << hol_1_min_length << " " << hol_1_flag << std::endl;
    std::cout << "hol_2: " << hol_2_min_angle << " " << hol_2_medium_dQ_dx << " " << hol_2_ncount << " " << hol_2_energy << " " << hol_2_flag << std::endl;
  }

  if (vis_print){
    std::cout << "vis_1: " << vis_1_filled << " " << vis_1_n_vtx_segs << " " << vis_1_energy << " " << vis_1_num_good_tracks << " " << vis_1_max_angle << " " << vis_1_max_shower_angle << " " << vis_1_tmp_length1 << " " << vis_1_tmp_length2 << " " << vis_1_particle_type << " " << vis_1_flag << std::endl;
    std::cout << "vis_2: " << vis_2_filled << " " << vis_2_n_vtx_segs << " " << vis_2_min_angle << " " << vis_2_min_weak_track << " " << vis_2_angle_beam << " " << vis_2_min_angle1 << " " << vis_2_iso_angle1 << " " << vis_2_min_medium_dQ_dx << " " << vis_2_min_length << " " << vis_2_sg_length << " " << vis_2_max_angle << " " << vis_2_max_weak_track << " " << vis_2_flag << std::endl;
  }

  if (stem_len_print){
    std::cout << "stem_len: " << stem_len_energy << " " << stem_len_length << " " << stem_len_flag_avoid_muon_check << " " << stem_len_num_daughters << " " << stem_len_daughter_length << " " << stem_len_flag << std::endl;
  }

  
}
