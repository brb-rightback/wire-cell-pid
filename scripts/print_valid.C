#include <vector>

void print_valid(TString run_info = "7019_5_264"){
  TChain *T = new TChain("T_tagger","T_tagger");
  TString filename = "nue_" + run_info + ".root";
  T->AddFile(filename);

  std::vector<int> flag_info_vec;
  
  // cosmic information
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
  
  // Now print out the information ...
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
