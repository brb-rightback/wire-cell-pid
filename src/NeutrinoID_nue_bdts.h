#include "TMVA/Reader.h"

float WCPPID::NeutrinoID::cal_bdts_xgboost(){
  double val = 0; // background like ...
  bool flag_print = true;

  // vector scores (15 of them)
  tagger_info.br3_3_score     = cal_br3_3_bdt(0.3);
  tagger_info.br3_5_score     = cal_br3_5_bdt(0.42);
  tagger_info.br3_6_score     = cal_br3_6_bdt(0.75);
  tagger_info.pio_2_score     = cal_pio_2_bdt(0.2);
  tagger_info.stw_2_score     = cal_stw_2_bdt(0.7);
  tagger_info.stw_3_score     = cal_stw_3_bdt(0.5);
  tagger_info.stw_4_score     = cal_stw_4_bdt(0.7);
  tagger_info.sig_1_score     = cal_sig_1_bdt(0.59);
  tagger_info.sig_2_score     = cal_sig_2_bdt(0.55);
  tagger_info.lol_1_score     = cal_lol_1_bdt(0.85);
  tagger_info.lol_2_score     = cal_lol_2_bdt(0.7);
  tagger_info.tro_1_score     = cal_tro_1_bdt(0.28);
  tagger_info.tro_2_score     = cal_tro_2_bdt(0.35);
  tagger_info.tro_4_score     = cal_tro_4_bdt(0.33);
  tagger_info.tro_5_score     = cal_tro_5_bdt(0.5);

  float default_val = -15;// background like

  TMVA::Reader reader;

  reader.AddVariable("cme_mu_energy",&tagger_info.cme_mu_energy);
  reader.AddVariable("cme_energy",&tagger_info.cme_energy);
  reader.AddVariable("cme_mu_length",&tagger_info.cme_mu_length);
  reader.AddVariable("cme_length",&tagger_info.cme_length);
  reader.AddVariable("cme_angle_beam",&tagger_info.cme_angle_beam);
  reader.AddVariable("anc_angle",&tagger_info.anc_angle);
  reader.AddVariable("anc_max_angle",&tagger_info.anc_max_angle);
  reader.AddVariable("anc_max_length",&tagger_info.anc_max_length);
  reader.AddVariable("anc_acc_forward_length",&tagger_info.anc_acc_forward_length);
  reader.AddVariable("anc_acc_backward_length",&tagger_info.anc_acc_backward_length);
  reader.AddVariable("anc_acc_forward_length1",&tagger_info.anc_acc_forward_length1);
  reader.AddVariable("anc_shower_main_length",&tagger_info.anc_shower_main_length);
  reader.AddVariable("anc_shower_total_length",&tagger_info.anc_shower_total_length);
  reader.AddVariable("anc_flag_main_outside",&tagger_info.anc_flag_main_outside);
  reader.AddVariable("gap_flag_prolong_u",&tagger_info.gap_flag_prolong_u);
  reader.AddVariable("gap_flag_prolong_v",&tagger_info.gap_flag_prolong_v);
  reader.AddVariable("gap_flag_prolong_w",&tagger_info.gap_flag_prolong_w);
  reader.AddVariable("gap_flag_parallel",&tagger_info.gap_flag_parallel);
  reader.AddVariable("gap_n_points",&tagger_info.gap_n_points);
  reader.AddVariable("gap_n_bad",&tagger_info.gap_n_bad);
  reader.AddVariable("gap_energy",&tagger_info.gap_energy);
  reader.AddVariable("gap_num_valid_tracks",&tagger_info.gap_num_valid_tracks);
  reader.AddVariable("gap_flag_single_shower",&tagger_info.gap_flag_single_shower);
  reader.AddVariable("hol_1_n_valid_tracks",&tagger_info.hol_1_n_valid_tracks);
  reader.AddVariable("hol_1_min_angle",&tagger_info.hol_1_min_angle);
  reader.AddVariable("hol_1_energy",&tagger_info.hol_1_energy);
  reader.AddVariable("hol_1_min_length",&tagger_info.hol_1_min_length);
  reader.AddVariable("hol_2_min_angle",&tagger_info.hol_2_min_angle);
  reader.AddVariable("hol_2_medium_dQ_dx",&tagger_info.hol_2_medium_dQ_dx);
  reader.AddVariable("hol_2_ncount",&tagger_info.hol_2_ncount);
  reader.AddVariable("lol_3_angle_beam",&tagger_info.lol_3_angle_beam);
  reader.AddVariable("lol_3_n_valid_tracks",&tagger_info.lol_3_n_valid_tracks);
  reader.AddVariable("lol_3_min_angle",&tagger_info.lol_3_min_angle);
  reader.AddVariable("lol_3_vtx_n_segs",&tagger_info.lol_3_vtx_n_segs);
  reader.AddVariable("lol_3_shower_main_length",&tagger_info.lol_3_shower_main_length);
  reader.AddVariable("lol_3_n_out",&tagger_info.lol_3_n_out);
  reader.AddVariable("lol_3_n_sum",&tagger_info.lol_3_n_sum);
  reader.AddVariable("hol_1_flag_all_shower",&tagger_info.hol_1_flag_all_shower); // naming issue
  reader.AddVariable("mgo_energy",&tagger_info.mgo_energy);
  reader.AddVariable("mgo_max_energy",&tagger_info.mgo_max_energy);
  reader.AddVariable("mgo_total_energy",&tagger_info.mgo_total_energy);
  reader.AddVariable("mgo_n_showers",&tagger_info.mgo_n_showers);
  reader.AddVariable("mgo_max_energy_1",&tagger_info.mgo_max_energy_1);
  reader.AddVariable("mgo_max_energy_2",&tagger_info.mgo_max_energy_2);
  reader.AddVariable("mgo_total_other_energy",&tagger_info.mgo_total_other_energy);
  reader.AddVariable("mgo_n_total_showers",&tagger_info.mgo_n_total_showers);
  reader.AddVariable("mgo_total_other_energy_1",&tagger_info.mgo_total_other_energy_1);
  reader.AddVariable("mgt_flag_single_shower",&tagger_info.mgt_flag_single_shower);
  reader.AddVariable("mgt_max_energy",&tagger_info.mgt_max_energy);
  reader.AddVariable("mgt_total_other_energy",&tagger_info.mgt_total_other_energy);
  reader.AddVariable("mgt_max_energy_1",&tagger_info.mgt_max_energy_1);
  reader.AddVariable("mgt_e_indirect_max_energy",&tagger_info.mgt_e_indirect_max_energy);
  reader.AddVariable("mgt_e_direct_max_energy",&tagger_info.mgt_e_direct_max_energy);
  reader.AddVariable("mgt_n_direct_showers",&tagger_info.mgt_n_direct_showers);
  reader.AddVariable("mgt_e_direct_total_energy",&tagger_info.mgt_e_direct_total_energy);
  reader.AddVariable("mgt_flag_indirect_max_pio",&tagger_info.mgt_flag_indirect_max_pio);
  reader.AddVariable("mgt_e_indirect_total_energy",&tagger_info.mgt_e_indirect_total_energy);
  reader.AddVariable("mip_quality_energy",&tagger_info.mip_quality_energy);
  reader.AddVariable("mip_quality_overlap",&tagger_info.mip_quality_overlap);
  reader.AddVariable("mip_quality_n_showers",&tagger_info.mip_quality_n_showers);
  reader.AddVariable("mip_quality_n_tracks",&tagger_info.mip_quality_n_tracks);
  reader.AddVariable("mip_quality_flag_inside_pi0",&tagger_info.mip_quality_flag_inside_pi0);
  reader.AddVariable("mip_quality_n_pi0_showers",&tagger_info.mip_quality_n_pi0_showers);
  reader.AddVariable("mip_quality_shortest_length",&tagger_info.mip_quality_shortest_length);
  reader.AddVariable("mip_quality_acc_length",&tagger_info.mip_quality_acc_length);
  reader.AddVariable("mip_quality_shortest_angle",&tagger_info.mip_quality_shortest_angle);
  reader.AddVariable("mip_quality_flag_proton",&tagger_info.mip_quality_flag_proton);
  reader.AddVariable("br1_1_shower_type",&tagger_info.br1_1_shower_type);
  reader.AddVariable("br1_1_vtx_n_segs",&tagger_info.br1_1_vtx_n_segs);
  reader.AddVariable("br1_1_energy",&tagger_info.br1_1_energy);
  reader.AddVariable("br1_1_n_segs",&tagger_info.br1_1_n_segs);
  reader.AddVariable("br1_1_flag_sg_topology",&tagger_info.br1_1_flag_sg_topology);
  reader.AddVariable("br1_1_flag_sg_trajectory",&tagger_info.br1_1_flag_sg_trajectory);
  reader.AddVariable("br1_1_sg_length",&tagger_info.br1_1_sg_length);
  reader.AddVariable("br1_2_n_connected",&tagger_info.br1_2_n_connected);
  reader.AddVariable("br1_2_max_length",&tagger_info.br1_2_max_length);
  reader.AddVariable("br1_2_n_connected_1",&tagger_info.br1_2_n_connected_1);
  reader.AddVariable("br1_2_n_shower_segs",&tagger_info.br1_2_n_shower_segs);
  reader.AddVariable("br1_2_max_length_ratio",&tagger_info.br1_2_max_length_ratio);
  reader.AddVariable("br1_2_shower_length",&tagger_info.br1_2_shower_length);
  reader.AddVariable("br1_3_n_connected_p",&tagger_info.br1_3_n_connected_p);
  reader.AddVariable("br1_3_max_length_p",&tagger_info.br1_3_max_length_p);
  reader.AddVariable("br1_3_n_shower_main_segs",&tagger_info.br1_3_n_shower_main_segs);
  reader.AddVariable("br3_1_energy",&tagger_info.br3_1_energy);
  reader.AddVariable("br3_1_n_shower_segments",&tagger_info.br3_1_n_shower_segments);
  reader.AddVariable("br3_1_sg_flag_trajectory",&tagger_info.br3_1_sg_flag_trajectory);
  reader.AddVariable("br3_1_sg_direct_length",&tagger_info.br3_1_sg_direct_length);
  reader.AddVariable("br3_1_sg_length",&tagger_info.br3_1_sg_length);
  reader.AddVariable("br3_1_total_main_length",&tagger_info.br3_1_total_main_length);
  reader.AddVariable("br3_1_total_length",&tagger_info.br3_1_total_length);
  reader.AddVariable("br3_1_iso_angle",&tagger_info.br3_1_iso_angle);
  reader.AddVariable("br3_1_sg_flag_topology",&tagger_info.br3_1_sg_flag_topology);
  reader.AddVariable("br3_2_n_ele",&tagger_info.br3_2_n_ele);
  reader.AddVariable("br3_2_n_other",&tagger_info.br3_2_n_other);
  reader.AddVariable("br3_2_other_fid",&tagger_info.br3_2_other_fid);
  reader.AddVariable("br3_4_acc_length",&tagger_info.br3_4_acc_length);
  reader.AddVariable("br3_4_total_length",&tagger_info.br3_4_total_length);
  reader.AddVariable("br3_7_min_angle",&tagger_info.br3_7_min_angle);
  reader.AddVariable("br3_8_max_dQ_dx",&tagger_info.br3_8_max_dQ_dx);
  reader.AddVariable("br3_8_n_main_segs",&tagger_info.br3_8_n_main_segs);
  reader.AddVariable("vis_1_n_vtx_segs",&tagger_info.vis_1_n_vtx_segs);
  reader.AddVariable("vis_1_energy",&tagger_info.vis_1_energy);
  reader.AddVariable("vis_1_num_good_tracks",&tagger_info.vis_1_num_good_tracks);
  reader.AddVariable("vis_1_max_angle",&tagger_info.vis_1_max_angle);
  reader.AddVariable("vis_1_max_shower_angle",&tagger_info.vis_1_max_shower_angle);
  reader.AddVariable("vis_1_tmp_length1",&tagger_info.vis_1_tmp_length1);
  reader.AddVariable("vis_1_tmp_length2",&tagger_info.vis_1_tmp_length2);
  reader.AddVariable("vis_2_n_vtx_segs",&tagger_info.vis_2_n_vtx_segs);
  reader.AddVariable("vis_2_min_angle",&tagger_info.vis_2_min_angle);
  reader.AddVariable("vis_2_min_weak_track",&tagger_info.vis_2_min_weak_track);
  reader.AddVariable("vis_2_angle_beam",&tagger_info.vis_2_angle_beam);
  reader.AddVariable("vis_2_min_angle1",&tagger_info.vis_2_min_angle1);
  reader.AddVariable("vis_2_iso_angle1",&tagger_info.vis_2_iso_angle1);
  reader.AddVariable("vis_2_min_medium_dQ_dx",&tagger_info.vis_2_min_medium_dQ_dx);
  reader.AddVariable("vis_2_min_length",&tagger_info.vis_2_min_length);
  reader.AddVariable("vis_2_sg_length",&tagger_info.vis_2_sg_length);
  reader.AddVariable("vis_2_max_angle",&tagger_info.vis_2_max_angle);
  reader.AddVariable("vis_2_max_weak_track",&tagger_info.vis_2_max_weak_track);
  reader.AddVariable("pio_1_mass",&tagger_info.pio_1_mass);
  reader.AddVariable("pio_1_pio_type",&tagger_info.pio_1_pio_type);
  reader.AddVariable("pio_1_energy_1",&tagger_info.pio_1_energy_1);
  reader.AddVariable("pio_1_energy_2",&tagger_info.pio_1_energy_2);
  reader.AddVariable("pio_1_dis_1",&tagger_info.pio_1_dis_1);
  reader.AddVariable("pio_1_dis_2",&tagger_info.pio_1_dis_2);
  reader.AddVariable("pio_mip_id",&tagger_info.pio_mip_id);
  reader.AddVariable("stem_dir_flag_single_shower",&tagger_info.stem_dir_flag_single_shower);
  reader.AddVariable("stem_dir_angle",&tagger_info.stem_dir_angle);
  reader.AddVariable("stem_dir_energy",&tagger_info.stem_dir_energy);
  reader.AddVariable("stem_dir_angle1",&tagger_info.stem_dir_angle1);
  reader.AddVariable("stem_dir_angle2",&tagger_info.stem_dir_angle2);
  reader.AddVariable("stem_dir_angle3",&tagger_info.stem_dir_angle3);
  reader.AddVariable("stem_dir_ratio",&tagger_info.stem_dir_ratio);
  reader.AddVariable("br2_num_valid_tracks",&tagger_info.br2_num_valid_tracks);
  reader.AddVariable("br2_n_shower_main_segs",&tagger_info.br2_n_shower_main_segs);
  reader.AddVariable("br2_max_angle",&tagger_info.br2_max_angle);
  reader.AddVariable("br2_sg_length",&tagger_info.br2_sg_length);
  reader.AddVariable("br2_flag_sg_trajectory",&tagger_info.br2_flag_sg_trajectory);
  reader.AddVariable("stem_len_energy",&tagger_info.stem_len_energy);
  reader.AddVariable("stem_len_length",&tagger_info.stem_len_length);
  reader.AddVariable("stem_len_flag_avoid_muon_check",&tagger_info.stem_len_flag_avoid_muon_check);
  reader.AddVariable("stem_len_num_daughters",&tagger_info.stem_len_num_daughters);
  reader.AddVariable("stem_len_daughter_length",&tagger_info.stem_len_daughter_length);
  reader.AddVariable("brm_n_mu_segs",&tagger_info.brm_n_mu_segs);
  reader.AddVariable("brm_Ep",&tagger_info.brm_Ep);
  reader.AddVariable("brm_acc_length",&tagger_info.brm_acc_length);
  reader.AddVariable("brm_shower_total_length",&tagger_info.brm_shower_total_length);
  reader.AddVariable("brm_connected_length",&tagger_info.brm_connected_length);
  reader.AddVariable("brm_n_size",&tagger_info.brm_n_size);
  reader.AddVariable("brm_n_shower_main_segs",&tagger_info.brm_n_shower_main_segs);
  reader.AddVariable("brm_n_mu_main",&tagger_info.brm_n_mu_main);
  reader.AddVariable("lem_shower_main_length",&tagger_info.lem_shower_main_length);
  reader.AddVariable("lem_n_3seg",&tagger_info.lem_n_3seg);
  reader.AddVariable("lem_e_charge",&tagger_info.lem_e_charge);
  reader.AddVariable("lem_e_dQdx",&tagger_info.lem_e_dQdx);
  reader.AddVariable("lem_shower_num_main_segs",&tagger_info.lem_shower_num_main_segs);
  reader.AddVariable("brm_acc_direct_length",&tagger_info.brm_acc_direct_length); // naming issue
  reader.AddVariable("stw_1_energy",&tagger_info.stw_1_energy);
  reader.AddVariable("stw_1_dis",&tagger_info.stw_1_dis);
  reader.AddVariable("stw_1_dQ_dx",&tagger_info.stw_1_dQ_dx);
  reader.AddVariable("stw_1_flag_single_shower",&tagger_info.stw_1_flag_single_shower);
  reader.AddVariable("stw_1_n_pi0",&tagger_info.stw_1_n_pi0);
  reader.AddVariable("stw_1_num_valid_tracks",&tagger_info.stw_1_num_valid_tracks);
  reader.AddVariable("spt_shower_main_length",&tagger_info.spt_shower_main_length);
  reader.AddVariable("spt_shower_total_length",&tagger_info.spt_shower_total_length);
  reader.AddVariable("spt_angle_beam",&tagger_info.spt_angle_beam);
  reader.AddVariable("spt_angle_vertical",&tagger_info.spt_angle_vertical);
  reader.AddVariable("spt_max_dQ_dx",&tagger_info.spt_max_dQ_dx);
  reader.AddVariable("spt_angle_beam_1",&tagger_info.spt_angle_beam_1);
  reader.AddVariable("spt_angle_drift",&tagger_info.spt_angle_drift);
  reader.AddVariable("spt_angle_drift_1",&tagger_info.spt_angle_drift_1);
  reader.AddVariable("spt_num_valid_tracks",&tagger_info.spt_num_valid_tracks);
  reader.AddVariable("spt_n_vtx_segs",&tagger_info.spt_n_vtx_segs);
  reader.AddVariable("spt_max_length",&tagger_info.spt_max_length);
  reader.AddVariable("mip_energy",&tagger_info.mip_energy);
  reader.AddVariable("mip_n_end_reduction",&tagger_info.mip_n_end_reduction);
  reader.AddVariable("mip_n_first_mip",&tagger_info.mip_n_first_mip);
  reader.AddVariable("mip_n_first_non_mip",&tagger_info.mip_n_first_non_mip);
  reader.AddVariable("mip_n_first_non_mip_1",&tagger_info.mip_n_first_non_mip_1);
  reader.AddVariable("mip_n_first_non_mip_2",&tagger_info.mip_n_first_non_mip_2);
  reader.AddVariable("mip_vec_dQ_dx_0",&tagger_info.mip_vec_dQ_dx_0);
  reader.AddVariable("mip_vec_dQ_dx_1",&tagger_info.mip_vec_dQ_dx_1);
  reader.AddVariable("mip_max_dQ_dx_sample",&tagger_info.mip_max_dQ_dx_sample);
  reader.AddVariable("mip_n_below_threshold",&tagger_info.mip_n_below_threshold);
  reader.AddVariable("mip_n_below_zero",&tagger_info.mip_n_below_zero);
  reader.AddVariable("mip_n_lowest",&tagger_info.mip_n_lowest);
  reader.AddVariable("mip_n_highest",&tagger_info.mip_n_highest);
  reader.AddVariable("mip_lowest_dQ_dx",&tagger_info.mip_lowest_dQ_dx);
  reader.AddVariable("mip_highest_dQ_dx",&tagger_info.mip_highest_dQ_dx);
  reader.AddVariable("mip_medium_dQ_dx",&tagger_info.mip_medium_dQ_dx);
  reader.AddVariable("mip_stem_length",&tagger_info.mip_stem_length);
  reader.AddVariable("mip_length_main",&tagger_info.mip_length_main);
  reader.AddVariable("mip_length_total",&tagger_info.mip_length_total);
  reader.AddVariable("mip_angle_beam",&tagger_info.mip_angle_beam);
  reader.AddVariable("mip_iso_angle",&tagger_info.mip_iso_angle);
  reader.AddVariable("mip_n_vertex",&tagger_info.mip_n_vertex);
  reader.AddVariable("mip_n_good_tracks",&tagger_info.mip_n_good_tracks);
  reader.AddVariable("mip_E_indirect_max_energy",&tagger_info.mip_E_indirect_max_energy);
  reader.AddVariable("mip_flag_all_above",&tagger_info.mip_flag_all_above);
  reader.AddVariable("mip_min_dQ_dx_5",&tagger_info.mip_min_dQ_dx_5);
  reader.AddVariable("mip_n_other_vertex",&tagger_info.mip_n_other_vertex);
  reader.AddVariable("mip_n_stem_size",&tagger_info.mip_n_stem_size);
  reader.AddVariable("mip_flag_stem_trajectory",&tagger_info.mip_flag_stem_trajectory);
  reader.AddVariable("mip_min_dis",&tagger_info.mip_min_dis);
  reader.AddVariable("mip_vec_dQ_dx_2",&tagger_info.mip_vec_dQ_dx_2);
  reader.AddVariable("mip_vec_dQ_dx_3",&tagger_info.mip_vec_dQ_dx_3);
  reader.AddVariable("mip_vec_dQ_dx_4",&tagger_info.mip_vec_dQ_dx_4);
  reader.AddVariable("mip_vec_dQ_dx_5",&tagger_info.mip_vec_dQ_dx_5);
  reader.AddVariable("mip_vec_dQ_dx_6",&tagger_info.mip_vec_dQ_dx_6);
  reader.AddVariable("mip_vec_dQ_dx_7",&tagger_info.mip_vec_dQ_dx_7);
  reader.AddVariable("mip_vec_dQ_dx_8",&tagger_info.mip_vec_dQ_dx_8);
  reader.AddVariable("mip_vec_dQ_dx_9",&tagger_info.mip_vec_dQ_dx_9);
  reader.AddVariable("mip_vec_dQ_dx_10",&tagger_info.mip_vec_dQ_dx_10);
  reader.AddVariable("mip_vec_dQ_dx_11",&tagger_info.mip_vec_dQ_dx_11);
  reader.AddVariable("mip_vec_dQ_dx_12",&tagger_info.mip_vec_dQ_dx_12);
  reader.AddVariable("mip_vec_dQ_dx_13",&tagger_info.mip_vec_dQ_dx_13);
  reader.AddVariable("mip_vec_dQ_dx_14",&tagger_info.mip_vec_dQ_dx_14);
  reader.AddVariable("mip_vec_dQ_dx_15",&tagger_info.mip_vec_dQ_dx_15);
  reader.AddVariable("mip_vec_dQ_dx_16",&tagger_info.mip_vec_dQ_dx_16);
  reader.AddVariable("mip_vec_dQ_dx_17",&tagger_info.mip_vec_dQ_dx_17);
  reader.AddVariable("mip_vec_dQ_dx_18",&tagger_info.mip_vec_dQ_dx_18);
  reader.AddVariable("mip_vec_dQ_dx_19",&tagger_info.mip_vec_dQ_dx_19);
  reader.AddVariable("br3_3_score",&tagger_info.br3_3_score);
  reader.AddVariable("br3_5_score",&tagger_info.br3_5_score);
  reader.AddVariable("br3_6_score",&tagger_info.br3_6_score);
  reader.AddVariable("pio_2_score",&tagger_info.pio_2_score);
  reader.AddVariable("stw_2_score",&tagger_info.stw_2_score);
  reader.AddVariable("stw_3_score",&tagger_info.stw_3_score);
  reader.AddVariable("stw_4_score",&tagger_info.stw_4_score);
  reader.AddVariable("sig_1_score",&tagger_info.sig_1_score);
  reader.AddVariable("sig_2_score",&tagger_info.sig_2_score);
  reader.AddVariable("lol_1_score",&tagger_info.lol_1_score);
  reader.AddVariable("lol_2_score",&tagger_info.lol_2_score);
  reader.AddVariable("tro_1_score",&tagger_info.tro_1_score);
  reader.AddVariable("tro_2_score",&tagger_info.tro_2_score);
  reader.AddVariable("tro_4_score",&tagger_info.tro_4_score);
  reader.AddVariable("tro_5_score",&tagger_info.tro_5_score);
  reader.AddVariable("br4_1_shower_main_length",&tagger_info.br4_1_shower_main_length);
  reader.AddVariable("br4_1_shower_total_length",&tagger_info.br4_1_shower_total_length);
  reader.AddVariable("br4_1_min_dis",&tagger_info.br4_1_min_dis);
  reader.AddVariable("br4_1_energy",&tagger_info.br4_1_energy);
  reader.AddVariable("br4_1_flag_avoid_muon_check",&tagger_info.br4_1_flag_avoid_muon_check);
  reader.AddVariable("br4_1_n_vtx_segs",&tagger_info.br4_1_n_vtx_segs);
  reader.AddVariable("br4_2_ratio_45",&tagger_info.br4_2_ratio_45);
  reader.AddVariable("br4_2_ratio_35",&tagger_info.br4_2_ratio_35);
  reader.AddVariable("br4_2_ratio_25",&tagger_info.br4_2_ratio_25);
  reader.AddVariable("br4_2_ratio_15",&tagger_info.br4_2_ratio_15);
  reader.AddVariable("br4_2_ratio1_45",&tagger_info.br4_2_ratio1_45);
  reader.AddVariable("br4_2_ratio1_35",&tagger_info.br4_2_ratio1_35);
  reader.AddVariable("br4_2_ratio1_25",&tagger_info.br4_2_ratio1_25);
  reader.AddVariable("br4_2_ratio1_15",&tagger_info.br4_2_ratio1_15);
  reader.AddVariable("br4_2_iso_angle",&tagger_info.br4_2_iso_angle);
  reader.AddVariable("br4_2_iso_angle1",&tagger_info.br4_2_iso_angle1);
  reader.AddVariable("br4_2_angle",&tagger_info.br4_2_angle);
  reader.AddVariable("tro_3_stem_length",&tagger_info.tro_3_stem_length);
  reader.AddVariable("tro_3_n_muon_segs",&tagger_info.tro_3_n_muon_segs);
  reader.AddVariable("br4_1_n_main_segs",&tagger_info.br4_1_n_main_segs); // naming issue
  
  reader.BookMVA( "MyBDT", "./input_data_files/xgboost_set8_kaicheng.xml");
  

  if (tagger_info.br_filled==1){
    // protection of variables ...
    if(tagger_info.mip_min_dis>1000) tagger_info.mip_min_dis = 1000.0;
    if(tagger_info.mip_quality_shortest_length>1000) tagger_info.mip_quality_shortest_length = 1000;
    if(std::isnan(tagger_info.mip_quality_shortest_angle)) tagger_info.mip_quality_shortest_angle = 0;
    if(std::isnan(tagger_info.stem_dir_ratio)) tagger_info.stem_dir_ratio = 1.0; 
    
    
    double val1 = reader.EvaluateMVA("MyBDT");
    val = TMath::Log10( (1+val1)/(1-val1) );
    
    if (flag_print){
      std::cout << "BDT_final: " << val << std::endl;
      std::cout << "BDT: " << tagger_info.tro_1_score << " " << tagger_info.tro_2_score << " " << tagger_info.tro_4_score << " " << tagger_info.tro_5_score << " " << tagger_info.pio_2_score << std::endl;
      std::cout << "BDT: " << tagger_info.lol_1_score << " " << tagger_info.lol_2_score << " " << tagger_info.sig_1_score << " " << tagger_info.sig_2_score << " " << tagger_info.br3_3_score  << std::endl;
      std::cout << "BDT: " <<  tagger_info.br3_5_score << " " << tagger_info.br3_6_score << " " <<  tagger_info.stw_2_score << " " << tagger_info.stw_3_score << " " << tagger_info.stw_4_score << std::endl;
    }
    
  }else{
    val = default_val;
  }

  
  
  return val;
}

float WCPPID::NeutrinoID::cal_bdts(){
  double val =0 ;

  bool flag_print = true;
  
  // PIDs ...
  tagger_info.mipid_score     = cal_mipid_bdt(0.26);
  tagger_info.gap_score       = cal_gap_bdt(0.32); 
  tagger_info.hol_lol_score   = cal_hol_lol_bdt();
  tagger_info.cme_anc_score   = cal_cme_anc_bdt();
  tagger_info.mgo_mgt_score   = cal_mgo_mgt_bdt();
  tagger_info.br1_score       = cal_br1_bdt();
  tagger_info.br3_score       = cal_br3_bdt();
  tagger_info.br3_3_score     = cal_br3_3_bdt(0.3);
  tagger_info.br3_5_score     = cal_br3_5_bdt(0.42);
  tagger_info.br3_6_score     = cal_br3_6_bdt(0.75);
  tagger_info.stemdir_br2_score=cal_stemdir_br2_bdt();
  tagger_info.trimuon_score   = cal_trimuon_bdt();
  tagger_info.br4_tro_score   = cal_br4_tro_bdt();
  tagger_info.mipquality_score= cal_mipquality_bdt();
  tagger_info.pio_1_score     = cal_pio_1_bdt(0.1);
  tagger_info.pio_2_score     = cal_pio_2_bdt(0.2);
  tagger_info.stw_spt_score   = cal_stw_spt_bdt();
  tagger_info.vis_1_score     = cal_vis_1_bdt(0.6);
  tagger_info.vis_2_score     = cal_vis_2_bdt(0.52);
  tagger_info.stw_2_score     = cal_stw_2_bdt(0.7);
  tagger_info.stw_3_score     = cal_stw_3_bdt(0.5);
  tagger_info.stw_4_score     = cal_stw_4_bdt(0.7);
  tagger_info.sig_1_score     = cal_sig_1_bdt(0.59);
  tagger_info.sig_2_score     = cal_sig_2_bdt(0.55);
  tagger_info.lol_1_score     = cal_lol_1_bdt(0.85);
  tagger_info.lol_2_score     = cal_lol_2_bdt(0.7);
  tagger_info.tro_1_score     = cal_tro_1_bdt(0.28);
  tagger_info.tro_2_score     = cal_tro_2_bdt(0.35);
  tagger_info.tro_4_score     = cal_tro_4_bdt(0.33);
  tagger_info.tro_5_score     = cal_tro_5_bdt(0.5);

  float default_val = -1;
  
  TMVA::Reader reader;
  
  reader.AddVariable("mip_energy",&tagger_info.mip_energy);
  reader.AddVariable("mip_vec_dQ_dx_0",&tagger_info.mip_vec_dQ_dx_0);
  reader.AddVariable("mip_vec_dQ_dx_1",&tagger_info.mip_vec_dQ_dx_1);
  reader.AddVariable("mip_vec_dQ_dx_2",&tagger_info.mip_vec_dQ_dx_2);
  reader.AddVariable("mip_vec_dQ_dx_3",&tagger_info.mip_vec_dQ_dx_3);
  reader.AddVariable("mip_vec_dQ_dx_4",&tagger_info.mip_vec_dQ_dx_4);
  reader.AddVariable("mipid_score",&tagger_info.mipid_score);
  reader.AddVariable("gap_score",&tagger_info.gap_score);
  reader.AddVariable("hol_lol_score",&tagger_info.hol_lol_score);
  reader.AddVariable("cme_anc_score",&tagger_info.cme_anc_score);
  reader.AddVariable("mgo_mgt_score",&tagger_info.mgo_mgt_score);
  reader.AddVariable("br1_score",&tagger_info.br1_score);
  reader.AddVariable("br3_score",&tagger_info.br3_score);
  reader.AddVariable("br3_3_score",&tagger_info.br3_3_score);
  reader.AddVariable("br3_5_score",&tagger_info.br3_5_score);
  reader.AddVariable("br3_6_score",&tagger_info.br3_6_score);
  reader.AddVariable("stemdir_br2_score",&tagger_info.stemdir_br2_score);
  reader.AddVariable("trimuon_score",&tagger_info.trimuon_score);
  reader.AddVariable("br4_tro_score",&tagger_info.br4_tro_score);
  reader.AddVariable("mipquality_score",&tagger_info.mipquality_score);
  reader.AddVariable("pio_1_score",&tagger_info.pio_1_score);
  reader.AddVariable("pio_2_score",&tagger_info.pio_2_score);
  reader.AddVariable("stw_spt_score",&tagger_info.stw_spt_score);
  reader.AddVariable("vis_1_score",&tagger_info.vis_1_score);
  reader.AddVariable("vis_2_score",&tagger_info.vis_2_score);
  reader.AddVariable("stw_2_score",&tagger_info.stw_2_score);
  reader.AddVariable("stw_3_score",&tagger_info.stw_3_score);
  reader.AddVariable("stw_4_score",&tagger_info.stw_4_score);
  reader.AddVariable("sig_1_score",&tagger_info.sig_1_score);
  reader.AddVariable("sig_2_score",&tagger_info.sig_2_score);
  reader.AddVariable("lol_1_score",&tagger_info.lol_1_score);
  reader.AddVariable("lol_2_score",&tagger_info.lol_2_score);
  reader.AddVariable("tro_1_score",&tagger_info.tro_1_score);
  reader.AddVariable("tro_2_score",&tagger_info.tro_2_score);
  reader.AddVariable("tro_4_score",&tagger_info.tro_4_score);
  reader.AddVariable("tro_5_score",&tagger_info.tro_5_score);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/BDTcombine_BDT800_3.weights.xml");

  if (tagger_info.br_filled==1){
    val = reader.EvaluateMVA("MyBDT");

    if (flag_print){
      std::cout << "BDT: " << val << std::endl;
      std::cout << "BDT: " << tagger_info.mipid_score << " " <<  tagger_info.tro_1_score << " " << tagger_info.tro_2_score << " " << tagger_info.tro_4_score << " " << tagger_info.tro_5_score << std::endl;
      std::cout << "BDT: " << tagger_info.br4_tro_score << " " << tagger_info.lol_1_score << " " << tagger_info.lol_2_score << " " << tagger_info.sig_1_score << " " << tagger_info.sig_2_score << std::endl;
      std::cout << "BDT: " << tagger_info.br3_score << " " << tagger_info.br3_3_score << " " << tagger_info.br3_5_score << " " << tagger_info.br3_6_score << " " << tagger_info.gap_score << std::endl;
      std::cout << "BDT: " << tagger_info.stw_spt_score << " " << tagger_info.stw_2_score << " " << tagger_info.stw_3_score << " " << tagger_info.stw_4_score << " " << tagger_info.mgo_mgt_score << std::endl;
      std::cout << "BDT: " << tagger_info.vis_1_score << " " << tagger_info.vis_2_score << " " << tagger_info.pio_1_score << " " << tagger_info.pio_2_score << " " << tagger_info.mipquality_score << std::endl;
      std::cout << "BDT: " << tagger_info.br1_score << " " << tagger_info.stemdir_br2_score << " " << tagger_info.trimuon_score  <<  " " <<  tagger_info.hol_lol_score << " " << tagger_info.cme_anc_score << std::endl;
    }
    
  }else{
    val = default_val;
  }
    
  
  return val;
}

float WCPPID::NeutrinoID::cal_mipid_bdt(float default_val){
  TMVA::Reader reader;
  float val = default_val ; // default 

  

  reader.AddVariable("mip_energy",&tagger_info.mip_energy);
  reader.AddVariable("mip_n_end_reduction",&tagger_info.mip_n_end_reduction);
  reader.AddVariable("mip_n_first_mip",&tagger_info.mip_n_first_mip);
  reader.AddVariable("mip_n_first_non_mip",&tagger_info.mip_n_first_non_mip);
  reader.AddVariable("mip_n_first_non_mip_1",&tagger_info.mip_n_first_non_mip_1);
  reader.AddVariable("mip_n_first_non_mip_2",&tagger_info.mip_n_first_non_mip_2);
  reader.AddVariable("mip_vec_dQ_dx_0",&tagger_info.mip_vec_dQ_dx_0);
  reader.AddVariable("mip_vec_dQ_dx_1",&tagger_info.mip_vec_dQ_dx_1);
  reader.AddVariable("mip_max_dQ_dx_sample",&tagger_info.mip_max_dQ_dx_sample);
  reader.AddVariable("mip_n_below_threshold",&tagger_info.mip_n_below_threshold);
  reader.AddVariable("mip_n_below_zero",&tagger_info.mip_n_below_zero);
  reader.AddVariable("mip_n_lowest",&tagger_info.mip_n_lowest);
  reader.AddVariable("mip_n_highest",&tagger_info.mip_n_highest);
  reader.AddVariable("mip_lowest_dQ_dx",&tagger_info.mip_lowest_dQ_dx);
  reader.AddVariable("mip_highest_dQ_dx",&tagger_info.mip_highest_dQ_dx);
  reader.AddVariable("mip_medium_dQ_dx",&tagger_info.mip_medium_dQ_dx);
  reader.AddVariable("mip_stem_length",&tagger_info.mip_stem_length);
  reader.AddVariable("mip_length_main",&tagger_info.mip_length_main);
  reader.AddVariable("mip_length_total",&tagger_info.mip_length_total);
  reader.AddVariable("mip_angle_beam",&tagger_info.mip_angle_beam);
  reader.AddVariable("mip_iso_angle",&tagger_info.mip_iso_angle);
  reader.AddVariable("mip_n_vertex",&tagger_info.mip_n_vertex);
  reader.AddVariable("mip_n_good_tracks",&tagger_info.mip_n_good_tracks);
  reader.AddVariable("mip_E_indirect_max_energy",&tagger_info.mip_E_indirect_max_energy);
  reader.AddVariable("mip_flag_all_above",&tagger_info.mip_flag_all_above);
  reader.AddVariable("mip_min_dQ_dx_5",&tagger_info.mip_min_dQ_dx_5);
  reader.AddVariable("mip_n_other_vertex",&tagger_info.mip_n_other_vertex);
  reader.AddVariable("mip_n_stem_size",&tagger_info.mip_n_stem_size);
  reader.AddVariable("mip_flag_stem_trajectory",&tagger_info.mip_flag_stem_trajectory);
  reader.AddVariable("mip_min_dis",&tagger_info.mip_min_dis);

  reader.BookMVA( "MyBDT", "input_data_files/weights/mipid_BDT.weights.xml");

  if (tagger_info.mip_filled == 1)
    val = reader.EvaluateMVA("MyBDT");

  return val;
  
}

float WCPPID::NeutrinoID::cal_gap_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;

  reader.AddVariable("gap_flag_prolong_u",&tagger_info.gap_flag_prolong_u);
  reader.AddVariable("gap_flag_prolong_v",&tagger_info.gap_flag_prolong_v);
  reader.AddVariable("gap_flag_prolong_w",&tagger_info.gap_flag_prolong_w);
  reader.AddVariable("gap_flag_parallel",&tagger_info.gap_flag_parallel);
  reader.AddVariable("gap_n_points",&tagger_info.gap_n_points);
  reader.AddVariable("gap_n_bad",&tagger_info.gap_n_bad);
  reader.AddVariable("gap_energy",&tagger_info.gap_energy);
  reader.AddVariable("gap_num_valid_tracks",&tagger_info.gap_num_valid_tracks);
  reader.AddVariable("gap_flag_single_shower",&tagger_info.gap_flag_single_shower);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/gap_BDT.weights.xml");

  if (tagger_info.gap_filled==1)
  val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_hol_lol_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;

  reader.AddVariable("hol_1_n_valid_tracks",&tagger_info.hol_1_n_valid_tracks);
  reader.AddVariable("hol_1_min_angle",&tagger_info.hol_1_min_angle);
  reader.AddVariable("hol_1_energy",&tagger_info.hol_1_energy);
  reader.AddVariable("hol_1_flag_all_shower",&tagger_info.hol_1_flag_all_shower);
  reader.AddVariable("hol_1_min_length",&tagger_info.hol_1_min_length);

  reader.AddVariable("hol_2_min_angle",&tagger_info.hol_2_min_angle);
  reader.AddVariable("hol_2_medium_dQ_dx",&tagger_info.hol_2_medium_dQ_dx);
  reader.AddVariable("hol_2_ncount",&tagger_info.hol_2_ncount);
  
  reader.AddVariable("lol_3_angle_beam",&tagger_info.lol_3_angle_beam);
  reader.AddVariable("lol_3_n_valid_tracks",&tagger_info.lol_3_n_valid_tracks);
  reader.AddVariable("lol_3_min_angle",&tagger_info.lol_3_min_angle);
  reader.AddVariable("lol_3_vtx_n_segs",&tagger_info.lol_3_vtx_n_segs);
  reader.AddVariable("lol_3_shower_main_length",&tagger_info.lol_3_shower_main_length);
  reader.AddVariable("lol_3_n_out",&tagger_info.lol_3_n_out);
  reader.AddVariable("lol_3_n_sum",&tagger_info.lol_3_n_sum);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/hol_lol_BDT.weights.xml");

  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}


float WCPPID::NeutrinoID::cal_cme_anc_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;

  reader.AddVariable("cme_mu_energy",&tagger_info.cme_mu_energy);
  reader.AddVariable("cme_energy",&tagger_info.cme_energy);
  reader.AddVariable("cme_mu_length",&tagger_info.cme_mu_length);
  reader.AddVariable("cme_length",&tagger_info.cme_length);
  reader.AddVariable("cme_angle_beam",&tagger_info.cme_angle_beam);

  reader.AddVariable("anc_angle",&tagger_info.anc_angle);
  reader.AddVariable("anc_max_angle",&tagger_info.anc_max_angle);
  reader.AddVariable("anc_max_length",&tagger_info.anc_max_length);
  reader.AddVariable("anc_acc_forward_length",&tagger_info.anc_acc_forward_length);
  reader.AddVariable("anc_acc_backward_length",&tagger_info.anc_acc_backward_length);
  reader.AddVariable("anc_acc_forward_length1",&tagger_info.anc_acc_forward_length1);
  reader.AddVariable("anc_shower_main_length",&tagger_info.anc_shower_main_length);
  reader.AddVariable("anc_shower_total_length",&tagger_info.anc_shower_total_length);
  reader.AddVariable("anc_flag_main_outside",&tagger_info.anc_flag_main_outside);

 
  reader.BookMVA( "MyBDT", "input_data_files/weights/cme_anc_BDT.weights.xml");
 if (tagger_info.br_filled==1)
   val = reader.EvaluateMVA("MyBDT");
  return val;
  
}

float WCPPID::NeutrinoID::cal_mgo_mgt_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  
  reader.AddVariable("mgo_energy",&tagger_info.mgo_energy);
  reader.AddVariable("mgo_max_energy",&tagger_info.mgo_max_energy);
  reader.AddVariable("mgo_total_energy",&tagger_info.mgo_total_energy);
  reader.AddVariable("mgo_n_showers",&tagger_info.mgo_n_showers);
  reader.AddVariable("mgo_max_energy_1",&tagger_info.mgo_max_energy_1);
  reader.AddVariable("mgo_max_energy_2",&tagger_info.mgo_max_energy_2);
  reader.AddVariable("mgo_total_other_energy",&tagger_info.mgo_total_other_energy);
  reader.AddVariable("mgo_n_total_showers",&tagger_info.mgo_n_total_showers);
  reader.AddVariable("mgo_total_other_energy_1",&tagger_info.mgo_total_other_energy_1);
  
  reader.AddVariable("mgt_flag_single_shower",&tagger_info.mgt_flag_single_shower);
  reader.AddVariable("mgt_max_energy",&tagger_info.mgt_max_energy);
  reader.AddVariable("mgt_total_other_energy",&tagger_info.mgt_total_other_energy);
  reader.AddVariable("mgt_max_energy_1",&tagger_info.mgt_max_energy_1);
  reader.AddVariable("mgt_e_indirect_max_energy",&tagger_info.mgt_e_indirect_max_energy);
  reader.AddVariable("mgt_e_direct_max_energy",&tagger_info.mgt_e_direct_max_energy);
  reader.AddVariable("mgt_n_direct_showers",&tagger_info.mgt_n_direct_showers);
  reader.AddVariable("mgt_e_direct_total_energy",&tagger_info.mgt_e_direct_total_energy);
  reader.AddVariable("mgt_flag_indirect_max_pio",&tagger_info.mgt_flag_indirect_max_pio);
  reader.AddVariable("mgt_e_indirect_total_energy",&tagger_info.mgt_e_indirect_total_energy);


  reader.BookMVA( "MyBDT", "input_data_files/weights/mgo_mgt_BDT.weights.xml");
  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_br1_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  
  reader.AddVariable("br1_1_shower_type",&tagger_info.br1_1_shower_type);
  reader.AddVariable("br1_1_vtx_n_segs",&tagger_info.br1_1_vtx_n_segs);
  reader.AddVariable("br1_1_energy",&tagger_info.br1_1_energy);
  reader.AddVariable("br1_1_n_segs",&tagger_info.br1_1_n_segs);
  reader.AddVariable("br1_1_flag_sg_topology",&tagger_info.br1_1_flag_sg_topology);
  reader.AddVariable("br1_1_flag_sg_trajectory",&tagger_info.br1_1_flag_sg_trajectory);
  reader.AddVariable("br1_1_sg_length",&tagger_info.br1_1_sg_length);

  reader.AddVariable("br1_2_n_connected",&tagger_info.br1_2_n_connected);
  reader.AddVariable("br1_2_max_length",&tagger_info.br1_2_max_length);
  reader.AddVariable("br1_2_n_connected_1",&tagger_info.br1_2_n_connected_1);
  reader.AddVariable("br1_2_n_shower_segs",&tagger_info.br1_2_n_shower_segs);
  reader.AddVariable("br1_2_max_length_ratio",&tagger_info.br1_2_max_length_ratio);
  reader.AddVariable("br1_2_shower_length",&tagger_info.br1_2_shower_length);

  reader.AddVariable("br1_3_n_connected_p",&tagger_info.br1_3_n_connected_p);
  reader.AddVariable("br1_3_max_length_p",&tagger_info.br1_3_max_length_p);
  reader.AddVariable("br1_3_n_shower_main_segs",&tagger_info.br1_3_n_shower_main_segs);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/br1_BDT.weights.xml");

  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_br3_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;

  reader.AddVariable("br3_1_energy",&tagger_info.br3_1_energy);
  reader.AddVariable("br3_1_n_shower_segments",&tagger_info.br3_1_n_shower_segments);
  reader.AddVariable("br3_1_sg_flag_trajectory",&tagger_info.br3_1_sg_flag_trajectory);
  reader.AddVariable("br3_1_sg_direct_length",&tagger_info.br3_1_sg_direct_length);
  reader.AddVariable("br3_1_sg_length",&tagger_info.br3_1_sg_length);
  reader.AddVariable("br3_1_total_main_length",&tagger_info.br3_1_total_main_length);
  reader.AddVariable("br3_1_total_length",&tagger_info.br3_1_total_length);
  reader.AddVariable("br3_1_iso_angle",&tagger_info.br3_1_iso_angle);
  reader.AddVariable("br3_1_sg_flag_topology",&tagger_info.br3_1_sg_flag_topology);
  reader.AddVariable("br3_2_n_ele",&tagger_info.br3_2_n_ele);
  reader.AddVariable("br3_2_n_other",&tagger_info.br3_2_n_other);
  reader.AddVariable("br3_2_other_fid",&tagger_info.br3_2_other_fid);
  reader.AddVariable("br3_4_acc_length",&tagger_info.br3_4_acc_length);
  reader.AddVariable("br3_4_total_length",&tagger_info.br3_4_total_length);
  reader.AddVariable("br3_7_min_angle",&tagger_info.br3_7_min_angle);
  reader.AddVariable("br3_8_max_dQ_dx",&tagger_info.br3_8_max_dQ_dx);
  reader.AddVariable("br3_8_n_main_segs",&tagger_info.br3_8_n_main_segs);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/br3_BDT.weights.xml");

  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_br3_3_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float br3_3_v_energy;
  float br3_3_v_angle;
  float br3_3_v_dir_length;
  float br3_3_v_length;
  
  reader.AddVariable("br3_3_v_energy",&br3_3_v_energy);
  reader.AddVariable("br3_3_v_angle",&br3_3_v_angle);
  reader.AddVariable("br3_3_v_dir_length",&br3_3_v_dir_length);
  reader.AddVariable("br3_3_v_length",&br3_3_v_length);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/br3_3_BDT.weights.xml");
  
  for (size_t i = 0; i!= tagger_info.br3_3_v_energy.size(); i++){
    br3_3_v_energy = tagger_info.br3_3_v_energy.at(i);
    br3_3_v_angle = tagger_info.br3_3_v_angle.at(i);
    br3_3_v_dir_length = tagger_info.br3_3_v_dir_length.at(i);
    br3_3_v_length = tagger_info.br3_3_v_length.at(i);

    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val)     val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_br3_5_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float br3_5_v_dir_length;
  float br3_5_v_total_length;
  float br3_5_v_flag_avoid_muon_check;
  float br3_5_v_n_seg;
  float br3_5_v_angle;
  float br3_5_v_sg_length;
  float br3_5_v_energy;
  float br3_5_v_n_main_segs;
  float br3_5_v_n_segs;
  float br3_5_v_shower_main_length;
  float br3_5_v_shower_total_length;
  
  reader.AddVariable("br3_5_v_dir_length",&br3_5_v_dir_length);
  reader.AddVariable("br3_5_v_total_length",&br3_5_v_total_length);
  reader.AddVariable("br3_5_v_flag_avoid_muon_check",&br3_5_v_flag_avoid_muon_check);
  reader.AddVariable("br3_5_v_n_seg",&br3_5_v_n_seg);
  reader.AddVariable("br3_5_v_angle",&br3_5_v_angle);
  reader.AddVariable("br3_5_v_sg_length",&br3_5_v_sg_length);
  reader.AddVariable("br3_5_v_energy",&br3_5_v_energy);
  //reader.AddVariable("br3_5_v_n_main_segs",&br3_5_v_n_main_segs);
  reader.AddVariable("br3_5_v_n_segs",&br3_5_v_n_segs);
  reader.AddVariable("br3_5_v_shower_main_length",&br3_5_v_shower_main_length);
  reader.AddVariable("br3_5_v_shower_total_length",&br3_5_v_shower_total_length);
    
  reader.BookMVA( "MyBDT", "input_data_files/weights/br3_5_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.br3_5_v_dir_length.size();i++){
    br3_5_v_dir_length = tagger_info.br3_5_v_dir_length.at(i);
    br3_5_v_total_length = tagger_info.br3_5_v_total_length.at(i);
    br3_5_v_flag_avoid_muon_check = tagger_info.br3_5_v_flag_avoid_muon_check.at(i);
    br3_5_v_n_seg = tagger_info.br3_5_v_n_seg.at(i);
    br3_5_v_angle = tagger_info.br3_5_v_angle.at(i);
    br3_5_v_sg_length = tagger_info.br3_5_v_sg_length.at(i);
    br3_5_v_energy = tagger_info.br3_5_v_energy.at(i);
    br3_5_v_n_segs = tagger_info.br3_5_v_n_segs.at(i);
    br3_5_v_shower_main_length = tagger_info.br3_5_v_shower_main_length.at(i);
    br3_5_v_shower_total_length = tagger_info.br3_5_v_shower_total_length.at(i);
    
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_br3_6_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;
  float br3_6_v_angle;
  float br3_6_v_angle1;
  float br3_6_v_flag_shower_trajectory;
  float br3_6_v_direct_length;
  float br3_6_v_length;
  float br3_6_v_n_other_vtx_segs;
  float br3_6_v_energy;
  
  reader.AddVariable("br3_6_v_angle",&br3_6_v_angle);
  reader.AddVariable("br3_6_v_angle1",&br3_6_v_angle1);
  reader.AddVariable("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
  reader.AddVariable("br3_6_v_direct_length",&br3_6_v_direct_length);
  reader.AddVariable("br3_6_v_length",&br3_6_v_length);
  reader.AddVariable("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
  reader.AddVariable("br3_6_v_energy",&br3_6_v_energy);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/br3_6_BDT.weights.xml");

  for (size_t i=0; i!= tagger_info.br3_6_v_angle.size();i++){
    br3_6_v_angle = tagger_info.br3_6_v_angle.at(i);
    br3_6_v_angle1 = tagger_info.br3_6_v_angle1.at(i);
    br3_6_v_flag_shower_trajectory = tagger_info.br3_6_v_flag_shower_trajectory.at(i);
    br3_6_v_direct_length = tagger_info.br3_6_v_direct_length.at(i);
    br3_6_v_length = tagger_info.br3_6_v_length.at(i);
    br3_6_v_n_other_vtx_segs = tagger_info.br3_6_v_n_other_vtx_segs.at(i);
    br3_6_v_energy = tagger_info.br3_6_v_energy.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_stemdir_br2_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  
  reader.AddVariable("stem_dir_flag_single_shower",&tagger_info.stem_dir_flag_single_shower);
  reader.AddVariable("stem_dir_angle",&tagger_info.stem_dir_angle);
  reader.AddVariable("stem_dir_energy",&tagger_info.stem_dir_energy);
  reader.AddVariable("stem_dir_angle1",&tagger_info.stem_dir_angle1);
  reader.AddVariable("stem_dir_angle2",&tagger_info.stem_dir_angle2);
  reader.AddVariable("stem_dir_angle3",&tagger_info.stem_dir_angle3);
  reader.AddVariable("stem_dir_ratio",&tagger_info.stem_dir_ratio);
  
  reader.AddVariable("br2_num_valid_tracks",&tagger_info.br2_num_valid_tracks);
  reader.AddVariable("br2_n_shower_main_segs",&tagger_info.br2_n_shower_main_segs);
  reader.AddVariable("br2_max_angle",&tagger_info.br2_max_angle);
  reader.AddVariable("br2_sg_length",&tagger_info.br2_sg_length);
  reader.AddVariable("br2_flag_sg_trajectory",&tagger_info.br2_flag_sg_trajectory);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/stem_dir_br2_BDT.weights.xml");

  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}


float WCPPID::NeutrinoID::cal_trimuon_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  reader.AddVariable("stem_len_energy",&tagger_info.stem_len_energy);
  reader.AddVariable("stem_len_length",&tagger_info.stem_len_length);
  reader.AddVariable("stem_len_flag_avoid_muon_check",&tagger_info.stem_len_flag_avoid_muon_check);
  reader.AddVariable("stem_len_num_daughters",&tagger_info.stem_len_num_daughters);
  reader.AddVariable("stem_len_daughter_length",&tagger_info.stem_len_daughter_length);
  
  reader.AddVariable("brm_n_mu_segs",&tagger_info.brm_n_mu_segs);
  reader.AddVariable("brm_Ep",&tagger_info.brm_Ep);
  reader.AddVariable("brm_acc_length",&tagger_info.brm_acc_length);
  reader.AddVariable("brm_shower_total_length",&tagger_info.brm_shower_total_length);
  reader.AddVariable("brm_connected_length",&tagger_info.brm_connected_length);
  reader.AddVariable("brm_n_size",&tagger_info.brm_n_size);
  reader.AddVariable("brm_acc_direct_length",&tagger_info.brm_acc_direct_length);
  reader.AddVariable("brm_n_shower_main_segs",&tagger_info.brm_n_shower_main_segs);
  reader.AddVariable("brm_n_mu_main",&tagger_info.brm_n_mu_main);
  
  reader.AddVariable("lem_shower_main_length",&tagger_info.lem_shower_main_length);
  reader.AddVariable("lem_n_3seg",&tagger_info.lem_n_3seg);
  reader.AddVariable("lem_e_charge",&tagger_info.lem_e_charge);
  reader.AddVariable("lem_e_dQdx",&tagger_info.lem_e_dQdx);
  reader.AddVariable("lem_shower_num_main_segs",&tagger_info.lem_shower_num_main_segs);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/stl_lem_brm_BDT.weights.xml");

  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}


float WCPPID::NeutrinoID::cal_br4_tro_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  
  reader.AddVariable("br4_1_shower_main_length",&tagger_info.br4_1_shower_main_length);
  reader.AddVariable("br4_1_shower_total_length",&tagger_info.br4_1_shower_total_length);
  reader.AddVariable("br4_1_min_dis",&tagger_info.br4_1_min_dis);
  reader.AddVariable("br4_1_energy",&tagger_info.br4_1_energy);
  reader.AddVariable("br4_1_flag_avoid_muon_check",&tagger_info.br4_1_flag_avoid_muon_check);
  reader.AddVariable("br4_1_n_vtx_segs",&tagger_info.br4_1_n_vtx_segs);
  reader.AddVariable("br4_1_n_main_segs",&tagger_info.br4_1_n_main_segs);
  
  reader.AddVariable("br4_2_ratio_45",&tagger_info.br4_2_ratio_45);
  reader.AddVariable("br4_2_ratio_35",&tagger_info.br4_2_ratio_35);
  reader.AddVariable("br4_2_ratio_25",&tagger_info.br4_2_ratio_25);
  reader.AddVariable("br4_2_ratio_15",&tagger_info.br4_2_ratio_15);
  reader.AddVariable("br4_2_ratio1_45",&tagger_info.br4_2_ratio1_45);
  reader.AddVariable("br4_2_ratio1_35",&tagger_info.br4_2_ratio1_35);
  reader.AddVariable("br4_2_ratio1_25",&tagger_info.br4_2_ratio1_25);
  reader.AddVariable("br4_2_ratio1_15",&tagger_info.br4_2_ratio1_15);
  reader.AddVariable("br4_2_iso_angle",&tagger_info.br4_2_iso_angle);
  reader.AddVariable("br4_2_iso_angle1",&tagger_info.br4_2_iso_angle1);
  reader.AddVariable("br4_2_angle",&tagger_info.br4_2_angle);

  reader.AddVariable("tro_3_stem_length",&tagger_info.tro_3_stem_length);
  reader.AddVariable("tro_3_n_muon_segs",&tagger_info.tro_3_n_muon_segs);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/br4_tro_BDT.weights.xml");

  if (tagger_info.br_filled==1) val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_mipquality_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  
  reader.AddVariable("mip_quality_energy",&tagger_info.mip_quality_energy);
  reader.AddVariable("mip_quality_overlap",&tagger_info.mip_quality_overlap);
  reader.AddVariable("mip_quality_n_showers",&tagger_info.mip_quality_n_showers);
  reader.AddVariable("mip_quality_n_tracks",&tagger_info.mip_quality_n_tracks);
  reader.AddVariable("mip_quality_flag_inside_pi0",&tagger_info.mip_quality_flag_inside_pi0);
  reader.AddVariable("mip_quality_n_pi0_showers",&tagger_info.mip_quality_n_pi0_showers);
  reader.AddVariable("mip_quality_shortest_length",&tagger_info.mip_quality_shortest_length);
  reader.AddVariable("mip_quality_acc_length",&tagger_info.mip_quality_acc_length);
  reader.AddVariable("mip_quality_shortest_angle",&tagger_info.mip_quality_shortest_angle);
  reader.AddVariable("mip_quality_flag_proton",&tagger_info.mip_quality_flag_proton);
  reader.BookMVA( "MyBDT", "input_data_files/weights/mipquality_BDT.weights.xml");

  if (tagger_info.mip_quality_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_pio_1_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;
  
  reader.AddVariable("pio_1_mass",&tagger_info.pio_1_mass);
  reader.AddVariable("pio_1_pio_type",&tagger_info.pio_1_pio_type);
  reader.AddVariable("pio_1_energy_1",&tagger_info.pio_1_energy_1);
  reader.AddVariable("pio_1_energy_2",&tagger_info.pio_1_energy_2);
  reader.AddVariable("pio_1_dis_1",&tagger_info.pio_1_dis_1);
  reader.AddVariable("pio_1_dis_2",&tagger_info.pio_1_dis_2);
  reader.AddVariable("pio_mip_id",&tagger_info.pio_mip_id);
  reader.BookMVA( "MyBDT", "input_data_files/weights/pio_1_BDT.weights.xml");

  if (tagger_info.pio_filled==1 && tagger_info.pio_flag_pio==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_pio_2_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float pio_2_v_dis2;
  float pio_2_v_angle2;
  float pio_2_v_acc_length;
  
  reader.AddVariable("pio_2_v_dis2",&pio_2_v_dis2);
  reader.AddVariable("pio_2_v_angle2",&pio_2_v_angle2);
  reader.AddVariable("pio_2_v_acc_length",&pio_2_v_acc_length);
  reader.AddVariable("pio_mip_id",&tagger_info.pio_mip_id);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/pio_2_BDT.weights.xml");
  
  if (tagger_info.pio_filled==1 && tagger_info.pio_flag_pio==0){
    for (size_t i=0;i!=tagger_info.pio_2_v_dis2.size();i++){
      pio_2_v_dis2 = tagger_info.pio_2_v_dis2.at(i);
      pio_2_v_angle2 = tagger_info.pio_2_v_angle2.at(i);
      pio_2_v_acc_length = tagger_info.pio_2_v_acc_length.at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  if (val > 1e8) val = default_val;

  
  return val;
  
}


float WCPPID::NeutrinoID::cal_stw_spt_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;

  reader.AddVariable("stw_1_energy",&tagger_info.stw_1_energy);
  reader.AddVariable("stw_1_dis",&tagger_info.stw_1_dis);
  reader.AddVariable("stw_1_dQ_dx",&tagger_info.stw_1_dQ_dx);
  reader.AddVariable("stw_1_flag_single_shower",&tagger_info.stw_1_flag_single_shower);
  reader.AddVariable("stw_1_n_pi0",&tagger_info.stw_1_n_pi0);
  reader.AddVariable("stw_1_num_valid_tracks",&tagger_info.stw_1_num_valid_tracks);
  
  reader.AddVariable("spt_shower_main_length",&tagger_info.spt_shower_main_length);
  reader.AddVariable("spt_shower_total_length",&tagger_info.spt_shower_total_length);
  reader.AddVariable("spt_angle_beam",&tagger_info.spt_angle_beam);
  reader.AddVariable("spt_angle_vertical",&tagger_info.spt_angle_vertical);
  reader.AddVariable("spt_max_dQ_dx",&tagger_info.spt_max_dQ_dx);
  reader.AddVariable("spt_angle_beam_1",&tagger_info.spt_angle_beam_1);
  reader.AddVariable("spt_angle_drift",&tagger_info.spt_angle_drift);
  reader.AddVariable("spt_angle_drift_1",&tagger_info.spt_angle_drift_1);
  reader.AddVariable("spt_num_valid_tracks",&tagger_info.spt_num_valid_tracks);
  reader.AddVariable("spt_n_vtx_segs",&tagger_info.spt_n_vtx_segs);
  reader.AddVariable("spt_max_length",&tagger_info.spt_max_length);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/stw_spt_BDT.weights.xml");

  if (tagger_info.br_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_vis_1_bdt(float default_val){
  
  float val = default_val;
  
  TMVA::Reader reader;
  reader.AddVariable("vis_1_n_vtx_segs",&tagger_info.vis_1_n_vtx_segs);
  reader.AddVariable("vis_1_energy",&tagger_info.vis_1_energy);
  reader.AddVariable("vis_1_num_good_tracks",&tagger_info.vis_1_num_good_tracks);
  reader.AddVariable("vis_1_max_angle",&tagger_info.vis_1_max_angle);
  reader.AddVariable("vis_1_max_shower_angle",&tagger_info.vis_1_max_shower_angle);
  reader.AddVariable("vis_1_tmp_length1",&tagger_info.vis_1_tmp_length1);
  reader.AddVariable("vis_1_tmp_length2",&tagger_info.vis_1_tmp_length2);
  //  reader.AddVariable("vis_1_particle_type",&tagger_info.vis_1_particle_type);
    
  reader.BookMVA( "MyBDT", "input_data_files/weights/vis_1_BDT.weights.xml");

  if (tagger_info.vis_1_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}


float WCPPID::NeutrinoID::cal_vis_2_bdt(float default_val){
  
  float val = default_val;
  TMVA::Reader reader;

  reader.AddVariable("vis_2_n_vtx_segs",&tagger_info.vis_2_n_vtx_segs);
  reader.AddVariable("vis_2_min_angle",&tagger_info.vis_2_min_angle);
  reader.AddVariable("vis_2_min_weak_track",&tagger_info.vis_2_min_weak_track);
  reader.AddVariable("vis_2_angle_beam", &tagger_info.vis_2_angle_beam);
  reader.AddVariable("vis_2_min_angle1",&tagger_info.vis_2_min_angle1);
  reader.AddVariable("vis_2_iso_angle1",&tagger_info.vis_2_iso_angle1);
  reader.AddVariable("vis_2_min_medium_dQ_dx",&tagger_info.vis_2_min_medium_dQ_dx);
  reader.AddVariable("vis_2_min_length",&tagger_info.vis_2_min_length);
  reader.AddVariable("vis_2_sg_length",&tagger_info.vis_2_sg_length);
  reader.AddVariable("vis_2_max_angle",&tagger_info.vis_2_max_angle);
  reader.AddVariable("vis_2_max_weak_track",&tagger_info.vis_2_max_weak_track);


  reader.BookMVA( "MyBDT", "input_data_files/weights/vis_2_BDT.weights.xml");

  if (tagger_info.vis_2_filled==1)
    val = reader.EvaluateMVA("MyBDT");
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_stw_2_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;
  
  float stw_2_v_medium_dQ_dx;
  float stw_2_v_energy;
  float stw_2_v_angle;
  float stw_2_v_dir_length;
  float stw_2_v_max_dQ_dx;
  
  reader.AddVariable("stw_2_v_medium_dQ_dx",&stw_2_v_medium_dQ_dx);
  reader.AddVariable("stw_2_v_energy",&stw_2_v_energy);
  reader.AddVariable("stw_2_v_angle",&stw_2_v_angle);
  reader.AddVariable("stw_2_v_dir_length",&stw_2_v_dir_length);
  reader.AddVariable("stw_2_v_max_dQ_dx",&stw_2_v_max_dQ_dx);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/stw_2_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.stw_2_v_medium_dQ_dx.size();i++){
    stw_2_v_medium_dQ_dx = tagger_info.stw_2_v_medium_dQ_dx.at(i);
    stw_2_v_energy = tagger_info.stw_2_v_energy.at(i);
    stw_2_v_angle = tagger_info.stw_2_v_angle.at(i);
    stw_2_v_dir_length = tagger_info.stw_2_v_dir_length.at(i);
    stw_2_v_max_dQ_dx = tagger_info.stw_2_v_max_dQ_dx.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  return val;
  
}

float WCPPID::NeutrinoID::cal_stw_3_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float stw_3_v_angle;
  float stw_3_v_dir_length;
  float stw_3_v_energy;
  float stw_3_v_medium_dQ_dx;
  
  reader.AddVariable("stw_3_v_angle",&stw_3_v_angle);
  reader.AddVariable("stw_3_v_dir_length",&stw_3_v_dir_length);
  reader.AddVariable("stw_3_v_energy",&stw_3_v_energy);
  reader.AddVariable("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/stw_3_BDT.weights.xml");

  for(size_t i=0;i!=tagger_info.stw_3_v_angle.size();i++){
    stw_3_v_angle = tagger_info.stw_3_v_angle.at(i);
    stw_3_v_dir_length = tagger_info.stw_3_v_dir_length.at(i);
    stw_3_v_energy = tagger_info.stw_3_v_energy.at(i);
    stw_3_v_medium_dQ_dx = tagger_info.stw_3_v_medium_dQ_dx.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_stw_4_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float stw_4_v_angle;
  float stw_4_v_dis;
  float stw_4_v_energy;
  
  reader.AddVariable("stw_4_v_angle",&stw_4_v_angle);
  reader.AddVariable("stw_4_v_dis",&stw_4_v_dis);
  reader.AddVariable("stw_4_v_energy",&stw_4_v_energy);
  reader.BookMVA( "MyBDT", "input_data_files/weights/stw_4_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.stw_4_v_angle.size();i++){
    stw_4_v_angle = tagger_info.stw_4_v_angle.at(i);
    stw_4_v_dis = tagger_info.stw_4_v_dis.at(i);
    stw_4_v_energy = tagger_info.stw_4_v_energy.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_sig_1_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float sig_1_v_angle;
  float sig_1_v_flag_single_shower;
  float sig_1_v_energy;
  float sig_1_v_energy_1;
  
  reader.AddVariable("sig_1_v_angle",&sig_1_v_angle);
  reader.AddVariable("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
  reader.AddVariable("sig_1_v_energy",&sig_1_v_energy);
  reader.AddVariable("sig_1_v_energy_1",&sig_1_v_energy_1);
  reader.BookMVA( "MyBDT", "input_data_files/weights/sig_1_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.sig_1_v_angle.size();i++){
    sig_1_v_angle = tagger_info.sig_1_v_angle.at(i);
    sig_1_v_flag_single_shower = tagger_info.sig_1_v_flag_single_shower.at(i);
    sig_1_v_energy = tagger_info.sig_1_v_energy.at(i);
    sig_1_v_energy_1 = tagger_info.sig_1_v_energy_1.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_sig_2_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;
  float sig_2_v_energy;
  float sig_2_v_shower_angle;
  float sig_2_v_flag_single_shower;
  float sig_2_v_medium_dQ_dx;
  float sig_2_v_start_dQ_dx;
  
  reader.AddVariable("sig_2_v_energy",&sig_2_v_energy);
  reader.AddVariable("sig_2_v_shower_angle",&sig_2_v_shower_angle);
  reader.AddVariable("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
  reader.AddVariable("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);  
  reader.AddVariable("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);

  reader.BookMVA( "MyBDT", "input_data_files/weights/sig_2_BDT.weights.xml");

  for (size_t i=0; i!= tagger_info.sig_2_v_energy.size();i++){
    sig_2_v_energy = tagger_info.sig_2_v_energy.at(i);
    sig_2_v_shower_angle = tagger_info.sig_2_v_shower_angle.at(i);
    sig_2_v_flag_single_shower = tagger_info.sig_2_v_flag_single_shower.at(i);
    sig_2_v_medium_dQ_dx = tagger_info.sig_2_v_medium_dQ_dx.at(i);
    sig_2_v_start_dQ_dx = tagger_info.sig_2_v_start_dQ_dx.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_lol_1_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float lol_1_v_energy;
  float lol_1_v_vtx_n_segs;
  float lol_1_v_nseg;
  float lol_1_v_angle;
  
  reader.AddVariable("lol_1_v_energy",&lol_1_v_energy);
  reader.AddVariable("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
  reader.AddVariable("lol_1_v_nseg",&lol_1_v_nseg);
  reader.AddVariable("lol_1_v_angle",&lol_1_v_angle);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/lol_1_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.lol_1_v_energy.size(); i++){
    lol_1_v_energy = tagger_info.lol_1_v_energy.at(i);
    lol_1_v_vtx_n_segs = tagger_info.lol_1_v_vtx_n_segs.at(i);
    lol_1_v_nseg = tagger_info.lol_1_v_nseg.at(i);
    lol_1_v_angle = tagger_info.lol_1_v_angle.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }

  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_lol_2_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float lol_2_v_length;
  float lol_2_v_angle;
  float lol_2_v_type;
  float lol_2_v_vtx_n_segs;
  float lol_2_v_energy;
  float lol_2_v_shower_main_length;
  float lol_2_v_flag_dir_weak;

  reader.AddVariable("lol_2_v_length",&lol_2_v_length);
  reader.AddVariable("lol_2_v_angle",&lol_2_v_angle);
  reader.AddVariable("lol_2_v_type",&lol_2_v_type);
  reader.AddVariable("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
  reader.AddVariable("lol_2_v_energy",&lol_2_v_energy);
  reader.AddVariable("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
  reader.AddVariable("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);

  reader.BookMVA( "MyBDT", "input_data_files/weights/lol_2_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.lol_2_v_length.size();i++){
    lol_2_v_length = tagger_info.lol_2_v_length.at(i);
    lol_2_v_angle = tagger_info.lol_2_v_angle.at(i);
    lol_2_v_type = tagger_info.lol_2_v_type.at(i);
    lol_2_v_vtx_n_segs = tagger_info.lol_2_v_vtx_n_segs.at(i);
    lol_2_v_energy = tagger_info.lol_2_v_energy.at(i);
    lol_2_v_shower_main_length = tagger_info.lol_2_v_shower_main_length.at(i);
    lol_2_v_flag_dir_weak = tagger_info.lol_2_v_flag_dir_weak.at(i);
      
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}


float WCPPID::NeutrinoID::cal_tro_1_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float tro_1_v_particle_type;
  float tro_1_v_flag_dir_weak;
  float tro_1_v_min_dis;
  float tro_1_v_sg1_length;
  float tro_1_v_shower_main_length;
  float tro_1_v_max_n_vtx_segs;
  float tro_1_v_tmp_length;
  float tro_1_v_medium_dQ_dx;
  float tro_1_v_dQ_dx_cut;
  float tro_1_v_flag_shower_topology;
  
  reader.AddVariable("tro_1_v_particle_type",&tro_1_v_particle_type);
  reader.AddVariable("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
  reader.AddVariable("tro_1_v_min_dis",&tro_1_v_min_dis);
  reader.AddVariable("tro_1_v_sg1_length",&tro_1_v_sg1_length);
  reader.AddVariable("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
  reader.AddVariable("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
  reader.AddVariable("tro_1_v_tmp_length",&tro_1_v_tmp_length);
  reader.AddVariable("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
  reader.AddVariable("tro_1_v_dQ_dx_cut", &tro_1_v_dQ_dx_cut);
  reader.AddVariable("tro_1_v_flag_shower_topology", &tro_1_v_flag_shower_topology);

  reader.BookMVA( "MyBDT", "input_data_files/weights/tro_1_BDT.weights.xml");
  
  for (size_t i=0;i!=tagger_info.tro_1_v_particle_type.size();i++){
    tro_1_v_particle_type = tagger_info.tro_1_v_particle_type.at(i);
    tro_1_v_flag_dir_weak = tagger_info.tro_1_v_flag_dir_weak.at(i);
    tro_1_v_min_dis = tagger_info.tro_1_v_min_dis.at(i);
    tro_1_v_sg1_length = tagger_info.tro_1_v_sg1_length.at(i);
    tro_1_v_shower_main_length = tagger_info.tro_1_v_shower_main_length.at(i);
    tro_1_v_max_n_vtx_segs = tagger_info.tro_1_v_max_n_vtx_segs.at(i);
    tro_1_v_tmp_length = tagger_info.tro_1_v_tmp_length.at(i);
    tro_1_v_medium_dQ_dx = tagger_info.tro_1_v_medium_dQ_dx.at(i);
    tro_1_v_dQ_dx_cut = tagger_info.tro_1_v_dQ_dx_cut.at(i);
    tro_1_v_flag_shower_topology = tagger_info.tro_1_v_flag_shower_topology.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_tro_2_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float tro_2_v_energy;
  float tro_2_v_stem_length;
  float tro_2_v_iso_angle;
  float tro_2_v_max_length;
  float tro_2_v_angle;
  
  reader.AddVariable("tro_2_v_energy",&tro_2_v_energy);
  reader.AddVariable("tro_2_v_stem_length",&tro_2_v_stem_length);
  reader.AddVariable("tro_2_v_iso_angle",&tro_2_v_iso_angle);
  reader.AddVariable("tro_2_v_max_length",&tro_2_v_max_length);
  reader.AddVariable("tro_2_v_angle",&tro_2_v_angle);

  reader.BookMVA( "MyBDT", "input_data_files/weights/tro_2_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.tro_2_v_energy.size();i++){

    tro_2_v_energy = tagger_info.tro_2_v_energy.at(i);
    tro_2_v_stem_length = tagger_info.tro_2_v_stem_length.at(i);
    tro_2_v_iso_angle = tagger_info.tro_2_v_iso_angle.at(i);
    tro_2_v_max_length = tagger_info.tro_2_v_max_length.at(i);
    tro_2_v_angle = tagger_info.tro_2_v_angle.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_tro_4_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float tro_4_v_dir2_mag;
  float tro_4_v_angle;
  float tro_4_v_angle1;
  float tro_4_v_angle2;
  float tro_4_v_length;
  float tro_4_v_length1;
  float tro_4_v_medium_dQ_dx;
  float tro_4_v_end_dQ_dx;
  float tro_4_v_energy;
  float tro_4_v_shower_main_length;
  float tro_4_v_flag_shower_trajectory;
  
  reader.AddVariable("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
  reader.AddVariable("tro_4_v_angle",&tro_4_v_angle);
  reader.AddVariable("tro_4_v_angle1",&tro_4_v_angle1);
  reader.AddVariable("tro_4_v_angle2",&tro_4_v_angle2);
  reader.AddVariable("tro_4_v_length",&tro_4_v_length);
  reader.AddVariable("tro_4_v_length1",&tro_4_v_length1);
  reader.AddVariable("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
  reader.AddVariable("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
  reader.AddVariable("tro_4_v_energy",&tro_4_v_energy);
  reader.AddVariable("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
  reader.AddVariable("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/tro_4_BDT.weights.xml");

  for (size_t i=0; i!= tagger_info.tro_4_v_dir2_mag.size(); i++){

    tro_4_v_dir2_mag = tagger_info.tro_4_v_dir2_mag.at(i);
    tro_4_v_angle = tagger_info.tro_4_v_angle.at(i);
    tro_4_v_angle1 = tagger_info.tro_4_v_angle1.at(i);
    tro_4_v_angle2 = tagger_info.tro_4_v_angle2.at(i);
    tro_4_v_length = tagger_info.tro_4_v_length.at(i);
    tro_4_v_length1 = tagger_info.tro_4_v_length1.at(i);
    tro_4_v_medium_dQ_dx = tagger_info.tro_4_v_medium_dQ_dx.at(i);
    tro_4_v_end_dQ_dx = tagger_info.tro_4_v_end_dQ_dx.at(i);
    tro_4_v_energy = tagger_info.tro_4_v_energy.at(i);
    tro_4_v_shower_main_length = tagger_info.tro_4_v_shower_main_length.at(i);
    tro_4_v_flag_shower_trajectory = tagger_info.tro_4_v_flag_shower_trajectory.at(i);
    
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}

float WCPPID::NeutrinoID::cal_tro_5_bdt(float default_val){
  
  float val = 1e9;
  TMVA::Reader reader;

  float tro_5_v_max_angle;
  float tro_5_v_min_angle;
  float tro_5_v_max_length;
  float tro_5_v_iso_angle;
  float tro_5_v_n_vtx_segs;
  float tro_5_v_min_count;
  float tro_5_v_max_count;
  float tro_5_v_energy;
  
  reader.AddVariable("tro_5_v_max_angle",&tro_5_v_max_angle);
  reader.AddVariable("tro_5_v_min_angle",&tro_5_v_min_angle);
  reader.AddVariable("tro_5_v_max_length",&tro_5_v_max_length);
  reader.AddVariable("tro_5_v_iso_angle",&tro_5_v_iso_angle);
  reader.AddVariable("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
  reader.AddVariable("tro_5_v_min_count",&tro_5_v_min_count);
  reader.AddVariable("tro_5_v_max_count",&tro_5_v_max_count);
  reader.AddVariable("tro_5_v_energy",&tro_5_v_energy);
  
  reader.BookMVA( "MyBDT", "input_data_files/weights/tro_5_BDT.weights.xml");

  for (size_t i=0;i!=tagger_info.tro_5_v_max_angle.size();i++){

    tro_5_v_max_angle = tagger_info.tro_5_v_max_angle.at(i);
    tro_5_v_min_angle = tagger_info.tro_5_v_min_angle.at(i);
    tro_5_v_max_length = tagger_info.tro_5_v_max_length.at(i);
    tro_5_v_iso_angle = tagger_info.tro_5_v_iso_angle.at(i);
    tro_5_v_n_vtx_segs = tagger_info.tro_5_v_n_vtx_segs.at(i);
    tro_5_v_min_count = tagger_info.tro_5_v_min_count.at(i);
    tro_5_v_max_count = tagger_info.tro_5_v_max_count.at(i);
    tro_5_v_energy = tagger_info.tro_5_v_energy.at(i);
        
    float tmp_val = reader.EvaluateMVA("MyBDT");
    if (tmp_val < val) val = tmp_val;
  }
  if (val > 1e8) val = default_val;
  
  return val;
  
}
