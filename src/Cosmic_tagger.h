 
std::tuple<int, WireCellPID::PR3DCluster*, WireCell::Opflash*> WireCellPID::ToyFiducial::cosmic_tagger(WireCell::OpflashSelection& flashes, WireCellPID::PR3DCluster* main_cluster, WireCell::Opflash* main_flash, std::tuple<int, double, double, int>& bundle_info, WireCell::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WireCell::ToyCTPointCloud& ct_point_cloud, int run_no, int subrun_no, int event_no, bool flag_data, bool debug_tagger){
  // 0 for original match results
  // 1 for TGM ...
  // 2 for STM ...
  int flag_background = 0;

  // bundle information ...
  int bundle_type = std::get<0>(bundle_info);
  double bundle_ks = std::get<1>(bundle_info);
  double bundle_chi2 = std::get<2>(bundle_info);
  int bundle_ndf = std::get<3>(bundle_info);
  //
  
  return std::make_tuple(0, main_cluster, main_flash);
}
