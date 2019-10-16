#include "WireCellPID/ToyFiducial.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

using namespace WireCell;

int pnpoly(std::vector<double>& vertx, std::vector<double>& verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = int(vertx.size())-1; i < int(vertx.size()); j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}


WireCellPID::ToyFiducial::ToyFiducial(int dead_region_ch_ext, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w, double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode, int flag_data)
  : dead_region_ch_ext(dead_region_ch_ext)
  , offset_t(offset_t)
  , offset_u(offset_u)
  , offset_v(offset_v)
  , offset_w(offset_w)
  , slope_t(slope_t)
  , slope_u(slope_u)
  , slope_v(slope_v)
  , slope_w(slope_w)
  , angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
  , m_top(top)
  , m_bottom(bottom)
  , m_upstream(upstream)
  , m_downstream(downstream)
  , m_anode(anode)
  , m_cathode(cathode)
{

  file = new TFile("./input_data_files/stopping_ave_dQ_dx.root");
  g_muon = (TGraph*)file->Get("muon");
  g_pion = (TGraph*)file->Get("pion");
  g_kaon = (TGraph*)file->Get("kaon");
  g_proton = (TGraph*)file->Get("proton");
  g_electron = (TGraph*)file->Get("electron");
  
  
  if (flag_data){
    std::cout << "Data Reco Fiducial Volume! " << std::endl;
    // data 
    m_sc_bottom_1_y=-116*units::cm;
    m_sc_bottom_1_x=80*units::cm;
    
    m_sc_bottom_2_y=-99*units::cm;
    m_sc_bottom_2_x=256*units::cm;
    
    m_sc_top_1_y = 116*units::cm; // used to be 118 cm
    m_sc_top_1_x = 100*units::cm;
    
    m_sc_top_2_y = 102*units::cm; // used to be 103 cm
    m_sc_top_2_x = 256*units::cm;
    
    m_sc_upstream_1_z = 0*units::cm;
    m_sc_upstream_1_x = 120*units::cm;
    
    m_sc_upstream_2_z = 11*units::cm;
    m_sc_upstream_2_x = 256*units::cm;
    
    m_sc_downstream_1_z=1037*units::cm;
    m_sc_downstream_1_x=120*units::cm;
    
    m_sc_downstream_2_z=1026*units::cm;
    m_sc_downstream_2_x=256*units::cm;
  }else{
    // MC
    std::cout << "MC Truth Fiducial Volume! " << std::endl;
    m_sc_bottom_1_y=-116*units::cm;
    m_sc_bottom_1_x=34*units::cm;
    
    m_sc_bottom_2_y=-98*units::cm;
    m_sc_bottom_2_x=256*units::cm;
    
    m_sc_top_1_y = 116*units::cm;
    m_sc_top_1_x = 70*units::cm;
    
    m_sc_top_2_y = 100*units::cm;
    m_sc_top_2_x = 256*units::cm;
    
    m_sc_upstream_1_z = 0*units::cm;
    m_sc_upstream_1_x = 50*units::cm;
    
    m_sc_upstream_2_z = 14*units::cm;
    m_sc_upstream_2_x = 256*units::cm;
    
    m_sc_downstream_1_z=1037*units::cm;
    m_sc_downstream_1_x=40*units::cm;
    
    m_sc_downstream_2_z=1023*units::cm;
    m_sc_downstream_2_x=256*units::cm;
  }
  

  //
  boundary_xy_x.clear(); boundary_xy_y.clear();
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_1_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_2_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_2_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_1_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_top - boundary_dis_cut);
  // boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  
  // for (size_t i=0;i!=boundary_xy_x.size();i++){
  //   std::cout << boundary_xy_x.at(i)/units::cm << " XY " << boundary_xy_y.at(i)/units::cm << std::endl;
  // }
  
  boundary_xz_x.clear(); boundary_xz_z.clear();
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_1_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_2_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_2_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_1_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_downstream - boundary_dis_cut-1*units::cm);
  // boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+2*units::cm);

  // for (size_t i=0;i!=boundary_xz_x.size();i++){
  //   std::cout << boundary_xz_x.at(i)/units::cm << " XZ " << boundary_xz_z.at(i)/units::cm << std::endl;
  // }
}


bool WireCellPID::ToyFiducial::check_stm(WireCellPID::PR3DCluster* main_cluster, double offset_x, double flash_time, WireCell::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map){

  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  Vector main_dir = main_cluster->get_PCA_axis(0);
  TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
  std::vector<WCPointCloud<double>::WCPoint> candidate_exit_wcps;
  std::set<int> temp_set;
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps;
  
  // first round check
  {
    // get two extreme points ...
    wcps = main_cluster->get_two_boundary_wcps(2,true);
    // figure out the end points ...
    std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();
    
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.first);
      out_vec_wcps.push_back(temp_wcps);
    }
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.second);
      out_vec_wcps.push_back(temp_wcps);
    }
    
    // boundary check
    for (size_t i=0;i!=out_vec_wcps.size();i++){
      bool flag_save = false;
      // check all the points ... 
      for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
	Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
	if (!inside_fiducial_volume(p1,offset_x)){
	  candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  flag_save = true;
	  break;
	}
      }
      
      if (!flag_save){
	// check direction
	Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
	TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
	dir *= (-1);
	
	// check U and V and W
	TVector3 dir_1(0,dir.Y(),dir.Z());
	double angle1 = dir_1.Angle(U_dir);
	TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;
	
	double angle2 = dir_1.Angle(V_dir);
	TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
	
	double angle3 = dir_1.Angle(W_dir);
	TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;
	
	
	if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	  if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
	    flag_save = true;
	    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  }
	}
	
	if (!flag_save){
	  if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
	    if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	      flag_save = true;
	      candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	    }
	  }
	}
      }
    }
    
    //std::cout << wcps.first.x << " " <<wcps.first.y << " " << wcps.first.z << std::endl;
    //std::cout << wcps.second.x << " " <<wcps.second.y << " " << wcps.second.z << std::endl;
    
    for (size_t i=0; i!=candidate_exit_wcps.size(); i++){
      double dis1 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.first.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.first.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.first.z,2));
      double dis2 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.second.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.second.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.second.z,2));
      
      //      std::cout << candidate_exit_wcps.at(i).x << " " << candidate_exit_wcps.at(i).y << " " << candidate_exit_wcps.at(i).z << " " << dis1 << " " << dis2 << std::endl;
      
      // essentially one of the extreme points ...
      if (dis1 < dis2){
	if (dis1 < 1.0*units::cm)  temp_set.insert(0);
      }else{
	if (dis2 < 1.0*units::cm)  temp_set.insert(1);
      }
    }

    //    std::cout << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;
    
    // protection against two end point situation
    if (temp_set.size()==2){
      WireCell::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WireCell::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
      temp_set.clear();
      
      if ((!inside_fiducial_volume(tp1,offset_x))) temp_set.insert(0);
      if ((!inside_fiducial_volume(tp2,offset_x))) temp_set.insert(1);
      if (temp_set.size()==0){
	temp_set.insert(0);
	temp_set.insert(1);
      }
    }
  }

  if (temp_set.size()==0){
    candidate_exit_wcps.clear();
       // get two extreme points ...
    wcps = main_cluster->get_two_boundary_wcps(2);
    // figure out the end points ...
    std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();
    
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.first);
      out_vec_wcps.push_back(temp_wcps);
    }
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.second);
      out_vec_wcps.push_back(temp_wcps);
    }
    
    // boundary check
    for (size_t i=0;i!=out_vec_wcps.size();i++){
      bool flag_save = false;
      // check all the points ... 
      for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
	Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
	if (!inside_fiducial_volume(p1,offset_x)){
	  candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  flag_save = true;
	  break;
	}
      }
      
      if (!flag_save){
	// check direction
	Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
	TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
	dir *= (-1);
	
	// check U and V and W
	TVector3 dir_1(0,dir.Y(),dir.Z());
	double angle1 = dir_1.Angle(U_dir);
	TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;
	
	double angle2 = dir_1.Angle(V_dir);
	TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
	
	double angle3 = dir_1.Angle(W_dir);
	TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;
	
	
	if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	  if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
	    flag_save = true;
	    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  }
	}
	
	if (!flag_save){
	  if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
	    if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	      flag_save = true;
	      candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	    }
	  }
	}
      }
    }
    
    //std::cout << wcps.first.x << " " <<wcps.first.y << " " << wcps.first.z << std::endl;
    //std::cout << wcps.second.x << " " <<wcps.second.y << " " << wcps.second.z << std::endl;
    
    for (size_t i=0; i!=candidate_exit_wcps.size(); i++){
      double dis1 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.first.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.first.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.first.z,2));
      double dis2 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.second.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.second.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.second.z,2));
      
      //std::cout << candidate_exit_wcps.at(i).x << " " << candidate_exit_wcps.at(i).y << " " << candidate_exit_wcps.at(i).z << " " << dis1 << " " << dis2 << std::endl;
      
      // essentially one of the extreme points ...
      if (dis1 < dis2){
	if (dis1 < 1.0*units::cm)  temp_set.insert(0);
      }else{
	if (dis2 < 1.0*units::cm)  temp_set.insert(1);
      }
    }

    //std::cout << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;
    
    
    // protection against two end point situation
    if (temp_set.size()==2){
      WireCell::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WireCell::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
      temp_set.clear();
      
      if ((!inside_fiducial_volume(tp1,offset_x))) temp_set.insert(0);
      if ((!inside_fiducial_volume(tp2,offset_x))) temp_set.insert(1);
      if (temp_set.size()==0){
	temp_set.insert(0);
	temp_set.insert(1);
      }
    }
  }
  

  
  // fully contained, so not a STM
  if (candidate_exit_wcps.size()==0) {
    std::cout << "Mid Point: " << std::endl;
    return false;
  }
  
  std::cout << "end_point: " << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;

  // It is possible that we have two points out of fiducial
  // Michel electron is outside the boundary ...

  
  
  // Crawl backward according to the graph ??? ...
  WCPointCloud<double>::WCPoint first_wcp;
  WCPointCloud<double>::WCPoint last_wcp;
  bool flag_double_end = false;

  if (temp_set.size()!=0){
    if (*temp_set.begin()==0){
      first_wcp = wcps.first;
      last_wcp = wcps.second;
    }else{
      first_wcp = wcps.second;
      last_wcp = wcps.first;
    }
    if (temp_set.size()==2) flag_double_end = true;
    
  }else{
    if (candidate_exit_wcps.size()==1){
      first_wcp = candidate_exit_wcps.at(0);
      TVector3 dir1(wcps.first.x - candidate_exit_wcps.at(0).x,
		    wcps.first.y - candidate_exit_wcps.at(0).y,
		    wcps.first.z - candidate_exit_wcps.at(0).z);
      TVector3 dir2(wcps.second.x - candidate_exit_wcps.at(0).x,
		    wcps.second.y - candidate_exit_wcps.at(0).y,
		    wcps.second.z - candidate_exit_wcps.at(0).z);
      double dis1 = dir1.Mag();
      double dis2 = dir2.Mag();

      if (dir1.Angle(dir2) > 120/180.*3.1415926 && dis1 > 20*units::cm &&
      	  dis2 > 20*units::cm){
	std::cout << "Mid Point: " << std::endl;
	return false;
      }else{
	if (dis1 < dis2){
	  last_wcp = wcps.second;	   
	}else{
	  last_wcp = wcps.first; 
	}
      }
      
    }else{
      std::cout << "Mid Point: " << std::endl;
      return false;
    }
  }
  
  // forward check ...
  {
    if (flag_double_end) std::cout << "Forward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path(first_wcp, last_wcp);
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    Point mid_p = main_cluster->adjust_rough_path(); 
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

    // check both end points for TGM ...
    WireCell::PointVector& pts = main_cluster->get_fine_tracking_path();
    std::vector<double>& dQ = main_cluster->get_dQ();
    std::vector<double>& dx = main_cluster->get_dx();
    
    int kink_num = find_first_kink(main_cluster);
    
    double left_L = 0;
    double left_Q = 0;
    double exit_L = 0;
    double exit_Q = 0;
    for (size_t i=0;i!=kink_num;i++){
      exit_L += dx.at(i);
      exit_Q += dQ.at(i);
    }
    for (size_t i = kink_num; i!=dx.size(); i++){
      left_L += dx.at(i);
      left_Q += dQ.at(i);
    }
    
    std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;
    
    if ((exit_L < 3*units::cm || left_L < 3*units::cm) && (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
      return true;
    }
  
    if (left_L > 40*units::cm || left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0){
      if (!flag_double_end){
	std::cout << "Mid Point A " << inside_dead_region(mid_p) << " " << mid_p << " " << left_L  << " " << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
	return false;
      }
    }else{
      bool flag_fix_end = false;
      if (exit_L < 35*units::cm || (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0&& left_L > 2*units::cm) flag_fix_end = true;
      
      if (left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5 ||
	  left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }
    
    
      bool flag_pass = false;
      if (left_L < 40*units::cm) {
	if (flag_fix_end){
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	    eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	}else{
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	    eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	}
	
	if (flag_pass)
	  return true;
	else{
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
	  }
	}
      }
      
      if (left_L < 20*units::cm){
	if (flag_pass)
	  return true;
	else{
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	}
	
	if (flag_pass)
	  return true;
	else{
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	  }
	}
      }
    }
  }
  
  // backward check ...
  if (flag_double_end){
    std::cout << "Backward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path( last_wcp, first_wcp);
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    Point mid_p = main_cluster->adjust_rough_path(); 
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

    // check both end points for TGM ...
    WireCell::PointVector& pts = main_cluster->get_fine_tracking_path();
    std::vector<double>& dQ = main_cluster->get_dQ();
    std::vector<double>& dx = main_cluster->get_dx();
    
    int kink_num = find_first_kink(main_cluster);
    
    double left_L = 0;
    double left_Q = 0;
    double exit_L = 0;
    double exit_Q = 0;
    for (size_t i=0;i!=kink_num;i++){
      exit_L += dx.at(i);
      exit_Q += dQ.at(i);
    }
    for (size_t i = kink_num; i!=dx.size(); i++){
      left_L += dx.at(i);
      left_Q += dQ.at(i);
    }
    
    std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;
    
    if ((exit_L < 3*units::cm || left_L < 3*units::cm) && (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
      return true;
    }
  
    if (left_L > 40*units::cm || left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0){
      std::cout << "Mid Point A " << inside_dead_region(mid_p) << " " << mid_p << " " << left_L  << " " << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
      return false;
    }else{
      bool flag_fix_end = false;
      if (exit_L < 35*units::cm || (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0 && left_L > 2*units::cm) flag_fix_end = true;
      
      if (left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5 ||
	  left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }

    
    
      bool flag_pass = false;
      if (left_L < 40*units::cm) {
	if (flag_fix_end){
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	    eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	}else{
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	    eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	}
	
	if (flag_pass)
	  return true;
	else{
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
	  }
	}
      }
      
      if (left_L < 20*units::cm){
	if (flag_pass)
	  return true;
	else{
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	}
	
	if (flag_pass)
	  return true;
	else{
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
	      eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	  }
	}
      }
    }
  }
  
  // check 5512-209-10491
  // if (inside_dead_region(mid_p)) return true;
  //  

  // end check ...
  std::cout << "Mid Point " << std::endl;  
  return false;
}

int WireCellPID::ToyFiducial::find_first_kink(WireCellPID::PR3DCluster* main_cluster){
  WireCell::PointVector& fine_tracking_path = main_cluster->get_fine_tracking_path();
  std::vector<double>& dQ = main_cluster->get_dQ();
  std::vector<double>& dx = main_cluster->get_dx();
  
  TVector3 drift_dir(1,0,0);

  std::vector<double> refl_angles(fine_tracking_path.size(),0);
  std::vector<double> para_angles(fine_tracking_path.size(),0);
  std::vector<double> ave_angles(fine_tracking_path.size(),0);
  std::vector<int> max_numbers(fine_tracking_path.size(),-1);
  
  for (size_t i=0;i!=fine_tracking_path.size(); i++){
    double angle1 = 0;
    double angle2 = 0;
    for (int j=0;j!=6;j++){    
      TVector3 v10(0,0,0);
      TVector3 v20(0,0,0);
      if (i>j)
	v10.SetXYZ(fine_tracking_path.at(i).x - fine_tracking_path.at(i-j-1).x,
		   fine_tracking_path.at(i).y - fine_tracking_path.at(i-j-1).y,
		   fine_tracking_path.at(i).z - fine_tracking_path.at(i-j-1).z);
      
      if (i+j+1<fine_tracking_path.size())
	v20.SetXYZ(fine_tracking_path.at(i+j+1).x - fine_tracking_path.at(i).x,
		   fine_tracking_path.at(i+j+1).y - fine_tracking_path.at(i).y,
		   fine_tracking_path.at(i+j+1).z - fine_tracking_path.at(i).z);
      
      if (j==0){
	angle1 = v10.Angle(v20)/3.1415926*180.;
	angle2 = std::max(fabs(v10.Angle(drift_dir)/3.1415926*180.-90.),
			  fabs(v20.Angle(drift_dir)/3.1415926*180.-90.));
      }else{
	if (v10.Mag()!=0 && v20.Mag()!=0){
	  angle1 = std::min(v10.Angle(v20)/3.1415926*180., angle1);
	  angle2 = std::min(std::max(fabs(v10.Angle(drift_dir)/3.1415926*180.-90.),
				     fabs(v20.Angle(drift_dir)/3.1415926*180.-90.)),angle2);
	}
      }
    }

    refl_angles.at(i) = angle1;
    para_angles.at(i) = angle2;
  }

  for (int i=0;i!=fine_tracking_path.size();i++){
    double sum_angles = 0;
    double nsum = 0;
    double max_angle = 0;
    int max_num = -1;
    
    for (int j = -2; j!=3;j++){

      if (i+j>=0 && i+j<fine_tracking_path.size()){
	if (para_angles.at(i+j)>12){
	  sum_angles += pow(refl_angles.at(i+j),2);
	  nsum ++;
	  if (refl_angles.at(i+j) > max_angle){
	    max_angle = refl_angles.at(i+j);
	    max_num = i+j;
	  }
	}
      }
    }

    if (nsum!=0) sum_angles=sqrt(sum_angles/nsum);
    ave_angles.at(i) = sum_angles;
    max_numbers.at(i) = max_num;

  }
    
  for (int i=0;i!=fine_tracking_path.size();i++){
    // std::cout << i << " " << refl_angles.at(i) << " " <<sum_angles << " " << inside_fiducial_volume(fine_tracking_path.at(i)) << std::endl;
    if ((refl_angles.at(i) > 25 && ave_angles.at(i) > 12.5 ) && inside_fiducial_volume(fine_tracking_path.at(i))){
      TVector3 v10(fine_tracking_path.at(i).x - fine_tracking_path.front().x,
		   fine_tracking_path.at(i).y - fine_tracking_path.front().y,
		   fine_tracking_path.at(i).z - fine_tracking_path.front().z);
      TVector3 v20(fine_tracking_path.back().x - fine_tracking_path.at(i).x,
		   fine_tracking_path.back().y - fine_tracking_path.at(i).y,
		   fine_tracking_path.back().z - fine_tracking_path.at(i).z);
      double angle3 = v10.Angle(v20)/3.1415926*180.;
      if (angle3 < 20 && ave_angles.at(i) < 20 || angle3 < 12.5 && inside_dead_region(fine_tracking_path.at(i)) || angle3 < 7.5 || i<=4) continue;
      if (angle3 > 30){
	double sum_fQ = 0;
	double sum_fx = 0;
	double sum_bQ = 0;
	double sum_bx = 0;
	for (int k=0;k!=10;k++){
	  if (i>=k+1){
	    sum_fQ += dQ.at(i-k-1);
	    sum_fx += dx.at(i-k-1);
	  }
	  if (i+k+1 < dQ.size()){
	    sum_bQ += dQ.at(i+k+1);
	    sum_bx += dx.at(i+k+1);
	  }
	}
	sum_fQ /= (sum_fx/units::cm+1e-9)*50e3;
	sum_bQ /= (sum_bx/units::cm+1e-9)*50e3;
	//std::cout << sum_fQ << " " << sum_bQ << std::endl;
	if (sum_fQ > 0.6 && sum_bQ > 0.6){
	  std::cout << "Kink: " << i << " " << refl_angles.at(i) << " " << para_angles.at(i) << " " << ave_angles.at(i) << " " << max_numbers.at(i) << " " << angle3 << " " << dQ.at(i)/dx.at(i)*units::cm/50e3 << std::endl;
	  return max_numbers.at(i);
	}
      }
    }
  }

  for (int i=0;i!=fine_tracking_path.size();i++){
    // std::cout << i << " " << refl_angles.at(i) << " " <<sum_angles << " " << inside_fiducial_volume(fine_tracking_path.at(i)) << std::endl;
    if ((refl_angles.at(i) > 20 && ave_angles.at(i) > 15 ) && inside_fiducial_volume(fine_tracking_path.at(i))){
      TVector3 v10(fine_tracking_path.at(i).x - fine_tracking_path.front().x,
  		   fine_tracking_path.at(i).y - fine_tracking_path.front().y,
  		   fine_tracking_path.at(i).z - fine_tracking_path.front().z);
      TVector3 v20(fine_tracking_path.back().x - fine_tracking_path.at(i).x,
  		   fine_tracking_path.back().y - fine_tracking_path.at(i).y,
  		   fine_tracking_path.back().z - fine_tracking_path.at(i).z);
      double angle3 = v10.Angle(v20)/3.1415926*180.;
      if (angle3 < 20 && ave_angles.at(i) < 20 || angle3 < 12.5 && inside_dead_region(fine_tracking_path.at(i)) || angle3 < 7.5 || i<=4) continue;
      if (angle3 > 30){
	double sum_fQ = 0;
	double sum_fx = 0;
	double sum_bQ = 0;
	double sum_bx = 0;
	for (int k=0;k!=10;k++){
	  if (i>=k+1){
	    sum_fQ += dQ.at(i-k-1);
	    sum_fx += dx.at(i-k-1);
	  }
	  if (i+k+1 < dQ.size()){
	    sum_bQ += dQ.at(i+k+1);
	    sum_bx += dx.at(i+k+1);
	  }
	}
	sum_fQ /= (sum_fx/units::cm+1e-9)*50e3;
	sum_bQ /= (sum_bx/units::cm+1e-9)*50e3;
	//std::cout << sum_fQ << " " << sum_bQ << std::endl;
	if (sum_fQ > 0.6 && sum_bQ > 0.6){
	  std::cout << "Kink: " << i << " " << refl_angles.at(i) << " " << para_angles.at(i) << " " << ave_angles.at(i) << " " << max_numbers.at(i) << " " << angle3 << " " << dQ.at(i)/dx.at(i)*units::cm/50e3 << std::endl;
	  return max_numbers.at(i);
	}
      }
    }
  }

  
  
  return fine_tracking_path.size();
}

bool WireCellPID::ToyFiducial::eval_stm(WireCellPID::PR3DCluster* main_cluster,int kink_num,double peak_range, double offset_length, double com_range){
  WireCell::PointVector& pts = main_cluster->get_fine_tracking_path();

  std::vector<double>& dQ = main_cluster->get_dQ();
  std::vector<double>& dx = main_cluster->get_dx();
  std::vector<double> L(pts.size(),0);
  std::vector<double> dQ_dx(pts.size(),0);
  double dis = 0;
  L.at(0) = dis;
  dQ_dx.at(0) = dQ.at(0)/(dx.at(0)/units::cm+1e-9);
  for (size_t i=1;i!=pts.size();i++){
    dis += sqrt(pow(pts.at(i).x-pts.at(i-1).x,2) + pow(pts.at(i).y-pts.at(i-1).y,2) + pow(pts.at(i).z - pts.at(i-1).z,2));
    L.at(i) = dis;
    dQ_dx.at(i) = dQ.at(i)/(dx.at(i)/units::cm+1e-9);
  }
  //std::cout << L.size() << " " << dQ.size() << " " << dx.size() << std::endl;
  double end_L; 
  double max_num; 
  if (kink_num == pts.size()){
    end_L = L.back();
    max_num = L.size();
  }else{
    end_L = L.at(kink_num)-0.5*units::cm;
    max_num = kink_num;
  }
  
  double max_bin = -1;
  double max_sum = 0;
  for (size_t i=0;i!=L.size();i++){
    double sum = 0;
    double nsum = 0;
    double temp_max_bin = i;
    double temp_max_val = dQ_dx.at(i);
    if (L.at(i) < end_L + 0.5*units::cm && L.at(i) > end_L - peak_range && i < max_num){
      sum += dQ_dx.at(i); nsum ++;
      if (i>=2){
	sum += dQ_dx.at(i-2); nsum++;
	if (dQ_dx.at(i-2) > temp_max_val && i-2 < max_num){
	  temp_max_val = dQ_dx.at(i-2);
	  temp_max_bin = i-2;
	}
      }
      if (i>=1){
	sum += dQ_dx.at(i-1); nsum++;
	if (dQ_dx.at(i-1) > temp_max_val && i-1 < max_num){
	  temp_max_val = dQ_dx.at(i-1);
	  temp_max_bin = i-1;
	}
      }
      if (i+1<L.size()){
	sum += dQ_dx.at(i+1); nsum++;
	if (dQ_dx.at(i+1) > temp_max_val && i+1 < max_num){
	  temp_max_val = dQ_dx.at(i+1);
	  temp_max_bin = i+1;
	}
      }
      if (i+2<L.size()){
	sum += dQ_dx.at(i+2); nsum++;
	if (dQ_dx.at(i+2) > temp_max_val && i+2 < max_num){
	  temp_max_val = dQ_dx.at(i+2);
	  temp_max_bin = i+2;
	}
      }
      sum /= nsum;
      if (sum>max_sum){
	max_sum = sum;
	max_bin = temp_max_bin;
      }
    }
  }
  //std::cout << max_bin << " " << max_sum << std::endl;
  end_L = L.at(max_bin)+0.2*units::cm;
  int ncount = 0;
  std::vector<double> vec_x;
  std::vector<double> vec_y;
  for (size_t i=0;i!=L.size(); i++){
    if (end_L - L.at(i) < com_range && end_L - L.at(i) > 0){
      vec_x.push_back(end_L-L.at(i));
      vec_y.push_back(dQ_dx.at(i));
      ncount ++;
    }
  }

  TH1F *h1 = new TH1F("h1","h1",ncount,0,ncount);
  TH1F *h2 = new TH1F("h2","h2",ncount,0,ncount);
  TH1F *h3 = new TH1F("h3","h3",ncount,0,ncount);

  for (size_t i=0;i!=ncount;i++){
    // std::cout << i << " " << vec_y.at(i) << std::endl;
    h1->SetBinContent(i+1,vec_y.at(i));
    h2->SetBinContent(i+1,g_muon->Eval((vec_x.at(i)+offset_length)/units::cm));
    h3->SetBinContent(i+1,50e3);
  }
  double ks1 = h2->KolmogorovTest(h1,"M");
  double ratio1 = h2->GetSum()/(h1->GetSum()+1e-9);
  double ks2 = h3->KolmogorovTest(h1,"M");
  double ratio2 = h3->GetSum()/(h1->GetSum()+1e-9);
  
  delete h1;
  delete h2;
  delete h3;

  std::cout << "KS value: " << ks1 << " " << ks2 << " " << ratio1 << " " << ratio2 << " " << ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << std::endl;
  
  // std::vector<double> results;
  // results.push_back(ks1);
  // results.push_back(ks2);
  // results.push_back(ratio1);
  // results.push_back(ratio2);
  // return results;

  if (ks1-ks2 >= 0.0) return false;
  if (sqrt(pow(ks2/0.06,2)+pow((ratio2-1)/0.06,2))< 1.4) return false;

  if (ks1 - ks2 < -0.02 && (ks2 > 0.09 || ratio2 > 1.5)) return true;
  if ( ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 < 0) return true;
 

  return false;
  
}

// bool WireCellPID::ToyFiducial::check_fully_contained(WireCell::FlashTPCBundle *bundle, double offset_x, WireCell::ToyCTPointCloud& ct_point_cloud,std::map<PR3DCluster*, PR3DCluster*>& old_new_cluster_map, unsigned int* fail_mode, int flag){
//   PR3DCluster *main_cluster; 
//   PR3DCluster *main_cluster1; 
//   if (flag==1){ // check the current
//     main_cluster = bundle->get_main_cluster();
//     main_cluster1  = main_cluster;
//     //replace it with the better one, which takes into account the dead channels ... 
//     if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
//       main_cluster1 = old_new_cluster_map[main_cluster];
//   }else{
//     main_cluster = bundle->get_orig_cluster();
//     main_cluster1  = main_cluster;
//     //replace it with the better one, which takes into account the dead channels ... 
//     if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
//       main_cluster1 = old_new_cluster_map[main_cluster];
//     if (main_cluster==0) return true;
//   }
  
  
//   Opflash *flash = bundle->get_flash();
//   std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();

//   TVector3 drift_dir(1,0,0);
//   // hard coded for U and V plane ... 
//   TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
//   TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
//   TVector3 W_dir(0,1,0);

//   Vector main_dir = main_cluster->get_PCA_axis(0);
//   TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
//   for (size_t i=0;i!=out_vec_wcps.size();i++){
//     // check all the points ... 
//     for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
//       Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
//       if (!inside_fiducial_volume(p1,offset_x)){
//         if(fail_mode) *fail_mode |= 1U<<2;
//         return false;
//       }
//     }


//     Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
//     TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
//     dir *= (-1);

//     // check U and V and W
//     TVector3 dir_1(0,dir.Y(),dir.Z());
//     double angle1 = dir_1.Angle(U_dir);
//     TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
//     double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;
    
//     double angle2 = dir_1.Angle(V_dir);
//     TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
//     double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
    
//     double angle3 = dir_1.Angle(W_dir);
//     TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
//     double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;
    
//     // not added for now, need to check to add this one in when more events are available ...
//     // XQ, 7/11/2018
//     // double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


//     //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

//     if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
//       if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
//         if(fail_mode) *fail_mode |= 1U<<1;
//         return false;
//       }
//     }
	    
//     if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
//       if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
//         if(fail_mode) *fail_mode |= 1U;
//         return false;
//       }
//     }
    
//   }



  
//   return true;
// }



WireCellPID::ToyFiducial::~ToyFiducial(){
  delete g_muon;
  delete g_pion;
  delete g_proton;
  delete g_kaon;
  delete g_electron;
  delete file;
		 
}

bool WireCellPID::ToyFiducial::check_signal_processing(WireCell::Point& p, TVector3& dir, WireCell::ToyCTPointCloud& ct_point_cloud, double step, double offset_x){

  if (dir.Mag()==0){
    return true;
  }else{
    dir *= 1./dir.Mag();
    Point temp_p = p;

    int num_points = 0;
    int num_points_dead = 0;

    //  std::cerr << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << std::endl;
    
    while(inside_fiducial_volume(temp_p,offset_x)){
      num_points ++;
      //if (inside_dead_region(temp_p))
      //	num_points_dead ++;

      //      std::cerr << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << " ";
      
      WireCell::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,0);
      WireCell::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,1);
      WireCell::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,2);
      
      //      std::cerr << cloud_u.pts.size() << " " << cloud_v.pts.size() << " " << cloud_w.pts.size() << std::endl;

      if (cloud_u.pts.size()>0 || cloud_v.pts.size()>0 || cloud_w.pts.size() > 0 || inside_dead_region(temp_p))
      	num_points_dead++;
      
      if (num_points - num_points_dead >=5) return true;
	
      temp_p.x += dir.X() * step;
      temp_p.y += dir.Y() * step;
      temp_p.z += dir.Z() * step;
    }

    //    std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
    
    if (num_points_dead > 0.8*num_points){
    	return false;
    }else{
    	return true;
    }
  }
  
  return true;
}

bool WireCellPID::ToyFiducial::check_dead_volume(WireCell::Point& p, TVector3& dir, double step, double offset_x){
  if (!inside_fiducial_volume(p,offset_x)){
    return false;
  }else{
    if (dir.Mag()==0){
      return true;
    }else{
      dir *= 1./dir.Mag();
      Point temp_p = p;
      int num_points = 0;
      int num_points_dead = 0;
      while(inside_fiducial_volume(temp_p,offset_x)){

	num_points ++;
	if (inside_dead_region(temp_p))
	  num_points_dead ++;

	if (num_points - num_points_dead >=4) return true;
	
	temp_p.x += dir.X() * step;
	temp_p.y += dir.Y() * step;
	temp_p.z += dir.Z() * step;

	//	std::cout << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << std::endl;
      }

      // std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
      
      if (num_points_dead > 0.81*num_points){
	return false;
      }else{
	return true;
      }
      
    }
  }
}


bool WireCellPID::ToyFiducial::inside_fiducial_volume(WireCell::Point& p, double offset_x){

  int c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x-offset_x, p.y);
  int c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x-offset_x, p.z);

  //  std::cout << (p.x-offset_x)/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
  //std::cout << c1 << " " << c2 << std::endl;
  
  if (c1 && c2){
    return true;
  }else{
    return false;
  }
}

bool WireCellPID::ToyFiducial::inside_dead_region(WireCell::Point& p){
  // convert the position into U, V, W, and T number ...
  int time_slice = p.x * slope_t + offset_t;
  double pos_u = cos(angle_u) * p.z - sin(angle_u) *p.y;
  double pos_v = cos(angle_v) * p.z - sin(angle_v) *p.y;
  double pos_w = cos(angle_w) * p.z - sin(angle_w) *p.y;
  int ch_u = pos_u * slope_u + offset_u;
  int ch_v = pos_v * slope_v + offset_v + 2400;
  int ch_w = pos_w * slope_w + offset_w + 4800;

  //std::cout << ch_u << " " << ch_v << " " << ch_w << " " << time_slice << std::endl;
  //  std::cout << slope_w << " " << offset_w << " " << pos_w << std::endl;
  
  if (time_slice <0 || time_slice >=2398) return false;
  if (ch_u <0 || ch_u>=2400)  return false;
  if (ch_v <2400 || ch_v>=4800)  return false;
  if (ch_w <4800 || ch_w>=8256)  return false;

  std::set<SlimMergeGeomCell*> dead_u_mcells;
  std::set<SlimMergeGeomCell*> dead_v_mcells;
  std::set<SlimMergeGeomCell*> dead_w_mcells;

  if (ch_mcell_set_map.find(ch_u)!=ch_mcell_set_map.end())
    dead_u_mcells = ch_mcell_set_map[ch_u];
  if (ch_mcell_set_map.find(ch_v)!=ch_mcell_set_map.end())
    dead_v_mcells = ch_mcell_set_map[ch_v];
  if (ch_mcell_set_map.find(ch_w)!=ch_mcell_set_map.end())
    dead_w_mcells = ch_mcell_set_map[ch_w];
  
  // std::cout << ch_u << " " << ch_v << " " << ch578d-_w << " " << dead_u_mcells.size() << " " << dead_v_mcells.size() << " " << dead_w_mcells.size() << std::endl;
  
  // find the dead region given the U, V, and W number
  std::set<SlimMergeGeomCell*> results;
  for (auto it = dead_u_mcells.begin(); it!=dead_u_mcells.end(); it++){
    if (dead_v_mcells.find(*it)!=dead_v_mcells.end()){
      // compare UV sets
      results.insert(*it);
    }else if (dead_w_mcells.find(*it)!=dead_w_mcells.end()){
      // compare UW sets
      results.insert(*it);
    }
  }
  // compare VW sets ...
  for (auto it = dead_v_mcells.begin(); it!=dead_v_mcells.end(); it++){
    if (dead_w_mcells.find(*it)!=dead_w_mcells.end()){
      // compare UW sets
      results.insert(*it);
    }
  }
  
  
  // Check the T number for the remaining T ... 
  for (auto it = results.begin(); it!=results.end(); it++){
    if (mcell_time_map[*it].first <= time_slice &&
	time_slice <= mcell_time_map[*it].second)
      return true;
  }
  
  return false;
}


void WireCellPID::ToyFiducial::AddDeadRegion(WireCell::SlimMergeGeomCell* mcell, std::vector<int>& time_slices){

  mcells.push_back(mcell);
  int start_time = time_slices.front() - dead_region_ch_ext ;
  int end_time = time_slices.back() + dead_region_ch_ext;
  mcell_time_map[mcell] = std::make_pair(start_time, end_time);

  GeomWireSelection& uwires = mcell->get_uwires();
  GeomWireSelection& vwires = mcell->get_vwires();
  GeomWireSelection& wwires = mcell->get_wwires();

  // std::cout << uwires.size() << " " << vwires.size() << " " << wwires.size() << " " << start_time << " " << end_time << std::endl;
  
  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();

  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(0))!=bad_planes.end()){
    int start_ch = uwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <0) start_ch = 0;
    int end_ch = uwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=2400) end_ch = 2399;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(1))!=bad_planes.end()){
    int start_ch = vwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <2400) start_ch = 2400;
    int end_ch = vwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=4800) end_ch = 4799;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(2))!=bad_planes.end()){
    int start_ch = wwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <4800) start_ch = 4800;
    int end_ch = wwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=8256) end_ch = 8255;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  
  
}


// bool WireCellPID::ToyFiducial::check_tgm(WireCell::FlashTPCBundle *bundle, double offset_x, WireCell::ToyCTPointCloud& ct_point_cloud,std::map<PR3DCluster*, PR3DCluster*>& old_new_cluster_map, int flag){

//   PR3DCluster *main_cluster; 
//   PR3DCluster *main_cluster1; 

//   if (flag==1){ // check the current main cluster 
//     main_cluster = bundle->get_main_cluster();
//     main_cluster1 = main_cluster;
//   //replace it with the better one, which takes into account the dead channels ... 
//     if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
//       main_cluster1 = old_new_cluster_map[main_cluster];
//   }else if (flag==2){
//     main_cluster = bundle->get_orig_cluster();
//     main_cluster1 = main_cluster;
//   //replace it with the better one, which takes into account the dead channels ... 
//     if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
//       main_cluster1 = old_new_cluster_map[main_cluster];
//     if (main_cluster==0) return false;
//   }

  
//   Opflash *flash = bundle->get_flash();

//   std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();

//   // if (main_cluster->get_cluster_id()==6)
//   //std::cout << main_cluster << " " << main_cluster1 << std::endl;
  

//   // int max_group = 0;
//   // int max_count = out_vec_wcps.at(max_group).size();

//   // for (size_t i=1; i!=out_vec_wcps.size();i++){
//   //   if (out_vec_wcps.at(i).size() > max_count){
//   //     max_group = i;
//   //     max_count = out_vec_wcps.at(max_group).size();
//   //   }
//   // }

//   TVector3 drift_dir(1,0,0);
//   // hard coded for U and V plane ... 
//   TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
//   TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
//   TVector3 W_dir(0,1,0);

//   double length_limit = sqrt(pow(out_vec_wcps.at(0).at(0).x-out_vec_wcps.at(1).at(0).x,2)+
// 			     pow(out_vec_wcps.at(0).at(0).y-out_vec_wcps.at(1).at(0).y,2)+
// 			     pow(out_vec_wcps.at(0).at(0).z-out_vec_wcps.at(1).at(0).z,2));
  
//   // std::cout << "Flash: " << flash->get_flash_id() << std::endl;

//   // take a look at the first point ...
//   for (size_t i=0;i!=out_vec_wcps.size();i++){
//     bool flag_p1_inside = true;
//     int p1_index = -1;
//     for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
//       Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
//       flag_p1_inside = flag_p1_inside && inside_fiducial_volume(p1,offset_x);
//       if (!flag_p1_inside){
// 	p1_index = j;
// 	break;
//       }
//     }
    
    
//     // loop through the remaining groups and check ...
//     for (size_t k=i+1;k!=out_vec_wcps.size();k++){
//       bool flag_p2_inside = true;
//       int p2_index = -1;
//       for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
// 	Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
// 	flag_p2_inside = flag_p2_inside && inside_fiducial_volume(p2,offset_x);
// 	if (!flag_p2_inside){
// 	  p2_index = j;
// 	  break;
// 	}
//       }
      

//       // if (main_cluster->get_cluster_id()==7){
//       //  	std::cout << main_cluster->get_cluster_id() << " " << i << " " <<
//       //  	  out_vec_wcps.at(i).at(0).x/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " << 
//       //  	  k << " " << out_vec_wcps.at(k).at(0).x/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm <<
//       //  	  " " << p1_index << " " << p2_index << " " << flag_p1_inside << " " << flag_p2_inside << std::endl;
//       // 	// for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
//       // 	//   Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
//       // 	//   std::cout << j << " A " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << std::endl;
//       // 	// }
//       // 	// for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
//       // 	//   Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
//       // 	//   std::cout << j << " B " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << std::endl;
//       // 	// }
//       // }
	
//       if ((!flag_p1_inside) && (!flag_p2_inside)){
// 	// if not a neutrino candidate ... to be worked out ...
// 	 // Point p1(out_vec_wcps.at(i).at(p1_index).x,out_vec_wcps.at(i).at(p1_index).y,out_vec_wcps.at(i).at(p1_index).z);
// 	 // Point p2(out_vec_wcps.at(k).at(p2_index).x,out_vec_wcps.at(k).at(p2_index).y,out_vec_wcps.at(k).at(p2_index).z);
// 	//	std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << inside_fiducial_volume(p1,offset_x) << " A " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << inside_fiducial_volume(p2,offset_x) << " " << offset_x/units::cm << std::endl;


// 	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(p1_index).x-offset_x)/units::cm <<3=r-t=-[\] " " << out_vec_wcps.at(i).at(p1_index).y/units::cm << " " << out_vec_wcps.at(i).at(p1_index).z/units::cm << " " ;
// 	// std::cout << (out_vec_wcps.at(k).at(p2_index).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(p2_index).y/units::cm << " " << out_vec_wcps.at(k).at(p2_index).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;

// 	//std::cout << p1_index << " " << p2_index << std::endl;
	
// 	// check two points in between
// 	bool flag_check = false;
// 	for (int kk=0;kk!=3;kk++){
// 	  Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
// 		   out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
// 		   out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
// 	  flag_check = flag_check || inside_fiducial_volume(p3,offset_x);
// 	}
	
// 	// if (main_cluster->get_cluster_id()==10)
// 	//   std::cout << flag_check << " " << out_vec_wcps.size() << std::endl;
	
// 	if (flag_check){
// 	  if (flash->get_type()==2){

// 	    double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
// 				      pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
// 				      pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));

// 	    //std::cout << temp_length/units::cm << " " << length_limit/units::cm << std::endl;

// 	    if (i==0&&k==1){
// 	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) &&
// 		   temp_length > 0.45*length_limit)
// 		return true;
// 	    }else{
// 	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) &&
// 		   temp_length > 0.45*length_limit)
// 		return true;
// 	    }
	    
// 	  }else{
// 	    return true; // through going muon ...
// 	  }
// 	}else{
	  
// 	  if (out_vec_wcps.size()==2){
// 	    return true;
// 	  }else{
// 	    // if (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x))
// 	    //   return true;

// 	    bool flag_check_again = false;
// 	    for (int kkk = 0;kkk!=out_vec_wcps.size(); kkk++){
// 	      if (kkk==i || kkk==k) continue;
// 	      for (int kk=0;kk!=4;kk++){
// 	    	Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(kkk).at(0).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
// 	    		 out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(kkk).at(0).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
// 	    		 out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(kkk).at(0).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
// 	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
// 	      }
	      
// 	      for (int kk=0;kk!=3;kk++){
// 	    	Point p3(out_vec_wcps.at(kkk).at(0).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(0).x)/4.*(kk+1),
// 	    		 out_vec_wcps.at(kkk).at(0).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(0).y)/4.*(kk+1),
// 	    		 out_vec_wcps.at(kkk).at(0).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(0).z)/4.*(kk+1));
// 	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
// 	      }
// 	    }
// 	    if (!flag_check_again){
// 	      //find the longest one ...
// 	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
// 					pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
// 					pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));
	      
// 	      if (i==0&&k==1){
// 		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) && temp_length > 0.45*length_limit)
// 		  return true;
// 	      }else{
// 		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) && temp_length > 0.45*length_limit)
// 		  return true;
// 	      }
// 	    }	    
// 	  }
// 	}
//       }else{
// 	Vector main_dir = main_cluster->get_PCA_axis(0);
// 	TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
// 	TVector3 dir_test(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,
// 			  out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,
// 			  out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z);

// 	//	std::cout << main_cluster->get_cluster_id() << " " << fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.) << std::endl;

// 	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " ;
// 	// std::cout << (out_vec_wcps.at(k).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;
	
// 	if (fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.)>75 || i==0 && k==1){
// 	  // check dead region ...
// 	  bool flag_p1_inside_p = flag_p1_inside;
// 	  if (flag_p1_inside_p){
// 	    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
// 	    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
// 	    dir *= (-1);

// 	    if (dir.Angle(dir_test) > 3.1415926*2./3.) continue;
	    
// 	    //	    std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;

// 	    // check U and V and W
// 	    TVector3 dir_1(0,dir.Y(),dir.Z());
// 	    double angle1 = dir_1.Angle(U_dir);
// 	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
// 	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

// 	    double angle2 = dir_1.Angle(V_dir);
// 	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
// 	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

// 	    double angle3 = dir_1.Angle(W_dir);
// 	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
// 	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

// 	    // not added for now, need to check to add this one in when more events are available ...
// 	    // XQ, 7/11/2018
// 	    double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


// 	    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

// 	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
// 	      flag_p1_inside_p = flag_p1_inside_p && check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x);
// 	    }
	    
// 	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
// 	      flag_p1_inside_p= flag_p1_inside_p && check_dead_volume(p1,dir,1*units::cm,offset_x);
// 	  }
	  
// 	  bool flag_p2_inside_p = flag_p2_inside;
// 	  if (flag_p2_inside_p){
// 	    Point p2(out_vec_wcps.at(k).at(0).x,out_vec_wcps.at(k).at(0).y,out_vec_wcps.at(k).at(0).z);
// 	    TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
// 	    dir *= (-1);

// 	    if (dir.Angle(dir_test) < 3.1415926/3.) continue;
	    
// 	    //	    std::cout << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;
// 	    TVector3 dir_1(0,dir.Y(),dir.Z());
// 	    double angle1 = dir_1.Angle(U_dir);
// 	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
// 	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

// 	    double angle2 = dir_1.Angle(V_dir);
// 	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
// 	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

// 	    double angle3 = dir_1.Angle(W_dir);
// 	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
// 	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

// 	    //	    std::cout << "B: " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " <<  angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

// 	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
// 	      flag_p2_inside_p = flag_p2_inside_p && check_signal_processing(p2,dir,ct_point_cloud,1*units::cm,offset_x);
// 	    }
	    
// 	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
// 	      flag_p2_inside_p=flag_p2_inside_p && check_dead_volume(p2,dir,1*units::cm,offset_x);
// 	  }
	  
// 	  if ((!flag_p1_inside_p) && (!flag_p2_inside_p)){
// 	    if (flash->get_type()==2){
// 	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,2)+
// 					pow(out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,2)+
// 					pow(out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z,2));
// 	      if (i==0&&k==1){
// 		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud,true)) &&
// 		    temp_length > 0.45*length_limit
// 		    )
// 		  return true;
// 	      }else{
// 		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud)) &&
// 		    temp_length > 0.45*length_limit
// 		    )
// 		  return true;
// 	      }
// 	    }else{
// 	      return true;
// 	    }
// 	  }
// 	  // check signal processing ...

// 	// {
// 	// 	if (flag_p1_inside)
// 	// 	  ;
	
// 	// 	if (flag_p2_inside)
// 	// 	  ;
	
// 	// }

// 	}
//       }
//     }
//   }
  

//   return false;

//   // also check against the dead channel ...  
//   // // check the fiducial volume ...
//   // std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_main_axis_wcps();
//   // Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
//   // Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
//   // //double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
//   // bool flag_inside_p1 = inside_fiducial_volume(p1,offset_x);
//   // bool flag_inside_p2 = inside_fiducial_volume(p2,offset_x);
//   // //std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;

  
  
//   // // check the dead region ... 
//   // if (flag_inside_p1){
//   //   // define a local direction ...
//   //   TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
//   //   dir *= (-1);
//   //   flag_inside_p1=check_dead_volume(p1,dir,1*units::cm,offset_x);
//   // }
//   // if (flag_inside_p2){
//   //   // define a  local direction ...
//   //   TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
//   //   dir *= (-1);
//   //   flag_inside_p2=check_dead_volume(p2,dir,1*units::cm,offset_x);
//   // }
//   // return (!flag_inside_p1)&&(!flag_inside_p2);


//   // bool flag_2nd = true;
//        // {
	 
	 
//        // 	 if ((!flag_inside_p1)&&(!flag_inside_p2)){
//        // 	   event_type |= 1UL << 3; // through going muon ... 
//        // 	   flag_2nd = false;
//        // 	 }
	 
//        // 	 if (flag_2nd && ((!flag_inside_p1)|| (!flag_inside_p2) )){
//        // 	   // check the fiducial volume ...
//        // 	   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_extreme_wcps();
	   
//        // 	   Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
//        // 	   Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
	   
//        // 	   flag_inside_p1 = fid->inside_fiducial_volume(p1,offset_x);
//        // 	   flag_inside_p2 = fid->inside_fiducial_volume(p2,offset_x);

//        // 	   std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;
	   
//        // 	   // check the dead region ...
//        // 	   if (flag_inside_p1){
//        // 	     // define a local direction ...
//        // 	     TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
//        // 	     dir *= (-1);
//        // 	     flag_inside_p1=fid->check_dead_volume(p1,dir,1*units::cm,offset_x);
//        // 	   }
//        // 	   if (flag_inside_p2){
//        // 	     // define a  local direction ...
//        // 	     TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
//        // 	     dir *= (-1);
//        // 	     flag_inside_p2=fid->check_dead_volume(p2,dir,1*units::cm,offset_x);
//        // 	   }
	   
//        // 	   if ((!flag_inside_p1)&&(!flag_inside_p2)){
//        // 	     event_type |= 1UL << 3; // through going muon ... 
//        // 	   }
//        // 	 }

  
//   return false;
// }

// bool WireCellPID::ToyFiducial::check_neutrino_candidate(WireCell::PR3DCluster *main_cluster,WCPointCloud<double>::WCPoint& wcp1 ,WCPointCloud<double>::WCPoint& wcp2, double offset_x, WireCell::ToyCTPointCloud& ct_point_cloud,bool flag_2view_check){
//   main_cluster->Create_graph(ct_point_cloud);
  
//   main_cluster->dijkstra_shortest_paths(wcp1);
//   main_cluster->cal_shortest_path(wcp2);

//   std::list<WCPointCloud<double>::WCPoint>& path_wcps = main_cluster->get_path_wcps();
  
  
  
//   PointVector path_wcps_vec;
//   PointVector path_wcps_vec1;  
//   double low_dis_limit = 0.5*units::cm;
//   for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
//     if (path_wcps_vec.size()==0){
//       Point p((*it).x,(*it).y,(*it).z);
//       path_wcps_vec.push_back(p);
//       path_wcps_vec1.push_back(p);
//     }else{
//       double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
// 			+pow((*it).y - path_wcps_vec.back().y,2)
// 			+pow((*it).z - path_wcps_vec.back().z,2));
//       if (dis > low_dis_limit){
// 	Point p((*it).x,(*it).y,(*it).z);
// 	path_wcps_vec.push_back(p);
//       }

//       dis = sqrt(pow((*it).x - path_wcps_vec1.back().x,2)
// 		 +pow((*it).y - path_wcps_vec1.back().y,2)
// 		 +pow((*it).z - path_wcps_vec1.back().z,2));
//       if (dis <= 2*low_dis_limit){
// 	Point p((*it).x,(*it).y,(*it).z);
// 	path_wcps_vec1.push_back(p);
//       }else{
// 	int nseg = dis/2./low_dis_limit+1;
// 	for (int i=0;i!=nseg;i++){
// 	  Point temp_p;
// 	  temp_p.x = path_wcps_vec1.back().x + (i+1.)*((*it).x - path_wcps_vec1.back().x)/nseg;
// 	  temp_p.y = path_wcps_vec1.back().y + (i+1.)*((*it).y - path_wcps_vec1.back().y)/nseg;
// 	  temp_p.z = path_wcps_vec1.back().z + (i+1.)*((*it).z - path_wcps_vec1.back().z)/nseg;
// 	  path_wcps_vec1.push_back(temp_p);
// 	}
//       }
//     }
//   }  

  
//   //check whether path is good ... 
//   {
//     int num_nth = 0;
//     double min_dis = 1e9;

//     // bool flag_2view_check = true;
    
//     // U and V induction view checks
//     if (flag_2view_check){
//       TVector3 drift_dir(1,0,0);
//       // hard coded for U and V plane ... 
//       TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
//       TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
//       TVector3 W_dir(0,1,0);

//       TVector3 dir(wcp2.x-wcp1.x,wcp2.y-wcp1.y,wcp2.z-wcp1.z);
      
//       TVector3 dir_1(0,dir.Y(),dir.Z());
//       double angle1 = dir_1.Angle(U_dir);
//       TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
//       double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

//       double angle2 = dir_1.Angle(V_dir);
//       TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
//       double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
      
//       double angle3 = dir_1.Angle(W_dir);
//       TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
//       double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

//       double angle4 = fabs(3.1415926/2.-drift_dir.Angle(dir))/3.1415926*180.;
      
      
//       if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5 || angle4 < 5.)){
// 	flag_2view_check = false;
//       }
      
//     }

//     //int num_total_dead = 0;
    
//     int num_bad = 0;
    
//     for (int i=0;i!=path_wcps_vec1.size();i++){
//       WireCell::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,0);
//       WireCell::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,1);
//       WireCell::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,2);

//       bool flag_reset = false;

//       if (flag_2view_check){
// 	if (cloud_u.pts.size() >0 && cloud_v.pts.size() > 0 || // require two planes to be good ...
// 	    cloud_u.pts.size() >0 && cloud_w.pts.size() > 0 ||
// 	    cloud_v.pts.size() >0 && cloud_w.pts.size() > 0){
// 	  flag_reset = true;
// 	}else{
// 	  if (inside_dead_region(path_wcps_vec1.at(i)))
// 	    flag_reset = true;
// 	}
//       }else{
// 	if (cloud_u.pts.size() >0  || // require one plane to be good ...
// 	    cloud_v.pts.size() >0  ||
// 	    cloud_w.pts.size() >0 ){
// 	  flag_reset = true;
// 	}else{
// 	  if (inside_dead_region(path_wcps_vec1.at(i)))
// 	    flag_reset = true;
// 	}
//       }

//       if (!ct_point_cloud.is_good_point(path_wcps_vec1.at(i)))
// 	num_bad ++;
      
//       // std::cout << "O: " << path_wcps_vec1.at(i).x/units::cm << " " 
//       // 		<< path_wcps_vec1.at(i).y/units::cm << " "
//       // 		<< path_wcps_vec1.at(i).z/units::cm << " " << flag_reset << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << std::endl;
      
//       if (flag_reset){
//   	num_nth =0;
// 	min_dis = 1e9;
// 	num_bad = 0;
//       }else{
// 	if (inside_fiducial_volume(path_wcps_vec1.at(i),offset_x)){
// 	  double dis1 = sqrt(pow(path_wcps_vec1.at(i).x-wcp1.x,2)+pow(path_wcps_vec1.at(i).y-wcp1.y,2)+pow(path_wcps_vec1.at(i).z-wcp1.z,2));
// 	  double dis2 = sqrt(pow(path_wcps_vec1.at(i).x-wcp2.x,2)+pow(path_wcps_vec1.at(i).y-wcp2.y,2)+pow(path_wcps_vec1.at(i).z-wcp2.z,2));
// 	  if (dis1 < min_dis) min_dis = dis1;
// 	  if (dis2 < min_dis) min_dis = dis2;
// 	  num_nth ++;
// 	  // num_total_dead ++;
	  
// 	  // if (main_cluster->get_cluster_id()==7)
// 	  //   std::cout << main_cluster->get_cluster_id() << " " << min_dis/units::cm << " " << flag_2view_check << " " << path_wcps_vec1.at(i).x/units::cm << " " 
// 	  // 	      << path_wcps_vec1.at(i).y/units::cm << " "
// 	  // 	      << path_wcps_vec1.at(i).z/units::cm << " "
// 	  // 	      << num_nth << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << " " << num_bad << std::endl;
	  
// 	}
// 	//std::cout << num_nth << std::endl;

// 	if (num_nth > 7 && min_dis < 25*units::cm && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
// 	//if (num_nth > 7 && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
//       }
//     }
//   }
  
  
//   // if (cluster_id == 13){
//   //   std::cout << wcp1.x/units::cm << " " << wcp1.y/units::cm << " " << wcp1.z/units::cm << " " << wcp2.x/units::cm << " " << wcp2.y/units::cm << " " << wcp2.z/units::cm << std::endl;
//   int count = 0;
//   double max_angle=0;
//   Point max_point(0,0,0);
//   TVector3 drift_dir(1,0,0);
//   for (size_t i=5;i+5<path_wcps_vec.size();i++){
//     TVector3 dir1(path_wcps_vec.at(i).x - path_wcps_vec.at(i-5).x,
// 		  path_wcps_vec.at(i).y - path_wcps_vec.at(i-5).y,
// 		  path_wcps_vec.at(i).z - path_wcps_vec.at(i-5).z);
//     TVector3 dir2(path_wcps_vec.at(i).x - path_wcps_vec.at(i+5).x,
// 		  path_wcps_vec.at(i).y - path_wcps_vec.at(i+5).y,
// 		  path_wcps_vec.at(i).z - path_wcps_vec.at(i+5).z);
    
//     TVector3 dir3, dir4, dir5, dir6;
//     {
//       PointVector pts;
//       double temp_x = 0;
//       double temp_y = 0;
//       double temp_z = 0;
//       double temp_count = 0;
//       for (size_t j=1;j!=15;j++){
//        	if (i>=j){
// 	  Point pt(path_wcps_vec.at(i-j).x,path_wcps_vec.at(i-j).y,path_wcps_vec.at(i-j).z);
	  
// 	  if (j<=12&&j>2){
// 	    temp_x += pt.x;
// 	    temp_y += pt.y;
// 	    temp_z += pt.z;
// 	    temp_count ++;
// 	  }
// 	  pts.push_back(pt);
// 	}
// 	Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
// 	dir3 = main_cluster->calc_PCA_dir(pt,pts);
// 	dir5.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
// 		    temp_y/temp_count - path_wcps_vec.at(i).y,
// 		    temp_z/temp_count - path_wcps_vec.at(i).z);
// 	if (dir3.Angle(dir1)>3.1415926/2.)
// 	  dir3 *= -1;
//       }
//     }
//     {
//       PointVector pts;
//       double temp_x = 0;
//       double temp_y = 0;
//       double temp_z = 0;
//       double temp_count = 0;
//       for (size_t j=1;j!=15;j++){
// 	if (i+j<path_wcps_vec.size()){
// 	  Point pt(path_wcps_vec.at(i+j).x,path_wcps_vec.at(i+j).y,path_wcps_vec.at(i+j).z);
// 	  if (j<=12&&j>2){
// 	    temp_x += pt.x;
// 	    temp_y += pt.y;
// 	    temp_z += pt.z;
// 	    temp_count ++;
// 	  }
// 	  pts.push_back(pt);
// 	}
//       }
//       Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
//       dir4 = main_cluster->calc_PCA_dir(pt,pts);
//       dir6.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
//       		  temp_y/temp_count - path_wcps_vec.at(i).y,
//       		  temp_z/temp_count - path_wcps_vec.at(i).z);
//       if (dir4.Angle(dir2)>3.1415926/2.)
// 	dir4 *= -1;
//     }

//     int cut1 = 0;
//     if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180.>25) cut1++;
//     if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180.>25) cut1++;
//     if ((3.1415926 - dir5.Angle(dir6))/3.1415926*180.>25) cut1++;
//     int cut2 = 0;
//     if (fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. > 5) cut2++;
//     if (fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. > 5) cut2++;
//     if (fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. > 5) cut2++;

    
//     // if (main_cluster->get_cluster_id()==7)
//     //   std::cout << i << " " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926 - dir1.Angle(dir2))/3.1415926*180. << " " << (3.1415926 - dir3.Angle(dir4))/3.1415926*180. << " " << (3.1415926 - dir5.Angle(dir6))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. << " " << cut1 << " " << cut2 << std::endl;
  
   
    
    
//     if (cut1>=3 && cut2>=2){
//       if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180. > max_angle){
// 	max_angle = (3.1415926 - dir3.Angle(dir4))/3.1415926*180.;
// 	max_point = path_wcps_vec.at(i);
//       }
      
//       count ++;
//       if (count >=3){
// 	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
// 		       path_wcps_vec.at(i).y-wcp1.y,
// 		       path_wcps_vec.at(i).z-wcp1.z);
// 	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
// 		       path_wcps_vec.at(i).y-wcp2.y,
// 		       path_wcps_vec.at(i).z-wcp2.z);

// 	// if (main_cluster->get_cluster_id()==7)
// 	//   std::cout << "A: " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180.<< " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;


// 	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 ) ||
// 	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >32 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )&& temp1.Mag()>10*units::cm && temp2.Mag()>10*units::cm ||
// 	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >25 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 || (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
// 	    && temp1.Mag()>15*units::cm && temp2.Mag()>15*units::cm){

// 	  // if (main_cluster->get_cluster_id()==7)
// 	  //   std::cout << "B: " <<  (!inside_fiducial_volume(max_point,offset_x)) << " " << inside_dead_region(max_point) << std::endl;
	  
// 	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
// 	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
// 	      ){ // should not too close to anode 
// 	  }else{
// 	    return true;
// 	  }
// 	}
//       } else if (count>=1){
//       	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
//       		       path_wcps_vec.at(i).y-wcp1.y,
//       		       path_wcps_vec.at(i).z-wcp1.z);
//       	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
//       		       path_wcps_vec.at(i).y-wcp2.y,
//       		       path_wcps_vec.at(i).z-wcp2.z);

// 	//	std::cout << "B: " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. << " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;
	
//       	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 ||
// 	     (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
// 	    && temp1.Mag()>5*units::cm && temp2.Mag()>5*units::cm 
// 	    ){
// 	  //	  std::cout << "AAA" << std::endl; 

//       	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
//       	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
//       	      ){ // should not too close to anode 
//       	  }else{
//       	    return true;
//       	  }
//       	}
//       }
//     }else{
//       count = 0 ;
//       max_angle = 0;
//       max_point.x = 0;
//       max_point.y = 0;
//       max_point.z = 0;
//     }
//   }
  
//   return false;
// }
