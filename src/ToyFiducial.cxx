#include "WCPPID/ToyFiducial.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

using namespace WCP;

#include "Cosmic_tagger.h"

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


WCPPID::ToyFiducial::ToyFiducial(int dead_region_ch_ext, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w, double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode, int flag_data)
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

bool WCPPID::ToyFiducial::check_full_detector_dead(){
  std::set<int> dead_chs;
  for (auto it = mcell_time_map.begin(); it!=mcell_time_map.end(); it++){
    SlimMergeGeomCell *mcell = it->first;
    double start_x = (it->second.first-offset_t)/slope_t;
    double end_x = (it->second.second-offset_t)/slope_t;
    
    GeomWireSelection& uwires = mcell->get_uwires();
    GeomWireSelection& vwires = mcell->get_vwires();
    GeomWireSelection& wwires = mcell->get_wwires();

    if (end_x < -10*units::cm || start_x > 266*units::cm) continue;
    
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(0))!=bad_planes.end()){
      for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	int ch = (*it1)->channel();
	dead_chs.insert(ch);
	//std::cout << ch << std::endl;
      }
    }
    if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(1))!=bad_planes.end()){
      for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	int ch = (*it1)->channel();
	dead_chs.insert(ch);
	//std::cout << ch << std::endl;
      }
    }
    if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(2))!=bad_planes.end()){
      for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	int ch = (*it1)->channel();
	dead_chs.insert(ch);
	//std::cout << ch << std::endl;
      }
    }
  }
  //std::cout << dead_chs.size() << std::endl;

  if (dead_chs.size()>8000) return true;
  return false;
  
}

bool WCPPID::ToyFiducial::check_other_clusters(WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*>& clusters){
  Int_t ncount = 0;
  WCP::ToyPointCloud *main_pcloud = main_cluster->get_point_cloud();

  Int_t number_clusters = 0;
  Double_t total_length = 0;
  
  for (auto it = clusters.begin(); it!=clusters.end(); it++){
    WCPPID::PR3DCluster *cluster = *it;
    std::pair<WCP::WCPointCloud<double>::WCPoint,WCP::WCPointCloud<double>::WCPoint> wcps = cluster->get_two_boundary_wcps();
    double coverage_x = wcps.first.x - wcps.second.x;
    double length = sqrt(pow(wcps.first.x-wcps.second.x,2) + pow(wcps.first.y-wcps.second.y,2) + pow(wcps.first.z-wcps.second.z,2));
    WCP::ToyPointCloud *pcloud = cluster->get_point_cloud();
    std::tuple<int,int,double> results = main_pcloud->get_closest_points(pcloud);
    // std::cout << "ABC: " << coverage_x/units::cm <<  " " << length/units::cm << " " << std::get<2>(results)/units::cm << std::endl;
    //    if (coverage_x > 0.75*units::cm && length > 3*units::cm)
    if (std::get<2>(results) < 25*units::cm && fabs(coverage_x)>0.75*units::cm && length > 5*units::cm){
      number_clusters ++;
      total_length += length;
    }
  }


  
  if (number_clusters >0 && (number_clusters/3. + total_length/(35*units::cm)/number_clusters) >=1)
    return true;
  
  return false;
}

bool WCPPID::ToyFiducial::check_stm(WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*>& additional_clusters, double offset_x, double flash_time, WCP::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, int& event_type){

  //  check_full_detector_dead();

  


  //std::cout << flag_other_clusters << std::endl;
  
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
  if (main_cluster->get_point_cloud_steiner()->get_cloud().pts.size()==0)
    return false;

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
      WCP::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WCP::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
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
      WCP::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WCP::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
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
    std::cout << "Mid Point: A" << std::endl;
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
	std::cout << "Mid Point: B" << std::endl;
	return false;
      }else{
	if (dis1 < dis2){
	  last_wcp = wcps.second;	   
	}else{
	  last_wcp = wcps.first; 
	}
      }
      
    }else{
      std::cout << "Mid Point: C" << std::endl;
      return false;
    }
  }

  bool flag_other_clusters = check_other_clusters(main_cluster, additional_clusters);
  
  // forward check ...
  {
    if (flag_double_end) std::cout << "Forward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path(first_wcp, last_wcp);
    //std::cout << "haha" << std::endl;
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    //std::cout << "haha" << std::endl;
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);

    if (main_cluster->get_fine_tracking_path().size()<=3) return false;
    
    //std::cout << "haha " << std::endl;
    Point mid_p = main_cluster->adjust_rough_path(); 
    //std::cout << "haha" << std::endl;
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
    
    // check both end points for TGM ...
    WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
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

     //     std::cout << pts.front() << " " << pts.back() << " " << (!inside_fiducial_volume(pts.front(),offset_x)) << " " << (!inside_fiducial_volume(pts.back(),offset_x)) << std::endl;
    if ( (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      // std::cout << exit_L_TGM/units::cm << " " << left_L_TGM/units::cm << std::endl;
      bool flag_TGM_anode = false;
      
      if ((pts.back().x < 2*units::cm || pts.front().x < 2*units::cm) && kink_num >=0 && kink_num < pts.size()){ // at Anode ...
	if (pts.at(kink_num).x < 6*units::cm){
	  TVector3 v10(pts.back().x-pts.at(kink_num).x,pts.back().y-pts.at(kink_num).y,pts.back().z-pts.at(kink_num).z);
	  TVector3 v20(pts.front().x-pts.at(kink_num).x,pts.front().y-pts.at(kink_num).y,pts.front().z-pts.at(kink_num).z);
	  if (fabs(v10.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.Mag()>15*units::cm || fabs(v20.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.Mag()>15*units::cm)
	      flag_TGM_anode = true;
	  //std::cout << v10.Angle(drift_dir)/3.1415926*180. << " " << v20.Angle(drift_dir)/3.1415926*180. << std::endl;
	}
      }
      if ((exit_L < 3*units::cm || left_L < 3*units::cm || flag_TGM_anode)){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
	event_type |= 1UL << 3;
	return true;
      }
    }else if ((!inside_fiducial_volume(pts.front(),offset_x)) && left_L < 3*units::cm){
      // check dead volume and SP ...
      Point p1 = pts.back();
      TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
      dir *= (-1);
      if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
	event_type |= 1UL << 3;
	return true;
      }
      //std::cout << "ABC: " << dir.X() << " " << dir.Y() << " " << dir.Z() << " " << check_dead_volume(p1,dir,1*units::cm,offset_x) << std::endl;
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
	  left_L < 5*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.8 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }
    
    
      bool flag_pass = false;

      if (!flag_other_clusters){
	if (left_L < 40*units::cm) {
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	  
	  if (!flag_pass){
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
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	    }
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
      }else{
	if (flag_fix_end)
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm, true);
	else
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm, 0., 35*units::cm, true);
      }
      
      if (flag_pass) {
	main_cluster->clear_fit_tracks();
	main_cluster->search_other_tracks(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

	std::cout << main_cluster->get_fit_tracks().size() << std::endl;
	
	if (!detect_proton(main_cluster, kink_num)) return true;
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
    WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
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
    
    if ( (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      bool flag_TGM_anode = false;
      
      if ((pts.back().x < 2*units::cm || pts.front().x < 2*units::cm) && kink_num >=0 && kink_num < pts.size()){ // at Anode ...
	if (pts.at(kink_num).x < 6*units::cm){
	  TVector3 v10(pts.back().x-pts.at(kink_num).x,pts.back().y-pts.at(kink_num).y,pts.back().z-pts.at(kink_num).z);
	  TVector3 v20(pts.front().x-pts.at(kink_num).x,pts.front().y-pts.at(kink_num).y,pts.front().z-pts.at(kink_num).z);
	  if (fabs(v10.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.Mag()>15*units::cm || fabs(v20.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.Mag()>15*units::cm)
	      flag_TGM_anode = true;
	  //std::cout << v10.Angle(drift_dir)/3.1415926*180. << " " << v20.Angle(drift_dir)/3.1415926*180. << std::endl;
	}
      }
      if ((exit_L < 3*units::cm || left_L < 3*units::cm) || flag_TGM_anode){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
	event_type |= 1UL << 3;
	return true;
      }
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
      if (!flag_other_clusters){
	if (left_L < 40*units::cm) {
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	  
	  if (!flag_pass){
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
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	    }
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
      }else{
	if (flag_fix_end)
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm, true);
	else
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm, 0., 35*units::cm, true);
      }

      if (flag_pass) {
	main_cluster->clear_fit_tracks();
	main_cluster->search_other_tracks(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
	if (!detect_proton(main_cluster, kink_num)) return true;
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

int WCPPID::ToyFiducial::find_first_kink(WCPPID::PR3DCluster* main_cluster){
  WCP::PointVector& fine_tracking_path = main_cluster->get_fine_tracking_path();
  std::vector<double>& dQ = main_cluster->get_dQ();
  std::vector<double>& dx = main_cluster->get_dx();
  std::vector<double>& pu = main_cluster->get_pu();
  std::vector<double>& pv = main_cluster->get_pv();
  std::vector<double>& pw = main_cluster->get_pw();
  std::vector<double>& pt = main_cluster->get_pt();

  //  std::cout << dQ.size() << " " << pu.size() << std::endl;
  
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
    //std::cout << i << " " << refl_angles.at(i) << " " << ave_angles.at(i)  << " " << inside_fiducial_volume(fine_tracking_path.at(i)) << std::endl;
    if ((refl_angles.at(i) > 25.5 && ave_angles.at(i) > 12.5 ) && inside_fiducial_volume(fine_tracking_path.at(i))){
      TVector3 v10(fine_tracking_path.at(i).x - fine_tracking_path.front().x,
		   fine_tracking_path.at(i).y - fine_tracking_path.front().y,
		   fine_tracking_path.at(i).z - fine_tracking_path.front().z);
      TVector3 v20(fine_tracking_path.back().x - fine_tracking_path.at(i).x,
		   fine_tracking_path.back().y - fine_tracking_path.at(i).y,
		   fine_tracking_path.back().z - fine_tracking_path.at(i).z);
      double angle3 = v10.Angle(v20)/3.1415926*180.;

      // std::cout << angle3 << std::endl;
      
      if (angle3 < 20 && ave_angles.at(i) < 20 || angle3 < 12.5 && inside_dead_region(fine_tracking_path.at(i)) || angle3 < 7.5 || i<=4) continue;
      if (angle3 > 30){

	// shorted Y ...
	if (pw.at(i) > 7135 && pw.at(i)< 7264){
	  bool flag_bad = false;
	  for (int k=-1;k!=2;k++){
	    auto it1 = ch_mcell_set_map.find(std::round(pv.at(i)+k));
	    if ( it1!=ch_mcell_set_map.end() && it1->second.size()>0){
	      flag_bad = true;
	      break;
	    }
	  }
	  if (flag_bad) continue;
	}
	
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
	
	//	std::cout << "Test: " << angle3 << " " << sum_fQ << " " << sum_bQ << std::endl;
	if (sum_fQ > 0.6 && sum_bQ > 0.6 || sum_fQ + sum_bQ > 1.4 && (sum_fQ > 0.8 || sum_bQ > 0.8) && v10.Mag() > 10*units::cm && v20.Mag() > 10*units::cm){
	  //std::cout << i << " " << max_numbers.at(i) << " " << dQ.size() << std::endl;
	  if (i+2<dQ.size()){
	    std::cout << "Kink: " << i << " " << refl_angles.at(i) << " " << para_angles.at(i) << " " << ave_angles.at(i) << " " << max_numbers.at(i) << " " << angle3 << " " << dQ.at(i)/dx.at(i)*units::cm/50e3 << " " << pu.at(i) << " " << pv.at(i) << " " << pw.at(i) << std::endl;
	    return max_numbers.at(i);
	  }
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
	// shorted Y ...
	if (pw.at(i) > 7135 && pw.at(i)< 7264){
	  bool flag_bad = false;
	  for (int k=-1;k!=2;k++){
	    auto it1 = ch_mcell_set_map.find(std::round(pv.at(i)+k));
	    if (it1!=ch_mcell_set_map.end() && it1->second.size()>0){
	      flag_bad = true;
	      break;
	    }
	  }
	  if (flag_bad) continue;
	}

	bool flag_bad_u = false;
	{
	  for (int k=-1;k!=2;k++){
	    auto it1 = ch_mcell_set_map.find(std::round(pu.at(i)+k));
	    if (it1!=ch_mcell_set_map.end() && it1->second.size()>0){
	      flag_bad_u = true;
	      break;
	    }
	  }
	}
	bool flag_bad_v = false;
	{
	  for (int k=-1;k!=2;k++){
	    auto it1 = ch_mcell_set_map.find(std::round(pv.at(i)+k));
	    if (it1!=ch_mcell_set_map.end() && it1->second.size()>0){
	      flag_bad_v = true;
	      break;
	    }
	  }
	}
	bool flag_bad_w = false;
	{
	  for (int k=-1;k!=2;k++){
	    auto it1 = ch_mcell_set_map.find(std::round(pw.at(i)+k));
	    if (it1!=ch_mcell_set_map.end() && it1->second.size()>0){
	      flag_bad_w = true;
	      break;
	    }
	  }
	}
	
	
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
	if (fabs(sum_fQ-sum_bQ) < 0.07*(sum_fQ+sum_bQ) && (flag_bad_u||flag_bad_v||flag_bad_w)) continue;
	
	if (sum_fQ > 0.6 && sum_bQ > 0.6 ){
	  if (i+2<dQ.size()){
	    std::cout << "Kink: " << i << " " << refl_angles.at(i) << " " << para_angles.at(i) << " " << ave_angles.at(i) << " " << max_numbers.at(i) << " " << angle3 << " " << dQ.at(i)/dx.at(i)*units::cm/50e3 << std::endl;
	    return max_numbers.at(i);
	  }
	}
      }
    }
  }

  
  
  return fine_tracking_path.size();
}

bool WCPPID::ToyFiducial::detect_proton(WCPPID::PR3DCluster* main_cluster,int kink_num){
  // common part ..
  WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
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
  double end_L; 
  double max_num; 
  if (kink_num == pts.size()){
    end_L = L.back();
    max_num = L.size();
  }else{
    end_L = L.at(kink_num)-0.5*units::cm;
    max_num = kink_num;
  }
  
  // find the maximum bin ...
  double max_bin = -1;
  double max_sum = 0;
  for (size_t i=0;i!=L.size();i++){
    double sum = 0;
    double nsum = 0;
    double temp_max_bin = i;
    double temp_max_val = dQ_dx.at(i);
    if (L.at(i) < end_L + 0.5*units::cm && L.at(i) > end_L - 40*units::cm && i < max_num){
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
  //  std::cout << max_bin << " " << max_sum << std::endl;
  end_L = L.at(max_bin)+0.2*units::cm;
  int ncount = 0;//, ncount_p = 0;
  std::vector<double> vec_x;//, vec_xp;
  std::vector<double> vec_y;//, vec_yp;

  for (size_t i=0;i!=L.size(); i++){
    if (end_L - L.at(i) < 35*units::cm && end_L - L.at(i) > 3*units::cm){
      vec_x.push_back(end_L-L.at(i));
      vec_y.push_back(dQ_dx.at(i));
      ncount ++;
    }

    // if (end_L - L.at(i) < 20*units::cm){
    //   vec_xp.push_back(end_L-L.at(i));
    //   vec_yp.push_back(dQ_dx.at(i));
    //   ncount_p ++;
    // }
  }
  
  if (ncount >=5){
    TH1F *h1 = new TH1F("h1","h1",ncount,0,ncount);
    TH1F *h2 = new TH1F("h2","h2",ncount,0,ncount);
    TH1F *h3 = new TH1F("h3","h3",ncount,0,ncount);

    // TH1F *h4 = new TH1F("h4","h4",ncount_p,0,ncount_p);
    // TH1F *h5 = new TH1F("h5","h5",ncount_p,0,ncount_p);
    
    for (size_t i=0;i!=ncount;i++){
      // std::cout << i << " " << vec_y.at(i) << std::endl;
      h1->SetBinContent(i+1,vec_y.at(i));
      h2->SetBinContent(i+1,g_muon->Eval((vec_x.at(i))/units::cm));
      h3->SetBinContent(i+1,50e3);
    }
    // for (size_t i=0;i!=ncount_p;i++){
    //   h4->SetBinContent(i+1,vec_yp.at(i));
    //   h5->SetBinContent(i+1,g_muon->Eval((vec_xp.at(i))/units::cm));
    // }
    
    double ks1 = h2->KolmogorovTest(h1,"M");
    double ratio1 = h2->GetSum()/(h1->GetSum()+1e-9);
    double ks2 = h3->KolmogorovTest(h1,"M");
    double ratio2 = h3->GetSum()/(h1->GetSum()+1e-9);
    //double ks3 = h4->KolmogorovTest(h5,"M");
    //double ratio3 = h4->GetSum()/(h5->GetSum()+1e-9);
    delete h1;
    delete h2;
    delete h3;
    // delete h4;
    // delete h5;
    
    std::cout << "End proton detection: " << ks1 << " " << ks2 << " " << ratio1 << " " << ratio2 << " " << ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << " " << dQ_dx.at(max_bin)/50e3 << " " << dQ_dx.size() - max_bin << " " << std::endl;
    
    if ( ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > 0.03 && dQ_dx.at(max_bin)/50e3 > 2.5 && (dQ_dx.size() - max_bin <= 3 || ks2 < 0.05 && dQ_dx.size() - max_bin <= 12) ) {

      if (dQ_dx.size()-max_bin<=1 &&  (dQ_dx.at(max_bin)/50e3 < 3.0 && (ks1 < 0.06 || ks1<0.065 && ks2>0.04) || ks1<0.035 && dQ_dx.at(max_bin)/50e3 < 4.0 )) return false;

      return true;
    }
  }
  
  return false;
}

bool WCPPID::ToyFiducial::eval_stm(WCPPID::PR3DCluster* main_cluster,int kink_num,  double peak_range, double offset_length, double com_range, bool flag_strong_check){
  WCP::PointVector& pts = main_cluster->get_fine_tracking_path();

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
  std::vector<double> vec_res_x;
  std::vector<double> vec_res_y;
  for (size_t i=0;i!=L.size(); i++){
    if (end_L - L.at(i) < com_range && end_L - L.at(i) > 0){
      vec_x.push_back(end_L-L.at(i));
      vec_y.push_back(dQ_dx.at(i));
      ncount ++;
    }else if (L.at(i)>end_L){
      vec_res_x.push_back(L.at(i)-end_L);
      vec_res_y.push_back(dQ_dx.at(i));
    }
  }

  double ave_res_dQ_dx = 0;
  double res_length = 0;
  
  for (size_t i=0;i!=vec_res_y.size();i++){
    ave_res_dQ_dx += vec_res_y.at(i);
  }
  if (vec_res_y.size()>0){
    res_length = vec_res_x.back();
    ave_res_dQ_dx /= 1.*vec_res_y.size();
  }
  
  //std::cout << "Test: " << res_length/units::cm << " " << ave_res_dQ_dx << std::endl;
  
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
  if (sqrt(pow(ks2/0.06,2)+pow((ratio2-1)/0.06,2))< 1.4 && ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > -0.02) return false;

  // if residual does not look like a michel electron
  if (res_length > 16 * units::cm && ave_res_dQ_dx > 72500 || 
      res_length > 10 * units::cm && ave_res_dQ_dx > 72500 && ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > -0.05 ||
      res_length > 10 * units::cm && ave_res_dQ_dx > 85000 ||
      res_length > 6 * units::cm && ave_res_dQ_dx > 92500 ||
      res_length > 6 * units::cm && ave_res_dQ_dx > 72500 && ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > -0.05)
    return false;
  
  if (!flag_strong_check){
    if (ks1 - ks2 < -0.02 && (ks2 > 0.09 || ratio2 > 1.5)) return true;
    if ( ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 < 0) return true;
  }else{
    if (ks1 - ks2 < -0.02 && (ks2 > 0.09 || ratio2 > 1.5) && ks1 < 0.05 && fabs(ratio1-1)<0.1) return true;
    if ( ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 < 0 && ks1 < 0.05 && fabs(ratio1-1)<0.1) return true;
  }

  return false;
  
}




WCPPID::ToyFiducial::~ToyFiducial(){
  delete g_muon;
  delete g_pion;
  delete g_proton;
  delete g_kaon;
  delete g_electron;
  delete file;
		 
}

bool WCPPID::ToyFiducial::check_signal_processing(WCP::Point& p, TVector3& dir, WCP::ToyCTPointCloud& ct_point_cloud, double step, double offset_x){

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
      
      WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,0);
      WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,1);
      WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,2);
      
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

bool WCPPID::ToyFiducial::check_dead_volume(WCP::Point& p, TVector3& dir, double step, double offset_x){
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

      //      std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
      
      if (num_points_dead > 0.81*num_points){
	return false;
      }else{
	return true;
      }
      
    }
  }
}


bool WCPPID::ToyFiducial::inside_fiducial_volume(WCP::Point& p, double offset_x){

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

bool WCPPID::ToyFiducial::inside_dead_region(WCP::Point& p){
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


void WCPPID::ToyFiducial::AddDeadRegion(WCP::SlimMergeGeomCell* mcell, std::vector<int>& time_slices){

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


