 
std::tuple<int, WCPPID::PR3DCluster*, WCP::Opflash*> WCPPID::ToyFiducial::cosmic_tagger(WCP::OpflashSelection& flashes, WCPPID::PR3DCluster* main_cluster, WCP::Opflash* main_flash, std::tuple<int, double, double, int>& bundle_info, WCP::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud, int run_no, int subrun_no, int event_no, bool flag_data, bool debug_tagger){
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

  // Note, you need to update both flash, and offset_x accordingly ...
  // example of the check fully contained ...
  // check_fully_contained(main_cluster, offset_x,  ct_point_cloud)

  // example of the TGM checking ...
  // check_tgm(main_cluster,flash, offset_x, ct_point_cloud)
  
  return std::make_tuple(0, main_cluster, main_flash);
}





//Enabled the fully contained check
bool WCPPID::ToyFiducial::check_fully_contained(WCPPID::PR3DCluster* main_cluster, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud){
  
  WCPPID::PR3DCluster *main_cluster1 = main_cluster;
  

  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();

  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  Vector main_dir = main_cluster->get_PCA_axis(0);
  TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    // check all the points ... 
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      if (!inside_fiducial_volume(p1,offset_x)){
        return false;
      }
    }


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
    
    // not added for now, need to check to add this one in when more events are available ...
    // XQ, 7/11/2018
    // double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
      if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
        return false;
      }
    }
	    
    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
      if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	return false;
      }
    }
    
  }

  
  return true;
}


// check TGM code ...

bool WCPPID::ToyFiducial::check_tgm(WCPPID::PR3DCluster* main_cluster,Opflash* main_flash, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud){

  WCPPID::PR3DCluster *main_cluster1 = main_cluster; 
  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();


  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  double length_limit = sqrt(pow(out_vec_wcps.at(0).at(0).x-out_vec_wcps.at(1).at(0).x,2)+
			     pow(out_vec_wcps.at(0).at(0).y-out_vec_wcps.at(1).at(0).y,2)+
			     pow(out_vec_wcps.at(0).at(0).z-out_vec_wcps.at(1).at(0).z,2));
  
  // std::cout << "Flash: " << flash->get_flash_id() << std::endl;

  // take a look at the first point ...
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    bool flag_p1_inside = true;
    int p1_index = -1;
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      flag_p1_inside = flag_p1_inside && inside_fiducial_volume(p1,offset_x);
      if (!flag_p1_inside){
	p1_index = j;
	break;
      }
    }
    
    
    // loop through the remaining groups and check ...
    for (size_t k=i+1;k!=out_vec_wcps.size();k++){
      bool flag_p2_inside = true;
      int p2_index = -1;
      for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
	Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
	flag_p2_inside = flag_p2_inside && inside_fiducial_volume(p2,offset_x);
	if (!flag_p2_inside){
	  p2_index = j;
	  break;
	}
      }
      

      // if (main_cluster->get_cluster_id()==7){
      //  	std::cout << main_cluster->get_cluster_id() << " " << i << " " <<
      //  	  out_vec_wcps.at(i).at(0).x/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " << 
      //  	  k << " " << out_vec_wcps.at(k).at(0).x/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm <<
      //  	  " " << p1_index << " " << p2_index << " " << flag_p1_inside << " " << flag_p2_inside << std::endl;
      // 	// for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      // 	//   Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      // 	//   std::cout << j << " A " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << std::endl;
      // 	// }
      // 	// for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
      // 	//   Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
      // 	//   std::cout << j << " B " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << std::endl;
      // 	// }
      // }
	
      if ((!flag_p1_inside) && (!flag_p2_inside)){
	// if not a neutrino candidate ... to be worked out ...
	 // Point p1(out_vec_wcps.at(i).at(p1_index).x,out_vec_wcps.at(i).at(p1_index).y,out_vec_wcps.at(i).at(p1_index).z);
	 // Point p2(out_vec_wcps.at(k).at(p2_index).x,out_vec_wcps.at(k).at(p2_index).y,out_vec_wcps.at(k).at(p2_index).z);
	//	std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << inside_fiducial_volume(p1,offset_x) << " A " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << inside_fiducial_volume(p2,offset_x) << " " << offset_x/units::cm << std::endl;


	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(p1_index).x-offset_x)/units::cm <<3=r-t=-[\] " " << out_vec_wcps.at(i).at(p1_index).y/units::cm << " " << out_vec_wcps.at(i).at(p1_index).z/units::cm << " " ;
	// std::cout << (out_vec_wcps.at(k).at(p2_index).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(p2_index).y/units::cm << " " << out_vec_wcps.at(k).at(p2_index).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;

	//std::cout << p1_index << " " << p2_index << std::endl;
	
	// check two points in between
	bool flag_check = false;
	for (int kk=0;kk!=3;kk++){
	  Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	  flag_check = flag_check || inside_fiducial_volume(p3,offset_x);
	}
	
	// if (main_cluster->get_cluster_id()==10)
	//   std::cout << flag_check << " " << out_vec_wcps.size() << std::endl;
	
	if (flag_check){
	  if (main_flash->get_type()==2){
		    
	    double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
				      pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
				      pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));

	    //std::cout << temp_length/units::cm << " " << length_limit/units::cm << std::endl;

	    if (i==0&&k==1){
	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) &&
		   temp_length > 0.45*length_limit)
		return true;
	    }else{
	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) &&
		   temp_length > 0.45*length_limit)
		return true;
	    }
	    
	  }else{
	    return true; // through going muon ...
	  }
	}else{
	  
	  if (out_vec_wcps.size()==2){
	    return true;
	  }else{
	    // if (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x))
	    //   return true;

	    bool flag_check_again = false;
	    for (int kkk = 0;kkk!=out_vec_wcps.size(); kkk++){
	      if (kkk==i || kkk==k) continue;
	      for (int kk=0;kk!=4;kk++){
	    	Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(kkk).at(0).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
	    		 out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(kkk).at(0).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
	    		 out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(kkk).at(0).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
	      }
	      
	      for (int kk=0;kk!=3;kk++){
	    	Point p3(out_vec_wcps.at(kkk).at(0).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(0).x)/4.*(kk+1),
	    		 out_vec_wcps.at(kkk).at(0).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(0).y)/4.*(kk+1),
	    		 out_vec_wcps.at(kkk).at(0).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(0).z)/4.*(kk+1));
	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
	      }
	    }
	    if (!flag_check_again){
	      //find the longest one ...
	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
					pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
					pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));
	      
	      if (i==0&&k==1){
		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) && temp_length > 0.45*length_limit)
		  return true;
	      }else{
		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) && temp_length > 0.45*length_limit)
		  return true;
	      }
	    }	    
	  }
	}
      }else{
	Vector main_dir = main_cluster->get_PCA_axis(0);
	TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
	TVector3 dir_test(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,
			  out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,
			  out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z);

	//	std::cout << main_cluster->get_cluster_id() << " " << fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.) << std::endl;

	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " ;
	// std::cout << (out_vec_wcps.at(k).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;
	
	if (fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.)>75 || i==0 && k==1){
	  // check dead region ...
	  bool flag_p1_inside_p = flag_p1_inside;
	  if (flag_p1_inside_p){
	    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
	    dir *= (-1);

	    if (dir.Angle(dir_test) > 3.1415926*2./3.) continue;
	    
	    //	    std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;

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

	    // not added for now, need to check to add this one in when more events are available ...
	    // XQ, 7/11/2018
	    double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


	    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p1_inside_p = flag_p1_inside_p && check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p1_inside_p= flag_p1_inside_p && check_dead_volume(p1,dir,1*units::cm,offset_x);
	  }
	  
	  bool flag_p2_inside_p = flag_p2_inside;
	  if (flag_p2_inside_p){
	    Point p2(out_vec_wcps.at(k).at(0).x,out_vec_wcps.at(k).at(0).y,out_vec_wcps.at(k).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
	    dir *= (-1);

	    if (dir.Angle(dir_test) < 3.1415926/3.) continue;
	    
	    //	    std::cout << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;
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

	    //	    std::cout << "B: " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " <<  angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p2_inside_p = flag_p2_inside_p && check_signal_processing(p2,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p2_inside_p=flag_p2_inside_p && check_dead_volume(p2,dir,1*units::cm,offset_x);
	  }
	  
	  if ((!flag_p1_inside_p) && (!flag_p2_inside_p)){
	    if (main_flash->get_type()==2){
	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,2)+
					pow(out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,2)+
					pow(out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z,2));
	      if (i==0&&k==1){
		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud,true)) &&
		    temp_length > 0.45*length_limit
		    )
		  return true;
	      }else{
		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud)) &&
		    temp_length > 0.45*length_limit
		    )
		  return true;
	      }
	    }else{
	      return true;
	    }
	  }
	  // check signal processing ...

	// {
	// 	if (flag_p1_inside)
	// 	  ;
	
	// 	if (flag_p2_inside)
	// 	  ;
	
	// }

	}
      }
    }
  }
  

  return false;

  // also check against the dead channel ...  
  // // check the fiducial volume ...
  // std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_main_axis_wcps();
  // Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
  // Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
  // //double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
  // bool flag_inside_p1 = inside_fiducial_volume(p1,offset_x);
  // bool flag_inside_p2 = inside_fiducial_volume(p2,offset_x);
  // //std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;

  
  
  // // check the dead region ... 
  // if (flag_inside_p1){
  //   // define a local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p1=check_dead_volume(p1,dir,1*units::cm,offset_x);
  // }
  // if (flag_inside_p2){
  //   // define a  local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p2=check_dead_volume(p2,dir,1*units::cm,offset_x);
  // }
  // return (!flag_inside_p1)&&(!flag_inside_p2);


  // bool flag_2nd = true;
       // {
	 
	 
       // 	 if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	   event_type |= 1UL << 3; // through going muon ... 
       // 	   flag_2nd = false;
       // 	 }
	 
       // 	 if (flag_2nd && ((!flag_inside_p1)|| (!flag_inside_p2) )){
       // 	   // check the fiducial volume ...
       // 	   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_extreme_wcps();
	   
       // 	   Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
       // 	   Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
	   
       // 	   flag_inside_p1 = fid->inside_fiducial_volume(p1,offset_x);
       // 	   flag_inside_p2 = fid->inside_fiducial_volume(p2,offset_x);

       // 	   std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;
	   
       // 	   // check the dead region ...
       // 	   if (flag_inside_p1){
       // 	     // define a local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p1=fid->check_dead_volume(p1,dir,1*units::cm,offset_x);
       // 	   }
       // 	   if (flag_inside_p2){
       // 	     // define a  local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p2=fid->check_dead_volume(p2,dir,1*units::cm,offset_x);
       // 	   }
	   
       // 	   if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	     event_type |= 1UL << 3; // through going muon ... 
       // 	   }
       // 	 }

  
  return false;
}



bool WCPPID::ToyFiducial::check_neutrino_candidate(WCPPID::PR3DCluster *main_cluster,WCPointCloud<double>::WCPoint& wcp1 ,WCPointCloud<double>::WCPoint& wcp2, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,bool flag_2view_check){
  main_cluster->Create_graph(ct_point_cloud);
  
  main_cluster->dijkstra_shortest_paths(wcp1);
  main_cluster->cal_shortest_path(wcp2);

  std::list<WCPointCloud<double>::WCPoint>& path_wcps = main_cluster->get_path_wcps();
  
  
  
  PointVector path_wcps_vec;
  PointVector path_wcps_vec1;  
  double low_dis_limit = 0.5*units::cm;
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_wcps_vec.size()==0){
      Point p((*it).x,(*it).y,(*it).z);
      path_wcps_vec.push_back(p);
      path_wcps_vec1.push_back(p);
    }else{
      double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
			+pow((*it).y - path_wcps_vec.back().y,2)
			+pow((*it).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec.push_back(p);
      }

      dis = sqrt(pow((*it).x - path_wcps_vec1.back().x,2)
		 +pow((*it).y - path_wcps_vec1.back().y,2)
		 +pow((*it).z - path_wcps_vec1.back().z,2));
      if (dis <= 2*low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec1.push_back(p);
      }else{
	int nseg = dis/2./low_dis_limit+1;
	for (int i=0;i!=nseg;i++){
	  Point temp_p;
	  temp_p.x = path_wcps_vec1.back().x + (i+1.)*((*it).x - path_wcps_vec1.back().x)/nseg;
	  temp_p.y = path_wcps_vec1.back().y + (i+1.)*((*it).y - path_wcps_vec1.back().y)/nseg;
	  temp_p.z = path_wcps_vec1.back().z + (i+1.)*((*it).z - path_wcps_vec1.back().z)/nseg;
	  path_wcps_vec1.push_back(temp_p);
	}
      }
    }
  }  

  
  //check whether path is good ... 
  {
    int num_nth = 0;
    double min_dis = 1e9;

    // bool flag_2view_check = true;
    
    // U and V induction view checks
    if (flag_2view_check){
      TVector3 drift_dir(1,0,0);
      // hard coded for U and V plane ... 
      TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
      TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
      TVector3 W_dir(0,1,0);

      TVector3 dir(wcp2.x-wcp1.x,wcp2.y-wcp1.y,wcp2.z-wcp1.z);
      
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

      double angle4 = fabs(3.1415926/2.-drift_dir.Angle(dir))/3.1415926*180.;
      
      
      if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5 || angle4 < 5.)){
	flag_2view_check = false;
      }
      
    }

    //int num_total_dead = 0;
    
    int num_bad = 0;
    
    for (int i=0;i!=path_wcps_vec1.size();i++){
      WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,0);
      WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,1);
      WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,2);

      bool flag_reset = false;

      if (flag_2view_check){
	if (cloud_u.pts.size() >0 && cloud_v.pts.size() > 0 || // require two planes to be good ...
	    cloud_u.pts.size() >0 && cloud_w.pts.size() > 0 ||
	    cloud_v.pts.size() >0 && cloud_w.pts.size() > 0){
	  flag_reset = true;
	}else{
	  if (inside_dead_region(path_wcps_vec1.at(i)))
	    flag_reset = true;
	}
      }else{
	if (cloud_u.pts.size() >0  || // require one plane to be good ...
	    cloud_v.pts.size() >0  ||
	    cloud_w.pts.size() >0 ){
	  flag_reset = true;
	}else{
	  if (inside_dead_region(path_wcps_vec1.at(i)))
	    flag_reset = true;
	}
      }

      if (!ct_point_cloud.is_good_point(path_wcps_vec1.at(i)))
	num_bad ++;
      
      // std::cout << "O: " << path_wcps_vec1.at(i).x/units::cm << " " 
      // 		<< path_wcps_vec1.at(i).y/units::cm << " "
      // 		<< path_wcps_vec1.at(i).z/units::cm << " " << flag_reset << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << std::endl;
      
      if (flag_reset){
  	num_nth =0;
	min_dis = 1e9;
	num_bad = 0;
      }else{
	if (inside_fiducial_volume(path_wcps_vec1.at(i),offset_x)){
	  double dis1 = sqrt(pow(path_wcps_vec1.at(i).x-wcp1.x,2)+pow(path_wcps_vec1.at(i).y-wcp1.y,2)+pow(path_wcps_vec1.at(i).z-wcp1.z,2));
	  double dis2 = sqrt(pow(path_wcps_vec1.at(i).x-wcp2.x,2)+pow(path_wcps_vec1.at(i).y-wcp2.y,2)+pow(path_wcps_vec1.at(i).z-wcp2.z,2));
	  if (dis1 < min_dis) min_dis = dis1;
	  if (dis2 < min_dis) min_dis = dis2;
	  num_nth ++;
	  // num_total_dead ++;
	  
	  // if (main_cluster->get_cluster_id()==7)
	  //   std::cout << main_cluster->get_cluster_id() << " " << min_dis/units::cm << " " << flag_2view_check << " " << path_wcps_vec1.at(i).x/units::cm << " " 
	  // 	      << path_wcps_vec1.at(i).y/units::cm << " "
	  // 	      << path_wcps_vec1.at(i).z/units::cm << " "
	  // 	      << num_nth << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << " " << num_bad << std::endl;
	  
	}
	//std::cout << num_nth << std::endl;

	if (num_nth > 7 && min_dis < 25*units::cm && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
	//if (num_nth > 7 && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
      }
    }
  }
  
  
  // if (cluster_id == 13){
  //   std::cout << wcp1.x/units::cm << " " << wcp1.y/units::cm << " " << wcp1.z/units::cm << " " << wcp2.x/units::cm << " " << wcp2.y/units::cm << " " << wcp2.z/units::cm << std::endl;
  int count = 0;
  double max_angle=0;
  Point max_point(0,0,0);
  TVector3 drift_dir(1,0,0);
  for (size_t i=5;i+5<path_wcps_vec.size();i++){
    TVector3 dir1(path_wcps_vec.at(i).x - path_wcps_vec.at(i-5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i-5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i-5).z);
    TVector3 dir2(path_wcps_vec.at(i).x - path_wcps_vec.at(i+5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i+5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i+5).z);
    
    TVector3 dir3, dir4, dir5, dir6;
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
       	if (i>=j){
	  Point pt(path_wcps_vec.at(i-j).x,path_wcps_vec.at(i-j).y,path_wcps_vec.at(i-j).z);
	  
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
	Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
	dir3 = main_cluster->calc_PCA_dir(pt,pts);
	dir5.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
		    temp_y/temp_count - path_wcps_vec.at(i).y,
		    temp_z/temp_count - path_wcps_vec.at(i).z);
	if (dir3.Angle(dir1)>3.1415926/2.)
	  dir3 *= -1;
      }
    }
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
	if (i+j<path_wcps_vec.size()){
	  Point pt(path_wcps_vec.at(i+j).x,path_wcps_vec.at(i+j).y,path_wcps_vec.at(i+j).z);
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
      }
      Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
      dir4 = main_cluster->calc_PCA_dir(pt,pts);
      dir6.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
      		  temp_y/temp_count - path_wcps_vec.at(i).y,
      		  temp_z/temp_count - path_wcps_vec.at(i).z);
      if (dir4.Angle(dir2)>3.1415926/2.)
	dir4 *= -1;
    }

    int cut1 = 0;
    if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180.>25) cut1++;
    if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180.>25) cut1++;
    if ((3.1415926 - dir5.Angle(dir6))/3.1415926*180.>25) cut1++;
    int cut2 = 0;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. > 5) cut2++;

    
    // if (main_cluster->get_cluster_id()==7)
    //   std::cout << i << " " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926 - dir1.Angle(dir2))/3.1415926*180. << " " << (3.1415926 - dir3.Angle(dir4))/3.1415926*180. << " " << (3.1415926 - dir5.Angle(dir6))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. << " " << cut1 << " " << cut2 << std::endl;
  
   
    
    
    if (cut1>=3 && cut2>=2){
      if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180. > max_angle){
	max_angle = (3.1415926 - dir3.Angle(dir4))/3.1415926*180.;
	max_point = path_wcps_vec.at(i);
      }
      
      count ++;
      if (count >=3){
	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
		       path_wcps_vec.at(i).y-wcp1.y,
		       path_wcps_vec.at(i).z-wcp1.z);
	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
		       path_wcps_vec.at(i).y-wcp2.y,
		       path_wcps_vec.at(i).z-wcp2.z);

	// if (main_cluster->get_cluster_id()==7)
	//   std::cout << "A: " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180.<< " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;


	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 ) ||
	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >32 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )&& temp1.Mag()>10*units::cm && temp2.Mag()>10*units::cm ||
	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >25 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 || (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
	    && temp1.Mag()>15*units::cm && temp2.Mag()>15*units::cm){

	  // if (main_cluster->get_cluster_id()==7)
	  //   std::cout << "B: " <<  (!inside_fiducial_volume(max_point,offset_x)) << " " << inside_dead_region(max_point) << std::endl;
	  
	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
	      ){ // should not too close to anode 
	  }else{
	    return true;
	  }
	}
      } else if (count>=1){
      	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
      		       path_wcps_vec.at(i).y-wcp1.y,
      		       path_wcps_vec.at(i).z-wcp1.z);
      	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
      		       path_wcps_vec.at(i).y-wcp2.y,
      		       path_wcps_vec.at(i).z-wcp2.z);

	//	std::cout << "B: " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. << " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;
	
      	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 ||
	     (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
	    && temp1.Mag()>5*units::cm && temp2.Mag()>5*units::cm 
	    ){
	  //	  std::cout << "AAA" << std::endl; 

      	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
      	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
      	      ){ // should not too close to anode 
      	  }else{
      	    return true;
      	  }
      	}
      }
    }else{
      count = 0 ;
      max_angle = 0;
      max_point.x = 0;
      max_point.y = 0;
      max_point.z = 0;
    }
  }
  
  return false;
}
