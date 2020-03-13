void WCPPID::PR3DCluster::dQ_dx_multi_fit(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,  std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, double flash_time , double dis_end_point_ext , bool flag_dQ_dx_fit_reg ){
  
   // Need to take into account the time, so one can properly calculate X value for diffusion ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_offset = mp.get_time_offset();
  int nrebin = mp.get_nrebin();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w + 0.5;
  double offset_u = -first_u_dis/pitch_u + 0.5;
  double offset_v = -first_v_dis/pitch_v + 0.5;
  
  double first_t_dis = point_cloud->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - point_cloud->get_cloud().pts[0].x;
  double slope_xt = 1./time_slice_width;
  double offset_t =  first_t_dis/time_slice_width + 0.5;
  
  // get the correct flash time matching TPC side
  flash_time = flash_time - time_offset*units::microsecond; // us

  // given an x position value, we want to convert this to a drift time 
  // pos_x/time_slice_width * nrebin * 0.5*units::us // us
  // difference between these two numbers are the time in us ... 
  
  // Now figure out the diffusion coefficients
  // these are current  numbers in WCT, not sure what would be the values for data ... 
  double DL = 6.4 * pow(units::cm,2)/units::second ;
  double DT = 9.8 * pow(units::cm,2)/units::second ;

  // these are the transverse broading due to software filters in the wire dimension
  // these should be quadrature added to the
  // there is some cancellation of this effect with the effect of using average field response in the deconvolution ... 0-
  double col_sigma_w_T = 0.188060 * pitch_w*0.2; // units::mm
  double ind_sigma_u_T = 0.402993 * pitch_u*0.3; // units::mm
  double ind_sigma_v_T = 0.402993 * pitch_v*0.5; // units::mm

  double rel_uncer_ind = 0.075; // original 0.1
  double rel_uncer_col = 0.05; // original 0.035

  double add_uncer_ind = 0.0;
  double add_uncer_col = 300.0; // if 800 can reach balance with induction planes ...
  
  // this is the longitudinal filters in the time dimension ...
  double add_sigma_L = 1.428249  * time_slice_width / nrebin / 0.5; // units::mm

  // Point-like case ... 
  // these should be the expected values:
  // U, V, W, T,  1252.01, 3819.54, 6799.62, 1485.81
  // reco position 1252.02, 3819.63, 6799.67, 1485.79
  // good ...
  update_data_dQ_dx_fit(global_wc_map, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

  int n_3D_pos = 0;
  int n_2D_u = map_2D_ut_charge.size();
  int n_2D_v = map_2D_vt_charge.size();
  int n_2D_w = map_2D_wt_charge.size();
  
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    WCPPID::ProtoVertex *start_v = 0, *end_v = 0;
    for (auto it1=it->second.begin(); it1!=it->second.end(); it1++){
      WCPPID::ProtoVertex *vt = *it1;
      if ( vt->get_wcpt().index == sg->get_wcpt_vec().front().index){
	start_v = vt;
      }else if ( vt->get_wcpt().index == sg->get_wcpt_vec().back().index){
	end_v = vt;
      }
    }
    std::vector<WCP::Point >& pts = sg->get_point_vec();
    std::vector<int>& indices = sg->get_fit_index_vec();
    std::cout << pts.size() << " " << indices.size() << std::endl;
    for (size_t i = 0;i!=indices.size(); i++){
      if (i==0){
	if (start_v->get_fit_index()==-1){
	  start_v->set_fit_index(n_3D_pos);
	  n_3D_pos++;
	}
      }else if (i+1 == pts.size()){
	if (end_v->get_fit_index()==-1){
	  end_v->set_fit_index(n_3D_pos);
	  n_3D_pos++;
	}
      }else{
	indices.at(i) = n_3D_pos;
	n_3D_pos++;
      }
    }
  }
  //  std::cout << n_3D_pos << std::endl;

  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w), pred_data_u_2D(n_2D_u), pred_data_v_2D(n_2D_v), pred_data_w_2D(n_2D_w);
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;

  // initial values ...
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  for (int i = 0;i!=n_3D_pos;i++){
    pos_3D_init(i) = 5000 * 6; // single MIP ...
  }

  // regular flag ...
  std::vector<int> reg_flag_u(n_3D_pos,0), reg_flag_v(n_3D_pos,0), reg_flag_w(n_3D_pos,0);

  // start to fill in the data
  {
    int n_u = 0;
    for (auto it = map_2D_ut_charge.begin(); it!= map_2D_ut_charge.end(); it++){
      if (std::get<0>(it->second) > 0)
	data_u_2D(n_u) = std::get<0>(it->second)/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*rel_uncer_ind,2) + pow(add_uncer_ind,2));
      else
	data_u_2D(n_u) = 0;
      if (std::isnan(data_u_2D(n_u)))
	std::cout << "U: " << data_u_2D(n_u) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*rel_uncer_ind << " " << std::endl;
      n_u ++;
    }
    int n_v = 0;
    for (auto it = map_2D_vt_charge.begin(); it!= map_2D_vt_charge.end(); it++){
      if (std::get<0>(it->second) > 0)
	data_v_2D(n_v) = std::get<0>(it->second)/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*rel_uncer_ind,2) + pow(add_uncer_ind,2));
      else
	data_v_2D(n_v) = 0;
      if (std::isnan(data_v_2D(n_v)))
	std::cout << "V: " << data_v_2D(n_v) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*rel_uncer_ind << " " << std::endl;
      n_v ++;
    }
    int n_w = 0;
    for (auto it = map_2D_wt_charge.begin(); it!= map_2D_wt_charge.end(); it++){
      if (std::get<0>(it->second)>0)
	data_w_2D(n_w) = std::get<0>(it->second)/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*rel_uncer_col,2) + pow(add_uncer_col,2));
      else
	data_w_2D(n_w) = 0;
      if (std::isnan(data_w_2D(n_w)))
	std::cout << "W: " << data_w_2D(n_w) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*rel_uncer_col << " " << std::endl;


      //      std::cout << "W: " << data_w_2D(n_w) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*rel_uncer_col << " " << std::endl;
      
      n_w ++;
    }
  }

  // loop over segments ...
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;

    PointVector& pts = sg->get_point_vec();
    std::vector<int>& indices = sg->get_fit_index_vec();
    std::vector<double>& dx_vec = sg->get_dx_vec();
    
    //    std::cout<< pts.size() << " " << dx_vec.size() << std::endl;

    for (size_t i1=1; i1+1<pts.size();i1++){
      Point prev_rec_pos, next_rec_pos;
      Point curr_rec_pos = pts.at(i1);
      prev_rec_pos.x = (pts.at(i1-1).x + pts.at(i1).x)/2.;
      prev_rec_pos.y = (pts.at(i1-1).y + pts.at(i1).y)/2.;
      prev_rec_pos.z = (pts.at(i1-1).z + pts.at(i1).z)/2.;
      next_rec_pos.x = (pts.at(i1+1).x + pts.at(i1).x)/2.;
      next_rec_pos.y = (pts.at(i1+1).y + pts.at(i1).y)/2.;
      next_rec_pos.z = (pts.at(i1+1).z + pts.at(i1).z)/2.;

      dx.at(indices.at(i1)) = sqrt(pow(curr_rec_pos.x-prev_rec_pos.x,2)+pow(curr_rec_pos.y-prev_rec_pos.y,2)+pow(curr_rec_pos.z-prev_rec_pos.z,2))
	+ sqrt(pow(curr_rec_pos.x-next_rec_pos.x,2)+pow(curr_rec_pos.y-next_rec_pos.y,2)+pow(curr_rec_pos.z-next_rec_pos.z,2));

      std::vector<double> centers_U ;
      std::vector<double> centers_V ;
      std::vector<double> centers_W ;
      std::vector<double> centers_T ;
      std::vector<double> sigmas_T;
      std::vector<double> sigmas_U;
      std::vector<double> sigmas_V;
      std::vector<double> sigmas_W;
      std::vector<double> weights;

      for (int j=0;j!=5;j++){
	Point reco_pos;
	reco_pos.x = prev_rec_pos.x + (curr_rec_pos.x-prev_rec_pos.x)/5.*(j+0.5);
	reco_pos.y = prev_rec_pos.y + (curr_rec_pos.y-prev_rec_pos.y)/5.*(j+0.5);
	reco_pos.z = prev_rec_pos.z + (curr_rec_pos.z-prev_rec_pos.z)/5.*(j+0.5);
	double central_T = offset_t + slope_xt * reco_pos.x ;
	double central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z);
	double central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400;
	double central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800;
	double weight = sqrt(pow(prev_rec_pos.x-curr_rec_pos.x,2)+
			     pow(prev_rec_pos.y-curr_rec_pos.y,2)+
			     pow(prev_rec_pos.z-curr_rec_pos.z,2));
	
	double drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ;
	if (drift_time <50 *units::microsecond ) drift_time = 50* units::microsecond;

	//  std::cout << drift_time/units::microsecond << std::endl;
	double diff_sigma_L = sqrt(2* DL * drift_time);
	double diff_sigma_T = sqrt(2* DT * drift_time);
	
	double sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width;
	double sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u;
	double sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v;
	double sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w;

	centers_U.push_back(central_U);
	centers_V.push_back(central_V);
	centers_W.push_back(central_W);
	centers_T.push_back(central_T);
	
	weights.push_back(weight);
	
	sigmas_U.push_back(sigma_T_u);
	sigmas_V.push_back(sigma_T_v);
	sigmas_W.push_back(sigma_T_w);
	sigmas_T.push_back(sigma_L);

	reco_pos.x = next_rec_pos.x + (curr_rec_pos.x-next_rec_pos.x)/5.*(j+0.5);
	reco_pos.y = next_rec_pos.y + (curr_rec_pos.y-next_rec_pos.y)/5.*(j+0.5);
	reco_pos.z = next_rec_pos.z + (curr_rec_pos.z-next_rec_pos.z)/5.*(j+0.5);
	central_T = offset_t + slope_xt * reco_pos.x ;
	central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z);
	central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400;
	central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800;
	weight = sqrt(pow(next_rec_pos.x-curr_rec_pos.x,2)+
		      pow(next_rec_pos.y-curr_rec_pos.y,2)+
		      pow(next_rec_pos.z-curr_rec_pos.z,2));

	drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ;
	if (drift_time <50 *units::microsecond ) drift_time = 50* units::microsecond;
	//  std::cout << drift_time/units::microsecond << std::endl;
	diff_sigma_L = sqrt(2* DL * drift_time);
	diff_sigma_T = sqrt(2* DT * drift_time);
	
	sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width;
	sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u;
	sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v;
	sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w;

	centers_U.push_back(central_U);
	centers_V.push_back(central_V);
	centers_W.push_back(central_W);
	centers_T.push_back(central_T);
	
	weights.push_back(weight);
	
	sigmas_U.push_back(sigma_T_u);
	sigmas_V.push_back(sigma_T_v);
	sigmas_W.push_back(sigma_T_w);
	sigmas_T.push_back(sigma_L);
      }

      int n_u = 0;
      double sum_u = 0;
      for (auto it = map_2D_ut_charge.begin(); it!= map_2D_ut_charge.end(); it++){
	if (fabs(it->first.first - centers_U.front()) <= 10 && fabs(it->first.second - centers_T.front()) <= 10 ){
	  double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_U, sigmas_U, weights , 0 , 4);
	  sum_u += value;
	  // near dead channels ...
	  if (std::get<2>(it->second)==0 && value > 0) reg_flag_u.at(indices.at(i1)) = 1;
	  if (value > 0 && std::get<0>(it->second) >0 && std::get<2>(it->second)!=0){
	    // if (i!=143)
	    RU.insert(n_u,indices.at(i1)) = value/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*rel_uncer_ind,2) + pow(add_uncer_ind,2));
	    // if (i==143) std::cout << "U: " << it->first.first << " " << it->first.second << " " << value << std::endl;
	  }
	}
	n_u ++;
      }
      
      int n_v = 0;
      double sum_v = 0;
      for (auto it = map_2D_vt_charge.begin(); it!= map_2D_vt_charge.end(); it++){
	if (fabs(it->first.first +  2400 - centers_V.front()) <= 10 && fabs(it->first.second - centers_T.front()) <= 10 ){
	  double value = cal_gaus_integral_seg(it->first.second, it->first.first + 2400, centers_T, sigmas_T, centers_V, sigmas_V, weights , 0 , 4);
	  sum_v += value;
	  // near dead channels
	  if (std::get<2>(it->second)==0 && value > 0) {
	    // if (i==136) std::cout << it->first.first << " " << it->first.second << std::endl;
	    reg_flag_v.at(indices.at(i1)) = 1;
	  }
	  
	  if (value > 0 && std::get<0>(it->second) >0 && std::get<2>(it->second)!=0){
	    //	  if (i==136) std::cout << "V: " << it->first.first << " " << it->first.second << " " << value << std::endl;
	    // if (i!=143)
	    RV.insert(n_v,indices.at(i1)) = value/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*rel_uncer_ind,2) + pow(add_uncer_ind,2));
	  }
	}
	n_v ++;
      }

      int n_w = 0;
      double sum_w = 0;
      for (auto it = map_2D_wt_charge.begin(); it!= map_2D_wt_charge.end(); it++){
	if (fabs(it->first.first + 4800 - centers_W.front()) <= 10 && fabs(it->first.second - centers_T.front()) <= 10 ){
	  double value = cal_gaus_integral_seg(it->first.second, it->first.first + 4800,centers_T, sigmas_T, centers_W, sigmas_W, weights , 0 , 4);
	  sum_w += value;
	  // near dead channels ...
	  if (std::get<2>(it->second)==0 && value > 0) reg_flag_w.at(indices.at(i1)) = 1;
	  
	  // relevant && charge > 0 && not dead channel ...
	  if (value > 0 && std::get<0>(it->second) >0 && std::get<2>(it->second)!=0){
	    // if (i!=143)
	    RW.insert(n_w,indices.at(i1)) = value/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*rel_uncer_col,2) + pow(add_uncer_col,2));
	    //  if (i==147) std::cout << "W: " << it->first.first << " " << it->first.second << " " << value << std::endl;
	  }
	}
	n_w ++;
      }
      
      
    } // loop over points ...
  } // loop over segments ...
  
  // loop over vertices ...
  
  

  
  
  
}
