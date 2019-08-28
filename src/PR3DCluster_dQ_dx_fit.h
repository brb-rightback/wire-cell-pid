#include <Eigen/SVD>

double WireCellPID::PR3DCluster::cal_gaus_integral_seg(int tbin, int wbin, std::vector<double>& t_centers, std::vector<double>& t_sigmas, std::vector<double>& w_centers, std::vector<double>& w_sigmas, std::vector<double>& weights, int flag, double nsigma){
  double result = 0;
  double result1 = 0;

  for (size_t i=0;i!=t_centers.size();i++){
    result += cal_gaus_integral(tbin,wbin,t_centers.at(i), t_sigmas.at(i), w_centers.at(i), w_sigmas.at(i),flag,nsigma) * weights.at(i);
    result1 += weights.at(i);
  }

  result /= result1;
  
  return result;
}


double WireCellPID::PR3DCluster::cal_gaus_integral(int tbin, int wbin, double t_center, double t_sigma, double w_center, double w_sigma, int flag, double nsigma){
  // flag  = 0 no boundary effect, pure Gaussian, time or collection plane
  // flag  = 1 taking into account boundary effect for induction plane
  double result = 0;
  // time domain ...
  if (fabs((tbin + 0.5)-t_center) <= nsigma * t_sigma &&
      fabs((wbin + 0.5)-w_center) <= nsigma * w_sigma){
    // time dimension ...
    result = 0.5*(std::erf((tbin+1-t_center)/sqrt(2.)/t_sigma)-std::erf((tbin-t_center)/sqrt(2.)/t_sigma));

    if (flag ==0){
      result *= 0.5*(std::erf((wbin+1-w_center)/sqrt(2.)/w_sigma)-std::erf((wbin-w_center)/sqrt(2.)/w_sigma));
    }else if (flag==1){

      double x2 = wbin + 1.5;
      double x1 = wbin + 0.5;
      double x0 = w_center;
      double content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      double w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);
      //w1 = 0.5;
      // std::cout << w1 << std::endl;
      /* double w1 = (1./2./sqrt(2*3.1415926)/pow(x1-x2,2) * */
      /* 	    (2*exp(-pow(x0-x1,2)/2./pow(w_sigma,2)) * w_sigma *(x0+x1-2*x2) - 2 * exp(-pow(x0-x2,2)/2./pow(w_sigma,2))*w_sigma*(x0-x2)- */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x1-x0)/sqrt(2.)/w_sigma) + */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x2-x0)/sqrt(2.)/w_sigma)))/content1; */

      x2 = wbin + 0.5;
      x1 = wbin - 0.5;
      double content2 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      double w2 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);
      // std::cout << w2 << std::endl
      /* double w2 = (1./2./sqrt(2*3.1415926)/pow(x1-x2,2) * */
      /* 	    (2*exp(-pow(x0-x1,2)/2./pow(w_sigma,2)) * w_sigma *(x0+x1-2*x2) - 2 * exp(-pow(x0-x2,2)/2./pow(w_sigma,2))*w_sigma*(x0-x2)- */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x1-x0)/sqrt(2.)/w_sigma) + */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x2-x0)/sqrt(2.)/w_sigma)))/content2; */

      /* if (w1 > 0.5) { */
      /* 	w1 = sqrt(w1); */
      /* }else if (w1<0.5){ */
      /* 	w1 = 1-sqrt(1-w1); */
      /* } */
      /*  if (w2 > 0.5) { */
      /* 	w2 = sqrt(w1); */
      /*  }else if (w2 < 0.5){ */
      /* 	w2 = 1-sqrt(1-w2); */
      /* } */
      
      result *= (content1 * w1+content2 * (1-w2));
    }else if (flag==2){
      double sum = 0;
      
      double x2 = wbin + 1.0;
      double x1 = wbin + 0.5;
      double x0 = w_center;
      double content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      double w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum += content1 * (0.545 + 0.697*w1);

      x2 = wbin + 1.5;
      x1 = wbin + 1.0;
      content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum+= content1 *(0.11364 + 0.1 * w1);

      x2 = wbin + 0.5;
      x1 = wbin;
      content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum+= content1 *(0.545 + 0.697*(1-w1));

      x2 = wbin;
      x1 = wbin - 0.5;
      content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum+= content1 *(0.11364 + 0.1 * (1-w1));
      
      result *= sum;
    }
  }
  return result;
}


void WireCellPID::PR3DCluster::dQ_dx_fit(std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,  std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, double flash_time){
  if (fine_tracking_path.size()<=1) return;
  pu.clear();
  pv.clear();
  pw.clear();
  pt.clear();
  dQ.clear();
  dx.clear();

  
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
  
  double first_t_dis = path_wcps.front().mcell->GetTimeSlice()*time_slice_width - path_wcps.front().x;
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
  
  // this is the longitudinal filters in the time dimension ...
  double add_sigma_L = 1.428249  * time_slice_width / nrebin / 0.5; // units::mm

  // Point-like case ... 
  // these should be the expected values:
  // U, V, W, T,  1252.01, 3819.54, 6799.62, 1485.81
  // reco position 1252.02, 3819.63, 6799.67, 1485.79
  // good ...

  //std::set<int> good_channels_set =
  update_data_dQ_dx_fit(global_wc_map, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

  int n_3D_pos = fine_tracking_path.size();
  int n_2D_u = map_2D_ut_charge.size();
  int n_2D_v = map_2D_vt_charge.size();
  int n_2D_w = map_2D_wt_charge.size();

  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w), pred_data_u_2D(n_2D_u), pred_data_v_2D(n_2D_v), pred_data_w_2D(n_2D_w);
  
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  std::vector<int> reg_flag_u(n_3D_pos,0), reg_flag_v(n_3D_pos,0), reg_flag_w(n_3D_pos,0);
  // initialize 
  for (int i = 0;i!=n_3D_pos;i++){
    pos_3D_init(i) = 50000;
  }

  // start to fill in the data ...
  {
    int n_u = 0;
    for (auto it = map_2D_ut_charge.begin(); it!= map_2D_ut_charge.end(); it++){
      if (std::get<0>(it->second) > 0)
	data_u_2D(n_u) = std::get<0>(it->second)/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2));
      else
	data_u_2D(n_u) = 0;
      if (std::isnan(data_u_2D(n_u)))
	std::cout << "U: " << data_u_2D(n_u) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*0.035 << " " << std::endl;
      n_u ++;
    }
    int n_v = 0;
    for (auto it = map_2D_vt_charge.begin(); it!= map_2D_vt_charge.end(); it++){
      if (std::get<0>(it->second) > 0)
	data_v_2D(n_v) = std::get<0>(it->second)/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2));
      else
	data_v_2D(n_v) = 0;
      if (std::isnan(data_v_2D(n_v)))
	std::cout << "V: " << data_v_2D(n_v) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*0.035 << " " << std::endl;
      n_v ++;
    }
    int n_w = 0;
    for (auto it = map_2D_wt_charge.begin(); it!= map_2D_wt_charge.end(); it++){
      if (std::get<0>(it->second)>0)
	data_w_2D(n_w) = std::get<0>(it->second)/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.035,2));
      else
	data_w_2D(n_w) = 0;
      if (std::isnan(data_w_2D(n_w)))
	std::cout << "W: " << data_w_2D(n_w) << " " << std::get<1>(it->second) << " " << std::get<0>(it->second)*0.035 << " " << std::endl;
      n_w ++;
    }
  }
  
  for (int i=0;i!=n_3D_pos;i++){
    
    Point prev_rec_pos, next_rec_pos;
    Point curr_rec_pos(fine_tracking_path.at(i).x, fine_tracking_path.at(i).y, fine_tracking_path.at(i).z);

    if (i==0){
      next_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i+1).x)/2.;
      next_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i+1).y)/2.;
      next_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i+1).z)/2.;
      double length = sqrt(pow(fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x,2)
			   +pow(fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y,2)
			   +pow(fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z,2) );
      prev_rec_pos.x = fine_tracking_path.at(i).x - (fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x)/length * 1.5*units::mm;
      prev_rec_pos.y = fine_tracking_path.at(i).y - (fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y)/length * 1.5*units::mm ;
      prev_rec_pos.z = fine_tracking_path.at(i).z - (fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z)/length * 1.5*units::mm;
    }else if (i==n_3D_pos-1){
      prev_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i-1).x)/2.;
      prev_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i-1).y)/2.;
      prev_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i-1).z)/2.;
      double length = sqrt(pow(fine_tracking_path.at(i-1).x-fine_tracking_path.at(i).x,2)
			   +pow(fine_tracking_path.at(i-1).y-fine_tracking_path.at(i).y,2)
			   +pow(fine_tracking_path.at(i-1).z-fine_tracking_path.at(i).z,2) );
      next_rec_pos.x = fine_tracking_path.at(i).x - (fine_tracking_path.at(i-1).x-fine_tracking_path.at(i).x);///length * 1.5*units::mm;
      next_rec_pos.y = fine_tracking_path.at(i).y - (fine_tracking_path.at(i-1).y-fine_tracking_path.at(i).y);///length * 1.5*units::mm ;
      next_rec_pos.z = fine_tracking_path.at(i).z - (fine_tracking_path.at(i-1).z-fine_tracking_path.at(i).z);///length * 1.5*units::mm;
      
    }else{
      prev_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i-1).x)/2.;
      prev_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i-1).y)/2.;
      prev_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i-1).z)/2.;

      next_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i+1).x)/2.;
      next_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i+1).y)/2.;
      next_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i+1).z)/2.;
    }

    // not needed ... 
    /* curr_rec_pos.x = (curr_rec_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
    /* prev_rec_pos.x = (prev_rec_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
    /* next_rec_pos.x = (next_rec_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
    
    dx.push_back(sqrt(pow(curr_rec_pos.x-prev_rec_pos.x,2)+pow(curr_rec_pos.y-prev_rec_pos.y,2)+pow(curr_rec_pos.z-prev_rec_pos.z,2))
		 +sqrt(pow(curr_rec_pos.x-next_rec_pos.x,2)+pow(curr_rec_pos.y-next_rec_pos.y,2)+pow(curr_rec_pos.z-next_rec_pos.z,2)));

    pu.push_back(offset_u + (slope_yu * curr_rec_pos.y + slope_zu * curr_rec_pos.z));
    pv.push_back(offset_v + (slope_yv * curr_rec_pos.y + slope_zv * curr_rec_pos.z)+2400);
    pw.push_back(offset_w + (slope_yw * curr_rec_pos.y + slope_zw * curr_rec_pos.z)+4800);
    pt.push_back(offset_t + slope_xt * curr_rec_pos.x );
    
    
    
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
      if (fabs(it->first.first - centers_U.front()) <= 10 &&
      	  fabs(it->first.second - centers_T.front()) <= 10 ){
	double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_U, sigmas_U, weights , 0 , 4);
	sum_u += value;
	if (std::get<2>(it->second)==0) reg_flag_u.at(i) = 1;
	if (value > 0 && std::get<0>(it->second) >0 && std::get<2>(it->second)!=0){
	  RU.insert(n_u,i) = value/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2));
	}
      }
      n_u ++;
    }

    int n_v = 0;
    double sum_v = 0;
    for (auto it = map_2D_vt_charge.begin(); it!= map_2D_vt_charge.end(); it++){
      if (fabs(it->first.first +  2400 - centers_V.front()) <= 10 &&
	  fabs(it->first.second - centers_T.front()) <= 10 ){
	double value = cal_gaus_integral_seg(it->first.second, it->first.first + 2400, centers_T, sigmas_T, centers_V, sigmas_V, weights , 0 , 4);
	sum_v += value;
	if (std::get<2>(it->second)==0) reg_flag_v.at(i) = 1;
	if (value > 0 && std::get<0>(it->second) >0 && std::get<2>(it->second)!=0)
	  RV.insert(n_v,i) = value/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2));
      }
      n_v ++;
    }

    
    int n_w = 0;
    double sum_w = 0;
    for (auto it = map_2D_wt_charge.begin(); it!= map_2D_wt_charge.end(); it++){
      if (fabs(it->first.first + 4800 - centers_W.front()) <= 10 &&
      	  fabs(it->first.second - centers_T.front()) <= 10 ){
	double value = cal_gaus_integral_seg(it->first.second, it->first.first + 4800,centers_T, sigmas_T, centers_W, sigmas_W, weights , 0 , 4);
	sum_w += value;
	if (std::get<2>(it->second)==0) reg_flag_w.at(i) = 1;
	if (value > 0 && std::get<0>(it->second) >0 && std::get<2>(it->second)!=0)
	  RW.insert(n_w,i) = value/sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.035,2));
      }
      n_w ++;
    }
  }


  Eigen::SparseMatrix<double> RUT = Eigen::SparseMatrix<double>(RU.transpose());
  Eigen::SparseMatrix<double> RVT = Eigen::SparseMatrix<double>(RV.transpose());
  Eigen::SparseMatrix<double> RWT = Eigen::SparseMatrix<double>(RW.transpose());

  Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos);

  double reg_cut = 4*units::mm;
  double ind_weight = 0.2;
  double col_weight = 1;
  for (size_t i=0;i!=n_3D_pos;i++){
    bool flag_u = reg_flag_u.at(i);
    bool flag_v = reg_flag_v.at(i);
    bool flag_w = reg_flag_w.at(i);
    
    if (n_3D_pos!=1){
      if (i==0){
	double weight = 0;
	if (flag_u) weight += ind_weight;
	if (flag_v) weight += ind_weight;
	if (flag_w) weight += col_weight;

	FMatrix.insert(0,0) = -weight; 
	FMatrix.insert(0,1) = weight;
      }else if (i==n_3D_pos-1){
	double weight = 0;
	if (flag_u) weight += ind_weight;
	if (flag_v) weight += ind_weight;
	if (flag_w) weight += col_weight;

	FMatrix.insert(i,i) = -weight; 
	FMatrix.insert(i,i-1) = weight;
      }else{
	double weight = 0;
	if (flag_u) weight += ind_weight;
	if (flag_v) weight += ind_weight;
	if (flag_w) weight += col_weight;

	FMatrix.insert(i,i)=-2.*weight; 
	FMatrix.insert(i,i+1)=weight; 
	FMatrix.insert(i,i-1)=weight;
      }
    }
  }
  
  double lambda = 0.001;
  FMatrix *= lambda;
  Eigen::SparseMatrix<double> FMatrixT = Eigen::SparseMatrix<double>(FMatrix.transpose());
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;
  Eigen::SparseMatrix<double> A =  RUT * RU + RVT * RV + RWT * RW + FMatrixT * FMatrix;//
  solver.compute(A);
  
  pos_3D = solver.solveWithGuess(b,pos_3D_init);

  if (std::isnan(solver.error())){
    pos_3D = solver.solve(b);
  }

  double sum = 0 ;
  for (int i=0;i!=n_3D_pos;i++){
    /* std::cout << i << " "<< pos_3D(i) << " " << dx.at(i)/units::cm << " " << pos_3D(i)/dx.at(i)*units::cm << std::endl; */
    double central_U = offset_u + (slope_yu * fine_tracking_path.at(i).y + slope_zu * fine_tracking_path.at(i).z);
    if (central_U >=296 && central_U <=327 ||
	central_U >=336 && central_U <=337 ||
	central_U >=343 && central_U <=351 ||
	central_U >=376 && central_U <=400 ||
	central_U >=410 && central_U <=484 ||
	central_U >=501 && central_U <=524 ||
	central_U >=536 && central_U <=671)
      dQ.push_back(pos_3D(i)/0.7);
    else
      dQ.push_back(pos_3D(i));
    
    sum += dQ.back();
  }
  std::cout << "total: " << sum << std::endl;


  // prediction ...
  pred_data_u_2D = RU * pos_3D;
  pred_data_v_2D = RV * pos_3D;
  pred_data_w_2D = RW * pos_3D;

  int n_u = 0;
  for (auto it = map_2D_ut_charge.begin(); it!=map_2D_ut_charge.end(); it++){
    proj_data_u_map[std::make_pair(it->first.first, it->first.second)] = std::make_tuple(std::get<0>(it->second), std::get<1>(it->second), pred_data_u_2D(n_u) * sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2)));
    n_u++;
  }
  int n_v = 0;
  for (auto it = map_2D_vt_charge.begin(); it!=map_2D_vt_charge.end(); it++){
    proj_data_v_map[std::make_pair(it->first.first+2400, it->first.second)] = std::make_tuple(std::get<0>(it->second), std::get<1>(it->second), pred_data_v_2D(n_v) * sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2)));
    n_v++;
  }
  int n_w = 0;
  for (auto it = map_2D_wt_charge.begin(); it!=map_2D_wt_charge.end(); it++){
    proj_data_w_map[std::make_pair(it->first.first+4800, it->first.second)] = std::make_tuple(std::get<0>(it->second), std::get<1>(it->second), pred_data_w_2D(n_w) * sqrt(pow(std::get<1>(it->second),2)+pow(std::get<0>(it->second)*0.1,2)));
    n_w++;
  }

  
  
  
}

void WireCellPID::PR3DCluster::update_data_dQ_dx_fit(std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_wt_charge){
  
  std::set<SlimMergeGeomCell*> cluster_mcells_set;
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = *it;
    cluster_mcells_set.insert(mcell);
  }

  for (auto it=mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    std::map<const GeomWire*, SMGCSelection >& timeslice_wc_map = global_wc_map[time_slice];

    WireChargeMap& wire_charge_map = mcell->get_wirecharge_map();
    WireChargeMap& wire_charge_err_map = mcell->get_wirechargeerr_map();
    for (auto it = wire_charge_map.begin(); it!= wire_charge_map.end(); it++){
      const GeomWire* wire = it->first;
      double charge = it->second;
      double charge_err = wire_charge_err_map[wire];
      
      if (charge > 0)
	if (timeslice_wc_map[wire].size()>1){
	  for (auto it1 = timeslice_wc_map[wire].begin(); it1!=timeslice_wc_map[wire].end(); it1++){
	    SlimMergeGeomCell *mcell1 = *it1;
	    if (cluster_mcells_set.find(mcell1)==cluster_mcells_set.end()) charge_err = 8000;
	  }
	}

      if (wire->iplane()==0){
	auto it = map_2D_ut_charge.find(std::make_pair(wire->index(),time_slice));
	if (it !=map_2D_ut_charge.end()){
	  it->second = std::make_tuple(charge, charge_err, std::get<2>(it->second));
	}else{
	  std::cout << "Mismatch U" << std::endl;
	}
      }else if (wire->iplane()==1){
	auto it = map_2D_vt_charge.find(std::make_pair(wire->index(),time_slice));
	if (it !=map_2D_vt_charge.end()){
	  it->second = std::make_tuple(charge, charge_err, std::get<2>(it->second));
	}else{
	  std::cout << "Mismatch V" << std::endl;
	}
      }else{
	auto it = map_2D_wt_charge.find(std::make_pair(wire->index(),time_slice));
	if (it !=map_2D_wt_charge.end()){
	  it->second = std::make_tuple(charge, charge_err, std::get<2>(it->second));
	}else{
	  std::cout << "Mismatch W" << std::endl;
	}
      }

      
    }
  }

  // fill in additional stuff ... 
  /* for (auto it = map_2D_ut_charge.begin(); it!=map_2D_ut_charge.end(); it++){  */
  /*   if (std::get<2>(it->second)==0) continue; */
  /*   if (proj_data_u_map.find(it->first)==proj_data_u_map.end()){ */
  /*     proj_data_u_map[it->first] = std::make_tuple(std::get<0>(it->second),std::get<1>(it->second),0); */
  /*   } */
  /* } */
  /* for (auto it = map_2D_vt_charge.begin(); it!=map_2D_vt_charge.end(); it++){  */
  /*   if (std::get<2>(it->second)==0) continue; */
  /*   if (proj_data_v_map.find(std::make_pair(it->first.first+2400,it->first.second))==proj_data_v_map.end()){ */
  /*     proj_data_v_map[std::make_pair(it->first.first+2400,it->first.second)] = std::make_tuple(std::get<0>(it->second),std::get<1>(it->second),0); */
  /*   } */
  /* } */
  /* for (auto it = map_2D_wt_charge.begin(); it!=map_2D_wt_charge.end(); it++){  */
  /*   if (std::get<2>(it->second)==0) continue; */
  /*   if (proj_data_w_map.find(std::make_pair(it->first.first+4800,it->first.second))==proj_data_w_map.end()){ */
  /*     proj_data_w_map[std::make_pair(it->first.first+4800,it->first.second)] = std::make_tuple(std::get<0>(it->second),std::get<1>(it->second),0); */
  /*   } */
  /* } */
 
  
}
