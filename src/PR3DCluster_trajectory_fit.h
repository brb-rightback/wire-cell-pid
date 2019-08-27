void WireCellPID::PR3DCluster::fill_data_map_trajectory(std::vector<int> indices, std::map<int, std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set, std::map<std::pair<int,int>,  std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge){
  proj_data_u_map.clear();
  proj_data_v_map.clear();
  proj_data_w_map.clear();

  if (indices.size()==0){
    for (auto it = map_2D_ut_charge.begin(); it!=map_2D_ut_charge.end(); it++){
      //good_channels_set.insert(it->first.first);
      proj_data_u_map[it->first] = std::make_tuple(std::get<0>(it->second),std::get<1>(it->second),0);
    }
    for (auto it = map_2D_vt_charge.begin(); it!=map_2D_vt_charge.end(); it++){
      //good_channels_set.insert(it->first.first+2400);
      proj_data_v_map[std::make_pair(it->first.first+2400,it->first.second)] = std::make_tuple(std::get<0>(it->second),std::get<1>(it->second),0);
    }
    for (auto it = map_2D_wt_charge.begin(); it!=map_2D_wt_charge.end(); it++){
      //good_channels_set.insert(it->first.first+4800);
      proj_data_w_map[std::make_pair(it->first.first+4800,it->first.second)] = std::make_tuple(std::get<0>(it->second),std::get<1>(it->second),0);
    }
  }else{
    for (auto it = indices.begin(); it!=indices.end(); it++){
      //      std::cout << pu.at(*it) << " " << pv.at(*it) << " " << pw.at(*it) << " " << pt.at(*it) << std::endl;
      for (auto it1 = map_3D_2DU_set[*it].first.begin(); it1!=map_3D_2DU_set[*it].first.end(); it1++){
	proj_data_u_map[*it1] = std::make_tuple(std::get<0>(map_2D_ut_charge[*it1]),std::get<1>(map_2D_ut_charge[*it1]),0);
	//std::cout << "U: " << it1->first << " " << it1->second << std::endl; 
      }
      for (auto it1 = map_3D_2DV_set[*it].first.begin(); it1!=map_3D_2DV_set[*it].first.end(); it1++){
	proj_data_v_map[std::make_pair(it1->first+2400,it1->second)] = std::make_tuple(std::get<0>(map_2D_vt_charge[*it1]),std::get<1>(map_2D_vt_charge[*it1]),0);
	//std::cout << "V: " << it1->first+2400 << " " << it1->second << std::endl; 
      }
      for (auto it1 = map_3D_2DW_set[*it].first.begin(); it1!=map_3D_2DW_set[*it].first.end(); it1++){
	proj_data_w_map[std::make_pair(it1->first+4800,it1->second)] = std::make_tuple(std::get<0>(map_2D_wt_charge[*it1]),std::get<1>(map_2D_wt_charge[*it1]),0);
	//std::cout << "W: " << it1->first+4800 << " " << it1->second << std::endl; 
      }
    }
    
  }
}

void WireCellPID::PR3DCluster::trajectory_fit(WireCell::PointVector& ps_vec, std::map<int, std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,std::map<std::pair<int,int>, std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, int charge_div_method, double div_sigma){
  
  // calculate the distance between points ...
  std::vector<double> distances;
  for (size_t i=0;i+1!=ps_vec.size();i++){
    distances.push_back(sqrt(pow(ps_vec.at(i+1).x-ps_vec.at(i).x,2) +
			     pow(ps_vec.at(i+1).y-ps_vec.at(i).y,2) +
			     pow(ps_vec.at(i+1).z-ps_vec.at(i).z,2)));
  }

  // map 2D points to its index
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DU_index;
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DV_index;
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DW_index;

  for (auto it = map_2DU_3D_set.begin(); it!= map_2DU_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DU_index.push_back(std::make_pair(it->first,*it1));
    }
  }
  for (auto it = map_2DV_3D_set.begin(); it!= map_2DV_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DV_index.push_back(std::make_pair(it->first,*it1));
    }
  }
  for (auto it = map_2DW_3D_set.begin(); it!= map_2DW_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DW_index.push_back(std::make_pair(it->first,*it1));
    }
  }

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();


  double slope_x = 1./time_slice_width;
  double first_t_dis = path_wcps.front().mcell->GetTimeSlice()*time_slice_width - path_wcps.front().x;
  double offset_t = first_t_dis / time_slice_width;

  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w;
  double offset_u = -first_u_dis/pitch_u;
  double offset_v = -first_v_dis/pitch_v;

   // form matrix ...
  int n_3D_pos = 3 * ps_vec.size();
  int n_2D_u = 2 * vec_2DU_index.size();
  int n_2D_v = 2 * vec_2DV_index.size();
  int n_2D_w = 2 * vec_2DW_index.size();

  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  for (size_t i=0;i!=ps_vec.size();i++){
    pos_3D_init(3*i) = ps_vec.at(i).x;
    pos_3D_init(3*i+1) = ps_vec.at(i).y;
    pos_3D_init(3*i+2) = ps_vec.at(i).z;
  }


  

  // 2D pixel and then the 3D index ...
  std::map<std::tuple<int,int, int>, double> map_Udiv_fac;
  std::map<std::tuple<int,int, int>, double> map_Vdiv_fac;
  std::map<std::tuple<int,int, int>, double> map_Wdiv_fac;
  
  // charge division method ...
  if (charge_div_method==1){
    // equal division
    for (auto it = map_2DU_3D_set.begin(); it!=map_2DU_3D_set.end(); it++){
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Udiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = 1./it->second.size();
  	//std::cout << it->first.first << " " << it->first.second << " " << *it1 << std::endl;
      }
    }
    for (auto it = map_2DV_3D_set.begin(); it!=map_2DV_3D_set.end(); it++){
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Vdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = 1./it->second.size();
  	//std::cout << it->first.first << " " << it->first.second << " " << *it1 << std::endl;
      }
    }
    for (auto it = map_2DW_3D_set.begin(); it!=map_2DW_3D_set.end(); it++){
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Wdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = 1./it->second.size();
  	//	std::cout << it->first.first << " " << it->first.second << " " << *it1 << std::endl;
      }
    }
  }else if (charge_div_method == 2){
    // use div_sigma ...
    for (auto it = map_2DU_3D_set.begin(); it!=map_2DU_3D_set.end(); it++){
      double sum = 0;
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	double central_t = slope_x * ps_vec.at(*it1).x + offset_t;
  	double central_ch = slope_yu * ps_vec.at(*it1).y + slope_zu * ps_vec.at(*it1).z + offset_u;
  	double factor = exp(-0.5 * (pow((central_t-it->first.second)*time_slice_width,2) + pow((central_ch-it->first.first)*pitch_u,2))/pow(div_sigma,2));
  	map_Udiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = factor;
  	sum += factor;
  	//std::cout << it->first.first << " " << it->first.second << " " << *it1 << std::endl;
      }
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Udiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] /= sum;
      }
    }

    for (auto it = map_2DV_3D_set.begin(); it!=map_2DV_3D_set.end(); it++){
      double sum = 0;
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	double central_t = slope_x * ps_vec.at(*it1).x + offset_t;
  	double central_ch = slope_yv * ps_vec.at(*it1).y + slope_zv * ps_vec.at(*it1).z + offset_v;
  	double factor = exp(-0.5 * (pow((central_t-it->first.second)*time_slice_width,2) + pow((central_ch-it->first.first)*pitch_v,2))/pow(div_sigma,2));
  	map_Vdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = factor;
  	sum += factor;
  	//std::cout << it->first.first << " " << it->first.second << " " << *it1 << std::endl;
      }
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Vdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] /= sum;
      }
    }

    for (auto it = map_2DW_3D_set.begin(); it!=map_2DW_3D_set.end(); it++){
      double sum = 0;
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	double central_t = slope_x * ps_vec.at(*it1).x + offset_t;
  	double central_ch = slope_yw * ps_vec.at(*it1).y + slope_zw * ps_vec.at(*it1).z + offset_w;
  	double factor = exp(-0.5 * (pow((central_t-it->first.second)*time_slice_width,2) + pow((central_ch-it->first.first)*pitch_w,2))/pow(div_sigma,2));
  	map_Wdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = factor;
  	sum += factor;
  	//std::cout << it->first.first << " " << it->first.second << " " << *it1 << " " << factor << std::endl;
      }
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Wdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] /= sum;
      }
    }
  }
  
  // additional error calculation based on mixture of dead channels
  for (size_t index = 0; index!=vec_2DU_index.size(); index++){
    double charge = std::get<0>(map_2D_ut_charge[vec_2DU_index.at(index).first]);
    double charge_err = std::get<1>(map_2D_ut_charge[vec_2DU_index.at(index).first]);
    int n_divide = map_2DU_3D_set[vec_2DU_index.at(index).first].size();
    // induction plane error is large, so they have less weights in the fit to decide X position
    double scaling = charge/charge_err * map_Udiv_fac[std::make_tuple(vec_2DU_index.at(index).first.first, vec_2DU_index.at(index).first.second, vec_2DU_index.at(index).second)];

    int index_3D = vec_2DU_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z

    if (map_3D_2DU_set[index_3D].second!=0 && map_3D_2DU_set[index_3D].first.size()!=0)
      scaling *= (map_3D_2DU_set[index_3D].first.size()*1.0 - map_3D_2DU_set[index_3D].second)/map_3D_2DU_set[index_3D].first.size();

    //    std::cout << scaling << std::endl;
    
    data_u_2D(2*index) =  scaling * (vec_2DU_index.at(index).first.first - offset_u);
    data_u_2D(2*index+1) = scaling * (vec_2DU_index.at(index).first.second - offset_t);

    RU.insert(2*index,3*index_3D+1) = scaling * slope_yu; // Y--> U
    RU.insert(2*index,3*index_3D+2) = scaling * slope_zu; // Z--> U
    RU.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
  }

  for (size_t index = 0; index!=vec_2DV_index.size(); index++){
    double charge = std::get<0>(map_2D_vt_charge[vec_2DV_index.at(index).first]);
    double charge_err = std::get<1>(map_2D_vt_charge[vec_2DV_index.at(index).first]);
    int n_divide = map_2DV_3D_set[vec_2DV_index.at(index).first].size();
    // induction plane error is large, so they have less weights in the fit to decide X position
    double scaling = charge/charge_err * map_Vdiv_fac[std::make_tuple(vec_2DV_index.at(index).first.first, vec_2DV_index.at(index).first.second, vec_2DV_index.at(index).second)];

    int index_3D = vec_2DV_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z

    if (map_3D_2DV_set[index_3D].second!=0 && map_3D_2DV_set[index_3D].first.size()!=0)
      scaling *= (map_3D_2DV_set[index_3D].first.size()*1.0 - map_3D_2DV_set[index_3D].second)/map_3D_2DV_set[index_3D].first.size();

    //    std::cout << scaling << std::endl;
    
    data_v_2D(2*index) =  scaling * (vec_2DV_index.at(index).first.first - offset_v);
    data_v_2D(2*index+1) = scaling * (vec_2DV_index.at(index).first.second - offset_t);

    RV.insert(2*index,3*index_3D+1) = scaling * slope_yv; // Y--> V
    RV.insert(2*index,3*index_3D+2) = scaling * slope_zv; // Z--> V
    RV.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
  }

  for (size_t index = 0; index!=vec_2DW_index.size(); index++){
    double charge = std::get<0>(map_2D_wt_charge[vec_2DW_index.at(index).first]);
    double charge_err = std::get<1>(map_2D_wt_charge[vec_2DW_index.at(index).first]);
    int n_divide = map_2DW_3D_set[vec_2DW_index.at(index).first].size();
    // induction plane error is large, so they have less weights in the fit to decide X position
    double scaling = charge/charge_err * map_Wdiv_fac[std::make_tuple(vec_2DW_index.at(index).first.first, vec_2DW_index.at(index).first.second, vec_2DW_index.at(index).second)];

    int index_3D = vec_2DW_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z

    if (map_3D_2DW_set[index_3D].second!=0 && map_3D_2DW_set[index_3D].first.size()!=0)
      scaling *= (map_3D_2DW_set[index_3D].first.size()*1.0 - map_3D_2DW_set[index_3D].second)/map_3D_2DW_set[index_3D].first.size();

    //    std::cout << scaling << std::endl;
    
    data_w_2D(2*index) =  scaling * (vec_2DW_index.at(index).first.first - offset_w);
    data_w_2D(2*index+1) = scaling * (vec_2DW_index.at(index).first.second - offset_t);

    RW.insert(2*index,3*index_3D+2) = scaling * slope_zw; // Z--> W
    RW.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
  }
  
  Eigen::SparseMatrix<double> RUT = Eigen::SparseMatrix<double>(RU.transpose());
  Eigen::SparseMatrix<double> RVT = Eigen::SparseMatrix<double>(RV.transpose());
  Eigen::SparseMatrix<double> RWT = Eigen::SparseMatrix<double>(RW.transpose());
  
  // when to apply regularization ... how???
  Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos) ;
   for (size_t i=0;i!=ps_vec.size();i++){
    
    if (i==0){
      if (map_3D_2DU_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ||
	  map_3D_2DU_set[i].first.size()==0 && map_3D_2DW_set[i].first.size()==0 ||
	  map_3D_2DW_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ){
	FMatrix.insert(0,0) = -1./distances.at(0); // X
	FMatrix.insert(0,3) = 1./distances.at(0);
	
	FMatrix.insert(1,1) = -1./distances.at(0); // Y
	FMatrix.insert(1,4) = 1./distances.at(0);
	
	FMatrix.insert(2,2) = -1./distances.at(0); // Z
	FMatrix.insert(2,5) = 1./distances.at(0);
      }
    }else if (i+1==ps_vec.size()){
      if (map_3D_2DU_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ||
	  map_3D_2DU_set[i].first.size()==0 && map_3D_2DW_set[i].first.size()==0 ||
	  map_3D_2DW_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ){
	FMatrix.insert(3*i,3*i) = -1./distances.at(ps_vec.size()-2); // X
	FMatrix.insert(3*i,3*i-3) = 1./distances.at(ps_vec.size()-2);

	FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(ps_vec.size()-2);
	FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(ps_vec.size()-2);

	FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(ps_vec.size()-2);
	FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(ps_vec.size()-2);
      }
    }else{
      if (map_3D_2DU_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ||
	  map_3D_2DU_set[i].first.size()==0 && map_3D_2DW_set[i].first.size()==0 ||
	  map_3D_2DW_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ){
	FMatrix.insert(3*i,3*i-3) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
	FMatrix.insert(3*i,3*i) = -1./distances.at(i-1)/distances.at(i);
	FMatrix.insert(3*i,3*i+3) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
	
	FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
	FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(i-1)/distances.at(i);
	FMatrix.insert(3*i+1,3*i+4) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
	
	FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
	FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(i-1)/distances.at(i);
	FMatrix.insert(3*i+2,3*i+5) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
      }
    }
  }

  double lambda = 0;
  lambda = 2 // strength 
    * sqrt(9. //  average chi2 guessted
 	   * (vec_2DU_index.size() + vec_2DV_index.size() + vec_2DW_index.size()) // how many of them
 	   * 6 * 6 // charge/charge_err estimation ... 
 	   /(ps_vec.size() * 1.)); //weighting
  double angle_range = 0.25;
  FMatrix *= lambda/angle_range ; // disable the angle cut ... 
  
  
  Eigen::SparseMatrix<double> FMatrixT = Eigen::SparseMatrix<double>(FMatrix.transpose());
  // Eigen::SparseMatrix<double> PMatrixT = Eigen::SparseMatrix<double>(PMatrix.transpose());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;// + PMatrixT * pos_3D_init * pow(lambda/dis_range,2);

  Eigen::SparseMatrix<double> A =   RUT * RU + RVT * RV + RWT * RW + FMatrixT * FMatrix;// + PMatrixT * PMatrix * pow(lambda/dis_range,2);

  //  std::cout << "Solve1 " << std::endl;
  
  solver.compute(A);

  //  std::cout << "Solve " << std::endl;
  
  pos_3D = solver.solveWithGuess(b,pos_3D_init);
  
  if (std::isnan(solver.error())){
    pos_3D = solver.solve(b);
    //std::cout << "#iterations: " << solver.iterations() << std::endl;
    //std::cout << "#estimated error: " << solver.error() << std::endl;
  }
  
  //std::cout << path_wcps_vec.size() << " " << map_2DU_index.size() << " " << map_2DV_index.size() << " " << map_2DW_index.size() << std::endl;

  if (std::isnan(solver.error())){
    //std::cout << "Bad fit" << std::endl;
    fine_tracking_path.clear();
    pu.clear();
    pv.clear();
    pw.clear();
    pt.clear();
    if (fine_tracking_path.size()==0){
      for (size_t i=0;i!=ps_vec.size();i++){
	Point p;
	p.x = ps_vec.at(i).x;
	p.y = ps_vec.at(i).y;
	p.z = ps_vec.at(i).z;
	fine_tracking_path.push_back(p);

	pu.push_back(offset_u + 0.5 + (slope_yu * p.y + slope_zu * p.z));
	pv.push_back(offset_v + 0.5 + (slope_yv * p.y + slope_zv * p.z)+2400);
	pw.push_back(offset_w + 0.5 + (slope_yw * p.y + slope_zw * p.z)+4800);
	pt.push_back(offset_t + 0.5 + slope_x * p.x );
      }
    }
  }else{
    // std::cout << "Good fit" << std::endl;
    flag_fine_tracking = true;
    fine_tracking_path.clear();
    pu.clear();
    pv.clear();
    pw.clear();
    pt.clear();
    for (size_t i=0;i!=ps_vec.size();i++){
      Point p;
      p.x = pos_3D(3*i);
      p.y = pos_3D(3*i+1);
      p.z = pos_3D(3*i+2);

      /* p.x = ps_vec.at(i).x; */
      /* p.y = ps_vec.at(i).y; */
      /* p.z = ps_vec.at(i).z; */
	
      if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)){
	//std::cout << "gaga " << std::endl;
      }else{
	fine_tracking_path.push_back(p);
	//fine_tracking_path.push_back(ps_vec.at(i));
	pu.push_back(offset_u + 0.5 + (slope_yu * p.y + slope_zu * p.z));
	pv.push_back(offset_v + 0.5 + (slope_yv * p.y + slope_zv * p.z)+2400);
	pw.push_back(offset_w + 0.5 + (slope_yw * p.y + slope_zw * p.z)+4800);
	pt.push_back(offset_t + 0.5 + slope_x * p.x );
	//	std::cout << ps_vec.at(i) << " " << p << " " << sqrt(pow(ps_vec.at(i).x-p.x,2)+pow(ps_vec.at(i).y-p.y,2)+pow(ps_vec.at(i).z-p.z,2))/units::mm <<std::endl;
      }
    }
  }
}


std::vector<int> WireCellPID::PR3DCluster::examine_point_association(std::set<std::pair<int,int> >& temp_2dut, std::set<std::pair<int,int> >& temp_2dvt, std::set<std::pair<int,int> >& temp_2dwt,
								     std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge, double charge_cut){

  std::set<int> temp_types_u;
  std::set<int> temp_types_v;
  std::set<int> temp_types_w;
    
  std::set<std::pair<int,int> > saved_2dut;
  std::set<std::pair<int,int> > saved_2dvt;
  std::set<std::pair<int,int> > saved_2dwt;

  std::vector<int> results;
  results.resize(3,0);
  
  for (auto it = temp_2dut.begin(); it!=temp_2dut.end(); it++){
    auto it1 = map_2D_ut_charge.find(*it);
    if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > charge_cut ){
      temp_types_u.insert(std::get<2>(it1->second));
      if (std::get<2>(it1->second)==0) results.at(0)++;
      saved_2dut.insert(*it);
    }
  }

  for (auto it = temp_2dvt.begin(); it!=temp_2dvt.end(); it++){
    auto it1 = map_2D_vt_charge.find(*it);
    if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > charge_cut ){
      temp_types_v.insert(std::get<2>(it1->second));
      if (std::get<2>(it1->second)==0) results.at(1)++;
      saved_2dvt.insert(*it);
    }
  }

  for (auto it = temp_2dwt.begin(); it!=temp_2dwt.end(); it++){
    auto it1 = map_2D_wt_charge.find(*it);
    if (it1!=map_2D_wt_charge.end() && std::get<0>(it1->second) > charge_cut ){
      temp_types_w.insert(std::get<2>(it1->second));
      if (std::get<2>(it1->second)==0) results.at(2)++;
      saved_2dwt.insert(*it);
    }
  }

  if (temp_types_u.find(0)!=temp_types_u.end() && temp_types_u.size()==1){
    saved_2dut.clear();
    results.at(0) = 0;
  }
  if (temp_types_v.find(0)!=temp_types_v.end() && temp_types_v.size()==1){
    saved_2dvt.clear();
    results.at(1) = 0;
  }
  if (temp_types_w.find(0)!=temp_types_w.end() && temp_types_w.size()==1){
    saved_2dwt.clear();
    results.at(2) = 0;
  }
  temp_2dut = saved_2dut;
  temp_2dvt = saved_2dvt;
  temp_2dwt = saved_2dwt;

  return results;
}

void WireCellPID::PR3DCluster::form_point_association(WireCell::Point &p, std::set<std::pair<int,int> >& temp_2dut, std::set<std::pair<int,int> >& temp_2dvt, std::set<std::pair<int,int> >& temp_2dwt, WireCell::ToyCTPointCloud& ct_point_cloud, double dis_cut, int nlevel, double time_cut ){
  // global information
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  
  float coef1 = 2 * pow(sin(angle_u),2);
  float coef2 = 2 * (pow(sin(angle_u),2) - pow(cos(angle_u),2));

  typedef boost::property_map<MCUGraph, boost::vertex_index_t>::type IndexMap;
  typedef boost::graph_traits<MCUGraph>::adjacency_iterator adjacency_iterator;
  
  // original point cloud ...
  if (point_cloud!=0 && graph!=0){
    WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
    IndexMap index = get(boost::vertex_index,*graph);
    if (cloud.pts.size()>0){
      WireCell::WCPointCloud<double>::WCPoint wcp = point_cloud->get_closest_wcpoint(p);
      double temp_dis = sqrt(pow(wcp.x-p.x,2)+pow(wcp.y-p.y,2)+pow(wcp.z-p.z,2));
      
      //std::cout << temp_dis/units::cm << " " << dis_cut/units::cm << std::endl;
      if (temp_dis < dis_cut){
	std::set<int> total_vertices_found;
	std::set<int> vertices_to_be_examined;
	std::set<int> vertices_saved_for_next;
	total_vertices_found.insert(wcp.index);
	vertices_to_be_examined.insert(wcp.index);

	//	std::cout << cloud.pts.size() << std::endl;
	
	for (int j=0;j!=nlevel;j++){
	  for (auto it = vertices_to_be_examined.begin(); it!=vertices_to_be_examined.end(); it++){
	    int temp_current_index = (*it);

	    //	    std::cout << temp_current_index << " " << cloud.pts.size() << std::endl;

	    std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(temp_current_index,*graph),*graph);
	    for (; neighbors.first!=neighbors.second; ++neighbors.first){
	      //std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
	      if (total_vertices_found.find(index(*neighbors.first))==total_vertices_found.end()){
		total_vertices_found.insert(index(*neighbors.first));
		vertices_saved_for_next.insert(index(*neighbors.first));
	      }
	    }
	  }
	  vertices_to_be_examined = vertices_saved_for_next;
	}
	
	SMGCSet nearby_mcells_set;
	for (auto it = total_vertices_found.begin(); it!=total_vertices_found.end(); it++){
	  SlimMergeGeomCell *mcell = cloud.pts[*it].mcell;
	  nearby_mcells_set.insert(mcell);
	}
	
	int cur_time_slice = cloud.pts[wcp.index].mcell->GetTimeSlice();
	int cur_wire_u = cloud.pts[wcp.index].index_u;
	int cur_wire_v = cloud.pts[wcp.index].index_v;
	int cur_wire_w = cloud.pts[wcp.index].index_w;

	//	std::cout << "A: " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_w << std::endl;
	
	double dis_cut_u = dis_cut;
	double dis_cut_v = dis_cut;
	double dis_cut_w = dis_cut;
	
	
	double max_time_slice_u = 0;
	double max_time_slice_v = 0;
	double max_time_slice_w = 0;
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_u)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
		max_time_slice_u = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_v)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
		max_time_slice_v = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_w)<=1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
		max_time_slice_w = fabs(this_time_slice-cur_time_slice);
	  }
	  
	}
	
	if (max_time_slice_u * time_slice_width*1.2 < dis_cut_u)
	  dis_cut_u = max_time_slice_u * time_slice_width*1.2;
	if (max_time_slice_v * time_slice_width*1.2 < dis_cut_v)
	  dis_cut_v = max_time_slice_v * time_slice_width*1.2;
	if (max_time_slice_w * time_slice_width*1.2 < dis_cut_w)
	  dis_cut_w = max_time_slice_w * time_slice_width*1.2;
	
	
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  
	  double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    float min_u_dis;
	    if (cur_wire_u < uwires.front()->index()){
	      min_u_dis = uwires.front()->index()-cur_wire_u;
	    }else if (cur_wire_u >= uwires.front()->index() &&
		      cur_wire_u <= uwires.back()->index()){
	      min_u_dis = 0;
	    }else{
	      min_u_dis = cur_wire_u-uwires.back()->index();
	    }
	    float min_v_dis;
	    if (cur_wire_v < vwires.front()->index()){
	      min_v_dis = vwires.front()->index()-cur_wire_v;
	    }else if (cur_wire_v >= vwires.front()->index() &&
		      cur_wire_v <= vwires.back()->index()){
	      min_v_dis = 0;
	    }else{
	      min_v_dis = cur_wire_v-vwires.back()->index();
	    }
	    float min_w_dis;
	    if (cur_wire_w < wwires.front()->index()){
	      min_w_dis = wwires.front()->index()-cur_wire_w;
	    }else if (cur_wire_w >= wwires.front()->index() &&
		      cur_wire_w <= wwires.back()->index()){
	      min_w_dis = 0;
	    }else{
	      min_w_dis = cur_wire_w-wwires.back()->index();
	    }
	    float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	    
	    if ( range_u > 0 && range_v >0 && range_w > 0){
	      float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
	      float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
	      float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
	      float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
	      float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
	      float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;
	      for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
		//auto it1 = map_2D_ut_charge.find(std::make_pair(j,this_time_slice));
		//if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > 0 )
		temp_2dut.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
		//auto it1 = map_2D_vt_charge.find(std::make_pair(j,this_time_slice));
		//if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > 0 )
		temp_2dvt.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
		//auto it1 = map_2D_wt_charge.find(std::make_pair(j,this_time_slice));
		//if (it1!=map_2D_wt_charge.end() && std::get<0>(it1->second) > 0 )
		temp_2dwt.insert(std::make_pair(j,this_time_slice));
	      }
	    } // cuts
	  } // cuts
	} // loop over all mcells
      } // find a good closest point?
    } // point cloud exist
  }
    
  // steiner tree point cloud ...
  if (point_cloud_steiner!=0 && graph_steiner!=0 ){
    WireCell::WCPointCloud<double>& cloud = point_cloud_steiner->get_cloud();
    IndexMap index = get(boost::vertex_index,*graph_steiner);

    if (cloud.pts.size()>0){
      WireCell::WCPointCloud<double>::WCPoint wcp = point_cloud_steiner->get_closest_wcpoint(p);
      double temp_dis = sqrt(pow(wcp.x-p.x,2)+pow(wcp.y-p.y,2)+pow(wcp.z-p.z,2));
      
      //      std::cout << temp_dis/units::cm << " " << dis_cut/units::cm << std::endl;
      if (temp_dis < dis_cut){
	std::set<int> total_vertices_found;
	std::set<int> vertices_to_be_examined;
	std::set<int> vertices_saved_for_next;
	total_vertices_found.insert(wcp.index);
	vertices_to_be_examined.insert(wcp.index);
	
	for (int j=0;j!=nlevel;j++){
	  for (auto it = vertices_to_be_examined.begin(); it!=vertices_to_be_examined.end(); it++){
	    int temp_current_index = (*it);
	    std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(temp_current_index,*graph_steiner),*graph_steiner);
	    for (; neighbors.first!=neighbors.second; ++neighbors.first){
	      //std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
	      if (total_vertices_found.find(index(*neighbors.first))==total_vertices_found.end()){
		total_vertices_found.insert(index(*neighbors.first));
		vertices_saved_for_next.insert(index(*neighbors.first));
	      }
	    }
	  }
	  vertices_to_be_examined = vertices_saved_for_next;
	}
	
	WireCell::Point temp_p(wcp.x, wcp.y, wcp.z);
	std::vector<int> temp_results = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
	// std::cout << cloud.pts[wcp.index].index_u << " "  << temp_results.at(0) << " " << temp_results.at(1) << std::endl;
	//	if (wcp.mcell!=0)
	// std::cout << temp_results.at(0) << " " << wcp.mcell->GetTimeSlice() << std::endl;
	int cur_time_slice = temp_results.at(0);
	int cur_wire_u = cloud.pts[wcp.index].index_u;
	int cur_wire_v = cloud.pts[wcp.index].index_v;
	int cur_wire_w = cloud.pts[wcp.index].index_w;

	//	std::cout << "B: " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_w << std::endl;
	
	SMGCSet nearby_mcells_set;
	std::vector<int> point_indices;
	std::vector<int> point_timeslices;
	
	for (auto it = total_vertices_found.begin(); it!=total_vertices_found.end(); it++){
	  SlimMergeGeomCell *mcell = cloud.pts[*it].mcell;
	  if (mcell!=0){
	    nearby_mcells_set.insert(mcell);
	  }else{
	    temp_p.x = cloud.pts[*it].x;
	    temp_p.y = cloud.pts[*it].y;
	    temp_p.z = cloud.pts[*it].z;
	    temp_results = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
	    
	    point_indices.push_back(*it);
	    point_timeslices.push_back(temp_results.at(0));
	  }
	}
	
	
	
	double dis_cut_u = dis_cut;
	double dis_cut_v = dis_cut;
	double dis_cut_w = dis_cut;
	
	double max_time_slice_u = 0;
	double max_time_slice_v = 0;
	double max_time_slice_w = 0;
	
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_u)<1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
		max_time_slice_u = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_v)<1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
		max_time_slice_v = fabs(this_time_slice-cur_time_slice);
	  }
	  for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	    if (fabs((*it1)->index()-cur_wire_w)<1)
	      if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
		max_time_slice_w = fabs(this_time_slice-cur_time_slice);
	  }
	}
	for (size_t i=0; i!=point_indices.size(); i++){
	  int this_time_slice = point_timeslices.at(i);
	  if (fabs(cloud.pts[point_indices.at(i)].index_u - cur_wire_u)<=1){
	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
	      max_time_slice_u = fabs(this_time_slice-cur_time_slice);
	  }
	  if (fabs(cloud.pts[point_indices.at(i)].index_v - cur_wire_v)<=1){
	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
	      max_time_slice_v = fabs(this_time_slice-cur_time_slice);
	  }
	  if (fabs(cloud.pts[point_indices.at(i)].index_w - cur_wire_w)<=1){
	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
	      max_time_slice_w = fabs(this_time_slice-cur_time_slice);
	  }
	}
	
	if (max_time_slice_u * time_slice_width*1.2 < dis_cut_u)
	  dis_cut_u = max_time_slice_u * time_slice_width*1.2;
	if (max_time_slice_v * time_slice_width*1.2 < dis_cut_v)
	  dis_cut_v = max_time_slice_v * time_slice_width*1.2;
	if (max_time_slice_w * time_slice_width*1.2 < dis_cut_w)
	  dis_cut_w = max_time_slice_w * time_slice_width*1.2;
	
	// actual cut ...
	for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
	  SlimMergeGeomCell *mcell = *it;
	  int this_time_slice = mcell->GetTimeSlice();
	  
	  double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    float min_u_dis;
	    if (cur_wire_u < uwires.front()->index()){
	      min_u_dis = uwires.front()->index()-cur_wire_u;
	    }else if (cur_wire_u >= uwires.front()->index() &&
		      cur_wire_u <= uwires.back()->index()){
	      min_u_dis = 0;
	    }else{
	      min_u_dis = cur_wire_u-uwires.back()->index();
	    }
	    float min_v_dis;
	    if (cur_wire_v < vwires.front()->index()){
	      min_v_dis = vwires.front()->index()-cur_wire_v;
	    }else if (cur_wire_v >= vwires.front()->index() &&
		      cur_wire_v <= vwires.back()->index()){
	      min_v_dis = 0;
	    }else{
	      min_v_dis = cur_wire_v-vwires.back()->index();
	    }
	    float min_w_dis;
	    if (cur_wire_w < wwires.front()->index()){
	      min_w_dis = wwires.front()->index()-cur_wire_w;
	    }else if (cur_wire_w >= wwires.front()->index() &&
		      cur_wire_w <= wwires.back()->index()){
	      min_w_dis = 0;
	    }else{
	      min_w_dis = cur_wire_w-wwires.back()->index();
	    }
	    float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	    
	    if ( range_u > 0 && range_v >0 && range_w > 0){
	      float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
	      float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
	      float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
	      float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
	      float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
	      float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;
	      for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
		temp_2dut.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
		temp_2dvt.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
		temp_2dwt.insert(std::make_pair(j,this_time_slice));
	      }
	    } // cuts
	  } // cuts
	} // loop over all mcells
	
	for (size_t i=0; i!=point_indices.size(); i++){

	  int this_time_slice = point_timeslices.at(i);
	  int this_index_u = cloud.pts[point_indices.at(i)].index_u;
	  int this_index_v = cloud.pts[point_indices.at(i)].index_v;
	  int this_index_w = cloud.pts[point_indices.at(i)].index_w;
	  
	  double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
	  if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
	    float min_u_dis;
	    if (cur_wire_u < this_index_u-1){
	      min_u_dis = this_index_u -1 - cur_wire_u;
	    }else if (cur_wire_u > this_index_u + 1){
	      min_u_dis = cur_wire_u - this_index_u - 1;
	    }else{
	      min_u_dis = 0;
	    }
	    float min_v_dis;
	    if (cur_wire_v < this_index_v-1){
	      min_v_dis = this_index_v -1 - cur_wire_v;
	    }else if (cur_wire_v > this_index_v + 1){
	      min_v_dis = cur_wire_v - this_index_v -1 ;
	    }else{
	      min_v_dis = 0;
	    }
	    float min_w_dis;
	    if (cur_wire_w < this_index_w-1){
	      min_w_dis = this_index_w  -1 - cur_wire_w;
	    }else if (cur_wire_w > this_index_w+1){
	      min_w_dis = cur_wire_w - this_index_w -1;
	    }else{
	      min_w_dis = 0;
	    }
	    float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
	    float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	    
	    if ( range_u > 0 && range_v >0 && range_w > 0){
	      float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
	      float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
	      float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
	      float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
	      float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
	      float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;
	      for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
		temp_2dut.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
		temp_2dvt.insert(std::make_pair(j,this_time_slice));
	      }
	      for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
		temp_2dwt.insert(std::make_pair(j,this_time_slice));
	      }
	    } // cuts
	  } // cuts
	  
	} // loop through points
	
      } // distance
    } // no steiner tree cloud
  }

  // just projection ...
  if (temp_2dut.size()==0 && temp_2dvt.size()==0 && temp_2dwt.size()==0){
    //std::cout << "haha " << std::endl;
    std::vector<int> temp_results = ct_point_cloud.convert_3Dpoint_time_ch(p);
    int cur_time_slice = temp_results.at(0);
    int cur_index_u = temp_results.at(1);
    int cur_index_v = temp_results.at(2);
    int cur_index_w = temp_results.at(3);

    for (int i=-time_cut; i!=time_cut+1;i++){
      for (int j=-time_cut; j!=time_cut+1; j++){
	if (abs(i)+abs(j) <= time_cut){
	  temp_2dut.insert(std::make_pair(cur_index_u+i, cur_time_slice+j));
	  temp_2dvt.insert(std::make_pair(cur_index_v+i, cur_time_slice+j));
	  temp_2dwt.insert(std::make_pair(cur_index_w+i, cur_time_slice+j));
	}
      }
    }
  }
  
  
}


void WireCellPID::PR3DCluster::form_map(WireCell::ToyCTPointCloud& ct_point_cloud, WireCell::PointVector& pts,
		  std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge,
		  std::map<int, std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set,
		  std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,
	       double end_point_factor, double mid_point_factor, int nlevel, double time_cut, double charge_cut){
  
  map_3D_2DU_set.clear();
  map_3D_2DV_set.clear();
  map_3D_2DW_set.clear();

  map_2DU_3D_set.clear();
  map_2DV_3D_set.clear();
  map_2DW_3D_set.clear();

 
  WireCell::PointVector saved_pts;
  int count = 0;
  
  // distance ...
  std::vector<double> distances;
  for (size_t i=0;i+1!=pts.size();i++){
    distances.push_back(sqrt(pow(pts.at(i+1).x-pts.at(i).x,2) +
			     pow(pts.at(i+1).y-pts.at(i).y,2) +
			     pow(pts.at(i+1).z-pts.at(i).z,2)));
  }

  

  // start to loop over the path ...
  for (size_t i=0;i!=pts.size();i++){
    double dis_cut;
    if (i==0){
      dis_cut = std::min(distances.at(i) * end_point_factor,4/3.*end_point_factor*units::cm);
    }else if (i+1==pts.size()){
      dis_cut = std::min(distances.back() * end_point_factor,4/3.*end_point_factor*units::cm);
    }else{
      dis_cut = std::min(std::max(distances.at(i-1)*mid_point_factor,distances.at(i)*mid_point_factor),4/3.*mid_point_factor*units::cm);
    }

    std::set<std::pair<int,int> > temp_2dut, temp_2dvt, temp_2dwt;
    form_point_association(pts.at(i), temp_2dut, temp_2dvt, temp_2dwt, ct_point_cloud, dis_cut, nlevel, time_cut);
    // examine ...
    std::vector<int> temp_flag = examine_point_association(temp_2dut, temp_2dvt, temp_2dwt, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,charge_cut);

    //    std::cout << temp_2dut.size() << " " << temp_2dvt.size() << " " << temp_2dwt.size() << " " << temp_flag.at(0) << " " << temp_flag.at(1) << " " << temp_flag.at(2) << std::endl;
    // just projection ...
    // fill the data ...
    if (temp_2dut.size() + temp_2dvt.size() + temp_2dwt.size() > 0){
      map_3D_2DU_set[count] = std::make_pair(temp_2dut,temp_flag.at(0));
      map_3D_2DV_set[count] = std::make_pair(temp_2dvt,temp_flag.at(1));
      map_3D_2DW_set[count] = std::make_pair(temp_2dwt,temp_flag.at(2));
      for (auto it = temp_2dut.begin(); it!=temp_2dut.end();it++){
	if (map_2DU_3D_set.find(*it)==map_2DU_3D_set.end()){
	  std::set<int>  temp_set;
	  temp_set.insert(count);
	  map_2DU_3D_set[*it] = temp_set;
	}else{
	  map_2DU_3D_set[*it].insert(count);
	}
      }
      for (auto it = temp_2dvt.begin(); it!=temp_2dvt.end();it++){
	if (map_2DV_3D_set.find(*it)==map_2DV_3D_set.end()){
	  std::set<int>  temp_set;
	  temp_set.insert(count);
	  map_2DV_3D_set[*it] = temp_set;
	}else{
	  map_2DV_3D_set[*it].insert(count);
	}
      }
      for (auto it = temp_2dwt.begin(); it!=temp_2dwt.end();it++){
	if (map_2DW_3D_set.find(*it)==map_2DW_3D_set.end()){
	  std::set<int>  temp_set;
	  temp_set.insert(count);
	  map_2DW_3D_set[*it] = temp_set;
	}else{
	  map_2DW_3D_set[*it].insert(count);
	}
      }

      
      saved_pts.push_back(pts.at(i));
      count ++;
    }
  }
  
  pts = saved_pts;
 }

void WireCellPID::PR3DCluster::prepare_data(WireCell::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge){
  
  std::vector<int> proj_channel;
  std::vector<int> proj_timeslice;
  std::vector<int> proj_charge;
  std::vector<int> proj_charge_err;
  std::vector<int> proj_flag;
  get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);

  int min_time = 1e9;
  int max_time = -1e9;
  int min_uch = 1e9;
  int max_uch = -1e9;
  int min_vch = 1e9;
  int max_vch = -1e9;
  int min_wch = 1e9;
  int max_wch = -1e9;
  
  // std::cout << proj_charge.size() << " " << proj_flag.size() << std::endl;
  for (size_t i=0;i!=proj_channel.size();i++){
    if (min_time > proj_timeslice.at(i)) min_time = proj_timeslice.at(i);
    if (max_time < proj_timeslice.at(i)) max_time = proj_timeslice.at(i);
	
    if (proj_channel.at(i)<2400){
      map_2D_ut_charge[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));

      if (min_uch > proj_channel.at(i)) min_uch = proj_channel.at(i);
      if (max_uch < proj_channel.at(i)) max_uch = proj_channel.at(i);
    }else if (proj_channel.at(i)<4800){
      map_2D_vt_charge[std::make_pair(proj_channel.at(i)-2400,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));

      if (min_vch > proj_channel.at(i)) min_vch = proj_channel.at(i);
      if (max_vch < proj_channel.at(i)) max_vch = proj_channel.at(i);
    }else{
      map_2D_wt_charge[std::make_pair(proj_channel.at(i)-4800,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));
      if (min_wch > proj_channel.at(i)) min_wch = proj_channel.at(i);
      if (max_wch < proj_channel.at(i)) max_wch = proj_channel.at(i);
    }
  }

  // flag 0 for dead or overlapping channels
  // flag 1 for good channels
  // flag 2 for isolated channels
  // flag 3 for additional channels ...
  
  // add the rest of live channels within range???
  std::map<std::pair<int,int>, std::pair<double,double> > map_u_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_uch, max_uch, 0);
  std::map<std::pair<int,int>, std::pair<double,double> > map_v_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_vch, max_vch, 1);
  std::map<std::pair<int,int>, std::pair<double,double> > map_w_tcc = ct_point_cloud.get_overlap_good_ch_charge(min_time, max_time, min_wch, max_wch, 2);

  for (auto it = map_u_tcc.begin(); it!=map_u_tcc.end(); it++){
    if (map_2D_ut_charge.find(std::make_pair(it->first.second, it->first.first))==map_2D_ut_charge.end()){
      map_2D_ut_charge[std::make_pair(it->first.second, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }

  for (auto it = map_v_tcc.begin(); it!=map_v_tcc.end(); it++){
    if (map_2D_vt_charge.find(std::make_pair(it->first.second-2400, it->first.first))==map_2D_vt_charge.end()){
      map_2D_vt_charge[std::make_pair(it->first.second-2400, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }

  for (auto it = map_w_tcc.begin(); it!=map_w_tcc.end(); it++){
    if (map_2D_wt_charge.find(std::make_pair(it->first.second-4800, it->first.first))==map_2D_wt_charge.end()){
      map_2D_wt_charge[std::make_pair(it->first.second, it->first.first)] = std::make_tuple(it->second.first, it->second.second, 3);
      //      std::cout << it->first.first << " " << it->first.second << std::endl;
    }else{
      //      std::cout << it->first.first << " " << it->first.second << std::endl;
    }
  }
  
  
  //  for (auto it = map_w_tcc.begin(); it!=map_w_tcc.end(); it++){
  // std::cout << it->first.first << " " << it->first.second << std::endl;
  // }
  
}

void WireCellPID::PR3DCluster::organize_ps_path(WireCell::PointVector& pts, double low_dis_limit, double end_point_limit){
  WireCell::PointVector ps_vec = pts;
  pts.clear();
  // fill in the beginning part
  {
    Point p1 = ps_vec.front();
    Point p2 = ps_vec.front();
    double dis1 = 0;
    for (auto it = ps_vec.begin(); it!=ps_vec.end(); it++){
      p2 = *it;
      dis1 = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis1 > low_dis_limit) break;
    }
    if (dis1!=0){
      p1.x += (p1.x - p2.x)/dis1 * end_point_limit;
      p1.y += (p1.y - p2.y)/dis1 * end_point_limit;
      p1.z += (p1.z - p2.z)/dis1 * end_point_limit;
      pts.push_back(p1);
    }
  }
  
  // fill in the middle part
  for (size_t i=0;i!=ps_vec.size(); i++){
    Point p1 = ps_vec.at(i);
    double dis = sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2));
    if (dis < low_dis_limit * 0.8 ){
      continue;
    }else if (dis < low_dis_limit * 1.6){
      pts.push_back(p1);
    }else{
      int npoints = std::round(dis/low_dis_limit);
      
      for (int j=0;j!=npoints;j++){
	Point p(pts.back().x + (p1.x-pts.back().x) / npoints * (j+1),
		pts.back().y + (p1.y-pts.back().y) / npoints * (j+1),
		pts.back().z + (p1.z-pts.back().z) / npoints * (j+1));
	pts.push_back(p);
      }
    }
  }
  

  // fill in the end part
  if (end_point_limit!=0){
    Point p1 = ps_vec.back();
    Point p2 = ps_vec.back();
    double dis1 = 0;
    for (auto it = ps_vec.rbegin(); it!=ps_vec.rend(); it++){
      p2 = *it;
      dis1 = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis1 > low_dis_limit) break;
    }
    if (dis1!=0){
      p1.x += (p1.x - p2.x)/dis1 * end_point_limit;
      p1.y += (p1.y - p2.y)/dis1 * end_point_limit;
      p1.z += (p1.z - p2.z)/dis1 * end_point_limit;
      pts.push_back(p1);
    }
  }
  
}

WireCell::PointVector WireCellPID::PR3DCluster::organize_wcps_path(std::list<WCPointCloud<double>::WCPoint>& path_wcps_list,  double low_dis_limit, double end_point_limit){

  PointVector pts;
  
  std::vector<WCPointCloud<double>::WCPoint> temp_wcps_vec(path_wcps_list.begin(), path_wcps_list.end());

  // fill in the beginning point ...
  {
    Point p1(temp_wcps_vec.front().x, temp_wcps_vec.front().y, temp_wcps_vec.front().z);
    Point p2(temp_wcps_vec.front().x, temp_wcps_vec.front().y, temp_wcps_vec.front().z);
    double dis1 = 0;
    for (auto it = temp_wcps_vec.begin(); it!=temp_wcps_vec.end(); it++){
      p2.x = (*it).x;
      p2.y = (*it).y;
      p2.z = (*it).z;
      dis1 = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis1 > low_dis_limit) break;
    }
    if (dis1!=0){
      p1.x += (p1.x - p2.x)/dis1 * end_point_limit;
      p1.y += (p1.y - p2.y)/dis1 * end_point_limit;
      p1.z += (p1.z - p2.z)/dis1 * end_point_limit;
      pts.push_back(p1);
    }
  }

  // fill in the middle part
  for (size_t i=0;i!=temp_wcps_vec.size(); i++){
    Point p1(temp_wcps_vec.at(i).x, temp_wcps_vec.at(i).y, temp_wcps_vec.at(i).z);

    double dis = sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2));
    
    if (dis < low_dis_limit * 0.8 ){
      continue;
    }else if (dis < low_dis_limit * 1.6){
      pts.push_back(p1);
    }else{
      int npoints = std::round(dis/low_dis_limit);
      
      for (int j=0;j!=npoints;j++){
	Point p(pts.back().x + (p1.x-pts.back().x) / npoints * (j+1),
		pts.back().y + (p1.y-pts.back().y) / npoints * (j+1),
		pts.back().z + (p1.z-pts.back().z) / npoints * (j+1));
	pts.push_back(p);
      }
    }
  }
  

  // fill in the end part
  {
    Point p1(temp_wcps_vec.back().x, temp_wcps_vec.back().y, temp_wcps_vec.back().z);
    Point p2(temp_wcps_vec.back().x, temp_wcps_vec.back().y, temp_wcps_vec.back().z);
    double dis1 = 0;
    for (auto it = temp_wcps_vec.rbegin(); it!=temp_wcps_vec.rend(); it++){
      p2.x = (*it).x;
      p2.y = (*it).y;
      p2.z = (*it).z;
      dis1 = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis1 > low_dis_limit) break;
    }
    if (dis1!=0){
      p1.x += (p1.x - p2.x)/dis1 * end_point_limit;
      p1.y += (p1.y - p2.y)/dis1 * end_point_limit;
      p1.z += (p1.z - p2.z)/dis1 * end_point_limit;
      pts.push_back(p1);
    }
  }

    
  
  
  return pts;
}
