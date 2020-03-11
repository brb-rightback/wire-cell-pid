
void WCPPID::PR3DCluster::do_multi_tracking(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, double time, bool flag_dQ_dx_fit_reg){
  // collect charge
  collect_charge_multi_trajectory(map_segment_vertices, ct_point_cloud);

  // clear things
  fine_tracking_path.clear();
  dQ.clear();
  dx.clear();
  pu.clear();
  pv.clear();
  pw.clear();
  pt.clear();
  reduced_chi2.clear();

  bool flag_1st_tracking = true;
  bool flag_2nd_tracking = true;
  bool flag_dQ_dx = true;

  // prepare the data for the fit, do not contain everything ...
  // form from the cluster ...
  std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_ut_charge;
  std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_vt_charge;
  std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_wt_charge;
  prepare_data(ct_point_cloud, global_wc_map, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

  // first round of organizing the path from the path_wcps (shortest path)
  double low_dis_limit = 1.2*units::cm;
  double end_point_limit = 0.6*units::cm;
  organize_segments_path(map_vertex_segments, map_segment_vertices, low_dis_limit, end_point_limit);
  
  // form association ...
  // map 3D index to set of 2D points
  std::map<int,std::pair<std::set<std::pair<int,int>>, float> > map_3D_2DU_set;
  std::map<int,std::pair<std::set<std::pair<int,int>>, float> > map_3D_2DV_set;
  std::map<int,std::pair<std::set<std::pair<int,int>>, float> > map_3D_2DW_set;
  std::map<int,std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, int> > map_3D_tuple;
  // map 2D points to 3D indices
  std::map<std::pair<int,int>,std::set<int>> map_2DU_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DV_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DW_3D_set;
  
  if (flag_1st_tracking){
    form_map_multi_segments(map_vertex_segments, map_segment_vertices, ct_point_cloud,
			    map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,
			    map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set, map_3D_tuple,
			    map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set);
    multi_trajectory_fit(map_vertex_segments, map_segment_vertices,
			 map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set, map_3D_tuple,
			 map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set,
			 map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);
  }
  
}

WCP::Point WCPPID::PR3DCluster::get_pos_multi(std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, int>& input){
  WCPPID::ProtoVertex* vtx = std::get<0>(input);
  WCPPID::ProtoSegment* seg = std::get<1>(input);
  int num = std::get<2>(input);
  Point p;
  if (vtx!=0){
    p = vtx->get_fit_pt();
  }else{
    p = seg->get_point_vec().at(num);
  }
  return p;
}

void WCPPID::PR3DCluster::multi_trajectory_fit(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<int, std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, int> >& map_3D_tuple, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,std::map<std::pair<int,int>, std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, int charge_div_method, double div_sigma){
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

  /* for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){ */
  /*   std::cout  << " " << it->mcell << std::endl; */
  /* } */
  
  double slope_x = 1./time_slice_width;
  //mcell
  double first_t_dis = point_cloud->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - point_cloud->get_cloud().pts[0].x;
    //  double first_t_dis = path_wcps.front().mcell->GetTimeSlice()*time_slice_width - path_wcps.front().x;
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

  // 2D pixel and then the 3D index ...
  std::map<std::tuple<int,int, int>, double> map_Udiv_fac;
  std::map<std::tuple<int,int, int>, double> map_Vdiv_fac;
  std::map<std::tuple<int,int, int>, double> map_Wdiv_fac;

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
	Point tmp_p = get_pos_multi(map_3D_tuple[*it1]);
  	double central_t = slope_x * tmp_p.x + offset_t;
  	double central_ch = slope_yu * tmp_p.y + slope_zu * tmp_p.z + offset_u;
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
	Point tmp_p = get_pos_multi(map_3D_tuple[*it1]);
  	double central_t = slope_x * tmp_p.x + offset_t;
  	double central_ch = slope_yv * tmp_p.y + slope_zv * tmp_p.z + offset_v;
  	double factor = exp(-0.5 * (pow((central_t-it->first.second)*time_slice_width,2) + pow((central_ch-it->first.first)*pitch_v,2))/pow(div_sigma,2));
  	map_Vdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] = factor;
  	sum += factor;
	//	if ( it->first.first == 2520 &&  it->first.second ==1)
	//  std::cout << it->first.first << " " << it->first.second << " " << central_t << " " << central_ch << " " << factor << std::endl;
      }
     
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
  	map_Vdiv_fac[std::make_tuple(it->first.first, it->first.second, *it1)] /= sum;
      }
    }

    
    for (auto it = map_2DW_3D_set.begin(); it!=map_2DW_3D_set.end(); it++){
      double sum = 0;
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
	Point tmp_p = get_pos_multi(map_3D_tuple[*it1]);
  	double central_t = slope_x * tmp_p.x + offset_t;
  	double central_ch = slope_yw * tmp_p.y + slope_zw * tmp_p.z + offset_w;
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

  // deal vertex
  for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    int i = vtx->get_fit_index();
    bool flag_fit_fix = vtx->get_flag_fit_fix();
    Point init_p = vtx->get_fit_pt();
    if (!flag_fit_fix){ // not fix the fit ...
      init_p = fit_point(init_p, i, map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
			 map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,
			 map_Udiv_fac, map_Vdiv_fac, map_Wdiv_fac, offset_t, slope_x,
			 offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv,
			 offset_w, slope_yw, slope_zw);
    } 
    vtx->set_fit(init_p, 0, -1, offset_u + 0.5 + (slope_yu * init_p.y + slope_zu * init_p.z), offset_v + 0.5 + (slope_yv * init_p.y + slope_zv * init_p.z)+2400, offset_w + 0.5 + (slope_yw * init_p.y + slope_zw * init_p.z)+4800, offset_t + 0.5 + slope_x * init_p.x, -1);
  }

  // deal tracks 
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
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
    PointVector init_ps = sg->get_point_vec();
    std::vector<int> init_indices = sg->get_fit_index_vec();
    std::vector<bool> init_fit_skip = sg->get_fit_flag_skip();
    //    std::cout << init_ps.size() << " " << init_indices.size() << " " << init_fit_skip.size() << " " << init_fit_skip.front() << " " << init_fit_skip.back() << std::endl;

    PointVector final_ps;
    
    for (size_t i=0;i!=init_ps.size();i++){
      if (i==0){
	final_ps.push_back(start_v->get_fit_pt());
      }else if (i+1==init_ps.size()){
	final_ps.push_back(end_v->get_fit_pt());
      }else{
	Point temp_p = init_ps.at(i);
	if (!init_fit_skip.at(i)){
	  temp_p = fit_point(init_ps.at(i), init_indices.at(i), map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
			     map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,
			     map_Udiv_fac, map_Vdiv_fac, map_Wdiv_fac, offset_t, slope_x,
			     offset_u, slope_yu, slope_zu, offset_v, slope_yv, slope_zv,
			     offset_w, slope_yw, slope_zw);
	}
	//	std::cout << i << " " << init_ps.at(i) << " " << temp_p << std::endl;
	final_ps.push_back(temp_p);
      }
    }
    // sort out the results ...
    
    
    // set results
    std::vector<double> tmp_dQ_vec(final_ps.size(),0);
    std::vector<double> tmp_dx_vec(final_ps.size(),-1);
    std::vector<double> tmp_pu_vec(final_ps.size(),0);
    std::vector<double> tmp_pv_vec(final_ps.size(),0);
    std::vector<double> tmp_pw_vec(final_ps.size(),0);
    std::vector<double> tmp_pt_vec(final_ps.size(),0);
    std::vector<double> tmp_reduced_chi2_vec(final_ps.size(),0);
    for (size_t i=0;i!=final_ps.size();i++){
      tmp_pu_vec.at(i) = offset_u + 0.5 + (slope_yu * final_ps.at(i).y + slope_zu * final_ps.at(i).z);
      tmp_pv_vec.at(i) = offset_v + 0.5 + (slope_yv * final_ps.at(i).y + slope_zv * final_ps.at(i).z)+2400;
      tmp_pw_vec.at(i) = offset_w + 0.5 + (slope_yw * final_ps.at(i).y + slope_zw * final_ps.at(i).z)+4800;
      tmp_pt_vec.at(i) = offset_t + 0.5 + slope_x * final_ps.at(i).x ;
    }
    sg->set_fit_vec(final_ps, tmp_dQ_vec, tmp_dx_vec, tmp_pu_vec, tmp_pv_vec, tmp_pw_vec, tmp_pt_vec, tmp_reduced_chi2_vec);
    
  }
  
}

Point WCPPID::PR3DCluster::fit_point(Point& init_p, int i, std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<std::pair<int,int>, std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, std::map<std::tuple<int,int, int>, double>& map_Udiv_fac, std::map<std::tuple<int,int, int>, double>& map_Vdiv_fac, std::map<std::tuple<int,int, int>, double>& map_Wdiv_fac, double offset_t, double slope_x, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw){
  int n_2D_u = 2* map_3D_2DU_set[i].first.size();
  int n_2D_v = 2* map_3D_2DV_set[i].first.size();
  int n_2D_w = 2* map_3D_2DW_set[i].first.size();
  Eigen::VectorXd temp_pos_3D(3),data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
  Eigen::VectorXd temp_pos_3D_init(3);
  Eigen::SparseMatrix<double> RU(n_2D_u, 3) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, 3) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, 3) ;
  // initialization
  temp_pos_3D_init(0) = init_p.x;
  temp_pos_3D_init(1) = init_p.y;
  temp_pos_3D_init(2) = init_p.z;

  for (size_t i=0;i!=n_2D_u;i++){
    data_u_2D(i) = 0;
  }
  for (size_t i=0;i!=n_2D_v;i++){
    data_v_2D(i) = 0;
  }
  for (size_t i=0;i!=n_2D_w;i++){
    data_w_2D(i) = 0;
  }
  
  int index = 0;
  for (auto it = map_3D_2DU_set[i].first.begin(); it!=map_3D_2DU_set[i].first.end(); it++){
    double charge, charge_err;
    if (map_2D_ut_charge.find(*it)!=map_2D_ut_charge.end()){
      charge = std::get<0>(map_2D_ut_charge[*it]);
      charge_err = std::get<1>(map_2D_ut_charge[*it]);
    }else{
      charge = 100; // serve as regularization ...
      charge_err = 1000;
    }
    if (charge<100){
      charge = 100;
      charge_err = 1000;
    }
    // induction plane error is large, so they have less weights in the fit to decide X position
    double scaling = charge/charge_err * map_Udiv_fac[std::make_tuple(it->first, it->second, i)];
    
    if (map_3D_2DU_set[i].second < 0.5){
      if (map_3D_2DU_set[i].second!=0)
	scaling *= pow(map_3D_2DU_set[i].second/0.5,1);
      else
	scaling *= 0.05;
    }
    
    if (scaling!=0){
      data_u_2D(2*index) =  scaling * (it->first - offset_u);
      data_u_2D(2*index+1) = scaling * (it->second - offset_t);
      
      RU.insert(2*index, 1) = scaling * slope_yu; // Y--> U
      RU.insert(2*index, 2) = scaling * slope_zu; // Z--> U
      RU.insert(2*index+1,0) = scaling * slope_x; // X --> T
    }
    index ++;
  }
  index = 0;
  for (auto it = map_3D_2DV_set[i].first.begin(); it!=map_3D_2DV_set[i].first.end(); it++){
    double charge, charge_err;
    if (map_2D_vt_charge.find(*it)!=map_2D_vt_charge.end()){
      charge = std::get<0>(map_2D_vt_charge[*it]);
      charge_err = std::get<1>(map_2D_vt_charge[*it]);
    }else{
      charge = 100; // serve as regularization ...
      charge_err = 1000;
    }
    if (charge < 100){
      charge = 100;
      charge_err = 1000;
    }
    // induction plane error is large, so they have less weights in the fit to decide X position
    double scaling = charge/charge_err * map_Vdiv_fac[std::make_tuple(it->first, it->second, i)];
    if (map_3D_2DV_set[i].second < 0.5){
      if (map_3D_2DV_set[i].second!=0)
	scaling *= pow(map_3D_2DV_set[i].second/0.5,1);
      else
	scaling *= 0.05;
    }
    
    if (scaling!=0){
      data_v_2D(2*index) =  scaling * (it->first - offset_v);
      data_v_2D(2*index+1) = scaling * (it->second - offset_t);
      
      RV.insert(2*index, 1) = scaling * slope_yv; // Y--> V
      RV.insert(2*index, 2) = scaling * slope_zv; // Z--> V
      RV.insert(2*index+1,0) = scaling * slope_x; // X --> T
    }
    index ++;
  }
  
  index = 0;
  for (auto it = map_3D_2DW_set[i].first.begin(); it!=map_3D_2DW_set[i].first.end(); it++){
    double charge, charge_err;
    if (map_2D_wt_charge.find(*it)!=map_2D_wt_charge.end()){
      charge = std::get<0>(map_2D_wt_charge[*it]);
      charge_err = std::get<1>(map_2D_wt_charge[*it]);
    }else{
      charge = 100; // serve as regularization ...
      charge_err = 1000;
    }
    if (charge<100){
      charge = 100;
      charge_err = 1000;
    }
    // induction plane error is large, so they have less weights in the fit to decide X position
    double scaling = charge/charge_err * map_Wdiv_fac[std::make_tuple(it->first, it->second, i)];
    
    if (map_3D_2DW_set[i].second < 0.5){
      if (map_3D_2DW_set[i].second!=0)
	scaling *= pow(map_3D_2DW_set[i].second/0.5,1);
      else
	scaling *= 0.05;
    }
    if (scaling!=0){
      data_w_2D(2*index) =  scaling * (it->first - offset_w);
      data_w_2D(2*index+1) = scaling * (it->second - offset_t);
      
      RW.insert(2*index, 2) = scaling * slope_zw; // Z--> W
      RW.insert(2*index+1,0) = scaling * slope_x; // X --> T
    }
    index ++;
  }
  
  Eigen::SparseMatrix<double> RUT = Eigen::SparseMatrix<double>(RU.transpose());
  Eigen::SparseMatrix<double> RVT = Eigen::SparseMatrix<double>(RV.transpose());
  Eigen::SparseMatrix<double> RWT = Eigen::SparseMatrix<double>(RW.transpose());
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;// + PMatrixT * pos_3D_init ;
  
  Eigen::SparseMatrix<double> A =   RUT * RU + RVT * RV + RWT * RW ;// + PMatrixT*PMatrix;
  solver.compute(A);
  
  temp_pos_3D = solver.solveWithGuess(b,temp_pos_3D_init);
  
  Point final_p;
  if (std::isnan(solver.error())){
    final_p.x = temp_pos_3D_init(0);
    final_p.y = temp_pos_3D_init(1);
    final_p.z = temp_pos_3D_init(2);
  }else{
    final_p.x = temp_pos_3D(0);
    final_p.y = temp_pos_3D(1);
    final_p.z = temp_pos_3D(2);
  }
  
  //std::cout << init_p << " " << final_p << std::endl;
  
  return final_p;
}


void WCPPID::PR3DCluster::form_map_multi_segments(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud,
						  std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge,
						  std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<int, std::tuple<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*, int> >& map_3D_tuple,
						  std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,
						  double end_point_factor, double mid_point_factor, int nlevel, double time_cut, double charge_cut){
  map_3D_2DU_set.clear();
  map_3D_2DV_set.clear();
  map_3D_2DW_set.clear();
  map_3D_tuple.clear();

  map_2DU_3D_set.clear();
  map_2DV_3D_set.clear();
  map_2DW_3D_set.clear();

  int count = 0;
  
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
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
    PointVector saved_pts;
    std::vector<int> saved_index;
    std::vector<bool> saved_skip;

    std::vector<double> distances;
    for (size_t i=0;i+1!=pts.size();i++){
      distances.push_back(sqrt(pow(pts.at(i+1).x-pts.at(i).x,2) +
			       pow(pts.at(i+1).y-pts.at(i).y,2) +
			       pow(pts.at(i+1).z-pts.at(i).z,2)));
    }
    
    for (size_t i=0;i!=pts.size();i++){
      double dis_cut;
      if (i==0){
	dis_cut = std::min(distances.at(i) * end_point_factor,4/3.*end_point_factor*units::cm);
	if (start_v->get_fit_range()<0){
	  start_v->set_fit_range(dis_cut);
	}else{
	  if (dis_cut < start_v->get_fit_range()) start_v->set_fit_range(dis_cut);
	}
      }else if (i+1==pts.size()){
	dis_cut = std::min(distances.back() * end_point_factor,4/3.*end_point_factor*units::cm);
	if (end_v->get_fit_range()<0){
	  end_v->set_fit_range(dis_cut);
	}else{
	  if (dis_cut < end_v->get_fit_range()) end_v->set_fit_range(dis_cut);
	}
      }else{
	dis_cut = std::min(std::max(distances.at(i-1)*mid_point_factor,distances.at(i)*mid_point_factor),4/3.*mid_point_factor*units::cm);
      }

      // not the first and last point
      if (i!=0 && i+1!=pts.size()){
	std::set<std::pair<int,int> > temp_2dut, temp_2dvt, temp_2dwt;
	form_point_association(pts.at(i), temp_2dut, temp_2dvt, temp_2dwt, ct_point_cloud, dis_cut, nlevel, time_cut);
	// examine ...
	std::vector<int> temp_results = ct_point_cloud.convert_3Dpoint_time_ch(pts.at(i));
	temp_results.at(2)-=2400;
	temp_results.at(3)-=4800;
	std::vector<float> temp_flag;
	if (i==0 || i==1 || i+1 ==pts.size() || i+2 == pts.size()){
	  temp_flag = examine_point_association(temp_results, temp_2dut, temp_2dvt, temp_2dwt, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,true,charge_cut);
	}else{
	  temp_flag = examine_point_association(temp_results, temp_2dut, temp_2dvt, temp_2dwt, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,false,charge_cut);
	}
	if (temp_flag.at(0) + temp_flag.at(1) + temp_flag.at(2) > 0){
	  map_3D_2DU_set[count] = std::make_pair(temp_2dut,temp_flag.at(0));
	  map_3D_2DV_set[count] = std::make_pair(temp_2dvt,temp_flag.at(1));
	  map_3D_2DW_set[count] = std::make_pair(temp_2dwt,temp_flag.at(2));
	  map_3D_tuple[count] = std::make_tuple((WCPPID::ProtoVertex*)0,sg,i);
	  
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
	  saved_index.push_back(count);
	  saved_skip.push_back(false);
	  count ++;
	}
      }else if (i==0){ // first point
	saved_pts.push_back(pts.at(i));
	saved_index.push_back(count);
	saved_skip.push_back(true);
	if (start_v->get_fit_index()==-1){
	  start_v->set_fit_index(count);
	  count ++;
	}
      }else if (i+1==pts.size()){ // last point
	saved_pts.push_back(pts.at(i));
	saved_index.push_back(count);
	saved_skip.push_back(true);
	if (end_v->get_fit_index()==-1){
	  end_v->set_fit_index(count);
	  count ++;
	}
      } 
    } // loop over points ...
    // fill the segments ...
    sg->set_fit_associate_vec(saved_pts, saved_index, saved_skip);
    //    std::cout << saved_pts.size() << " " << saved_index.size() << " " << saved_skip.size() << std::endl;
  } // loop over segment

  // deal with all the vertex again ...
  for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
    WCPPID::ProtoVertex *vtx = it->first;
    double dis_cut = vtx->get_fit_range();
    count = vtx->get_fit_index();
    Point pt = vtx->get_fit_pt();
    
    std::set<std::pair<int,int> > temp_2dut, temp_2dvt, temp_2dwt; 
    form_point_association(pt, temp_2dut, temp_2dvt, temp_2dwt, ct_point_cloud, dis_cut, nlevel, time_cut); 
    // examine ... 
    std::vector<int> temp_results = ct_point_cloud.convert_3Dpoint_time_ch(pt); 
    temp_results.at(2)-=2400; 
    temp_results.at(3)-=4800; 
    std::vector<float> temp_flag;

    temp_flag = examine_point_association(temp_results, temp_2dut, temp_2dvt, temp_2dwt, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,true,charge_cut); 
  
    
    map_3D_2DU_set[count] = std::make_pair(temp_2dut,temp_flag.at(0)); 
    map_3D_2DV_set[count] = std::make_pair(temp_2dvt,temp_flag.at(1)); 
    map_3D_2DW_set[count] = std::make_pair(temp_2dwt,temp_flag.at(2)); 
    map_3D_tuple[count] = std::make_tuple(vtx,(WCPPID::ProtoSegment*)0,0); 
	  
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
   
  }
}

void WCPPID::PR3DCluster::organize_segments_path(std::map<WCPPID::ProtoVertex*, WCPPID::ProtoSegmentSet >& map_vertex_segments, std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet >& map_segment_vertices, double low_dis_limit, double end_point_limit){
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
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
    bool flag_startv_end = true;
    bool flag_endv_end = true;
    if (map_vertex_segments[start_v].size()>1) flag_startv_end = false;
    if (map_vertex_segments[end_v].size()>1) flag_endv_end = false;

    PointVector pts;
    std::vector<WCPointCloud<double>::WCPoint> temp_wcps_vec = sg->get_wcpt_vec();
    Point start_p, end_p;
    
    // start_ point
    start_p.x = temp_wcps_vec.front().x;
    start_p.y = temp_wcps_vec.front().y;
    start_p.z = temp_wcps_vec.front().z;

    if (flag_startv_end){
      Point p2(temp_wcps_vec.front().x, temp_wcps_vec.front().y, temp_wcps_vec.front().z);
      double dis1 = 0;
      for (auto it = temp_wcps_vec.begin(); it!=temp_wcps_vec.end(); it++){
	p2.x = (*it).x;
	p2.y = (*it).y;
	p2.z = (*it).z;
	dis1 = sqrt(pow(start_p.x-p2.x,2)+pow(start_p.y-p2.y,2)+pow(start_p.z-p2.z,2));
	if (dis1 > low_dis_limit) break;
      }
      if (dis1!=0){
	start_p.x += (start_p.x-p2.x)/dis1*end_point_limit;
	start_p.y += (start_p.y-p2.y)/dis1*end_point_limit;
	start_p.z += (start_p.z-p2.z)/dis1*end_point_limit;
      }
    }
    
    // end_point
    end_p.x = temp_wcps_vec.back().x;
    end_p.y = temp_wcps_vec.back().y;
    end_p.z = temp_wcps_vec.back().z;
    if (flag_endv_end){
      Point p2(temp_wcps_vec.back().x, temp_wcps_vec.back().y, temp_wcps_vec.back().z);
      double dis1 = 0;
      for (auto it = temp_wcps_vec.rbegin(); it!=temp_wcps_vec.rend(); it++){
	p2.x = (*it).x;
	p2.y = (*it).y;
	p2.z = (*it).z;
	dis1 = sqrt(pow(end_p.x-p2.x,2)+pow(end_p.y-p2.y,2)+pow(end_p.z-p2.z,2));
	if (dis1 > low_dis_limit) break;
      }
      if (dis1!=0){
	end_p.x += (end_p.x - p2.x)/dis1 * end_point_limit;
	end_p.y += (end_p.y - p2.y)/dis1 * end_point_limit;
	end_p.z += (end_p.z - p2.z)/dis1 * end_point_limit;
      }
    }
    start_v->set_fit_pt(start_p);
    end_v->set_fit_pt(end_p);

    // middle points
    pts.push_back(start_p);
    for (size_t i=0;i!=temp_wcps_vec.size(); i++){
      Point p1(temp_wcps_vec.at(i).x, temp_wcps_vec.at(i).y, temp_wcps_vec.at(i).z);
      double dis = low_dis_limit;
      double dis1 = sqrt(pow(p1.x-end_p.x,2) + pow(p1.y-end_p.y,2) + pow(p1.z-end_p.z,2));
      if (pts.size()>0)
	dis = sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2));

      if (dis1 < low_dis_limit * 0.8){
	continue;
      }else if (dis < low_dis_limit * 0.8 ){
	continue;
      }else if (dis < low_dis_limit * 1.6){
	pts.push_back(p1);
	//std::cout << p1 << " " << sqrt(pow(p1.x-pts.back().x,2)+pow(p1.y-pts.back().y,2)+pow(p1.z-pts.back().z,2))/units::cm << std::endl;
      }else{
	int npoints = std::round(dis/low_dis_limit);
	
	Point p_save = pts.back();
	for (int j=0;j!=npoints;j++){
	  Point p(p_save.x + (p1.x-p_save.x) / npoints * (j+1),
		  p_save.y + (p1.y-p_save.y) / npoints * (j+1),
		  p_save.z + (p1.z-p_save.z) / npoints * (j+1));
	  pts.push_back(p);
	  //	std::cout << p << " " << sqrt(pow(p.x-pts.back().x,2)+pow(p.y-pts.back().y,2)+pow(p.z-pts.back().z,2))/units::cm << std::endl;
	}
	
      }
    }
    pts.push_back(end_p);

    sg->set_point_vec(pts);
    //std::cout << sg << " " << start_p << " " << end_p << " " << temp_wcps_vec.size() << " " << pts.size() << std::endl;
    //    std::cout << sg << flag_startv_end << " " << flag_endv_end << std::endl;
    //    std::cout << sg << " " << start_v << " " << end_v << std::endl;
  }
}
    


void WCPPID::PR3DCluster::collect_charge_multi_trajectory(std::map<WCPPID::ProtoSegment*, WCPPID::ProtoVertexSet>& map_segment_vertices, WCP::ToyCTPointCloud& ct_point_cloud, double dis_cut, double range_cut){
  //clear up ...
  collected_charge_map.clear();
  
  std::set<std::pair<int,int>> existing_tcs;

  // form a set cotaining everything inside the cluster
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    GeomWireSelection& uwires = mcell->get_uwires();
    GeomWireSelection& vwires = mcell->get_vwires();
    GeomWireSelection& wwires = mcell->get_wwires();
    for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
      const GeomWire *wire = (*it1);
      existing_tcs.insert(std::make_pair(time_slice,wire->channel()));
    }
    for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
      const GeomWire *wire = (*it1);
      existing_tcs.insert(std::make_pair(time_slice,wire->channel()));
    }
    for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
      const GeomWire *wire = (*it1);
      existing_tcs.insert(std::make_pair(time_slice,wire->channel()));
    }
  }


  for (auto it2 = map_segment_vertices.begin(); it2 != map_segment_vertices.end(); it2++){
    WCPPID::ProtoSegment* sg1 = it2->first;
    
    // form a trajectory according to dis and fine tracking?
    PointVector traj_pts;
    //  PointVector& pts = get_fine_tracking_path();
    
    // not using fine tracking results ... 
    PointVector pts;
    std::vector<WCPointCloud<double>::WCPoint>& path_wcps = sg1->get_wcpt_vec();
    for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
      Point p((*it).x, (*it).y, (*it).z);
      pts.push_back(p);
    }

     //std::cout << "trajectory points " << pts.size() << std::endl;

    // if (cluster_id==13){
    //   for (int i=0; i!=pts.size(); i++){
    //     std::cout << i << " " << pts.at(i).x/units::cm << " " << pts.at(i).y/units::cm << " "<< pts.at(i).z/units::cm << std::endl;
    //   }
    // }

    for (int i=0; i!=pts.size(); i++){
      if (pts.at(i).y <-120*units::cm || pts.at(i).y > 120*units::cm ||
	  pts.at(i).z < -5*units::cm || pts.at(i).z > 1070*units::cm) continue;
      
      if (traj_pts.size()==0){
	traj_pts.push_back(pts.at(i));
      }else{
	double dis = sqrt(pow(pts.at(i).x-pts.at(i-1).x,2) +
			  pow(pts.at(i).y-pts.at(i-1).y,2) +
			  pow(pts.at(i).z-pts.at(i-1).z,2));
	if (dis <= dis_cut){
	  traj_pts.push_back(pts.at(i));
	}else{
	  int nseg = dis / dis_cut + 1;
	  //	std::cout << i << " " << pts.at(i).x/units::cm << " " << pts.at(i).y/units::cm << " "<< pts.at(i).z/units::cm << " " << dis << " " << dis_cut << " " << nseg << std::endl;
	  for (int j=0; j!=nseg;j++){
	    Point temp_pt;
	    temp_pt.x = pts.at(i-1).x + (pts.at(i).x-pts.at(i-1).x) *(j+1.)/nseg;
	    temp_pt.y = pts.at(i-1).y + (pts.at(i).y-pts.at(i-1).y) *(j+1.)/nseg;
	    temp_pt.z = pts.at(i-1).z + (pts.at(i).z-pts.at(i-1).z) *(j+1.)/nseg;
	    traj_pts.push_back(temp_pt);
	  }
	}
      }
    }

    // collect the nearby points, and compare with existing maps
    for (size_t i=0;i!=traj_pts.size();i++){
      //std::cout << i << " " << traj_pts.at(i).x/units::cm << " " << traj_pts.at(i).y/units::cm << " " << traj_pts.at(i).z/units::cm << " " << range_cut << std::endl;
      WCP::CTPointCloud<double> nearby_points = ct_point_cloud.get_closest_points(traj_pts.at(i),range_cut,0);
      
      //    std::cout << "0 " << nearby_points.pts.size() << std::endl;
      
      for (size_t j=0;j!=nearby_points.pts.size();j++){
	if (existing_tcs.find(std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)) == existing_tcs.end()){
	  collected_charge_map[std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)] = std::make_pair(nearby_points.pts.at(j).charge,nearby_points.pts.at(j).charge_err);
	}
      }
      nearby_points = ct_point_cloud.get_closest_points(traj_pts.at(i),range_cut,1);
      //std::cout << "1 " << nearby_points.pts.size() << std::endl;
      
      for (size_t j=0;j!=nearby_points.pts.size();j++){
	if (existing_tcs.find(std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)) == existing_tcs.end()){
	  collected_charge_map[std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)] = std::make_pair(nearby_points.pts.at(j).charge,nearby_points.pts.at(j).charge_err);
	}
      }
      nearby_points = ct_point_cloud.get_closest_points(traj_pts.at(i),range_cut,2);
      
      //    std::cout << "2 " << nearby_points.pts.size() << std::endl;
      
      for (size_t j=0;j!=nearby_points.pts.size();j++){
	if (existing_tcs.find(std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)) == existing_tcs.end()){
	  collected_charge_map[std::make_pair(nearby_points.pts.at(j).time_slice,nearby_points.pts.at(j).channel)] = std::make_pair(nearby_points.pts.at(j).charge,nearby_points.pts.at(j).charge_err);
	}
      }
    }
    
    //    std::cout << collected_charge_map.size() << std::endl;
    
  }

  
}
