

std::vector<SMGCSelection> WCPPID::PR3DCluster::Examine_graph(WCP::ToyCTPointCloud& ct_point_cloud){
  if (graph!=(MCUGraph*)0)
    delete graph;
  if (point_cloud==(ToyPointCloud*)0)
    Create_point_cloud();

  // form connected_pieces ...
  const int N = point_cloud->get_num_points();
  graph = new MCUGraph(N);
  Establish_close_connected_graph();

  Connect_graph_overclustering_protection(ct_point_cloud);
  std::vector<int> component(num_vertices(*graph));
  const int num = connected_components(*graph,&component[0]);

  std::vector<SMGCSelection> sep_mcells;
  std::set<SlimMergeGeomCell*> used_mcells;
  for (int i=0;i!=num;i++){
    SMGCSelection mcells;
    sep_mcells.push_back(mcells);
  }

  
  
  WCP::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  std::vector<int>::size_type i;
  for (i=0;i!=component.size(); ++i){
    SlimMergeGeomCell *mcell = cloud.pts[i].mcell;
    if (used_mcells.find(mcell)==used_mcells.end()){
      used_mcells.insert(mcell);
      sep_mcells[component[i]].push_back(mcell);
    }
    //pt_clouds.at(component[i])->AddPoint(cloud.pts[i],cloud_u.pts[i],cloud_v.pts[i],cloud_w.pts[i]);
  }

  // std::cout << num << std::endl;
  // for (int i=0;i!=num;i++){
  //   std::cout << i << " " << sep_mcells.at(i).size() << std::endl;
  // }
  
  
  return sep_mcells;
}


void WCPPID::PR3DCluster::Connect_graph_overclustering_protection(WCP::ToyCTPointCloud& ct_point_cloud){
  WCP::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCP::WC2DPointCloud<double>& cloud_u = point_cloud->get_cloud_u();
  WCP::WC2DPointCloud<double>& cloud_v = point_cloud->get_cloud_v();
  WCP::WC2DPointCloud<double>& cloud_w = point_cloud->get_cloud_w();

  
  // now form the connected components
  std::vector<int> component(num_vertices(*graph));
  const int num = connected_components(*graph,&component[0]);

  if (num > 1){

    int max_cluster = -1; int max_size = -1;
    std::vector<ToyPointCloud*> pt_clouds;
    for (int j=0;j!=num;j++){
      ToyPointCloud *pt_cloud = new ToyPointCloud();
      pt_clouds.push_back(pt_cloud);
    }
    
    std::vector<int>::size_type i;
    for (i=0;i!=component.size(); ++i){
      pt_clouds.at(component[i])->AddPoint(cloud.pts[i],cloud_u.pts[i],cloud_v.pts[i],cloud_w.pts[i]);
    }
    for (int j=0;j!=num;j++){
      pt_clouds.at(j)->build_kdtree_index();

      if (pt_clouds.at(j)->get_num_points() > max_size){
	max_cluster = j;
	max_size = pt_clouds.at(j)->get_num_points();
      }
    }
    // actual content ...
    
    // closest distance approach ...
    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis(num, std::vector< std::tuple<int,int,double> >(num));
    //    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_mst(num, std::vector< std::tuple<int,int,double> >(num));

    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_dir1(num, std::vector< std::tuple<int,int,double> >(num));
    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_dir2(num, std::vector< std::tuple<int,int,double> >(num));
    //std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_dir_mst(num, std::vector< std::tuple<int,int,double> >(num));
     
    // initialization ...
    for (int j=0;j!=num;j++){
      for (int k=0;k!=num;k++){
	index_index_dis[j][k] = std::make_tuple(-1,-1,1e9);
	//	index_index_dis_mst[j][k] = std::make_tuple(-1,-1,1e9);

	index_index_dis_dir1[j][k] = std::make_tuple(-1,-1,1e9);
  	index_index_dis_dir2[j][k] = std::make_tuple(-1,-1,1e9);
	//index_index_dis_dir_mst[j][k] = std::make_tuple(-1,-1,1e9);
      }
    }
    for (int j=0;j!=num;j++){
      for (int k=j+1;k!=num;k++){
	// closest distance ...
	std::tuple<int, int, double>  temp_index_index_dis = pt_clouds.at(j)->get_closest_points(pt_clouds.at(k));
	
	// close distance ...
	if (std::get<0>(temp_index_index_dis) != -1){
	  index_index_dis[j][k] = temp_index_index_dis; 
	  
	  bool flag = check_connectivity(index_index_dis[j][k], cloud, ct_point_cloud, pt_clouds.at(j), pt_clouds.at(k));
	  if (!flag) index_index_dis[j][k] = std::make_tuple(-1,-1,1e9);
	  
	  index_index_dis[k][j] = index_index_dis[j][k];

	  
	  // direction ...
	  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(temp_index_index_dis));
  	  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(temp_index_index_dis));
  	  Point p1(wp1.x,wp1.y,wp1.z);
  	  Point p2(wp2.x,wp2.y,wp2.z);
	  
  	  TVector3 dir1 = VHoughTrans(p1, 30*units::cm, pt_clouds.at(j));
  	  TVector3 dir2 = VHoughTrans(p2, 30*units::cm, pt_clouds.at(k));
  	  dir1 *= -1;
  	  dir2 *= -1;
	  
  	  std::pair<int,double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(p1, dir1, 80*units::cm, 5*units::cm, 7.5, 3*units::cm);
	  if (result1.first >=0){
	    index_index_dis_dir1[j][k] = std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
	    bool flag = check_connectivity(index_index_dis_dir1[j][k], cloud, ct_point_cloud, pt_clouds.at(j), pt_clouds.at(k));
	    if (!flag) index_index_dis_dir1[j][k] = std::make_tuple(-1,-1,1e9);
	    index_index_dis_dir1[k][j] = index_index_dis_dir1[j][k];
	  }
  	  
	  std::pair<int,double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(p2, dir2, 80*units::cm, 5*units::cm, 7.5, 3*units::cm);
	  if (result2.first >=0){
	    index_index_dis_dir2[j][k] = std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
	    bool flag = check_connectivity(index_index_dis_dir2[j][k], cloud, ct_point_cloud, pt_clouds.at(j), pt_clouds.at(k));
	    if (!flag) index_index_dis_dir2[j][k] = std::make_tuple(-1,-1,1e9);
	    index_index_dis_dir2[k][j] = index_index_dis_dir2[j][k];
	  }
	}
	  
      } // loop over separated pieces ...
    }
    
    for (int j=0; j!=num; j++){
      for (int k=j+1; k!=num; k++){
	// adding edges ...
	if (std::get<0>(index_index_dis[j][k])>=0){
	  auto edge = add_edge(std::get<0>(index_index_dis[j][k]),std::get<1>(index_index_dis[j][k]), WCPPID::EdgeProp(std::get<2>(index_index_dis[j][k])), *graph);
	}

	if (std::get<0>(index_index_dis_dir1[j][k])>=0){
	  if (std::get<2>(index_index_dis_dir1[j][k])>5*units::cm){
	    auto edge = add_edge(std::get<0>(index_index_dis_dir1[j][k]),std::get<1>(index_index_dis_dir1[j][k]), WCPPID::EdgeProp(std::get<2>(index_index_dis_dir1[j][k])*1.2), *graph);
	  }else{
	    auto edge = add_edge(std::get<0>(index_index_dis_dir1[j][k]),std::get<1>(index_index_dis_dir1[j][k]), WCPPID::EdgeProp(std::get<2>(index_index_dis_dir1[j][k])), *graph);
	  }
	}

	if (std::get<0>(index_index_dis_dir2[j][k])>=0){
	  if (std::get<2>(index_index_dis_dir2[j][k])>5*units::cm){
	    auto edge = add_edge(std::get<0>(index_index_dis_dir2[j][k]),std::get<1>(index_index_dis_dir2[j][k]), WCPPID::EdgeProp(std::get<2>(index_index_dis_dir2[j][k])*1.2), *graph);
	  }else{
	    auto edge = add_edge(std::get<0>(index_index_dis_dir2[j][k]),std::get<1>(index_index_dis_dir2[j][k]), WCPPID::EdgeProp(std::get<2>(index_index_dis_dir2[j][k])), *graph);
	  }
	}
	
      }
    }

    // test check the main
    /* for (int i=0;i!=num;i++){ */
    /*   for (int j=i+1;j!=num;j++){ */
    /* 	if (fabs(pt_clouds.at(j)->get_num_points()-425)<10 || pt_clouds.at(j)->get_num_points()>1000){ */
    /* 	  std::cout << "number points: " << pt_clouds.at(j)->get_num_points() << " " << pt_clouds.at(i)->get_num_points() << " " << std::get<2>(index_index_dis[j][i])/units::cm << " " << std::get<2>(index_index_dis_dir1[j][i])/units::cm << " " << std::get<2>(index_index_dis_dir2[j][i])/units::cm << std::endl; */
      
    /* 	  if (pt_clouds.at(i)->get_num_points() ==9){ */
    /* 	    bool flag1 = check_connectivity(index_index_dis[j][i], cloud, ct_point_cloud, pt_clouds.at(max_cluster), pt_clouds.at(i)); */
    /* 	    bool flag2 = check_connectivity(index_index_dis_dir1[j][i], cloud, ct_point_cloud, pt_clouds.at(max_cluster), pt_clouds.at(i)); */
    /* 	    bool flag3 = check_connectivity(index_index_dis_dir2[j][i], cloud, ct_point_cloud, pt_clouds.at(max_cluster), pt_clouds.at(i)); */
    /* 	    std::cout << flag1 << " " << flag2 << " " << flag3 << std::endl; */
    /* 	  } */
    /* 	} */
    /*   } */
    /* } */




    
    for (int i=0;i!=num;i++){
      delete pt_clouds.at(i);
    }
  }
  
}


bool WCPPID::PR3DCluster::check_connectivity(std::tuple<int, int, double>& index_index_dis, WCP::WCPointCloud<double>& cloud, WCP::ToyCTPointCloud& ct_point_cloud, WCP::ToyPointCloud* pc1, WCP::ToyPointCloud* pc2){
  if (std::get<0>(index_index_dis)==-1 || std::get<1>(index_index_dis)==-1) return false;
  
  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis));
  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis));
  Point p1(wp1.x,wp1.y,wp1.z);
  Point p2(wp2.x,wp2.y,wp2.z);

  // std::cout << wp1.index_u << " " << wp1.index_v << " " << wp1.index_w << " " << p1 << std::endl;
  //std::cout << wp2.index_u << " " << wp2.index_v << " " << wp2.index_w << " " << p2 << std::endl;
  
  // directions
  TVector3 dir1 = VHoughTrans(p1, 15*units::cm, pc1);  dir1 *= -1;
  TVector3 dir2 = VHoughTrans(p2, 15*units::cm, pc2);   dir2 *= -1;
  TVector3 dir3(p1.x - p2.x, p1.y - p2.y, p1.z- p2.z);

  std::vector<bool> flag_1 = check_direction(dir1);
  std::vector<bool> flag_2 = check_direction(dir2);
  std::vector<bool> flag_3 = check_direction(dir3);

  bool flag_prolonged_u = false;
  bool flag_prolonged_v = false;
  bool flag_prolonged_w = false;
  bool flag_parallel = false;

  if (flag_3.at(0) && (flag_1.at(0) || flag_2.at(0))) flag_prolonged_u = true;
  if (flag_3.at(1) && (flag_1.at(1) || flag_2.at(1))) flag_prolonged_v = true;
  if (flag_3.at(2) && (flag_1.at(2) || flag_2.at(2))) flag_prolonged_w = true;
  if (flag_3.at(3) && (flag_1.at(3) && flag_2.at(3))) flag_parallel = true;
  
  //std::cout << flag_prolonged_u << " " << flag_prolonged_v << " " << flag_prolonged_w << " " << flag_parallel << std::endl;
  
  
  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
  // check 2 pitch ...
  double step_size = 0.6*units::cm;
  int num_steps = std::round(dis/step_size);

  int num_bad[5]={0,0,0,0,0};
 
  
  for (int i=0;i!=num_steps;i++){
    Point test_p;
    test_p.x = p1.x + (p2.x-p1.x)/(num_steps+1.)*(i+1);
    test_p.y = p1.y + (p2.y-p1.y)/(num_steps+1.)*(i+1);
    test_p.z = p1.z + (p2.z-p1.z)/(num_steps+1.)*(i+1);

    std::vector<int> scores;
    if (i==0 || i+1 == num_steps)
      scores = ct_point_cloud.test_good_point(test_p, dis/(num_steps+1.)*0.98);
    else
      scores = ct_point_cloud.test_good_point(test_p);

    int num_bad_details= 0;
    if (scores.at(0) + scores.at(3)==0){
      if (!flag_prolonged_u) num_bad[0] ++;
      num_bad_details ++;
    }
    if (scores.at(1) + scores.at(4)==0){
      if (!flag_prolonged_v) num_bad[1] ++;
      num_bad_details ++;
    }
    if (scores.at(2) + scores.at(5)==0){
      if (!flag_prolonged_w) num_bad[2] ++;
      num_bad_details ++;
    }

    // parallel case ...
    if (flag_parallel){
      if (num_bad_details>1) num_bad[3] ++;
    }else{
      if (num_bad_details>0) num_bad[3] ++;
    }
    //    std::cout << i << " " << num_bad_details << std::endl;
    //    std::cout << i << " " << scores.at(0) << " " << scores.at(1) << " " << scores.at(2) << " " << scores.at(3) << " " << scores.at(4) << " " << scores.at(5) << std::endl;
  }
  //  std::cout << num_bad[0] << " " << num_bad[1] << " " << num_bad[2] << " " << num_bad[3] << " " << num_steps << std::endl;

  
  
  // prolonged case ...
  if (num_bad[0] <=2 && num_bad[1] <= 2 && num_bad[2] <=2 &&
      (num_bad[0] + num_bad[1] + num_bad[2] <=3) && 
      num_bad[0] < 0.1 * num_steps && num_bad[1] < 0.1 * num_steps && num_bad[2] < 0.1 * num_steps &&
      (num_bad[0] + num_bad[1] + num_bad[2]) < 0.15 * num_steps){
    if (flag_prolonged_u && flag_prolonged_v && flag_prolonged_w)
      if (num_bad[3] >0.6 * num_steps) return false;
    return true;
  }else if (num_bad[3] <=2 && num_bad[3] < 0.1*num_steps){
    return true;
  }
  

  return false;
}


std::vector<bool> WCPPID::PR3DCluster::check_direction(TVector3& v1){
   // parallel case 1 and perpendicular case 2 
  TVector3 drift_dir(1,0,0);
  // pronlonged case for U 3 and V 4 ...
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);


  TVector3 tempV1(0, v1.Y(), v1.Z());
  TVector3 tempV5;
  // prolonged U ...
  double angle1 = tempV1.Angle(U_dir);
  tempV5.SetXYZ(fabs(v1.X()),sqrt(pow(v1.Y(),2)+pow(v1.Z(),2))*sin(angle1),0);
  angle1 = tempV5.Angle(drift_dir);
  // prolonged V ...
  double angle2 = tempV1.Angle(V_dir);
  tempV5.SetXYZ(fabs(v1.X()),sqrt(pow(v1.Y(),2)+pow(v1.Z(),2))*sin(angle2),0);
  angle2 = tempV5.Angle(drift_dir);
  // prolonged W ...
  double angle3 = tempV1.Angle(W_dir);
  tempV5.SetXYZ(fabs(v1.X()),sqrt(pow(v1.Y(),2)+pow(v1.Z(),2))*sin(angle3),0);
  angle3 = tempV5.Angle(drift_dir);
  
  double angle4 = v1.Angle(drift_dir); // parallel ...
  std::vector<bool> results(4, false);

  // prolonged U
  if (angle1<12.5/180.*3.1415926) results.at(0) = true;
  if (angle2<12.5/180.*3.1415926) results.at(1) = true;
  if (angle3<12.5/180.*3.1415926) results.at(2) = true;
  if (fabs(angle4-3.1415926/2.)<10./180.*3.1415925) results.at(3) = true;
  
  return results;
}
