

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

    for (int j=0;j!=num;j++){
      for (int k=0;k!=num;k++){
	index_index_dis[j][k] = std::make_tuple(-1,-1,1e9);
      }
    }

    for (int j=0;j!=num;j++){
      for (int k=j+1;k!=num;k++){
	index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(pt_clouds.at(k));
	index_index_dis[k][j] = index_index_dis[j][k];
      }
    }


    // test check the main
    for (int i=0;i!=num;i++){
      if (i==max_cluster) continue;
      std::cout << "number points: " << pt_clouds.at(max_cluster)->get_num_points() << " " << pt_clouds.at(i)->get_num_points() << " " << std::get<2>(index_index_dis[max_cluster][i])/units::cm << std::endl;
      
      if (fabs(pt_clouds.at(i)->get_num_points()-310)>10) continue;
      bool flag = check_connectivity(index_index_dis[max_cluster][i], cloud, ct_point_cloud, pt_clouds.at(max_cluster), pt_clouds.at(i) );
      
      
      
      
      
    }




    
    for (int i=0;i!=num;i++){
      delete pt_clouds.at(i);
    }
  }
  
}


bool WCPPID::PR3DCluster::check_connectivity(std::tuple<int, int, double>& index_index_dis, WCP::WCPointCloud<double>& cloud, WCP::ToyCTPointCloud& ct_point_cloud, WCP::ToyPointCloud* pc1, WCP::ToyPointCloud* pc2){


  // parallel case 1 and perpendicular case 2 
  TVector3 drift_dir(1,0,0);
  // pronlonged case for U 3 and V 4 ...
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);
  
  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis));
  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis));
  Point p1(wp1.x,wp1.y,wp1.z);
  Point p2(wp2.x,wp2.y,wp2.z);

  // directions
  TVector3 dir1 = VHoughTrans(p1, 15*units::cm, pc1);
  TVector3 dir2 = VHoughTrans(p2, 15*units::cm, pc2);
  dir1 *= -1;
  dir2 *= -1;

  TVector3 tempV1(0, p2.y - p1.y, p2.z - p1.z);
  TVector3 tempV5;
  // prolonged U ...
  double angle1 = tempV1.Angle(U_dir);
  tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle1),0);
  angle1 = tempV5.Angle(drift_dir);
  // prolonged V ...
  double angle2 = tempV1.Angle(V_dir);
  tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle2),0);
  angle2 = tempV5.Angle(drift_dir);
  tempV5.SetXYZ(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
  // prolonged W ...
  double angle3 = tempV1.Angle(W_dir);
  tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle3),0);
  angle3 = tempV5.Angle(drift_dir);
  


  
  
  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
  // check 2 pitch ...
  double step_size = 0.6*units::cm;
  int num_steps = std::round(dis/step_size);
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




    
    std::cout << i << " " << scores.at(0) << " " << scores.at(1) << " " << scores.at(2) << " " << scores.at(3) << " " << scores.at(4) << " " << scores.at(5) << std::endl;
  }
  
}


std::vector<bool> WCPPID::PR3DCluster::check_direction(TVector3& v1){
  
}
