//#include "Minuit2/FunctionMinimum.h"
//#include "Minuit2/MnUserParameterState.h"
//#include "Minuit2/MnMigrad.h"
//#include "WCPData/Line.h"

#include "TMatrixDEigen.h"

#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

void WCPPID::NeutrinoID::improve_vertex(WCPPID::PR3DCluster* temp_cluster){
  

   /* for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){ */
   /*  WCPPID::ProtoSegment *sg = it->first; */
   /*  if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue; */
   /*  std::cout << sg->get_associate_points().size() << std::endl; */
   /* } */


  // find the vertex
  for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster->get_cluster_id() || it->second.size()<=1) continue;
    fit_vertex(vtx, it->second, temp_cluster);
  }
}

WCPPID::MyFCN::MyFCN(WCPPID::ProtoVertex* vtx, bool flag_vtx_constraint, double vtx_constraint_range, double vertex_protect_dis, double point_track_dis, double fit_dis) 
: vtx(vtx)
, flag_vtx_constraint(flag_vtx_constraint)
  , vtx_constraint_range(vtx_constraint_range) 
  , vertex_protect_dis(vertex_protect_dis)
  , point_track_dis(point_track_dis)
  , fit_dis(fit_dis)
{
  segments.clear();
  vec_points.clear();
}

WCPPID::MyFCN::~MyFCN(){
}

void WCPPID::MyFCN::print_points(){
  for (size_t i=0;i!=vec_points.size();i++){
    for (size_t j=0;j!=vec_points.at(i).size();j++){
      std::cout <<i << " " << j << " " << vec_points.at(i).at(j).x/units::cm << " " << vec_points.at(i).at(j).y/units::cm  << " " << vec_points.at(i).at(j).z/units::cm << std::endl;
    }
  }
}

void WCPPID::MyFCN::AddSegment(ProtoSegment *sg){
  // push in ...
  segments.push_back(sg);
  {
    WCP::PointVector pts;
    vec_points.push_back(pts);
  }

  WCP::Point center;
  double min_dis = 1e9;
  WCP::PointVector& pts = sg->get_point_vec();
  for (size_t i=0;i!=pts.size();i++){
    double dis_to_vertex = sqrt(pow(pts.at(i).x - vtx->get_fit_pt().x,2) + pow(pts.at(i).y - vtx->get_fit_pt().y,2) + pow(pts.at(i).z - vtx->get_fit_pt().z,2));
  if (dis_to_vertex < vertex_protect_dis || dis_to_vertex > fit_dis) continue;
  if (sg->get_closest_point(pts.at(i)).first > point_track_dis) continue;
    vec_points.back().push_back(pts.at(i));
    if (dis_to_vertex < min_dis){
      center = pts.at(i);
      min_dis = dis_to_vertex;
    }
  }

  // calculate the PCA ...
  if (vec_points.back().size()>1){

    // calculate center ...
    
    int nsum = vec_points.back().size();
    
    //    for (auto it = vec_points.back().begin(); it!=vec_points.back().end();it++){
    // }
    /*   center.x += (*it).x; */
    /*   center.y += (*it).y; */
    /*   center.z += (*it).z; */
    /*   nsum ++; */
    /* center.x /= nsum; center.y /= nsum; center.z /= nsum; */

    // Eigen vectors ...
    PointVector PCA_axis(3);
    for (int i=0;i!=3;i++){
      PCA_axis[i].x = 0;       PCA_axis[i].y = 0;       PCA_axis[i].z = 0;
    }
    double PCA_values[3]={0,0,0};
    
    TMatrixD cov_matrix(3,3);
    for (int i=0;i!=3;i++){
      for (int j=0;j!=3;j++){
	for (size_t k=0;k!=vec_points.back().size();k++){
	  if (i==0 && j==0){
	    cov_matrix(i,j) += (vec_points.back().at(k).x - center.x) * (vec_points.back().at(k).x - center.x);
	  }else if (i==0 && j==1){
	    cov_matrix(i,j) += (vec_points.back().at(k).x - center.x) * (vec_points.back().at(k).y - center.y);
	  }else if (i==0 && j==2){
	    cov_matrix(i,j) += (vec_points.back().at(k).x - center.x) * (vec_points.back().at(k).z - center.z);
	  }else if (i==1 && j==1){
	    cov_matrix(i,j) += (vec_points.back().at(k).y - center.y) * (vec_points.back().at(k).y - center.y);
	  }else if (i==1 && j==2){
	    cov_matrix(i,j) += (vec_points.back().at(k).y - center.y) * (vec_points.back().at(k).z - center.z);
	  }else if (i==2 && j==2){
	    cov_matrix(i,j) += (vec_points.back().at(k).z - center.z) * (vec_points.back().at(k).z - center.z);
	  }
	}
      }
    }
    cov_matrix(1,0) = cov_matrix(0,1);
    cov_matrix(2,0) = cov_matrix(0,2);
    cov_matrix(2,1) = cov_matrix(1,2);
    for (int i=0;i!=3;i++){
      for (int j=0;j!=3;j++){
	cov_matrix(i,j) = cov_matrix(i,j)/nsum;
      }
    }
    
    TMatrixDEigen eigen(cov_matrix);
    TMatrixD eigen_values = eigen.GetEigenValues();
    TMatrixD eigen_vectors = eigen.GetEigenVectors();

    PCA_values[0] = eigen_values(0,0) + pow(0.1*units::cm,2);
    PCA_values[1] = eigen_values(1,1) + pow(0.1*units::cm,2);
    PCA_values[2] = eigen_values(2,2) + pow(0.1*units::cm,2);

    for (int i=0;i!=3;i++){
      PCA_axis[i].x = eigen_vectors(0,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
      PCA_axis[i].y = eigen_vectors(1,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
      PCA_axis[i].z = eigen_vectors(2,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    }

    //    TVector3 v1(vec_points.back().front().x - vec_points.back().back().x, vec_points.back().front().y - vec_points.back().back().y, vec_points.back().front().z - vec_points.back().back().z);
    //TVector3 v2(PCA_axis[0].x, PCA_axis[0].y, PCA_axis[0].z);
    //std::cout << vec_points.back().size() << " " << sqrt(PCA_values[0])/units::cm << " " << sqrt(PCA_values[1])/units::cm << " " << sqrt(PCA_values[2])/units::cm << " " << v1.Angle(v2)/3.1415926*180. << std::endl;

    vec_PCA_dirs.push_back(std::make_tuple(PCA_axis[0], PCA_axis[1], PCA_axis[2]));
    vec_PCA_vals.push_back(std::make_tuple(PCA_values[0], PCA_values[1], PCA_values[2]));
    vec_centers.push_back(center);
    
  }else{
    Point a(0,0,0);
    vec_PCA_dirs.push_back(std::make_tuple(a,a,a));
    vec_PCA_vals.push_back(std::make_tuple(0, 0, 0));
    if (vec_points.back().size()==1){
      vec_centers.push_back(vec_points.back().back());
    }else{
      vec_centers.push_back(a);
    }
  }
  
  std::cout << vec_points.size() << " " << vec_points.back().size() << std::endl;
}

void WCPPID::MyFCN::update_fit_range(double tmp_vertex_protect_dis, double tmp_point_track_dis, double tmp_fit_dis){
  vertex_protect_dis = tmp_vertex_protect_dis;
  point_track_dis = tmp_point_track_dis;
  fit_dis = tmp_fit_dis;

  std::vector<WCPPID::ProtoSegment* > tmp_segments = segments;
  segments.clear();
  vec_points.clear();
  for (auto it = tmp_segments.begin(); it != tmp_segments.end(); it++){
    AddSegment(*it);
  }  
}

int WCPPID::MyFCN::get_fittable_tracks(){
  int ncount = 0;
  for (size_t i=0;i!=vec_points.size();i++){
    if (vec_points.at(i).size()>1) ncount++;
  }
  return ncount;
}

std::pair<WCPPID::ProtoSegment*, int> WCPPID::MyFCN::get_seg_info(int i){
  if (i < segments.size()){
    return std::make_pair(segments.at(i), vec_points.at(i).size());
  }
  return std::make_pair((WCPPID::ProtoSegment*)0, 0);
}


std::pair<bool, WCP::Point> WCPPID::MyFCN::FitVertex(){
  Point fit_pos = vtx->get_fit_pt();
  bool fit_flag = false;

  int ntracks = get_fittable_tracks();

  int npoints = 0;
  
  if (ntracks >1){
    // start the fit ...
    Eigen::VectorXd temp_pos_3D_init(3), temp_pos_3D(3); // to be fitted
    temp_pos_3D_init(0) = fit_pos.x;
    temp_pos_3D_init(1) = fit_pos.y;
    temp_pos_3D_init(2) = fit_pos.z;

    Eigen::Vector3d b(0,0,0);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3);
    
    
    for (size_t i=0;i!=vec_PCA_vals.size();i++){
      if (std::get<0>(vec_PCA_vals.at(i))>0){
	npoints += vec_points.at(i).size();
	// fill the matrix ... first row, second column
	Eigen::MatrixXd R(3,3);
	R(0,0) = 0; R(0,1) = 0; R(0,2) = 0;
	double val1 = sqrt(std::get<0>(vec_PCA_vals.at(i))/std::get<1>(vec_PCA_vals.at(i)));
	R(1,0) = val1 * std::get<1>(vec_PCA_dirs.at(i)).x;
	R(1,1) = val1 * std::get<1>(vec_PCA_dirs.at(i)).y;
	R(1,2) = val1 * std::get<1>(vec_PCA_dirs.at(i)).z;
	val1 = sqrt(std::get<0>(vec_PCA_vals.at(i))/std::get<2>(vec_PCA_vals.at(i)));
	R(1,0) = val1 * std::get<2>(vec_PCA_dirs.at(i)).x;
	R(1,1) = val1 * std::get<2>(vec_PCA_dirs.at(i)).y;
	R(1,2) = val1 * std::get<2>(vec_PCA_dirs.at(i)).z;

	Eigen::Vector3d data(vec_centers.at(i).x, vec_centers.at(i).y, vec_centers.at(i).z );
	data = R * data;
	Eigen::MatrixXd RT = R.transpose();

	b += RT * data;
	A += RT * R;
      }
    }
    // 1.69 is a factor to be tuned ...
    
    //add constraint ...
    if (flag_vtx_constraint){
      Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3,3);
      R(0,0) = 1./vtx_constraint_range*sqrt(npoints*1.69);
      R(1,1) = 1./vtx_constraint_range*sqrt(npoints*1.69);
      R(2,2) = 1./vtx_constraint_range*sqrt(npoints*1.69);
      
      Eigen::Vector3d data(fit_pos.x/vtx_constraint_range*sqrt(npoints*1.69), fit_pos.y/vtx_constraint_range*sqrt(npoints*1.69), fit_pos.z/vtx_constraint_range*sqrt(npoints*1.69));
      Eigen::MatrixXd RT = R.transpose();
      b += RT * data;
      A += RT * R;
    }
    
    Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
    solver.compute(A);
    temp_pos_3D = solver.solveWithGuess(b,temp_pos_3D_init);
    if (!std::isnan(solver.error())){
      fit_pos.x = temp_pos_3D(0);
      fit_pos.y = temp_pos_3D(1);
      fit_pos.z = temp_pos_3D(2);
      fit_flag = true;
    }
    
    std::cout << temp_pos_3D(0)/units::cm << " " << temp_pos_3D(1)/units::cm << " " << temp_pos_3D(2)/units::cm << " " << temp_pos_3D_init(0)/units::cm << " " << temp_pos_3D_init(1)/units::cm << " " << temp_pos_3D_init(2)/units::cm << " " << sqrt(pow(temp_pos_3D(0) - temp_pos_3D_init(0),2) + pow(temp_pos_3D(1) - temp_pos_3D_init(1),2) + pow(temp_pos_3D(2) - temp_pos_3D_init(2),2))/units::cm << std::endl;
  }
  
  return std::make_pair(fit_flag, fit_pos);
}

void WCPPID::MyFCN::UpdateInfo(WCP::Point fit_pos, WCPPID::PR3DCluster* temp_cluster){
  // convertion to the u, v, w, z ...
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
  //mcell
  double first_t_dis = temp_cluster->get_point_cloud()->get_cloud().pts[0].mcell->GetTimeSlice()*time_slice_width - temp_cluster->get_point_cloud()->get_cloud().pts[0].x;
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
  
  
  // quick hack ...
  for (auto it = segments.begin(); it!=segments.end(); it++){
    PointVector& pts = (*it)->get_point_vec();
    std::vector<double>& pu = (*it)->get_pu_vec();
    std::vector<double>& pv = (*it)->get_pv_vec();
    std::vector<double>& pw = (*it)->get_pw_vec();
    std::vector<double>& pt = (*it)->get_pt_vec();
    
    double dis1 = pow(pts.front().x - vtx->get_fit_pt().x,2) + pow(pts.front().y - vtx->get_fit_pt().y,2) + pow(pts.front().z - vtx->get_fit_pt().z,2);
    double dis2 = pow(pts.back().x - vtx->get_fit_pt().x,2) + pow(pts.back().y - vtx->get_fit_pt().y,2) + pow(pts.back().z - vtx->get_fit_pt().z,2);

    if (dis1 < dis2){
      pts.front() = fit_pos;
      pu.front() = offset_u + 0.5 + (slope_yu * fit_pos.y + slope_zu * fit_pos.z);
      pv.front() = offset_v + 0.5 + (slope_yv * fit_pos.y + slope_zv * fit_pos.z)+2400;
      pw.front() = offset_w + 0.5 + (slope_yw * fit_pos.y + slope_zw * fit_pos.z)+4800;
      pt.front() = offset_t + 0.5 + slope_x * fit_pos.x;
    }else{
      pts.back() = fit_pos;
      pu.back() = offset_u + 0.5 + (slope_yu * fit_pos.y + slope_zu * fit_pos.z);
      pv.back() = offset_v + 0.5 + (slope_yv * fit_pos.y + slope_zv * fit_pos.z)+2400;
      pw.back() = offset_w + 0.5 + (slope_yw * fit_pos.y + slope_zw * fit_pos.z)+4800;
      pt.back() = offset_t + 0.5 + slope_x * fit_pos.x;
    }
  }
  vtx->set_fit(fit_pos, vtx->get_dQ(), vtx->get_dx(), offset_u + 0.5 + (slope_yu * fit_pos.y + slope_zu * fit_pos.z), offset_v + 0.5 + (slope_yv * fit_pos.y + slope_zv * fit_pos.z)+2400, offset_w + 0.5 + (slope_yw * fit_pos.y + slope_zw * fit_pos.z)+4800, offset_t + 0.5 + slope_x * fit_pos.x, vtx->get_reduced_chi2());
}

void WCPPID::NeutrinoID::fit_vertex(WCPPID::ProtoVertex *vtx, WCPPID::ProtoSegmentSet& sg_set, WCPPID::PR3DCluster* temp_cluster){

  WCPPID::MyFCN fcn(vtx, false, 0.58*units::cm, 1.0*units::cm, 2*units::cm, 6*units::cm);
  for (auto it = sg_set.begin(); it!=sg_set.end(); it++){
    fcn.AddSegment(*it);
  }
  
  std::pair<bool, Point> results = fcn.FitVertex();

  if (results.first) fcn.UpdateInfo(results.second, temp_cluster);
  

  
   //  fcn.print_points();
}






//double WCPPID::MyFCN::operator() (const std::vector<double> & xx) const{
//  return get_chi2(xx);
//}

/*
double WCPPID::MyFCN::get_chi2(const std::vector<double> & xx) const{
  double chi2 = 0;
  const int ntracks = segments.size();
  Point p(vtx->get_fit_pt().x + xx[0] ,vtx->get_fit_pt().y + xx[1] , vtx->get_fit_pt().z + xx[2] );

  std::vector<TVector3 > dir(ntracks);
  for (size_t i=0;i!=ntracks;i++){
    double theta = xx[3+2*i];
    double phi = xx[3+2*i+1];
    dir.at(i).SetXYZ(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
  }


  double fit_res = 0.1*units::cm;
  for (size_t i=0;i!=ntracks;i++){
    WCP::Line l1(p, dir.at(i));
    //    double tmp_chi2 = 0;
    for (size_t j=0;j!=vec_points.at(i).size();j++){
      WCP::Point p1 =vec_points[i].at(j); 
      chi2 += pow(l1.closest_dis(p1),2)/pow(fit_res,2);
    }
    // if (vec_points.at(i).size()>0)
    // chi2 += tmp_chi2/vec_points.at(i).size();
     
  }
  
  return chi2;
}
*/
