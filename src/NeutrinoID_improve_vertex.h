//#include "Minuit2/FunctionMinimum.h"
//#include "Minuit2/MnUserParameterState.h"
//#include "Minuit2/MnMigrad.h"
//#include "WCPData/Line.h"

#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

void WCPPID::NeutrinoID::improve_vertex(WCPPID::PR3DCluster* temp_cluster){
  std::map<int, WCPPID::ProtoSegment*> map_id_seg;
  std::map<WCPPID::ProtoSegment*, int> map_seg_id;
  
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;

    map_id_seg[sg->get_id()] = sg;
    map_seg_id[sg] = sg->get_id();
    sg->clear_associate_points();
  }
  
  // find the relevant point clouds ...
  WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud()->get_cloud();
  std::vector<int>& point_sub_cluster_ids = temp_cluster->get_point_sub_cluster_ids();
  for (size_t i=0;i!=point_sub_cluster_ids.size();i++){
    if (point_sub_cluster_ids.at(i) == -1) continue;
    WCP::Point p(cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z);
    map_id_seg[point_sub_cluster_ids.at(i)]->add_associate_point(p);
  }

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

  //WCP::PointVector& pts = sg->get_associate_points();
  WCP::PointVector& pts = sg->get_point_vec();
  for (size_t i=0;i!=pts.size();i++){
    double dis_to_vertex = sqrt(pow(pts.at(i).x - vtx->get_fit_pt().x,2) + pow(pts.at(i).y - vtx->get_fit_pt().y,2) + pow(pts.at(i).z - vtx->get_fit_pt().z,2));
  if (dis_to_vertex < vertex_protect_dis || dis_to_vertex > fit_dis) continue;
  if (sg->get_closest_point(pts.at(i)).first > point_track_dis) continue;
    vec_points.back().push_back(pts.at(i));
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
    if (vec_points.at(i).size()>0) ncount++;
  }
  return ncount;
}

std::pair<WCPPID::ProtoSegment*, int> WCPPID::MyFCN::get_seg_info(int i){
  if (i < segments.size()){
    return std::make_pair(segments.at(i), vec_points.at(i).size());
  }
  return std::make_pair((WCPPID::ProtoSegment*)0, 0);
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



void WCPPID::NeutrinoID::fit_vertex(WCPPID::ProtoVertex *vtx, WCPPID::ProtoSegmentSet& sg_set, WCPPID::PR3DCluster* temp_cluster){

  WCPPID::MyFCN fcn(vtx, false, 1*units::cm, 1.0*units::cm, 2*units::cm, 6*units::cm);
  for (auto it = sg_set.begin(); it!=sg_set.end(); it++){
    fcn.AddSegment(*it);
  }
  std::vector<WCP::PointVector>& vec_points = fcn.get_vec_points();

  int npar = 3 + sg_set.size()*2;
  std::vector<double> xx(npar);
  
  
  xx[0] = 0;
  xx[1] = 0;
  xx[2] = 0;
  for (size_t i=0;i!=vec_points.size();i++){
    double theta = 0, phi = 0;
    if (vec_points.at(i).size()>0){
      Point p(0,0,0);
      for (size_t j=0;j!=vec_points.at(i).size();j++){
	p.x += vec_points.at(i).at(j).x - vtx->get_fit_pt().x;
	p.y += vec_points.at(i).at(j).y - vtx->get_fit_pt().y;
	p.z += vec_points.at(i).at(j).z - vtx->get_fit_pt().z;
      }
      TVector3 dir(p.x, p.y, p.z);
      theta = dir.Theta();
      phi = dir.Phi();
      xx[3 + 2*i] = theta;
      xx[3 + 2*i + 1] = phi;
    }
  }

  //  ROOT::Minuit2::MnUserParameters upar;
  //for (int i=0;i!=3;i++){
  // upar.Add(Form("x_%d",i),xx.at(i),0.01*units::cm);
  //}
  //for (int i=3;i!=npar;i++){
  // upar.Add(Form("x_%d",i),xx.at(i),0.01);
  //}
  // create MIGRAD minimizer
  //ROOT::Minuit2::MnMigrad migrad(fcn, upar);

  // minimize
  //ROOT::Minuit2::FunctionMinimum min = migrad();
  //const ROOT::Minuit2::MnUserParameterState& state = min.UserState();

  //  Point fit_vtx_p(vtx->get_fit_pt().x + state.Value(0), vtx->get_fit_pt().y + state.Value(1), vtx->get_fit_pt().z + state.Value(2));
  //  std::cout << vtx->get_fit_pt() << " " << fit_vtx_p << std::endl;


  /*
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
  for (auto it = sg_set.begin(); it!=sg_set.end(); it++){
    PointVector& pts = (*it)->get_point_vec();
    std::vector<double>& pu = (*it)->get_pu_vec();
    std::vector<double>& pv = (*it)->get_pv_vec();
    std::vector<double>& pw = (*it)->get_pw_vec();
    std::vector<double>& pt = (*it)->get_pt_vec();
    
    double dis1 = pow(pts.front().x - vtx->get_fit_pt().x,2) + pow(pts.front().y - vtx->get_fit_pt().y,2) + pow(pts.front().z - vtx->get_fit_pt().z,2);
    double dis2 = pow(pts.back().x - vtx->get_fit_pt().x,2) + pow(pts.back().y - vtx->get_fit_pt().y,2) + pow(pts.back().z - vtx->get_fit_pt().z,2);
    if (dis1 < dis2){
      pts.front() = fit_vtx_p;
      pu.front() = offset_u + 0.5 + (slope_yu * fit_vtx_p.y + slope_zu * fit_vtx_p.z);
      pv.front() = offset_v + 0.5 + (slope_yv * fit_vtx_p.y + slope_zv * fit_vtx_p.z)+2400;
      pw.front() = offset_w + 0.5 + (slope_yw * fit_vtx_p.y + slope_zw * fit_vtx_p.z)+4800;
      pt.front() = offset_t + 0.5 + slope_x * fit_vtx_p.x;
    }else{
      pts.back() = fit_vtx_p;
      pu.back() = offset_u + 0.5 + (slope_yu * fit_vtx_p.y + slope_zu * fit_vtx_p.z);
      pv.back() = offset_v + 0.5 + (slope_yv * fit_vtx_p.y + slope_zv * fit_vtx_p.z)+2400;
      pw.back() = offset_w + 0.5 + (slope_yw * fit_vtx_p.y + slope_zw * fit_vtx_p.z)+4800;
      pt.back() = offset_t + 0.5 + slope_x * fit_vtx_p.x;
    }
  }
  vtx->set_fit(fit_vtx_p, vtx->get_dQ(), vtx->get_dx(), offset_u + 0.5 + (slope_yu * fit_vtx_p.y + slope_zu * fit_vtx_p.z), offset_v + 0.5 + (slope_yv * fit_vtx_p.y + slope_zv * fit_vtx_p.z)+2400, offset_w + 0.5 + (slope_yw * fit_vtx_p.y + slope_zw * fit_vtx_p.z)+4800, offset_t + 0.5 + slope_x * fit_vtx_p.x, vtx->get_reduced_chi2());

  */
  
  //  std::cout << min.IsValid() << " " << min.Fval() << " " << fcn.get_chi2(xx) << std::endl;

  // std::cout << std::endl;

  //  fcn.print_points();
}



