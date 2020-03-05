#include "WCPPID/ProtoSegment.h"

using namespace WCP;

WCPPID::ProtoSegment::ProtoSegment(std::list<WCP::WCPointCloud<double>::WCPoint >& path_wcps )
  : flag_fit(false)
  , pcloud_fit(0)
{
  
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    wcpt_vec.push_back(*it);
  }
  
}

WCPPID::ProtoSegment::~ProtoSegment(){
  if (pcloud_fit != (ToyPointCloud*)0)
    delete pcloud_fit;
}



void WCPPID::ProtoSegment::set_fit_vec(std::vector<WCP::Point >& tmp_fit_pt_vec, std::vector<double>& tmp_dQ_vec, std::vector<double>& tmp_dx_vec, std::vector<double>& tmp_pu_vec, std::vector<double>& tmp_pv_vec, std::vector<double>& tmp_pw_vec, std::vector<double>& tmp_pt_vec, std::vector<double>& tmp_reduced_chi2_vec){
  flag_fit = true;
  
  fit_pt_vec = tmp_fit_pt_vec;
  dQ_vec = tmp_dQ_vec;
  dx_vec = tmp_dx_vec;
  dQ_dx_vec.resize(dQ_vec.size(),0);
  for (size_t i=0;i!=dQ_vec.size();i++){
    dQ_dx_vec.at(i) = dQ_vec.at(i)/(dx_vec.at(i)+1e-9);
  }

  pu_vec = tmp_pu_vec;
  pv_vec = tmp_pv_vec;
  pw_vec = tmp_pw_vec;
  pt_vec = tmp_pt_vec;
  reduced_chi2_vec = tmp_reduced_chi2_vec;

  if (pcloud_fit != (ToyPointCloud*)0) delete pcloud_fit;
  pcloud_fit = new ToyPointCloud();
  pcloud_fit->AddPoints(fit_pt_vec);
  pcloud_fit->build_kdtree_index();
}

void WCPPID::ProtoSegment::clear_fit(){
  flag_fit = false;
  fit_pt_vec.clear();
  dQ_vec.clear();
  dx_vec.clear();
  dQ_dx_vec.clear();

  pu_vec.clear();
  pv_vec.clear();
  pw_vec.clear();
  pt_vec.clear();
  reduced_chi2_vec.clear();
  
}

void WCPPID::ProtoSegment::print_dis(){
  for (size_t i=0;i!=wcpt_vec.size(); i++){
    Point p(wcpt_vec.at(i).x, wcpt_vec.at(i).y, wcpt_vec.at(i).z);
    std::pair<double, WCP::Point> results = get_closest_point(p);
    if (results.first > 0.6*units::cm)
      std::cout << i << " " << results.first/units::cm << std::endl;
  }
}

std::pair<double, WCP::Point> WCPPID::ProtoSegment::get_closest_point(WCP::Point &p){
  if (pcloud_fit != (ToyPointCloud*)0) 
    return pcloud_fit->get_closest_point(p);
  else{
    WCP::Point p1(0,0,0);
    return std::make_pair(-1,p1);
  }
}

std::tuple<double, double, double> WCPPID::ProtoSegment::get_closest_2d_dis(WCP::Point &p){
  if (pcloud_fit != (ToyPointCloud*)0) {
    std::pair<int, double> results_u = pcloud_fit->get_closest_2d_dis(p, 0);
    std::pair<int, double> results_v = pcloud_fit->get_closest_2d_dis(p, 1);
    std::pair<int, double> results_w = pcloud_fit->get_closest_2d_dis(p, 2);
    
    return std::make_tuple(results_u.second, results_v.second, results_w.second);
  }else{
    return std::make_tuple(-1,-1,-1);
  }
}
