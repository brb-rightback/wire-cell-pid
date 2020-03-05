#include "WCPPID/ProtoSegment.h"

using namespace WCP;

WCPPID::ProtoSegment::ProtoSegment(std::list<WCP::WCPointCloud<double>::WCPoint >& path_wcps )
  : flag_fit(false)
{
  
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    wcpt_vec.push_back(*it);
  }
  
}

WCPPID::ProtoSegment::~ProtoSegment(){
}



void WCPPID::ProtoSegment::set_fit_vec(std::vector<WCP::Point >& tmp_fit_pt_vec, std::vector<double>& tmp_dQ_vec, std::vector<double>& tmp_dx_vec){
  flag_fit = true;
  
  fit_pt_vec = tmp_fit_pt_vec;
  dQ_vec = tmp_dQ_vec;
  dx_vec = tmp_dx_vec;
  dQ_dx_vec.resize(dQ_vec.size(),0);
  for (size_t i=0;i!=dQ_vec.size();i++){
    dQ_dx_vec.at(i) = dQ_vec.at(i)/(dx_vec.at(i)+1e-9);
  }
}

void WCPPID::ProtoSegment::clear_fit(){
  flag_fit = false;
  fit_pt_vec.clear();
  dQ_vec.clear();
  dx_vec.clear();
  dQ_dx_vec.clear();
}
