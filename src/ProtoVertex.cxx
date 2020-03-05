#include "WCPPID/ProtoVertex.h"

using namespace WCP;

WCPPID::ProtoVertex::ProtoVertex(WCP::WCPointCloud<double>::WCPoint& wcpt)
  : wcpt(wcpt)
  , flag_fit(false)
  , dQ(0)
  , dx(-1)
{
  fit_pt.x = wcpt.x;
  fit_pt.y = wcpt.y;
  fit_pt.z = wcpt.z;
}

WCPPID::ProtoVertex::~ProtoVertex(){
}

void WCPPID::ProtoVertex::set_wcpt(WCP::WCPointCloud<double>::WCPoint& tmp_pt){
  wcpt = tmp_pt;
}

void WCPPID::ProtoVertex::set_fit(WCP::Point& tmp_fit_pt, double tmp_dQ, double tmp_dx, double tmp_pu, double tmp_pv, double tmp_pw, double tmp_pt, double tmp_reduced_chi2){
  fit_pt = tmp_fit_pt;
  flag_fit = true;

  dQ = tmp_dQ;
  dx = tmp_dx;
  pu = tmp_pu;
  pv = tmp_pv;
  pw = tmp_pw;
  pt = tmp_pt;
  reduced_chi2 = tmp_reduced_chi2;
  
}


double WCPPID::ProtoVertex::get_dis(WCP::Point& p){
  double dis = sqrt(pow(p.x - fit_pt.x,2)
		    + pow(p.y - fit_pt.y,2)
		    + pow(p.z - fit_pt.z,2) );
  return dis;
}

double WCPPID::ProtoVertex::get_fit_init_dis(){
  double dis = sqrt(pow(wcpt.x - fit_pt.x,2) +
		    pow(wcpt.y - fit_pt.y,2) +
		    pow(wcpt.z - fit_pt.z,2));
  return dis;
}
