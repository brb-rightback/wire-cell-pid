#include "WCPPID/ProtoVertex.h"

using namespace WCP;

WCPPID::ProtoVertex::ProtoVertex(WCP::WCPointCloud<double>::WCPoint& wcpt)
  : wcpt(wcpt)
  , flag_fit(false)
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

void WCPPID::ProtoVertex::set_fit_pt(WCP::Point& tmp_fit_pt){
  fit_pt = tmp_fit_pt;
  flag_fit = true;
}
