#include "WCPPID/ProtoSegment.h"

using namespace WCP;

WCPPID::ProtoSegment::ProtoSegment(std::list<WCP::WCPointCloud<double>::WCPoint >& path_wcps )
  : flag_fit(false)
{
  v1 = new ProtoVertex(path_wcps.front() );
  v2 = new ProtoVertex(path_wcps.back() );

  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    wcpt_vec.push_back(*it);
    Point p(it->x, it->y, it->z);
    fit_pt_vec.push_back(p);
  }
}

WCPPID::ProtoSegment::~ProtoSegment(){
  delete v1;
  delete v2;
}
