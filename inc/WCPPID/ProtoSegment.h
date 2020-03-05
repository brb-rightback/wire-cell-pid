#ifndef WCPPID_PROTOSEGMENT_H
#define WCPPID_PROTOSEGMENT_H

#include "WCPPID/ProtoVertex.h"

namespace WCPPID{
  class ProtoSegment{
  public:
    ProtoSegment(std::list<WCP::WCPointCloud<double>::WCPoint >& path_wcps );
    ~ProtoSegment();

    std::vector<WCP::WCPointCloud<double>::WCPoint >& get_wcpt_vec(){return wcpt_vec;};
    std::vector<WCP::Point >& get_pt_vec(){return fit_pt_vec;};
    void clear_pt_vec(){fit_pt_vec.clear();};

    void set_fit_pt_vec(std::vector<WCP::Point >& tmp_fit_pt_vec){fit_pt_vec = tmp_fit_pt_vec;};
    
    bool get_fit_flag(){return flag_fit;};
    void set_fit_flag(bool flag){flag_fit = flag;};
    
  protected:
    std::vector<WCP::WCPointCloud<double>::WCPoint > wcpt_vec;
    std::vector<WCP::Point > fit_pt_vec;

    ProtoVertex *v1;
    ProtoVertex *v2;

    bool flag_fit;
  };
  typedef std::vector<ProtoSegment*> ProtoSegmentSelection;
}

#endif 
