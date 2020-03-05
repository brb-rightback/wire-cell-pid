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
    std::vector<double>& get_dQ_vec(){return dQ_vec;};
    std::vector<double>& get_dx_vec(){return dx_vec;};

    std::vector<double> get_dQ_dx_vec(){return dQ_dx_vec;};
    
    // after fit
    void clear_fit();
    void set_fit_vec(std::vector<WCP::Point >& tmp_fit_pt_vec, std::vector<double>& tmp_dQ_vec, std::vector<double>& tmp_dx_vec);
    

    // if fit_flag is set, the fit_pts are useful ...
    bool get_fit_flag(){return flag_fit;};
    void set_fit_flag(bool flag){flag_fit = flag;};
    
  protected:
    std::vector<WCP::WCPointCloud<double>::WCPoint > wcpt_vec;

    std::vector<WCP::Point > fit_pt_vec;
    std::vector<double> dQ_vec;
    std::vector<double> dx_vec;
    std::vector<double> dQ_dx_vec;
    
    bool flag_fit;
  };
  typedef std::vector<ProtoSegment*> ProtoSegmentSelection;
}

#endif 
