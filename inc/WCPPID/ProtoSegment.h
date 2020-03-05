#ifndef WCPPID_PROTOSEGMENT_H
#define WCPPID_PROTOSEGMENT_H

#include "WCPPID/ProtoVertex.h"
#include "WCPData/ToyPointCloud.h"

namespace WCPPID{
  class ProtoSegment{
  public:
    // initial creation
    ProtoSegment(std::list<WCP::WCPointCloud<double>::WCPoint >& path_wcps );
    ~ProtoSegment();

    std::vector<WCP::WCPointCloud<double>::WCPoint >& get_wcpt_vec(){return wcpt_vec;};
    std::vector<WCP::Point >& get_point_vec(){return fit_pt_vec;};
    std::vector<double>& get_dQ_vec(){return dQ_vec;};
    std::vector<double>& get_dx_vec(){return dx_vec;};
    std::vector<double>& get_dQ_dx_vec(){return dQ_dx_vec;};
    std::vector<double>& get_pu_vec(){return pu_vec;};
    std::vector<double>& get_pv_vec(){return pv_vec;};
    std::vector<double>& get_pw_vec(){return pw_vec;};
    std::vector<double>& get_pt_vec(){return pt_vec;};
    std::vector<double>& get_reduced_chi2_vec(){return reduced_chi2_vec;};
    
    // after fit
    void clear_fit();
    // set the data ...
    void set_fit_vec(std::vector<WCP::Point >& tmp_fit_pt_vec, std::vector<double>& tmp_dQ_vec, std::vector<double>& tmp_dx_vec, std::vector<double>& tmp_pu_vec, std::vector<double>& tmp_pv_vec, std::vector<double>& tmp_pw_vec, std::vector<double>& tmp_pt_vec, std::vector<double>& tmp_reduced_chi2_vec);
    

    // if fit_flag is set, the fit_pts are useful ...
    bool get_fit_flag(){return flag_fit;};
    void set_fit_flag(bool flag){flag_fit = flag;};
    
    // get point
    std::pair<double, WCP::Point> get_closest_point(WCP::Point &p);
    std::tuple<double, double, double> get_closest_2d_dis(WCP::Point &p);
    
  protected:
    std::vector<WCP::WCPointCloud<double>::WCPoint > wcpt_vec;

    std::vector<WCP::Point > fit_pt_vec;
    std::vector<double> dQ_vec;
    std::vector<double> dx_vec;
    std::vector<double> dQ_dx_vec;
    std::vector<double> pu_vec;
    std::vector<double> pv_vec;
    std::vector<double> pw_vec;
    std::vector<double> pt_vec;
    std::vector<double> reduced_chi2_vec;

    // point cloud ...
    WCP::ToyPointCloud* pcloud_fit;
    
    bool flag_fit;
  };
  typedef std::vector<ProtoSegment*> ProtoSegmentSelection;
}

#endif 
