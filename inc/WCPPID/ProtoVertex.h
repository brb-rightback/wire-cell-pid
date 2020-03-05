#ifndef WCPPID_PROTOVERTEX_H
#define WCPPID_PROTOVERTEX_H

#include "WCPData/WCPointCloud.h"

namespace WCPPID{
  class ProtoVertex{
  public:
    ProtoVertex(WCP::WCPointCloud<double>::WCPoint& wcpt);
    ~ProtoVertex();

    void set_wcpt(WCP::WCPointCloud<double>::WCPoint& tmp_pt);
    WCP::WCPointCloud<double>::WCPoint& get_wcpt() {return wcpt;};

    //get fit flag ...
    bool get_fit_flag(){return flag_fit;};
    void set_fit_flag(bool flag){flag_fit = flag;};

    // after trajectory fit ...
    void set_fit(WCP::Point& tmp_fit_pt, double tmp_dQ, double tmp_dx, double tmp_pu, double tmp_pv, double tmp_pw, double tmp_pt, double tmp_reduced_chi2);
    WCP::Point& get_fit_pt() {return fit_pt;};

    // after dQ/dx fit ...
    double get_dQ(){return dQ;};
    double get_dx(){return dx;};
    double get_dQ_dx(){return dQ/(dx+1e-9);};
    double get_pu(){return pu;};
    double get_pv(){return pv;};
    double get_pw(){return pw;};
    double get_pt(){return pt;};
    double get_reduced_chi2(){return reduced_chi2;};
    
    // get_distance ...  
    double get_dis(WCP::Point& p);
    
  protected:
    WCP::WCPointCloud<double>::WCPoint wcpt; // initial WCP point from the graph ...
    WCP::Point fit_pt; // best fit points ...

    double dQ;
    double dx;
    double pu, pv, pw, pt;
    double reduced_chi2;

    bool flag_fit;
  };
  typedef std::vector<ProtoVertex*> ProtoVertexSelection;
  typedef std::set<ProtoVertex*> ProtoVertexSet;
}

#endif
