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

    bool get_fit_flag(){return flag_fit;};
    void set_fit_flag(bool flag){flag_fit = flag;};

    // after trajectory fit ...
    void set_fit_pt(WCP::Point& tmp_fit_pt);
    WCP::Point& get_fit_pt() {return fit_pt;};

    // after dQ/dx fit ...
    void set_dQ(double val){dQ = val;};
    double get_dQ(){return dQ;};
    void set_dx(double val){dx = val;};
    double get_dx(){return dx;};

    double get_dQ_dx(){return dQ/(dx+1e-9);};

    
    
  protected:
    WCP::WCPointCloud<double>::WCPoint wcpt; // initial WCP point from the graph ...
    WCP::Point fit_pt; // best fit points ...

    double dQ;
    double dx;

    bool flag_fit;
  };
  typedef std::vector<ProtoVertex*> ProtoVertexSelection;
}

#endif
