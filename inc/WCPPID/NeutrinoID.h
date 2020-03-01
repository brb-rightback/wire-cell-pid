#ifndef WIRECELLPID_NEUTRINOID_H
#define WIRECELLPID_NEUTRINOID_H

#include "WCPPID/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"

#include "WCPPID/WCVertex.h"
#include "WCPPID/WCParticle.h"


namespace WCPPID{
  class NeutrinoID{
  public:
    NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, WCP::ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, double flash_time);
    ~NeutrinoID();

    // actual functions ...
    void process_main_cluster();
    void process_other_clusters();
    
  protected:
    // input ...
    WCPPID::PR3DCluster *main_cluster;
    std::vector<WCPPID::PR3DCluster*> other_clusters;
    WCP::ToyCTPointCloud* ct_point_cloud;
    std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > > global_wc_map;
    double flash_time;

    // output ...
    int type; // nue, numu, NC for 1,2,3,    0 for no ID

    WCVertexSelection vertices;
    WCParticleSelection particles;
    
  };
  
}

#endif
