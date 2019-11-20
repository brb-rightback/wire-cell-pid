#ifndef WIRECELLPID_IMPROVEPR3DCluster_H
#define WIRECELLPID_IMPROVEPR3DCluster_H

#include "WCPPID/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCPSst/GeomDataSource.h"
#include "WCP2dToy/WCPHolder.h"

namespace WCPPID{
  WCPPID::PR3DCluster* Improve_PR3DCluster(WCPPID::PR3DCluster* cluster, WCP::ToyCTPointCloud& ct_point_cloud,WCPSst::GeomDataSource& gds, WCP2dToy::WCPHolder *holder);
  
  WCPPID::PR3DCluster* Improve_PR3DCluster(WCPPID::PR3DCluster* orig_cluster, WCPPID::PR3DCluster* new_cluster, WCP::ToyCTPointCloud& ct_point_cloud,WCPSst::GeomDataSource& gds, WCP2dToy::WCPHolder *holder);

  WCPPID::PR3DCluster* Improve_PR3DCluster_1(WCPPID::PR3DCluster* cluster, WCP::ToyCTPointCloud& ct_point_cloud,WCPSst::GeomDataSource& gds, WCP2dToy::WCPHolder *holder);

  // final
  WCPPID::PR3DCluster* Improve_PR3DCluster_2(WCPPID::PR3DCluster* cluster, WCP::ToyCTPointCloud& ct_point_cloud,WCPSst::GeomDataSource& gds, WCP2dToy::WCPHolder *holder, int nrebin, int frame_length, double unit_dis);
  
}

#endif 
