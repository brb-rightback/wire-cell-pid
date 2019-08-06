#ifndef WIRECELLPID_IMPROVEPR3DCluster_H
#define WIRECELLPID_IMPROVEPR3DCluster_H

#include "WireCellPID/PR3DCluster.h"
#include "WireCellData/ToyCTPointCloud.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCellPID{
  WireCellPID::PR3DCluster* Improve_PR3DCluster(WireCellPID::PR3DCluster* cluster, WireCell::ToyCTPointCloud& ct_point_cloud,WireCellSst::GeomDataSource& gds);
  
  WireCellPID::PR3DCluster* Improve_PR3DCluster(WireCellPID::PR3DCluster* orig_cluster, WireCellPID::PR3DCluster* new_cluster, WireCell::ToyCTPointCloud& ct_point_cloud,WireCellSst::GeomDataSource& gds);

  WireCellPID::PR3DCluster* Improve_PR3DCluster_1(WireCellPID::PR3DCluster* cluster, WireCell::ToyCTPointCloud& ct_point_cloud,WireCellSst::GeomDataSource& gds);

  WireCellPID::PR3DCluster* Improve_PR3DCluster_2(WireCellPID::PR3DCluster* cluster, WireCell::ToyCTPointCloud& ct_point_cloud,WireCellSst::GeomDataSource& gds, int nrebin, int frame_length, double unit_dis);
  
}

#endif 
