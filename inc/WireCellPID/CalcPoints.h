#ifndef WIRECELLPID_CALC_POINTS_H
#define WIRECELLPID_CALC_POINTS_H

#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/SlimMergeGeomCell.h"

namespace WireCellPID{
  void calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell, int nrebin, int frame_length, double unit_dis);
}

#endif 
