#include "WireCellPID/PR3DCluster.h"

using namespace WireCell;

WireCellPID::PR3DCluster::PR3DCluster(int cluster_id)
  : cluster_id(cluster_id)
{
  // point_cloud = 0;
  // graph = 0;
  // source_wcp_index = -1;
  // flag_fine_tracking = false;
  // flag_PCA = false;
}

WireCellPID::PR3DCluster::~PR3DCluster(){
  // if (point_cloud!=(ToyPointCloud*)0)
  //   delete point_cloud;
  // if (graph!=(MCUGraph*)0)
  //   delete graph;
}

void WireCellPID::PR3DCluster::AddCell(SlimMergeGeomCell* mcell, int time_slice){
  if (cell_times_set_map.find(mcell)==cell_times_set_map.end()){
    std::set<int> times;
    times.insert(time_slice);
    cell_times_set_map[mcell]=times;
    mcells.push_back(mcell);
  }else{
    std::set<int>& times = cell_times_set_map[mcell];
    //if (find(times.begin(),times.end(),time_slice)==times.end()){
    times.insert(time_slice);
    // }
  }
  
  if (time_cells_set_map.find(time_slice)==time_cells_set_map.end()){
    SMGCSet mcells_1;
    mcells_1.insert(mcell);
    time_cells_set_map[time_slice] = mcells_1;
  }else{
    SMGCSet& mcells_1 = time_cells_set_map[time_slice];
    //if (find(mcells_1.begin(),mcells_1.end(), mcell) == mcells_1.end()){
    mcells_1.insert(mcell);
    //}
  }
}
