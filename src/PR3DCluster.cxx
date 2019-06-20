#include "WireCellPID/PR3DCluster.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

using namespace WireCell;

WireCellPID::PR3DCluster::PR3DCluster(int cluster_id)
  : cluster_id(cluster_id)
{
  point_cloud = 0;
  // graph = 0;
  // source_wcp_index = -1;
  // flag_fine_tracking = false;
  // flag_PCA = false;
}

WireCellPID::PR3DCluster::~PR3DCluster(){
  if (point_cloud!=(ToyPointCloud*)0)
    delete point_cloud;
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


void WireCellPID::PR3DCluster::Create_point_cloud(WireCell::ToyPointCloud *global_point_cloud){
  if (point_cloud!=(ToyPointCloud*)0)
    return;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  
  
  point_cloud = new ToyPointCloud(angle_u, angle_v, angle_w);
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    PointVector pts = mcell->get_sampling_points();

    if (global_point_cloud!=(ToyPointCloud*)0)
      global_point_cloud->AddPoints(pts,mcell->get_sampling_points_wires(),mcell);  
    
    point_cloud->AddPoints(pts,mcell->get_sampling_points_wires(),mcell);
  }
  point_cloud->build_kdtree_index();
  //  std::cout << point_cloud->get_num_points() << std::endl;

}
