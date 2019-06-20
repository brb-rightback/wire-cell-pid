#ifndef WIRECELLPID_PR3DCLUSTER_H
#define WIRECELLPID_PR3DCLUSTER_H

#include "WireCellData/SlimMergeGeomCell.h"

namespace WireCellPID{
  class PR3DCluster{
  public:
    PR3DCluster(int cluster_id);
    ~PR3DCluster();
    void AddCell(WireCell::SlimMergeGeomCell* mcell, int time_slice);
    
    int get_cluster_id(){return cluster_id;};
    void set_cluster_id(int value){cluster_id=value;};
    int get_num_mcells(){return mcells.size();};
    WireCell::SMGCSelection& get_mcells(){return mcells;};
  protected:
    int cluster_id;
    
    WireCell::SMGCSelection mcells;
    std::map<int,WireCell::SMGCSet> time_cells_set_map;
    std::map<WireCell::SlimMergeGeomCell*, std::set<int>> cell_times_set_map;   
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}

#endif
