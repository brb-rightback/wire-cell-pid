#ifndef WIRECELLPID_PROTECTOVERCLUSTERING_H
#define WIRECELLPID_PROTECTOVERCLUSTERING_H

#include "WCPPID/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCPData/Opflash.h"
#include "WCPData/PhotonLibrary.h"

namespace WCPPID{
  void Protect_Over_Clustering(double eventTime, std::vector<std::pair<int, WCP::Opflash*> >& to_be_checked,
			       WCPPID::PR3DClusterSelection& live_clusters,
			       std::map<WCPPID::PR3DCluster*, int>& map_cluster_parent_id, // cluster to main cluster
			       std::map<int, std::vector<WCPPID::PR3DCluster*> >& map_parentid_clusters, // main cluster to clusters
			       WCP::ToyCTPointCloud& ct_point_cloud, bool flag_match_data, int run_no, double time_offset, double nrebin, double time_slice_width, bool flag_timestamp
			       );
  int get_next_cluster_id(int acc_cluster_id, std::set<int>& used_cluster_ids);

  std::pair<double, double> compare_pe_pattern(double eventTime, int run_no, double offset_x, WCP::Photon_Library *pl, WCP::SMGCSelection& mcells, WCP::Opflash* flash, bool flag_match_data, bool flag_timestamp = false);

  int convert_xyz_voxel_id(WCP::Point& p);
}

#endif 
