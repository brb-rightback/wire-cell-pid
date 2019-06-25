#ifndef WIRECELLPID_PR3DCLUSTER_H
#define WIRECELLPID_PR3DCLUSTER_H

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/ToyPointCloud.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

namespace WireCellPID{

  struct VertexProp {
    int index;
    //WCPointCloud<double>::WCPoint wcpoint;
    // add pointer to merged cell
  };
  struct EdgeProp {
    float dist; // edge distance
  };
  
  typedef adjacency_list<vecS, vecS, undirectedS, VertexProp, EdgeProp> MCUGraph;
  typedef graph_traits<MCUGraph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<MCUGraph>::edge_descriptor edge_descriptor;


  
  class PR3DCluster{
  public:
    PR3DCluster(int cluster_id);
    ~PR3DCluster();
    void AddCell(WireCell::SlimMergeGeomCell* mcell, int time_slice);
    
    int get_cluster_id(){return cluster_id;};
    void set_cluster_id(int value){cluster_id=value;};
    int get_num_mcells(){return mcells.size();};
    int get_num_points(){return point_cloud->get_num_points(); };
    int get_num_time_slices(){return time_cells_set_map.size();};
    WireCell::SMGCSelection& get_mcells(){return mcells;};

    // Point cloud related ...
    void Create_point_cloud(WireCell::ToyPointCloud *global_point_cloud = 0);
    WireCell::ToyPointCloud* get_point_cloud(){return point_cloud;};

    WireCell::Point calc_ave_pos(WireCell::Point& p, double dis);
    WireCell::Point calc_ave_pos(WireCell::Point& p, int N);
    int get_num_points(WireCell::Point& p_test, double dis);


    //PCA related ...
    void Calc_PCA();
    bool flag_PCA;
    void Calc_PCA(WireCell::PointVector& points);
    WireCell::Vector get_center(){return center;};
    WireCell::Vector get_PCA_axis(int axis){return PCA_axis[axis];};
    double get_PCA_value(int axis){return PCA_values[axis];};

    //Hough Transformation ...
    std::pair<double,double> HoughTrans(WireCell::Point& p, double dis);
    TVector3 VHoughTrans(WireCell::Point& p, double dis);
    std::pair<double,double> HoughTrans(WireCell::Point& p, double dis, WireCell::ToyPointCloud *point_cloud1);
    TVector3 VHoughTrans(WireCell::Point& p, double dis, WireCell::ToyPointCloud *point_cloud1);

    // graph related ...
    void Create_graph();
    //void Create_graph(WireCell::ToyCTPointCloud& ct_point_cloud);
    void Establish_close_connected_graph();
    void Connect_graph();
    //void Connect_graph(WireCell::ToyCTPointCloud& ct_point_cloud);
    void Del_graph();

    
    
  protected:
    int cluster_id;
    
    WireCell::SMGCSelection mcells;
    std::map<int,WireCell::SMGCSet> time_cells_set_map;
    std::map<WireCell::SlimMergeGeomCell*, std::set<int>> cell_times_set_map;

    WireCell::ToyPointCloud *point_cloud;
    WireCell::Vector center;
    
    WireCell::Vector PCA_axis[3];
    double PCA_values[3];

    // graph
    MCUGraph *graph;
    
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}

#endif
