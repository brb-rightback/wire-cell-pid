#ifndef WIRECELLPID_PR3DCLUSTER_H
#define WIRECELLPID_PR3DCLUSTER_H

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/ToyPointCloud.h"
#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/ToyCTPointCloud.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

namespace WireCellPID{

  struct VertexProp {
    int index;
    // WCPointCloud<double>::WCPoint wcpoint;
    // add pointer to merged cell
  };
  /* struct EdgeProp { */
  /*   float dist; // edge distance */
  /* }; */

  using EdgeProp = boost::property<boost::edge_weight_t, float>; 
  // using VertexProp = boost::property(boost::vertex_color_t, int>; 
  
  typedef adjacency_list<setS, vecS, undirectedS, VertexProp, EdgeProp> MCUGraph;
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
    void Del_point_cloud();
    WireCell::ToyPointCloud* get_point_cloud(){return point_cloud;};
    WireCell::ToyPointCloud* get_point_cloud_steiner(){return point_cloud_steiner;};
    WireCell::ToyPointCloud* get_point_cloud_steiner_terminal(){return point_cloud_steiner_terminal;};
    
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
    void Create_graph(WireCell::ToyPointCloud* ref_point_cloud = 0);
    void Create_graph(WireCell::ToyCTPointCloud& ct_point_cloud, WireCell::ToyPointCloud* ref_point_cloud = 0);
    void Establish_close_connected_graph();
    void Connect_graph(WireCell::ToyPointCloud* ref_point_cloud = 0);
    void Connect_graph(WireCell::ToyCTPointCloud& ct_point_cloud, WireCell::ToyPointCloud* ref_point_cloud = 0);
    void Del_graph();

    // Steiner tree
    void create_steiner_graph(WireCell::ToyCTPointCloud& ct_point_cloud, WireCellSst::GeomDataSource& gds, int nrebin, int frame_length, double unit_dis);
    MCUGraph* Create_steiner_tree(WireCell::ToyPointCloud *point_cloud_steiner, std::vector<bool>& flag_steiner_terminal, WireCell::GeomDataSource& gds, WireCell::SMGCSelection& mcells, bool flag_path = false, bool disable_dead_mix_cell = true);
    // example for more sophisticated algorithm later ...
    // use the full steiner terminal graph ...
    void recover_steiner_graph();

    void establish_same_mcell_steiner_edges(WireCell::GeomDataSource& gds, bool disable_dead_mix_cell = true, int flag=1);
    void remove_same_mcell_steiner_edges(int flag=1);
    
    
    void find_steiner_terminals(WireCell::GeomDataSource& gds, bool disable_dead_mix_cell = true);
    void form_cell_points_map();
    std::set<int> get_steiner_terminals(){return steiner_terminal_indices;};
    std::set<int> get_selected_terminals(){return selected_terminal_indices;};

    std::set<int> get_steiner_graph_terminals(){return steiner_graph_terminal_indices;};
    std::set<int> get_steiner_graph_selected_terminals(){return steiner_graph_selected_terminal_indices;};
    
    // find peak points within the mcells ...
    std::set<int> find_peak_point_indices(WireCell::SMGCSelection mcells,WireCell::GeomDataSource& gds, bool disable_dead_mix_cell = true, int nlevel = 1);

    std::pair<bool,double> calc_charge_wcp(WireCell::WCPointCloud<double>::WCPoint& wcp, WireCell::GeomDataSource& gds, bool disable_dead_mix_cell = true, double charge_cut = 4000);


    // path related
     // create things for Dijkstra
    std::vector<vertex_descriptor> parents;
    std::vector<int> distances;
    int source_wcp_index;
    int dest_wcp_index;
    
    std::list<WireCell::WCPointCloud<double>::WCPoint>& get_path_wcps(){return path_wcps;};
    std::list<WireCell::SlimMergeGeomCell*>& get_path_mcells(){return path_mcells;};

    // get exterme points ...
    std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> get_highest_lowest_wcps(int flag = 1);
    std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> get_front_back_wcps(int flag = 1);
    std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> get_earliest_latest_wcps(int flag = 1);
    std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> get_main_axis_wcps(int flag = 1);

    std::vector<std::vector<WireCell::WCPointCloud<double>::WCPoint>> get_extreme_wcps(int flag = 1, std::map<int,WireCell::SMGCSelection>* old_time_mcells_map=0);
    std::pair<WireCell::WCPointCloud<double>::WCPoint,WireCell::WCPointCloud<double>::WCPoint> get_two_boundary_wcps(int flag = 1);
    
    std::pair<WireCell::Point,WireCell::Point> get_two_extreme_points(int flag = 1);


    
    void dijkstra_shortest_paths(WireCell::WCPointCloud<double>::WCPoint& wcp_source, int flag = 1);
    void cal_shortest_path(WireCell::WCPointCloud<double>::WCPoint& wcp_target, int flag = 1);


    
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


    std::vector<edge_descriptor> same_mcell_steiner_edges;
    std::set<int> excluded_points;
    // saved steiner tree related products
    
    
    // prepare for steiner tree
    std::map<WireCell::SlimMergeGeomCell*, std::set<int>> cell_point_indices_map;
    std::set<int> steiner_terminal_indices;
    std::set<int> selected_terminal_indices;

    // data product for Steinter Tree Results
    MCUGraph *graph_steiner;
    WireCell::ToyPointCloud *point_cloud_steiner;
    WireCell::ToyPointCloud *point_cloud_steiner_terminal;
    std::vector<bool> flag_steiner_terminal;
    // more derived quantities to come ...
    std::set<int> steiner_graph_terminal_indices;
    std::set<int> steiner_graph_selected_terminal_indices;
    
    //path related
    std::list<WireCell::WCPointCloud<double>::WCPoint> path_wcps;
    std::list<WireCell::SlimMergeGeomCell*> path_mcells;
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}

#endif
