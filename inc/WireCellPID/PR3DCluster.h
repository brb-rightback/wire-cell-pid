#ifndef WIRECELLPID_PR3DCLUSTER_H
#define WIRECELLPID_PR3DCLUSTER_H

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/ToyPointCloud.h"
#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/ToyCTPointCloud.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <Eigen/Sparse>

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

    // projection related
    void get_projection(std::vector<int>& proj_channel, std::vector<int>& proj_timeslice, std::vector<int>& proj_charge, std::vector<int>& proj_charge_err , std::vector<int>& proj_flag, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map);
    void collect_charge_trajectory(WireCell::ToyCTPointCloud& ct_point_cloud, double dis_cut = 0.6*units::cm, double range_cut = 1.0*units::cm);


    //fine tracking related ...
    WireCell::PointVector& get_fine_tracking_path(){return fine_tracking_path;};
    std::vector<double>& get_dQ(){return dQ;};
    std::vector<double>& get_dx(){return dx;};
    std::vector<double>& get_pu(){return pu;};
    std::vector<double>& get_pv(){return pv;};
    std::vector<double>& get_pw(){return pw;};
    std::vector<double>& get_pt(){return pt;};
    std::vector<double>& get_reduced_chi2(){return reduced_chi2;};
    
    // main function to do the overall tracking, 
    void do_tracking(WireCell::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, double time = 4*units::microsecond, bool flag_dQ_dx_fit_reg = true);
    // organize path from the shortest path 
    WireCell::PointVector organize_wcps_path(std::list<WireCell::WCPointCloud<double>::WCPoint>& path_wcps_list,  double low_dis_limit, double end_point_limit);
    void organize_ps_path(WireCell::PointVector& ps_vec, double low_dis_limit, double end_point_limit);
    // prepare data
    void prepare_data(WireCell::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge);
    // form map ...
    void form_map(WireCell::ToyCTPointCloud& ct_point_cloud, WireCell::PointVector& pts,
		  std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge,
		  std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set,
		  std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,
		  double end_point_factor=0.6, double mid_point_factor=0.9, int nlevel = 3, double time_cut = 5, double charge_cut = 2000);

    void form_point_association(WireCell::Point &p, std::set<std::pair<int,int> >& temp_2dut, std::set<std::pair<int,int> >& temp_2dvt, std::set<std::pair<int,int> >& temp_2dwt, WireCell::ToyCTPointCloud& ct_point_cloud, double dis_cut, int nlevel = 3, double time_cut = 5);
    
    std::vector<float> examine_point_association(std::vector<int>& temp_results, std::set<std::pair<int,int> >& temp_2dut, std::set<std::pair<int,int> >& temp_2dvt, std::set<std::pair<int,int> >& temp_2dwt,
						 std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge, bool flag_end_point = false, double charge_cut = 2000);
    
    void trajectory_fit(WireCell::PointVector& ps_vec, std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,std::map<std::pair<int,int>, std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, int charge_div_method = 1, double div_sigma = 0.6*units::cm);

    bool skip_trajectory_point(WireCell::Point& p_fit, int i, WireCell::PointVector& ps_vec, std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<std::pair<int,int>, std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, WireCell::PointVector& fine_tracking_path,
				 double offset_t, double slope_x, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw);
    
    void fill_data_map_trajectory(std::vector<int> indices, std::map<int, std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, float> >& map_3D_2DW_set, std::map<std::pair<int,int>,  std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge);

    // dQ/dx fit part
    double cal_gaus_integral(int tbin, int wbin, double t_center, double t_sigma, double w_center, double w_sigma, int flag=0, double nsigma = 3);
    double cal_gaus_integral_seg(int tbin, int wbin, std::vector<double>& t_centers, std::vector<double>& t_sigmas, std::vector<double>& w_centers, std::vector<double>& w_sigmas, std::vector<double>& weights, int flag=0, double nsigma = 3);

    void dQ_dx_fit(std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,  std::tuple<double, double, int > >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_wt_charge, double flash_time = 4*units::microsecond, double dis_end_point_ext = 4.5*units::mm, bool flag_dQ_dx_fit_reg = true);
    void update_data_dQ_dx_fit(std::map<int,std::map<const WireCell::GeomWire*, WireCell::SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>,std::tuple<double, double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_wt_charge);

    std::vector<std::pair<double,double> > cal_compact_matrix(Eigen::SparseMatrix<double>& MW, Eigen::SparseMatrix<double>& RWT, int n_2D_w, int n_3D_pos, double cut_pos = 2);
    
    
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & get_proj_data_u_map(){return proj_data_u_map;};
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & get_proj_data_v_map(){return proj_data_v_map;};
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & get_proj_data_w_map(){return proj_data_w_map;};
    // crawl alg for stm
    void do_rough_path(WireCell::WCPointCloud<double>::WCPoint& first_wcp, WireCell::WCPointCloud<double>::WCPoint& last_wcp);
    WireCell::Point adjust_rough_path();
    
    WireCell::Point do_stm_crawl(WireCell::WCPointCloud<double>::WCPoint& first_wcp, WireCell::WCPointCloud<double>::WCPoint& last_wcp, int flag_end = 1);

    
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

    // projection related
    std::map<std::pair<int,int>,std::pair<int,int>> collected_charge_map;


     // fine tracking related ...
    WireCell::PointVector fine_tracking_path;
    std::vector<double> dQ;
    std::vector<double> dx;
    std::vector<double> pu;
    std::vector<double> pv;
    std::vector<double> pw;
    std::vector<double> pt;
    std::vector<double> reduced_chi2;

    std::map<std::pair<int,int>, std::tuple<double,double,double> > proj_data_u_map;
    std::map<std::pair<int,int>, std::tuple<double,double,double> > proj_data_v_map;
    std::map<std::pair<int,int>, std::tuple<double,double,double> > proj_data_w_map;
    
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}

#endif
