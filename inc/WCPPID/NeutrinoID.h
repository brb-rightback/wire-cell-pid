#ifndef WIRECELLPID_NEUTRINOID_H
#define WIRECELLPID_NEUTRINOID_H

#include "WCPSst/GeomDataSource.h"

#include "WCPPID/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"

#include "WCPPID/WCVertex.h"
#include "WCPPID/WCParticle.h"

#include "WCPPID/Map_Proto_Vertex_Segment.h"

//#include "Minuit2/FCNBase.h"

namespace WCPPID{
  struct WCRecoTree
  {
    int mc_Ntrack;  // number of tracks in MC
    int mc_id[1000];  // track id; size == mc_Ntrack
    int mc_pdg[1000];  // track particle pdg; size == mc_Ntrack
    int mc_process[1000];  // track generation process code; size == mc_Ntrack
    int mc_mother[1000];  // mother id of this track; size == mc_Ntrack
    float mc_startXYZT[1000][4];  // start position of this track; size == mc_Ntrack
    float mc_endXYZT[1000][4];  // end position of this track; size == mc_Ntrack
    float mc_startMomentum[1000][4];  // start momentum of this track; size == mc_Ntrack
    float mc_endMomentum[1000][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > *mc_daughters;
  };

  
  //  class MyFCN : public ROOT::Minuit2::FCNBase {
  class MyFCN {
  
  public: 
    //    double Up() const { return 1.;}

    MyFCN(ProtoVertex* vtx, bool flag_vtx_constraint = false, double vtx_constraint_range = 1*units::cm, double vertex_protect_dis = 1.5*units::cm, double vertex_protect_dis_short_track = 0.9*units::cm, double fit_dis = 6 * units::cm);    
    ~MyFCN();

    void update_fit_range(double tmp_vertex_protect_dis = 1.5*units::cm, double  tmp_vertex_protect_dis_short_track = 0.9*units::cm, double tmp_fit_dis = 6 * units::cm);
    void AddSegment(ProtoSegment *sg);
    std::pair<bool, WCP::Point> FitVertex();
    void UpdateInfo(WCP::Point fit_pos, WCPPID::PR3DCluster* temp_cluster, double default_dis_cut = 4.0*units::cm);
    
    std::pair<ProtoSegment*, int> get_seg_info(int i);
    int get_fittable_tracks();
    bool get_flag_vtx_constraint(){return flag_vtx_constraint;};
    void set_flag_vtx_constraint(bool val){flag_vtx_constraint = val;};

    void set_vtx_constraint_range(double val){vtx_constraint_range = val;};

    std::vector<ProtoSegment*>& get_segments(){return segments;};
    std::vector<WCP::PointVector>& get_vec_points(){return vec_points;};

    void print_points();
    
    //  double operator() (const std::vector<double> & xx) const;
    //  double get_chi2(const std::vector<double> & xx) const;
    
  private:
    ProtoVertex *vtx;
    bool flag_vtx_constraint;
    double vtx_constraint_range;
    
    double vertex_protect_dis;
    double vertex_protect_dis_short_track;
    double fit_dis;
    
    std::vector<ProtoSegment* > segments;
    std::vector<WCP::PointVector> vec_points;

    std::vector<std::tuple<WCP::Point, WCP::Point, WCP::Point> > vec_PCA_dirs;
    std::vector<std::tuple<double, double, double> > vec_PCA_vals;
    std::vector<WCP::Point> vec_centers;
  };
  
  
  class NeutrinoID{
  public:
    NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, std::vector<WCPPID::PR3DCluster*>& all_clusters, WCPSst::GeomDataSource& gds, int nrebin, int frame_length, float unit_dis,	WCP::ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, double flash_time);
    ~NeutrinoID();

    // deal with the map ...
    bool del_proto_vertex(ProtoVertex *pv);
    bool del_proto_segment(ProtoSegment *ps);
    bool add_proto_connection(ProtoVertex *pv, ProtoSegment *ps, WCPPID::PR3DCluster* cluster);
    void organize_vertices_segments();
    
    std::tuple<ProtoVertex*, ProtoSegment*, WCP::Point> check_end_point(WCP::PointVector& tracking_path, bool flag_front = true, double vtx_cut1 = 0.9*units::cm, double vtx_cut2 = 2.0*units::cm, double sg_cut1 = 2.0 * units::cm, double sg_cut2 = 1.2*units::cm);
    
    // get segments
    int get_num_segments(ProtoVertex *pv);
    
    // actual functions ...
    void process_main_cluster();
    void process_other_clusters();

    // fill reco information
    void fill_reco_simple_tree(WCRecoTree& rtree);
    
    
    // proto-vertex finder
    bool find_proto_vertex(WCPPID::PR3DCluster *cluster, bool flag_break_trak = true, int nrounds_find_other_tracks = 2);
    WCPPID::ProtoSegment* init_first_segment(WCPPID::PR3DCluster *cluster);
    void init_point_segment(WCPPID::PR3DCluster *cluster);
    void break_segments(std::vector<WCPPID::ProtoSegment*>& remaining_segments, WCPPID::PR3DCluster* temp_cluster);
    void examine_vertices(WCPPID::PR3DCluster* temp_cluster);
    void examine_vertices_1(WCPPID::PR3DCluster* temp_cluster);
    void examine_vertices_2(WCPPID::PR3DCluster* temp_cluster);
    bool examine_vertices(WCPPID::ProtoVertex* v1, WCPPID::ProtoVertex *v2, double offset_t, double slope_xt, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw);
    
    void find_other_segments(WCPPID::PR3DCluster* temp_cluster, bool flag_break_track = true, double search_range = 1.5*units::cm, double scaling_2d = 0.8);
    ProtoVertex* find_vertex_other_segment(WCPPID::PR3DCluster *temp_cluster, bool flag_forward, WCP::WCPointCloud<double>::WCPoint& wcp);


    // improve vertex ...
    void improve_vertex(WCPPID::PR3DCluster* temp_cluster);
    bool fit_vertex(WCPPID::ProtoVertex *vtx, WCPPID::ProtoSegmentSet& sg_set, WCPPID::PR3DCluster* temp_cluster);
    bool search_for_vertex_activities(WCPPID::ProtoVertex *vtx, WCPPID::ProtoSegmentSet& sg_set, WCPPID::PR3DCluster* temp_cluster);

    // get direction 
    TVector3 get_dir(WCPPID::ProtoVertex *vtx, WCPPID::ProtoSegment *sg, double dis = 2*units::cm);
    void determine_direction(WCPPID::PR3DCluster* temp_cluster);
    
    // clustering points
    void clustering_points(WCPPID::PR3DCluster* temp_cluster);

    WCPPID::ProtoSegment* find_segment(WCPPID::ProtoVertex *v1, WCPPID::ProtoVertex *v2);
    std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*> find_vertices(WCPPID::ProtoSegment* sg);
    WCPPID::ProtoVertex* find_other_vertex(WCPPID::ProtoSegment *sg, WCPPID::ProtoVertex* v1); 
    WCPPID::ProtoVertexSelection find_vertices(WCPPID::PR3DCluster* temp_cluster);
    WCPPID::ProtoSegmentSelection find_segmentss(WCPPID::PR3DCluster* temp_cluster);
    
    Map_Proto_Vertex_Segments& get_map_vertex_segments(){return map_vertex_segments;};
    Map_Proto_Segment_Vertices& get_map_segment_verteices(){return map_segment_vertices;};

    // deghost
    void deghost_clusters();
    void order_clusters(WCPPID::PR3DClusterSelection& ordered_clusters, std::map<int, WCPPID::ProtoSegmentSelection>& map_cluster_id_segments, std::map<WCPPID::PR3DCluster*, double>& map_cluster_total_length);
    
    // track shower separation
    void separate_track_shower();
    
    // fill_fit_parameters();
    void fill_fit_parameters();
    
  protected:
    int acc_vertex_id;
    int acc_segment_id;
    // input ...
    WCPPID::PR3DCluster *main_cluster;
    std::vector<WCPPID::PR3DCluster*> other_clusters;
    std::vector<WCPPID::PR3DCluster*> all_clusters;

    
    std::map<int, WCPPID::PR3DCluster*> map_id_cluster;
    
    WCP::ToyCTPointCloud* ct_point_cloud;
    std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > > global_wc_map;
    double flash_time;

    // output ...
    int type; // nue, numu, NC for 1,2,3,    0 for no ID

    // graph ...
    ProtoVertexSelection proto_vertices;
    ProtoSegmentSelection proto_segments;
    Map_Proto_Vertex_Segments map_vertex_segments;
    Map_Proto_Segment_Vertices map_segment_vertices;
    // map the cluster to the vertices/segments
    std::map<PR3DCluster*, ProtoVertexSet> map_cluster_vertices;
    std::map<ProtoVertex*, PR3DCluster*> map_vertex_cluster;
    std::map<PR3DCluster*, ProtoSegmentSet> map_cluster_segments;
    std::map<ProtoSegment*, PR3DCluster*> map_segment_cluster;

    std::vector<std::tuple<PR3DCluster*, int, int> > residual_segment_candidates;
    
    
    // after fit, for alter direction
    WCVertexSelection vertices;
    WCParticleSelection particles;
    
  };

  struct Res_proto_segment
  {
    int group_num;
    int number_points;
    int special_A, special_B;
    double length;
    int number_not_faked;
    double max_dis_u;
    double max_dis_v;
    double max_dis_w;
  };
  
}

#endif
