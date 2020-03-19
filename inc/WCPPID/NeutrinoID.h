#ifndef WIRECELLPID_NEUTRINOID_H
#define WIRECELLPID_NEUTRINOID_H

#include "WCPSst/GeomDataSource.h"

#include "WCPPID/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"

#include "WCPPID/WCVertex.h"
#include "WCPPID/WCParticle.h"

#include "WCPPID/Map_Proto_Vertex_Segment.h"

namespace WCPPID{
  class NeutrinoID{
  public:
    NeutrinoID(WCPPID::PR3DCluster *main_cluster, std::vector<WCPPID::PR3DCluster*>& other_clusters, WCPSst::GeomDataSource& gds, int nrebin, int frame_length, float unit_dis,	WCP::ToyCTPointCloud* ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, double flash_time);
    ~NeutrinoID();

    // deal with the map ...
    bool del_proto_vertex(ProtoVertex *pv);
    bool del_proto_segment(ProtoSegment *ps);
    bool add_proto_connection(ProtoVertex *pv, ProtoSegment *ps, WCPPID::PR3DCluster* cluster);
    void organize_vertices_segments();
    
    std::tuple<ProtoVertex*, ProtoSegment*, WCP::Point> check_end_point(WCP::PointVector& tracking_path, bool flag_front = true);
    
    // get segments
    int get_num_segments(ProtoVertex *pv);
    
    // actual functions ...
    void process_main_cluster();
    void process_other_clusters();

    
    // proto-vertex finder
    void find_proto_vertex(WCPPID::PR3DCluster *cluster, bool flag_break_trak = true, int nrounds_find_other_tracks = 2);
    WCPPID::ProtoSegment* init_first_segment(WCPPID::PR3DCluster *cluster);
    void break_segments(std::vector<WCPPID::ProtoSegment*>& remaining_segments, WCPPID::PR3DCluster* temp_cluster);

    void find_other_segments(WCPPID::PR3DCluster* temp_cluster, bool flag_break_track = true, double search_range = 1.5*units::cm, double scaling_2d = 0.8);
    ProtoVertex* find_vertex_other_segment(WCPPID::PR3DCluster *temp_cluster, bool flag_forward, WCP::WCPointCloud<double>::WCPoint& wcp);
    
    Map_Proto_Vertex_Segments& get_map_vertex_segments(){return map_vertex_segments;};
    Map_Proto_Segment_Vertices& get_map_segment_verteices(){return map_segment_vertices;};
    
  protected:
    int acc_vertex_id;
    int acc_segment_id;
    // input ...
    WCPPID::PR3DCluster *main_cluster;
    std::vector<WCPPID::PR3DCluster*> other_clusters;
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
