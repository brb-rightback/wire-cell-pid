#ifndef WIRECELLPID_NEUTRINOID_H
#define WIRECELLPID_NEUTRINOID_H

#include "WCPSst/GeomDataSource.h"

#include "WCPPID/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"

#include "WCPPID/WCVertex.h"
#include "WCPPID/WCParticle.h"

#include "WCPPID/ProtoVertex.h"
#include "WCPPID/ProtoSegment.h"


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
    
    // get segments
    int get_num_segments(ProtoVertex *pv);
    
    // actual functions ...
    void process_main_cluster();
    void process_other_clusters();

    
    // proto-vertex finder
    void find_proto_vertex(WCPPID::PR3DCluster *cluster);

    
    
  protected:
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
    std::map<ProtoVertex*, ProtoSegmentSet> map_vertex_segments;
    std::map<ProtoSegment*, ProtoVertexSet> map_segment_vertices;

    // map the cluster to the vertices/segments
    std::map<PR3DCluster*, ProtoVertexSet> map_cluster_vertices;
    std::map<ProtoVertex*, PR3DCluster*> map_vertex_cluster;
    std::map<PR3DCluster*, ProtoSegmentSet> map_cluster_segments;
    std::map<ProtoSegment*, PR3DCluster*> map_segment_cluster;


    
    // after fit, for alter direction
    WCVertexSelection vertices;
    WCParticleSelection particles;
    
  };
  
}

#endif
