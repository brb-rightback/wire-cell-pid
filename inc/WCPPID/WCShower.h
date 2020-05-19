#ifndef WCPPID_WCSHOWER_H
#define WCPPID_WCSHOWER_H

#include "WCPPID/Map_Proto_Vertex_Segment.h"
#include <vector>

namespace WCPPID{
  class WCShower{
  public:
    WCShower();
    ~WCShower();
    
    // define the initial condition
    void set_start_vertex(ProtoVertex* vertex, int type);
    void set_start_segment(ProtoSegment* seg);
    void set_start_point(WCP::Point p);
    void set_start_segment(ProtoSegment* seg, Map_Proto_Segment_Vertices& map_segment_vertices);
    
    std::pair<ProtoVertex*, int> get_start_vertex(){return std::make_pair(start_vertex, start_connection_type);};
    ProtoSegment* get_start_segment(){return start_segment;}
    
    std::pair<double, WCP::Point> get_closest_point(WCP::Point& p);
    double get_closest_dis(WCPPID::ProtoSegment* seg);
    
    // provide internal information to outside structure ...
    void fill_maps(std::map<WCPPID::ProtoVertex*, WCPPID::WCShower* >& map_vertex_in_shower, std::map<WCPPID::ProtoSegment*, WCPPID::WCShower*>& map_segment_in_shower);
    void fill_sets( std::set<WCPPID::ProtoVertex* >& used_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments, bool flag_exclude_start_segment = true);

    // complete structure ...
    void complete_structure_with_start_segment(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments);
    void add_segment(WCPPID::ProtoSegment* seg, Map_Proto_Segment_Vertices& map_segment_vertices);
    void add_shower(WCPPID::WCShower* temp_shower);
    double get_dis(WCPPID::ProtoSegment* seg);

    void calculate_kinematics();
    void calculate_kinematics_long_muon(std::set<WCPPID::ProtoSegment*> segments_in_muons);
    std::pair<WCPPID::ProtoSegment*, WCPPID::ProtoVertex*> get_last_segment_vertex_long_muon(std::set<WCPPID::ProtoSegment*> segments_in_muons);
    
    TVector3 cal_dir_3vector(WCP::Point p, double dis_cut = 15*units::cm);
    
    void set_flag_kinematics(bool val){flag_kinematics = val;};
    bool get_flag_kinematics(){return flag_kinematics;};
    int get_particle_type(){return particle_type;};
    void set_particle_type(int val){particle_type = val;};
    bool get_flag_shower(){return flag_shower;};
    double get_kine_range(){return kenergy_range;};
    double get_kine_dQdx(){return kenergy_dQdx;};
    void set_kine_charge(double val){kenergy_charge = val;};
    double get_kine_charge(){return kenergy_charge;};
    double get_kine_best(){return kenergy_best;};
    WCP::Point& get_start_point(){return start_point;};
    WCP::Point& get_end_point(){return end_point;};
    TVector3& get_init_dir(){return init_dir;};

 

    double get_total_length();
    double get_total_length(int tmp_cluster_id);
    double get_total_track_length();
    
    int get_num_segments(){return  map_seg_vtxs.size();};
    void update_particle_type();
    
    void rebuild_point_clouds();
    void build_point_clouds();
    WCP::ToyPointCloud* get_fit_pcloud(){return pcloud_fit;};
    WCP::ToyPointCloud* get_associated_pcloud(){return pcloud_associated;};

    //
    std::pair<std::set<WCPPID::ProtoSegment*>, std::set<WCPPID::ProtoVertex*> > get_connected_pieces(WCPPID::ProtoSegment* seg);
    
    Map_Proto_Segment_Vertices& get_map_seg_vtxs(){return map_seg_vtxs; };
    Map_Proto_Vertex_Segments& get_map_vtx_segs(){return map_vtx_segs; };
  protected:
    int particle_type;
    bool flag_shower;
    bool flag_kinematics;
    
    double kenergy_range;
    double kenergy_dQdx;
    double kenergy_charge;
    double kenergy_best;
    
    WCP::Point start_point;
    WCP::Point end_point;
    TVector3 init_dir;

    WCP::ToyPointCloud* pcloud_fit;
    WCP::ToyPointCloud* pcloud_associated;
    
    
    ProtoVertex* start_vertex;
    int start_connection_type;
    // 1 for direct connection, 2 for indirection connection with a gap, 3 for associations (not clear if this should be connected or not

    ProtoSegment* start_segment;

    Map_Proto_Vertex_Segments map_vtx_segs;
    Map_Proto_Segment_Vertices map_seg_vtxs; 
    
  };
  
  typedef std::vector<WCShower*> WCShowerSelection;
}

#endif
