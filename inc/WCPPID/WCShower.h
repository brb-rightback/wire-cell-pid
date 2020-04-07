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
    
    std::pair<ProtoVertex*, int> get_start_vertex(){return std::make_pair(start_vertex, start_connection_type);};
    ProtoSegment* get_start_segment(){return start_segment;}
    bool get_flag_shower(){return flag_shower;};
    
    // provide internal information to outside structure ...
    void fill_maps(std::map<WCPPID::ProtoVertex*, WCPPID::WCShower* >& map_vertex_in_shower, std::map<WCPPID::ProtoSegment*, WCPPID::WCShower*>& map_segment_in_shower);
    void fill_sets( std::set<WCPPID::ProtoVertex* >& used_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments, bool flag_exclude_start_segment = true);

    // complete structure ...
    void complete_structure_with_start_segment(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments);

    
    
  protected:
    int particle_type;
    bool flag_shower;

    double kenergy_range;
    double kenergy_dQdx;
    double kenergy_charge;
    WCP::Point start_point;
    WCP::Point end_point;
    TVector3 init_dir;
    
    
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
