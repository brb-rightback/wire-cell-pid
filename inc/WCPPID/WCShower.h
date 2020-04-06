#ifndef WCPPID_WCSHOWER_H
#define WCPPID_WCSHOWER_H

#include "WCPPID/ProtoVertex.h"
#include "WCPPID/ProtoSegment.h"
#include <vector>

namespace WCPPID{
  class WCShower{
  public:
    WCShower();
    ~WCShower();
    void set_start_vertex(ProtoVertex* vertex, int type);
    void set_end_vertex(ProtoVertex* vertex, int type);
    void set_start_segment(ProtoSegment* seg);
    
    std::pair<ProtoVertex*, int> get_start_vertex(){return std::make_pair(start_vertex, start_connection_type);};
    std::pair<ProtoVertex*, int> get_end_vertex(){return std::make_pair(end_vertex, end_connection_type);};
    ProtoSegment* get_start_segment(){return start_segment;}
    
    void add_contained_vertices(ProtoVertex *vertex);
    void add_contained_segments(ProtoSegment *segment);

    ProtoVertexSelection& get_contained_vertices(){return contained_vertices;};
    ProtoSegmentSelection& get_contained_segments(){return contained_segments;};
    
  protected:
    int particle_type;
    
    ProtoVertexSelection contained_vertices;
    ProtoSegmentSelection contained_segments;
    
    ProtoVertex* start_vertex;
    int start_connection_type;
    // 1 for direct connection, 2 for indirection connection with a gap, 3 for associations (not clear if this should be connected or not

    ProtoSegment* start_segment;
    
    ProtoVertex* end_vertex;
    int end_connection_type;
  };
  
  typedef std::vector<WCShower*> WCShowerSelection;
}

#endif
