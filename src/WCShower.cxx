#include "WCPPID/WCShower.h"

using namespace WCP;

WCPPID::WCShower::WCShower()
  : particle_type(0)
  , start_vertex(0)
  , start_connection_type(0)
  , start_segment(0)
  , end_vertex(0)
  , end_connection_type(0)
{
}

void WCPPID::WCShower::add_contained_vertices(WCPPID::ProtoVertex *vertex){
  contained_vertices.push_back(vertex);
}

void WCPPID::WCShower::add_contained_segments(WCPPID::ProtoSegment *segment){
  contained_segments.push_back(segment);
}

WCPPID::WCShower::~WCShower(){
  
}
