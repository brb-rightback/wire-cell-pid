#include "WCPPID/WCParticle.h"

using namespace WCP;

WCPPID::WCParticle::WCParticle()
  : particle_type(0)
  , start_vertex(0)
  , start_connection_type(0)
  , start_segment(0)
  , end_vertex(0)
  , end_connection_type(0)
{
}

void WCPPID::WCParticle::add_contained_vertices(WCPPID::ProtoVertex *vertex){
  contained_vertices.push_back(vertex);
}

void WCPPID::WCParticle::add_contained_segments(WCPPID::ProtoSegment *segment){
  contained_segments.push_back(segment);
}

WCPPID::WCParticle::~WCParticle(){
  
}
