#include "WCPPID/WCShower.h"

using namespace WCP;

WCPPID::WCShower::WCShower()
  : particle_type(0)
  , start_vertex(0)
  , start_connection_type(0)
  , start_segment(0)
{
}

WCPPID::WCShower::~WCShower(){
  
}

void WCPPID::WCShower::set_start_vertex(ProtoVertex* vertex, int type){
  start_vertex = vertex;
  start_connection_type = type;
}

void WCPPID::WCShower::set_start_segment(ProtoSegment* seg){
  start_segment = seg;
}
