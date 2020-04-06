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

void WCPPID::WCShower::fill_maps(std::map<WCPPID::ProtoVertex*, WCPPID::WCShower* >& map_vertex_in_shower, std::map<WCPPID::ProtoSegment*, WCPPID::WCShower*>& map_segment_in_shower){
  for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
    map_vertex_in_shower[it->first] = this;
  }
  for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
    map_segment_in_shower[it->first] = this;
  }
}





void WCPPID::WCShower::complete_structure_with_start_segment(Map_Proto_Vertex_Segments& map_vertex_segments, Map_Proto_Segment_Vertices& map_segment_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments){
  // fill the start segment ...
  std::vector<ProtoSegment* > new_segments;
  std::vector<ProtoVertex* > new_vertices;
  
  for (auto it = map_segment_vertices[start_segment].begin(); it!= map_segment_vertices[start_segment].end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (vtx == start_vertex) continue;
    map_vtx_segs[vtx].insert(start_segment);
    map_seg_vtxs[start_segment].insert(vtx);

    new_vertices.push_back(vtx);
  }

  while( new_vertices.size()>0 || new_segments.size()>0 ){
    if (new_vertices.size()>0){
      ProtoVertex *vtx = new_vertices.back();
      new_vertices.pop_back();
      for (auto it = map_vertex_segments[vtx].begin(); it != map_vertex_segments[vtx].end(); it++){
    	ProtoSegment *seg = *it;
    	if (used_segments.find(seg)!=used_segments.end()) continue;
    	map_vtx_segs[vtx].insert(seg);
    	map_seg_vtxs[seg].insert(vtx);
    	new_segments.push_back(seg);
	used_segments.insert(seg);
      }
    }

    if (new_segments.size()>0){
      ProtoSegment *seg = new_segments.back();
      new_segments.pop_back();
      for (auto it = map_segment_vertices[seg].begin(); it!= map_segment_vertices[seg].end(); it++){
  	ProtoVertex *vtx = *it;
  	if (map_vtx_segs.find(vtx)!= map_vtx_segs.end() || vtx == start_vertex) continue;
  	map_vtx_segs[vtx].insert(seg);
  	map_seg_vtxs[seg].insert(vtx);
  	new_vertices.push_back(vtx);
      }
    }
  }
  
  //std::cout << map_vtx_segs.size() << " " << map_seg_vtxs.size() << std::endl;
}
