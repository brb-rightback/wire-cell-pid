#ifndef WCPPID_MAP_PROTO_VERTEX_SEGMENT_H
#define WCPPID_MAP_PROTO_VERTEX_SEGMENT_H

#include "WCPPID/ProtoVertex.h"
#include "WCPPID/ProtoSegment.h"

namespace WCPPID{
  typedef std::map<ProtoSegment*, ProtoVertexSet, ProtoSegmentCompare> Map_Proto_Segment_Vertices;
  typedef std::map<ProtoVertex*, ProtoSegmentSet, ProtoVertexCompare> Map_Proto_Vertex_Segments;
}

#endif
