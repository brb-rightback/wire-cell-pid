void WCPPID::NeutrinoID::separate_track_shower(){
  for (auto it = map_segment_vertices.begin(); it != map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    sg->is_shower_trajectory();
  }
}
