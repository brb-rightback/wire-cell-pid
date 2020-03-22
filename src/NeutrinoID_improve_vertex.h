void WCPPID::NeutrinoID::improve_vertex(WCPPID::PR3DCluster* temp_cluster){
  std::map<int, WCPPID::ProtoSegment*> map_id_seg;
  std::map<WCPPID::ProtoSegment*, int> map_seg_id;
  
  for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue;

    map_id_seg[sg->get_id()] = sg;
    map_seg_id[sg] = sg->get_id();
    sg->clear_associate_points();
  }
  
  // find the relevant point clouds ...
  WCP::WCPointCloud<double>& cloud = temp_cluster->get_point_cloud()->get_cloud();
  std::vector<int>& point_sub_cluster_ids = temp_cluster->get_point_sub_cluster_ids();
  for (size_t i=0;i!=point_sub_cluster_ids.size();i++){
    if (point_sub_cluster_ids.at(i) == -1) continue;
    WCP::Point p(cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z);
    map_id_seg[point_sub_cluster_ids.at(i)]->add_associate_point(p);
  }

   /* for (auto it = map_segment_vertices.begin(); it!= map_segment_vertices.end(); it++){ */
   /*  WCPPID::ProtoSegment *sg = it->first; */
   /*  if (sg->get_cluster_id() != temp_cluster->get_cluster_id()) continue; */
   /*  std::cout << sg->get_associate_points().size() << std::endl; */
   /* } */


  // find the vertex
  for (auto it = map_vertex_segments.begin(); it!= map_vertex_segments.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster->get_cluster_id() || it->second.size()<=1) continue;
    //    std::cout << it->second.size() << std::endl;
    
  }
  
  
}
