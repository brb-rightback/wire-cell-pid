#include "WCPPID/WCShower.h"

using namespace WCP;

WCPPID::WCShower::WCShower()
  : particle_type(0)
  , flag_shower(true)
  , kenergy_range(0)
  , kenergy_dQdx(0)
  , kenergy_charge(0)
  , start_vertex(0)
  , start_connection_type(0)
  , start_segment(0)
  , pcloud_fit(0)
  , pcloud_associated(0)
{
  start_point.x = 0;
  start_point.y = 0;
  start_point.z = 0;
  
  end_point.x = 0;
  end_point.y = 0;
  end_point.z = 0;

  init_dir.SetXYZ(0,0,0);
}

WCPPID::WCShower::~WCShower(){
  if (pcloud_fit != (ToyPointCloud*)0) delete pcloud_fit;
  if (pcloud_associated != (ToyPointCloud*)0) delete pcloud_associated;
}

void WCPPID::WCShower::rebuild_point_clouds(){
  if (pcloud_fit != (ToyPointCloud*)0) delete pcloud_fit;
  if (pcloud_associated != (ToyPointCloud*)0) delete pcloud_associated;
  pcloud_fit = 0;
  pcloud_associated = 0;
  build_point_clouds();
}

void WCPPID::WCShower::build_point_clouds(){
  if (pcloud_fit == 0){
    pcloud_fit = new ToyPointCloud();
    for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      WCP::PointVector& pts = sg->get_point_vec();
      pcloud_fit->AddPoints(pts);
    }
    pcloud_fit->build_kdtree_index();
  }
  if (pcloud_associated == 0){
    pcloud_associated = new ToyPointCloud();
    for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *sg = it->first;
      ToyPointCloud* sg_pcloud = sg->get_associated_pcloud();
      WCP::WCPointCloud<double>& sg_cloud = sg_pcloud->get_cloud();
      WCP::WC2DPointCloud<double>& sg_cloud_u = sg_pcloud->get_cloud_u();
      WCP::WC2DPointCloud<double>& sg_cloud_v = sg_pcloud->get_cloud_v();
      WCP::WC2DPointCloud<double>& sg_cloud_w = sg_pcloud->get_cloud_w();
      for (size_t i=0;i!=sg_cloud.pts.size();i++){
	pcloud_associated->AddPoint(sg_cloud.pts.at(i), sg_cloud_u.pts.at(i), sg_cloud_v.pts.at(i), sg_cloud_w.pts.at(i));
      }
    }
    pcloud_associated->build_kdtree_index();
  }
}



void WCPPID::WCShower::calculate_kinematics(){
  if (map_seg_vtxs.size()==1){
    particle_type = start_segment->get_particle_type();
    flag_shower = start_segment->get_flag_shower();
    kenergy_range = start_segment->cal_kine_range();
    kenergy_dQdx = start_segment->cal_kine_dQdx();
    if (start_segment->get_flag_dir()==1){
      start_point = start_segment->get_point_vec().front();
      end_point = start_segment->get_point_vec().back();
    }else if (start_segment->get_flag_dir()==-1){
      start_point = start_segment->get_point_vec().back();
      end_point = start_segment->get_point_vec().front();
    }
    // initial direction ...
    if (start_connection_type == 1){
      init_dir = start_segment->cal_dir_3vector();
    }else if (start_connection_type == 2){
      init_dir.SetXYZ(start_point.x - start_vertex->get_fit_pt().x, start_point.y - start_vertex->get_fit_pt().y, start_point.z - start_vertex->get_fit_pt().z);
    }else if (start_connection_type == 3){
      init_dir.SetXYZ(start_point.x - start_vertex->get_fit_pt().x, start_point.y - start_vertex->get_fit_pt().y, start_point.z - start_vertex->get_fit_pt().z);
    }
    init_dir = init_dir.Unit();
  }else{
    // two different types ...
  }
}

void WCPPID::WCShower::set_start_vertex(ProtoVertex* vertex, int type){
  start_vertex = vertex;
  start_connection_type = type;
}

void WCPPID::WCShower::set_start_segment(ProtoSegment* seg){
  start_segment = seg;
}

void WCPPID::WCShower::fill_sets( std::set<WCPPID::ProtoVertex* >& used_vertices,  std::set<WCPPID::ProtoSegment* >& used_segments, bool flag_exclude_start_segment){
  for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
    used_vertices.insert(it->first);
  }
  for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
    if (flag_exclude_start_segment){
      if (it->first == start_segment) continue;
    }
    used_segments.insert(it->first);
  }
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
