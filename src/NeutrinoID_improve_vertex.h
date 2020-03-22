#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "TVector3.h"

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
    fit_vertex(vtx, it->second);
  }
}

WCPPID::MyFCN::MyFCN(WCPPID::ProtoVertex* vtx, bool flag_vtx_constraint, double vertex_protect_dis, double point_track_dis, double fit_dis) 
: vtx(vtx)
, flag_vtx_constraint(flag_vtx_constraint)
, vertex_protect_dis(vertex_protect_dis)
  , point_track_dis(point_track_dis)
  , fit_dis(fit_dis)
{
  segments.clear();
  vec_points.clear();
}

WCPPID::MyFCN::~MyFCN(){
}

void WCPPID::MyFCN::AddSegment(ProtoSegment *sg){
  // push in ...
  segments.push_back(sg);
  {
    WCP::PointVector pts;
    vec_points.push_back(pts);
  }

  WCP::PointVector& pts = sg->get_associate_points();
  for (size_t i=0;i!=pts.size();i++){
    double dis_to_vertex = sqrt(pow(pts.at(i).x - vtx->get_fit_pt().x,2) + pow(pts.at(i).y - vtx->get_fit_pt().y,2) + pow(pts.at(i).z - vtx->get_fit_pt().z,2));
    if (dis_to_vertex < vertex_protect_dis || dis_to_vertex > fit_dis) continue;
    if (sg->get_closest_point(pts.at(i)).first > point_track_dis) continue;
    vec_points.back().push_back(pts.at(i));
  }
}

void WCPPID::MyFCN::update_fit_range(double tmp_vertex_protect_dis, double tmp_point_track_dis, double tmp_fit_dis){
  vertex_protect_dis = tmp_vertex_protect_dis;
  point_track_dis = tmp_point_track_dis;
  fit_dis = tmp_fit_dis;

  std::vector<WCPPID::ProtoSegment* > tmp_segments = segments;
  segments.clear();
  vec_points.clear();
  for (auto it = tmp_segments.begin(); it != tmp_segments.end(); it++){
    AddSegment(*it);
  }  
}

int WCPPID::MyFCN::get_fittable_tracks(){
  int ncount = 0;
  for (size_t i=0;i!=vec_points.size();i++){
    if (vec_points.at(i).size()>0) ncount++;
  }
  return ncount;
}


std::pair<WCPPID::ProtoSegment*, int> WCPPID::MyFCN::get_seg_info(int i){
  if (i < segments.size()){
    return std::make_pair(segments.at(i), vec_points.at(i).size());
  }
  return std::make_pair((WCPPID::ProtoSegment*)0, 0);
}

double WCPPID::MyFCN::operator() (const std::vector<double> & xx) const{
  return get_chi2(xx);
}

double WCPPID::MyFCN::get_chi2(const std::vector<double> & xx) const{
  return 1;
}



void WCPPID::NeutrinoID::fit_vertex(WCPPID::ProtoVertex *vtx, WCPPID::ProtoSegmentSet& sg_set){
  WCPPID::MyFCN fcn(vtx);
  for (auto it = sg_set.begin(); it!=sg_set.end(); it++){
    fcn.AddSegment(*it);
  }
}



