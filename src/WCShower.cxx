#include "WCPPID/WCShower.h"

#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"
using namespace WCP;

WCPPID::WCShower::WCShower()
  : particle_type(0)
  , flag_shower(true)
  , flag_kinematics(false)
  , kenergy_range(0)
  , kenergy_dQdx(0)
  , kenergy_charge(0)
  , kenergy_best(0)
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


std::vector<double> WCPPID::WCShower::get_stem_dQ_dx(WCPPID::ProtoVertex *vertex, WCPPID::ProtoSegment *sg, int limit){
  std::vector<double> vec_dQ_dx;
  {
    std::vector<double>& vec_dQ = sg->get_dQ_vec();
    std::vector<double>& vec_dx = sg->get_dx_vec();    
    if (vertex->get_wcpt().index == sg->get_wcpt_vec().front().index){
      for (size_t i=0;i!=vec_dQ.size();i++){
	vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
	if (vec_dQ_dx.size() >= limit) break;
      }
    }else{
      for (int i=int(vec_dQ.size())-1;i>=0;i--){
	vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
	if (vec_dQ_dx.size() >= limit) break;
      }
    }
  }
  
  if (sg == start_segment && vec_dQ_dx.size() < limit){
    WCPPID::ProtoVertex *curr_vertex = vertex;
    WCPPID::ProtoSegment *curr_segment = sg;
    WCPPID::ProtoVertex *next_vertex = 0;
    WCPPID::ProtoSegment *next_segment = 0;
    int count = 0;
    while(vec_dQ_dx.size()< limit && count < 3){
      // find next vertex
      for (auto it = map_seg_vtxs[curr_segment].begin(); it != map_seg_vtxs[curr_segment].end(); it++){
	if ( (*it) == curr_vertex ) continue;
	next_vertex = *it;
      }
      //      std::cout << "V: " << next_vertex << std::endl;
      if (next_vertex ==0) break;

      // find the next segment
      TVector3 dir1(curr_vertex->get_fit_pt().x - next_vertex->get_fit_pt().x,
		    curr_vertex->get_fit_pt().y - next_vertex->get_fit_pt().y,
		    curr_vertex->get_fit_pt().z - next_vertex->get_fit_pt().z);
      TVector3 dir2;
      double max_angle = 0;
      for (auto it = map_vtx_segs[next_vertex].begin(); it != map_vtx_segs[next_vertex].end(); it++){
	if ((*it) == curr_segment) continue;
	TVector3 tmp_dir = (*it)->cal_dir_3vector(next_vertex->get_fit_pt(), 10*units::cm);
	double angle = dir1.Angle(tmp_dir)/3.1415926*180.;
	if (angle > max_angle){
	  max_angle = angle;
	  next_segment = *it;
	  dir2 = tmp_dir;
	}
      }
      //std::cout << "S: " << next_segment << std::endl;
      
      if (next_segment ==0) break;
      // check the next segment
      bool flag_bad = false;
      for (auto it = map_vtx_segs[next_vertex].begin(); it != map_vtx_segs[next_vertex].end(); it++){
	if ((*it) == curr_segment) continue;
	if ((*it) == next_segment) continue;
	TVector3 tmp_dir = (*it)->cal_dir_3vector(next_vertex->get_fit_pt(), 10*units::cm);
	if ((*it)->get_length() > 3*units::cm && dir2.Angle(tmp_dir)/3.1415926*180. < 25) flag_bad = true;
      }
      //std::cout << "C: " << flag_bad << std::endl;
      if (flag_bad) break;
      // keep add dQ/dx

      vec_dQ_dx.pop_back();
      std::vector<double>& vec_dQ = next_segment->get_dQ_vec();
      std::vector<double>& vec_dx = next_segment->get_dx_vec();    
      if (next_vertex->get_wcpt().index == next_segment->get_wcpt_vec().front().index){
	for (size_t i=0;i!=vec_dQ.size();i++){
	  vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
	  if (vec_dQ_dx.size() >= limit) break;
	}
      }else{
	for (int i=int(vec_dQ.size())-1;i>=0;i--){
	  vec_dQ_dx.push_back(vec_dQ.at(i)/(vec_dx.at(i)+1e-9)/(43e3/units::cm));
	  if (vec_dQ_dx.size() >= limit) break;
	}
      }
      
      if (vec_dQ_dx.size() >= limit) break;
      // prepare the next ...
      curr_vertex = next_vertex;
      next_vertex = 0;
      curr_segment = next_segment;
      next_segment = 0;
      count ++;
    }
  }
  
  return vec_dQ_dx;
}


int WCPPID::WCShower::get_num_main_segments(){
  int num = 0;
  for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() == start_segment->get_cluster_id()) num ++;
  }
  return num;
}


double WCPPID::WCShower::get_closest_dis(WCPPID::ProtoSegment *seg){
  Point test_p = get_closest_point(seg->get_point_vec().front()).second;
  test_p = seg->get_closest_point(test_p).second;
  test_p = get_closest_point(test_p).second;
  return seg->get_closest_point(test_p).first;
}

std::pair<double, WCP::Point> WCPPID::WCShower::get_closest_point(WCP::Point& p){
  if (pcloud_fit == (ToyPointCloud*)0)
    rebuild_point_clouds();
  

  return pcloud_fit->get_closest_point(p);
}


void WCPPID::WCShower::rebuild_point_clouds(){
  if (pcloud_fit != (ToyPointCloud*)0) delete pcloud_fit;
  if (pcloud_associated != (ToyPointCloud*)0) delete pcloud_associated;
  pcloud_fit = 0;
  pcloud_associated = 0;
  build_point_clouds();

}

double WCPPID::WCShower::get_total_track_length(){
  double total_length = 0;
  for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (!sg->get_flag_shower())     total_length += sg->get_length();
  }
  return total_length;
}

double WCPPID::WCShower::get_total_length(int tmp_cluster_id){
  double total_length = 0;
  for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (sg->get_cluster_id() == tmp_cluster_id)
      total_length += sg->get_length();
  }
  return total_length;
}

double WCPPID::WCShower::get_total_length(){
  double total_length = 0;
  for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    total_length += sg->get_length();
  }
  return total_length;
}

void WCPPID::WCShower::fill_point_vec(PointVector& tmp_pts, bool flag_main){
  for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    if (flag_main && sg->get_cluster_id() != start_segment->get_cluster_id()) continue;
    WCP::PointVector& pts = sg->get_point_vec();
    for (size_t i=1;i+1<pts.size();i++){
      tmp_pts.push_back(pts.at(i));
    }
  }
  for (auto it = map_vtx_segs.begin(); it!=map_vtx_segs.end();it++){
    WCPPID::ProtoVertex *vtx = it->first;
    if (flag_main && vtx->get_cluster_id() != start_segment->get_cluster_id()) continue;
    tmp_pts.push_back(vtx->get_fit_pt());
  }

  
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
      if (sg_pcloud!=0){
	WCP::WCPointCloud<double>& sg_cloud = sg_pcloud->get_cloud();
	WCP::WC2DPointCloud<double>& sg_cloud_u = sg_pcloud->get_cloud_u();
	WCP::WC2DPointCloud<double>& sg_cloud_v = sg_pcloud->get_cloud_v();
	WCP::WC2DPointCloud<double>& sg_cloud_w = sg_pcloud->get_cloud_w();
	for (size_t i=0;i!=sg_cloud.pts.size();i++){
	  pcloud_associated->AddPoint(sg_cloud.pts.at(i), sg_cloud_u.pts.at(i), sg_cloud_v.pts.at(i), sg_cloud_w.pts.at(i));
	}
      }
    }
    pcloud_associated->build_kdtree_index();
  }
}

std::pair<WCPPID::ProtoSegment*, WCPPID::ProtoVertex*> WCPPID::WCShower::get_last_segment_vertex_long_muon(std::set<WCPPID::ProtoSegment*> segments_in_muons){
  WCPPID::ProtoVertex *s_vtx = start_vertex;
  WCPPID::ProtoSegment *s_seg = start_segment;

  // std::cout << map_vtx_segs.size() << " " << map_seg_vtxs.size() << std::endl;
  // for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
  //   std::cout << it->first->get_id() << std::endl;
  // }
  //  std::cout << map_vtx_segs[s].size() << " " << map_seg_vtxs.size() << std::endl;
  std::set<ProtoSegment*> used_segments;
  used_segments.insert(s_seg);
  
  bool flag_continue = true;
  while(flag_continue){
    flag_continue = false;
    if (s_vtx == start_vertex){
      flag_continue = true;
    }else{
      for (auto it = map_vtx_segs[s_vtx].begin(); it!= map_vtx_segs[s_vtx].end(); it++){
	WCPPID::ProtoSegment *sg = *it;
	if (segments_in_muons.find(sg) != segments_in_muons.end() &&
	    used_segments.find(sg) == used_segments.end()){
	  //	    sg != s_seg){
	  //	  std::cout << sg->get_id() << " " << s_seg->get_id() << std::endl;
	  s_seg = sg;
	  used_segments.insert(s_seg);
	  flag_continue = true;
	  break;
	}
      }
    }
    
    if (flag_continue){
      for (auto it = map_seg_vtxs[s_seg].begin(); it != map_seg_vtxs[s_seg].end(); it++){
	//	std::cout <<s_seg->get_id() <<  " " << map_seg_vtxs[s_seg].size() << std::endl;
	WCPPID::ProtoVertex *vtx = *it;
	if (vtx != s_vtx){
	 
	  s_vtx = vtx;
	  break;
	}
      }
    }
  }
  return std::make_pair(s_seg, s_vtx);
}

void WCPPID::WCShower::calculate_kinematics_long_muon(std::set<WCPPID::ProtoSegment*> segments_in_muons){
   particle_type = start_segment->get_particle_type();
   flag_shower = false;

   double L=0;
   std::vector<double> vec_dQ, vec_dx;
   for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
     WCPPID::ProtoSegment *sg = it->first;
     if (segments_in_muons.find(sg) != segments_in_muons.end())   L += sg->get_length();
     vec_dQ.insert(vec_dQ.end(), sg->get_dQ_vec().begin(), sg->get_dQ_vec().end());
     vec_dx.insert(vec_dx.end(), sg->get_dx_vec().begin(), sg->get_dx_vec().end());
   }
   kenergy_range = start_segment->cal_kine_range(L);
   kenergy_dQdx = start_segment->cal_kine_dQdx(vec_dQ, vec_dx);
   // long muon ... // should be improve in the future using range ...
   kenergy_best = kenergy_dQdx;

   //std::cout << kenergy_range <<  " " << kenergy_dQdx << " " << kenergy_best << std::endl;
   
   init_dir = start_segment->cal_dir_3vector();
   if (start_segment->get_flag_dir()==1){
     start_point = start_segment->get_point_vec().front();
   }else if (start_segment->get_flag_dir()==-1){
     start_point = start_segment->get_point_vec().back();
   }
   double max_dis = 0; Point max_point;
   for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
     WCPPID::ProtoVertex *vtx = it->first;

     bool flag_contain = false;
     
     for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){ 
       if (segments_in_muons.find(*it1) != segments_in_muons.end()) {
	 flag_contain = true;
	 break;
       }
     }

     if (flag_contain){
       
       double dis = sqrt(pow(start_point.x - vtx->get_fit_pt().x,2) + pow(start_point.y - vtx->get_fit_pt().y,2) + pow(start_point.z - vtx->get_fit_pt().z,2));
       if (dis > max_dis){
	 max_dis = dis;
	   max_point = vtx->get_fit_pt();
       }
     }
   } // far of the vertex ...
   end_point = max_point;
}


void WCPPID::WCShower::calculate_kinematics(){

  //  std::cout << map_seg_vtxs.size() << " " << start_segment->get_point_vec().size() << " " << start_segment->get_particle_type() << std::endl;
  
  if (map_seg_vtxs.size()==1){
    particle_type = start_segment->get_particle_type();
    flag_shower = start_segment->get_flag_shower();
    kenergy_range = start_segment->cal_kine_range();
    kenergy_dQdx = start_segment->cal_kine_dQdx();

    // single track
    if (start_connection_type == 1){
      // whether it is shower or not
      if (start_segment->get_length() < 4*units::cm){
	kenergy_best = kenergy_dQdx;
      }else{
	kenergy_best = kenergy_range;
      }
    }else{
      if (flag_shower){
	kenergy_best = 0;
      }else{
	if (start_segment->get_length() < 4*units::cm){
	  kenergy_best = kenergy_dQdx;
	}else{
	  kenergy_best = kenergy_range;
	}
      }
    }
    if (start_connection_type ==1|| pcloud_fit == 0 ){
      if (start_segment->get_flag_dir()==1 ){
	start_point = start_segment->get_point_vec().front();
	end_point = start_segment->get_point_vec().back();
      }else if (start_segment->get_flag_dir()==-1){
	start_point = start_segment->get_point_vec().back();
	end_point = start_segment->get_point_vec().front();
      }
    }else{
      start_point = get_closest_point(start_vertex->get_fit_pt()).second;
      double max_dis = 0; Point max_point;
      for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
	WCPPID::ProtoVertex *vtx = it->first;
	double dis = sqrt(pow(start_point.x - vtx->get_fit_pt().x,2) + pow(start_point.y - vtx->get_fit_pt().y,2) + pow(start_point.z - vtx->get_fit_pt().z,2));
	if (dis > max_dis){
	  max_dis = dis;
	  max_point = vtx->get_fit_pt();
	}
      } // far of the vertex ...
      end_point = max_point;
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
    // more than one segment ...
    int nsegments = map_seg_vtxs.size();
    int nconnected_segs = 0;
    {
      std::pair<std::set<WCPPID::ProtoSegment*>, std::set<WCPPID::ProtoVertex*> > results = get_connected_pieces(start_segment);
      nconnected_segs = results.first.size();
    }
    // two different types ...
    if (nsegments == nconnected_segs){
      // single track
      particle_type = start_segment->get_particle_type();
      flag_shower = start_segment->get_flag_shower();

      if (start_connection_type==1 || pcloud_fit == 0){
	if (start_segment->get_flag_dir()==1){
	  start_point = start_segment->get_point_vec().front();
	}else if (start_segment->get_flag_dir()==-1){
	  start_point = start_segment->get_point_vec().back();
	}
      }else{
	start_point = get_closest_point(start_vertex->get_fit_pt()).second;
      }
      // initial direction ...
      if (start_connection_type == 1){
	if (start_segment->get_length() > 8*units::cm)
	  init_dir = start_segment->cal_dir_3vector();
	else
	  init_dir = cal_dir_3vector(start_vertex->get_fit_pt(),12*units::cm);
      }else if (start_connection_type == 2){
	//std::cout << start_segment->get_id() << " " << start_point << std::endl;
	init_dir.SetXYZ(start_point.x - start_vertex->get_fit_pt().x, start_point.y - start_vertex->get_fit_pt().y, start_point.z - start_vertex->get_fit_pt().z);
      }else if (start_connection_type == 3){
	init_dir.SetXYZ(start_point.x - start_vertex->get_fit_pt().x, start_point.y - start_vertex->get_fit_pt().y, start_point.z - start_vertex->get_fit_pt().z);
      }
      init_dir = init_dir.Unit();
      
      double max_dis = 0; Point max_point;
      for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
	WCPPID::ProtoVertex *vtx = it->first;
	double dis = sqrt(pow(start_point.x - vtx->get_fit_pt().x,2) + pow(start_point.y - vtx->get_fit_pt().y,2) + pow(start_point.z - vtx->get_fit_pt().z,2));
	if (dis > max_dis){
	  max_dis = dis;
	  max_point = vtx->get_fit_pt();
	}
      } // far of the vertex ...
      end_point = max_point;

      double L=0;
      std::vector<double> vec_dQ, vec_dx;
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg = it->first;
	L += sg->get_length();
	vec_dQ.insert(vec_dQ.end(), sg->get_dQ_vec().begin(), sg->get_dQ_vec().end());
	vec_dx.insert(vec_dx.end(), sg->get_dx_vec().begin(), sg->get_dx_vec().end());
      }
      kenergy_range = start_segment->cal_kine_range(L);
      kenergy_dQdx = start_segment->cal_kine_dQdx(vec_dQ, vec_dx);

      if (start_connection_type == 1){
	// whether it is shower or not
	if (start_segment->get_length() < 4*units::cm){
	  kenergy_best = kenergy_dQdx;
	}else{
	  kenergy_best = kenergy_range;
	}
      }else{
	if (flag_shower){
	  kenergy_best = 0;
	}else{
	  if (start_segment->get_length() < 4*units::cm){
	    kenergy_best = kenergy_dQdx;
	  }else{
	    kenergy_best = kenergy_range;
	  }
	}
      }
    }else{
      // multiple tracks ...
      particle_type = start_segment->get_particle_type();
      flag_shower = start_segment->get_flag_shower();
      if (start_connection_type==1|| pcloud_fit == 0){
	if (start_segment->get_flag_dir()==1){
	  start_point = start_segment->get_point_vec().front();
	}else if (start_segment->get_flag_dir()==-1){
	  start_point = start_segment->get_point_vec().back();
	}
      }else{
	start_point = get_closest_point(start_vertex->get_fit_pt()).second;
      }
       
      // initial direction ...
      if (start_connection_type == 1){
	if (start_segment->get_length() > 8*units::cm)
	  init_dir = start_segment->cal_dir_3vector();
	else
	  init_dir = cal_dir_3vector(start_vertex->get_fit_pt(),12*units::cm);
      }else if (start_connection_type == 2){
	init_dir.SetXYZ(start_point.x - start_vertex->get_fit_pt().x, start_point.y - start_vertex->get_fit_pt().y, start_point.z - start_vertex->get_fit_pt().z);
      }else if (start_connection_type == 3){
	init_dir.SetXYZ(start_point.x - start_vertex->get_fit_pt().x, start_point.y - start_vertex->get_fit_pt().y, start_point.z - start_vertex->get_fit_pt().z);
      }
      init_dir = init_dir.Unit();
      
     
      double max_dis = 0; Point max_point;
      for (auto it = map_vtx_segs.begin(); it != map_vtx_segs.end(); it++){
	WCPPID::ProtoVertex *vtx = it->first;
	double dis = sqrt(pow(start_point.x - vtx->get_fit_pt().x,2) + pow(start_point.y - vtx->get_fit_pt().y,2) + pow(start_point.z - vtx->get_fit_pt().z,2));
	if (dis > max_dis){
	  max_dis = dis;
	  max_point = vtx->get_fit_pt();
	}
      } // far of the vertex ...
      end_point = max_point;
      
      std::vector<double> vec_dQ, vec_dx;
      for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
	WCPPID::ProtoSegment *sg = it->first;
	vec_dQ.insert(vec_dQ.end(), sg->get_dQ_vec().begin(), sg->get_dQ_vec().end());
	vec_dx.insert(vec_dx.end(), sg->get_dx_vec().begin(), sg->get_dx_vec().end());
      }
      kenergy_range = 0;
      kenergy_dQdx = start_segment->cal_kine_dQdx(vec_dQ, vec_dx);
      kenergy_best = 0;
    }
  }

  // std::cout << map_seg_vtxs.size() << " " << start_connection_type  << " " << init_dir.X() << " " <<  init_dir.Y() << " " << init_dir.Z() << " " << kenergy_dQdx/units::MeV << " " << start_point << " " << start_vertex->get_fit_pt() << std::endl;
}

void WCPPID::WCShower::set_start_vertex(ProtoVertex* vertex, int type){
  start_vertex = vertex;
  start_connection_type = type;
}

void WCPPID::WCShower::set_start_segment(ProtoSegment* seg){
  start_segment = seg;
}

void WCPPID::WCShower::set_start_point(Point p){
  start_point = p;
}

void WCPPID::WCShower::set_start_segment(ProtoSegment* seg, Map_Proto_Segment_Vertices& map_segment_vertices){
  start_segment = seg;
  for (auto it = map_segment_vertices[start_segment].begin(); it!= map_segment_vertices[start_segment].end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    if (vtx == start_vertex) continue;
    map_vtx_segs[vtx].insert(start_segment);
    map_seg_vtxs[start_segment].insert(vtx);
  }
}

void WCPPID::WCShower::update_particle_type(){
  double track_length = 0; double shower_length = 0;
  if (map_seg_vtxs.size()>1){
    for (auto it = map_seg_vtxs.begin(); it!=map_seg_vtxs.end(); it++){
      WCPPID::ProtoSegment *seg = it->first;
      if (seg->get_flag_shower() || fabs(seg->get_particle_type())!=2212 ){ // not proton ...
 	shower_length += seg->get_length();
      }else{
	track_length += seg->get_length();
      }
    }
    if (shower_length > track_length){
      start_segment->set_particle_type(11);
      TPCParams& mp = Singleton<TPCParams>::Instance();
      start_segment->set_particle_mass(mp.get_mass_electron());
      start_segment->cal_4mom();
    }
    //    std::cout << start_segment->get_cluster_id() << " " << shower_length << " " << track_length << std::endl;
  }
}

TVector3 WCPPID::WCShower::cal_dir_3vector(WCP::Point p, double dis_cut ){
  WCP::Point p1(0,0,0);
  Int_t ncount = 0;
  
  for (auto it = map_seg_vtxs.begin(); it != map_seg_vtxs.end(); it++){
    PointVector& pts = it->first->get_point_vec();
    for (size_t i=0;i!=pts.size();i++){
      double dis = sqrt(pow(pts.at(i).x - p.x,2)+pow(pts.at(i).y - p.y,2)+pow(pts.at(i).z - p.z,2));
      //    std::cout << dis/units::cm << std::endl;
      if (dis < dis_cut){
	p1.x += pts.at(i).x;
	p1.y += pts.at(i).y;
	p1.z += pts.at(i).z;
	ncount ++;
      }
    }
  }
  
  TVector3 v1(p1.x/ncount - p.x, p1.y/ncount - p.y, p1.z/ncount - p.z);
  v1 = v1.Unit();
  return v1;
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


std::pair<std::set<WCPPID::ProtoSegment*>, std::set<WCPPID::ProtoVertex*> > WCPPID::WCShower::get_connected_pieces(WCPPID::ProtoSegment* tseg){
  std::set<WCPPID::ProtoSegment* > used_segments;
  std::set<WCPPID::ProtoVertex* > used_vertices;

  std::vector<ProtoSegment* > new_segments;
  std::vector<ProtoVertex* > new_vertices;

  new_segments.push_back(tseg);
  used_segments.insert(tseg);

  while(new_vertices.size()>0 || new_segments.size()>0 ){
    if (new_vertices.size()>0){
      ProtoVertex *vtx = new_vertices.back();
      new_vertices.pop_back();
      for (auto it = map_vtx_segs[vtx].begin(); it != map_vtx_segs[vtx].end(); it++){
    	ProtoSegment *seg1 = *it;
    	if (used_segments.find(seg1)!=used_segments.end()) continue;
    	new_segments.push_back(seg1);
	used_segments.insert(seg1);
      }
    }

    if (new_segments.size()>0){
      ProtoSegment *seg1 = new_segments.back();
      new_segments.pop_back();
      for (auto it = map_seg_vtxs[seg1].begin(); it!= map_seg_vtxs[seg1].end(); it++){
  	ProtoVertex *vtx = *it;
  	if (used_vertices.find(vtx)!=used_vertices.end()) continue;
	new_vertices.push_back(vtx);
	used_vertices.insert(vtx);
      }
    }
  }
  
  
  return std::make_pair(used_segments, used_vertices);
}

double WCPPID::WCShower::get_dis(WCPPID::ProtoSegment* seg){
  double min_dis = 1e9;
  Point min_point;
  Point test_p = seg->get_point_vec().front();
  for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    std::pair<double, WCP::Point> results = sg->get_closest_point(test_p);
    if (results.first < min_dis){
      min_dis = results.first;
      min_point = results.second;
    }
  }
  std::pair<double, WCP::Point> results1 = seg->get_closest_point(min_point);
  test_p = results1.second;
  for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *sg = it->first;
    std::pair<double, WCP::Point> results = sg->get_closest_point(test_p);
    if (results.first < min_dis){
      min_dis = results.first;
      min_point = results.second;
    }
  }
  return min_dis;
}

void WCPPID::WCShower::add_shower(WCPPID::WCShower* temp_shower){
  Map_Proto_Segment_Vertices& tmp_map_seg_vtxs = temp_shower->get_map_seg_vtxs();
  for (auto it = tmp_map_seg_vtxs.begin(); it!= tmp_map_seg_vtxs.end(); it++){
    WCPPID::ProtoSegment *seg = it->first;
    for (auto it1 = it->second.begin(); it1!= it->second.end(); it1++){
      WCPPID::ProtoVertex *vtx = *it1;
      map_seg_vtxs[seg].insert(vtx);
      map_vtx_segs[vtx].insert(seg);
    }
  }
  rebuild_point_clouds();
}

void WCPPID::WCShower::add_segment(WCPPID::ProtoSegment* seg, Map_Proto_Segment_Vertices& map_segment_vertices){
  for (auto it = map_segment_vertices[seg].begin(); it != map_segment_vertices[seg].end(); it++){
    WCPPID::ProtoVertex *vtx = *it;
    map_seg_vtxs[seg].insert(vtx);
    map_vtx_segs[vtx].insert(seg);
  }
  rebuild_point_clouds();
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
	if (vtx == start_vertex) continue;
	if (map_vtx_segs.find(vtx)!= map_vtx_segs.end()){
	  map_vtx_segs[vtx].insert(seg);
	  map_seg_vtxs[seg].insert(vtx);
	}else{
	  map_vtx_segs[vtx].insert(seg);
	  map_seg_vtxs[seg].insert(vtx);
	  new_vertices.push_back(vtx);
	}
  	
      }
    }
  }

  // for (auto it = map_seg_vtxs.begin(); it!= map_seg_vtxs.end(); it++){
  //   std::cout << it->first->get_id() << " " << it->second.size() << std::endl;
  // }
  
  rebuild_point_clouds();
  //std::cout << map_vtx_segs.size() << " " << map_seg_vtxs.size() << std::endl;
}
