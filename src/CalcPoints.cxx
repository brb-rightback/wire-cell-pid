#include "WireCellPID/CalcPoints.h"

using namespace WireCell;

void WireCellPID::calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell, int nrebin, int frame_length, double unit_dis){
  GeomWireSelection wires_u = mcell->get_uwires();
  GeomWireSelection wires_v = mcell->get_vwires();
  GeomWireSelection wires_w = mcell->get_wwires();

  GeomWireSelection max_wires = wires_u;
  WirePlaneType_t max_wire_plane_type = WirePlaneType_t(0);
  GeomWireSelection min_wires = wires_v;
  WirePlaneType_t min_wire_plane_type = WirePlaneType_t(1);
  GeomWireSelection other_wires;
  WirePlaneType_t other_wire_plane_type;
  
  if (wires_v.size() > max_wires.size()){
    max_wires = wires_v;
    max_wire_plane_type = WirePlaneType_t(1);
  }
  if (wires_w.size() > max_wires.size()){
    max_wires = wires_w;
    max_wire_plane_type = WirePlaneType_t(2);
  }
  
  if (wires_u.size() < min_wires.size()){
    min_wires = wires_u;
    min_wire_plane_type = WirePlaneType_t(0);
  }
  if (wires_w.size() < min_wires.size()){
    min_wires = wires_w;
    min_wire_plane_type = WirePlaneType_t(2);
  }

  if (max_wire_plane_type==WirePlaneType_t(0) &&
      min_wire_plane_type==WirePlaneType_t(1) ||
      max_wire_plane_type==WirePlaneType_t(1) &&
      min_wire_plane_type==WirePlaneType_t(0)){
    other_wire_plane_type = WirePlaneType_t(2);
    other_wires = wires_w;
  }else if (max_wire_plane_type==WirePlaneType_t(0) &&
	    min_wire_plane_type==WirePlaneType_t(2) ||
	    max_wire_plane_type==WirePlaneType_t(2) &&
	    min_wire_plane_type==WirePlaneType_t(0)){
    other_wire_plane_type = WirePlaneType_t(1);
    other_wires = wires_v;
  }else if (max_wire_plane_type==WirePlaneType_t(2) &&
	    min_wire_plane_type==WirePlaneType_t(1) ||
	    max_wire_plane_type==WirePlaneType_t(1) &&
	    min_wire_plane_type==WirePlaneType_t(2)){
    other_wire_plane_type = WirePlaneType_t(0);
    other_wires = wires_u;
  }

  PointVector sampling_points;
  std::vector<std::tuple<int,int,int>> sampling_points_wires;
  
  float other_pitch = gds.pitch(other_wire_plane_type);
  float dis_limit[2];
  float tolerance = 0.1 * units::mm;
  const GeomWire *other_wire_1 = other_wires.front();
  const GeomWire *other_wire_2 = other_wires.back();
  dis_limit[0] = gds.wire_dist(*other_wire_1)-other_pitch/2.-tolerance;
  dis_limit[1] = gds.wire_dist(*other_wire_2)+other_pitch/2.+tolerance;

  //std::cout << dis_limit[0] << " " << dis_limit[1] << std::endl;

  
  
  int max_step = std::max(3.0,max_wires.size()/12.);
  int min_step = std::max(3.0,min_wires.size()/12.);
  //std::cout << min_step << " " << max_step << std::endl;
  GeomWireSetp max_wires_set;
  GeomWireSetp min_wires_set;
  max_wires_set.insert(max_wires.front());
  for (size_t i=0;i<max_wires.size();i+=max_step){
    max_wires_set.insert(max_wires.at(i));
  }
  max_wires_set.insert(max_wires.back());

  min_wires_set.insert(min_wires.front());
  for (size_t i=0;i<min_wires.size();i+=min_step){
    min_wires_set.insert(min_wires.at(i));
  }
  min_wires_set.insert(min_wires.back());

  for (auto it = max_wires_set.begin(); it!=max_wires_set.end(); it++){
    for (auto it1 = min_wires_set.begin(); it1!=min_wires_set.end();it1++){
      Vector point;
      gds.crossing_point(*(*it),*(*it1),point);
      float dis = gds.wire_dist(point,other_wire_plane_type);
      if (dis>=dis_limit[0]&&dis<=dis_limit[1]){
	point.x = (mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;

	int index_u;
	int index_v;
	int index_w;
	const GeomWire* wire_other = gds.closest(point,other_wire_plane_type);
	if (wire_other == 0 ) continue;
	if (max_wire_plane_type==WirePlaneType_t(0)){
	  index_u = (*it)->index();
	  if (min_wire_plane_type==WirePlaneType_t(1)){
	    index_v = (*it1)->index();
	    index_w = wire_other->index();
	  }else{
	    index_w = (*it1)->index();
	    index_v = wire_other->index();
	  }
	}else if (max_wire_plane_type==WirePlaneType_t(1)){
	  index_v = (*it)->index();
	  if (min_wire_plane_type==WirePlaneType_t(0)){
	    index_u = (*it1)->index();
	    index_w = wire_other->index();
	  }else{
	    index_w = (*it1)->index();
	    index_u = wire_other->index();
	  }
	}else{
	  index_w = (*it)->index();
	  if (min_wire_plane_type==WirePlaneType_t(0)){
	    index_u = (*it1)->index();
	    index_v = wire_other->index();
	  }else if (min_wire_plane_type==WirePlaneType_t(1)){
	    index_v = (*it1)->index();
	    index_u = wire_other->index();
	  }
	}
	sampling_points.push_back(point);
	sampling_points_wires.push_back(std::make_tuple(index_u, index_v, index_w));
	//std::cout << index_u << " " << index_v << " " <<index_w << std::endl;
      }
      //std::cout << "A: " <<dis_limit[0] << " " << dis << " " << dis_limit[1] << std::endl; 
    }
  }
  
  //std::cout << sampling_points.size() << " " << wires_u.size() <<  " " << wires_v.size() << " " << wires_w.size() << " " << max_wires.size() << " " << min_wires.size() << std::endl;
  
  if (sampling_points.size()>0){
    mcell->AddSamplingPoints(sampling_points);
    mcell->AddSamplingPointsWires(sampling_points_wires);
    mcell->SetMaxWireInterval(max_wire_plane_type,max_step);
    mcell->SetMinWireInterval(min_wire_plane_type,min_step);
    // std::cout << max_wire_plane_type << " " << max_step << " " << min_wire_plane_type << " " << min_step << std::endl;
  }
  // std::cout << min_wires.size() << " " << max_wires.size() << " " << wires_u.size() << " " << wires_v.size() << " " << wires_w.size() << std::endl;

  
  
}
