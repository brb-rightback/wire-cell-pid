#include "TString.h"
#include "WCPData/Line.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
  TString filename = argv[1];
    std::ifstream infile(filename);
  Int_t type;
  infile >> type;
  Double_t vx, vy, vz;
  Double_t px, py, pz, x, y, z;
  if (type==1){
    infile >> vx >> vy >> vz;
    WCP::Point center(vx, vy, vz);
    infile >> px >> py >> pz >> x >> y >> z;
    Double_t mom1 = sqrt(px*px + py*py + pz*pz);
    TVector3 dir11(x-vx,y-vy, z-vz);
    infile >> px >> py >> pz >> x >> y >> z;
    Double_t mom2 = sqrt(px*px + py*py + pz*pz);
    TVector3 dir22(x-vx,y-vy,z-vz);
    double angle = dir11.Angle(dir22);
    double mass_pio = sqrt(4*mom1*mom2*pow(sin(angle/2.),2))*0.95;
    std::cout << mass_pio << std::endl;
    
  }else if (type==2){
    infile >> px >> py >> pz >> x >> y >> z;
    Double_t mom1 = sqrt(px*px + py*py + pz*pz);
    TVector3 dir1(px/mom1, py/mom1, pz/mom1);
    WCP::Point p1(x,y,z);
    WCP::Line l1(p1, dir1);
    
    infile >> px >> py >> pz >> x >> y >> z;
    Double_t mom2 = sqrt(px*px + py*py + pz*pz);
    TVector3 dir2(px/mom1, py/mom1, pz/mom1);
    WCP::Point p2(x,y,z);
    WCP::Line l2(p2, dir2);
    
    auto pair_points = l1.closest_dis_points(l2);
    
    WCP::Point center(0,0,0);
    center.x = (pair_points.first.x + pair_points.second.x)/2.;
    center.y = (pair_points.first.y + pair_points.second.y)/2.;
    center.z = (pair_points.first.z + pair_points.second.z)/2.;
    //calculate mass ...
    TVector3 dir11(l1.get_p1().x - center.x, l1.get_p1().y - center.y, l1.get_p1().z - center.z);
    TVector3 dir22(l2.get_p1().x - center.x, l2.get_p1().y - center.y, l2.get_p1().z - center.z);
    double angle = dir11.Angle(dir22);
    double mass_pio = sqrt(4*mom1*mom2*pow(sin(angle/2.),2))*0.95;
    std::cout << mass_pio << std::endl;
  }
}
