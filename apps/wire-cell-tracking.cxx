#include "WireCellSst/GeomDataSource.h"

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/ToyCTPointCloud.h"

#include "WireCellPID/CalcPoints.h"
#include "WireCellPID/PR3DCluster.h"

#include "WireCellPID/ExecMon.h"
#include "WireCellPID/ImprovePR3DCluster.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/matching.root -t[0,1? in time flash only] -c[0,1? main cluster only] -b[0,1? bee output]" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  
  bool flag_in_time_only = false; // default to run all code
  bool flag_main_cluster_only = true; // default to run only on the main cluster
  bool flag_bee_output = true; // default save 
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 't':
      flag_in_time_only = atoi(&argv[i][2]); 
      break;
    case 'c':
      flag_main_cluster_only = atoi(&argv[i][2]);
      break;
    case 'b':
      flag_bee_output = atoi(&argv[i][2]);
      break;
    }
  }
    
  WireCellPID::ExecMon em("starting");
  cout << em("load geometry") << endl;
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cout << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;

  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");

  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  float unit_dis;

  std::vector<int> *timesliceId = new std::vector<int>;
  std::vector<std::vector<int>> *timesliceChannel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge_err = new std::vector<std::vector<int>>;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  unsigned int triggerbits;
  Trun->SetBranchAddress("triggerBits",&triggerbits);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  
  Trun->SetBranchAddress("timesliceId",&timesliceId);
  Trun->SetBranchAddress("timesliceChannel",&timesliceChannel);
  Trun->SetBranchAddress("raw_charge",&raw_charge);
  Trun->SetBranchAddress("raw_charge_err",&raw_charge_err);
    
  Trun->GetEntry(0);
  double lowerwindow = 0;
  double upperwindow = 0;
  if(triggerbits==2048) { lowerwindow = 3.0; upperwindow = 5.0; }// bnb  
  if(triggerbits==512) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb
  
   // define singleton ... 
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  double angle_u = gds.angle(WirePlaneType_t(0));
  double angle_v = gds.angle(WirePlaneType_t(1));
  double angle_w = gds.angle(WirePlaneType_t(2));
  
  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_angle_u(angle_u);
  mp.set_angle_v(angle_v);
  mp.set_angle_w(angle_w);
  mp.set_ts_width(time_slice_width);


   // test geometry ...
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),0);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),0);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),0);
  double first_u_dis = gds.wire_dist(*uwire) ; // first U wire center ...
  double first_v_dis = gds.wire_dist(*vwire) ; // first V wire center ...
  double first_w_dis = gds.wire_dist(*wwire) ; // first W wire center ... 
  
  mp.set_first_u_dis(first_u_dis);
  mp.set_first_v_dis(first_v_dis);
  mp.set_first_w_dis(first_w_dis);

  std::map<int,std::pair<double,double>> dead_u_index;
  std::map<int,std::pair<double,double>> dead_v_index;
  std::map<int,std::pair<double,double>> dead_w_index;

  TTree *T_flash = (TTree*)file->Get("T_flash");
  Double_t time;
  Int_t type;
  Int_t flash_id;
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("flash_id",&flash_id);
  
  TTree *T_match = (TTree*)file->Get("T_match");
  Int_t tpc_cluster_id;
  Int_t event_type;
  T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id); // parent cluster id
  T_match->SetBranchAddress("flash_id",&flash_id);  // flash id 
  T_match->SetBranchAddress("event_type",&event_type);

  std::map<int, std::pair<int, double> > map_flash_info;
  std::map<int, int> map_flash_tpc_ids;
  std::map<int, int> map_tpc_flash_ids;
  
  
  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    map_flash_info[flash_id] = std::make_pair(type, time);
  }
  for (int i=0;i!=T_match->GetEntries();i++){
    T_match->GetEntry(i);
    if (flash_id==-1) continue;
    map_flash_tpc_ids[flash_id] = tpc_cluster_id;
    map_tpc_flash_ids[tpc_cluster_id] = flash_id;
  }
  
 

  // load mcell
  TTree *TC = (TTree*)file->Get("TC");
  std::vector<int> *cluster_id_vec = new std::vector<int>;
  std::vector<int> *parent_cluster_id = new std::vector<int>;
  std::vector<int> *time_slice_vec = new std::vector<int>;
  std::vector<double> *q_vec = new std::vector<double>;
  std::vector<double> *uq_vec = new std::vector<double>;
  std::vector<double> *vq_vec = new std::vector<double>;
  std::vector<double> *wq_vec = new std::vector<double>;
  std::vector<double> *udq_vec = new std::vector<double>;
  std::vector<double> *vdq_vec = new std::vector<double>;
  std::vector<double> *wdq_vec = new std::vector<double>;

  std::vector<int> *nwire_u_vec = new  std::vector<int>;
  std::vector<int> *nwire_v_vec = new  std::vector<int>;
  std::vector<int> *nwire_w_vec = new  std::vector<int>;
  std::vector<int> *flag_u_vec = new  std::vector<int>;
  std::vector<int> *flag_v_vec = new  std::vector<int>;
  std::vector<int> *flag_w_vec = new  std::vector<int>;

  std::vector<std::vector<int>> *wire_index_u_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_v_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_w_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<double>> *wire_charge_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_w_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_w_vec = new std::vector<std::vector<double>>;
  
  TC->SetBranchAddress("cluster_id",&cluster_id_vec);
  TC->SetBranchAddress("parent_cluster_id",&parent_cluster_id);
  TC->SetBranchAddress("time_slice",&time_slice_vec);
  TC->SetBranchAddress("q",&q_vec);
  TC->SetBranchAddress("uq",&uq_vec);
  TC->SetBranchAddress("vq",&vq_vec);
  TC->SetBranchAddress("wq",&wq_vec);
  TC->SetBranchAddress("udq",&udq_vec);
  TC->SetBranchAddress("vdq",&vdq_vec);
  TC->SetBranchAddress("wdq",&wdq_vec);
  TC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TC->SetBranchAddress("flag_u",&flag_u_vec);
  TC->SetBranchAddress("flag_v",&flag_v_vec);
  TC->SetBranchAddress("flag_w",&flag_w_vec);
  TC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TC->SetBranchAddress("wire_index_w",&wire_index_w_vec);
  TC->SetBranchAddress("wire_charge_u",&wire_charge_u_vec);
  TC->SetBranchAddress("wire_charge_v",&wire_charge_v_vec);
  TC->SetBranchAddress("wire_charge_w",&wire_charge_w_vec);
  TC->SetBranchAddress("wire_charge_err_u",&wire_charge_err_u_vec);
  TC->SetBranchAddress("wire_charge_err_v",&wire_charge_err_v_vec);
  TC->SetBranchAddress("wire_charge_err_w",&wire_charge_err_w_vec);

  //load mcell
  TTree *TDC = (TTree*)file->Get("TDC");
  std::vector<int> *ntime_slice_vec = new std::vector<int>;
  std::vector<std::vector<int>> *time_slices_vec = new std::vector<std::vector<int>>;
  TDC->SetBranchAddress("cluster_id",&cluster_id_vec);
  TDC->SetBranchAddress("ntime_slice",&ntime_slice_vec);
  TDC->SetBranchAddress("time_slice",&time_slices_vec);

  TDC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TDC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TDC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TDC->SetBranchAddress("flag_u",&flag_u_vec);
  TDC->SetBranchAddress("flag_v",&flag_v_vec);
  TDC->SetBranchAddress("flag_w",&flag_w_vec);
  TDC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TDC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TDC->SetBranchAddress("wire_index_w",&wire_index_w_vec);

  // load cells ... 
  GeomCellSelection mcells;
  //  CellIndexMap map_mcell_cluster_id;
  WireCellPID::PR3DClusterSelection live_clusters;
  WireCellPID::PR3DCluster *cluster;
  std::map<WireCellPID::PR3DCluster*, int> map_cluster_parent_id;
  std::map<int, std::vector<WireCellPID::PR3DCluster*> > map_parentid_clusters;
  int prev_cluster_id=-1;
  int ident = 0;
  TC->GetEntry(0);
  for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    int time_slice = time_slice_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);
    std::vector<double> wire_charge_u = wire_charge_u_vec->at(i);
    std::vector<double> wire_charge_v = wire_charge_v_vec->at(i);
    std::vector<double> wire_charge_w = wire_charge_w_vec->at(i);
    std::vector<double> wire_charge_err_u = wire_charge_err_u_vec->at(i);
    std::vector<double> wire_charge_err_v = wire_charge_err_v_vec->at(i);
    std::vector<double> wire_charge_err_w = wire_charge_err_w_vec->at(i);

    
    
    mcell->SetTimeSlice(time_slice);
    
    mcell->set_uq(uq_vec->at(i));
    mcell->set_vq(vq_vec->at(i));
    mcell->set_wq(wq_vec->at(i));
    
    mcell->set_udq(udq_vec->at(i));
    mcell->set_vdq(vdq_vec->at(i));
    mcell->set_wdq(wdq_vec->at(i));
    
    mcell->set_q(q_vec->at(i));
    
    double temp_x = (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    
    if (flag_u_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      for (int j=0;j!=nwire_u_vec->at(i);j++){
	if (dead_u_index.find(wire_index_u.at(j))==dead_u_index.end()){
	  dead_u_index[wire_index_u.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_u_index[wire_index_u.at(j)].first){
	    dead_u_index[wire_index_u.at(j)].first = temp_x - 0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_u_index[wire_index_u.at(j)].second){
	    dead_u_index[wire_index_u.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
	
      }
    }
    if (flag_v_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      for (int j=0;j!=nwire_v_vec->at(i);j++){
	if (dead_v_index.find(wire_index_v.at(j))==dead_v_index.end()){
	  dead_v_index[wire_index_v.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_v_index[wire_index_v.at(j)].first){
	    dead_v_index[wire_index_v.at(j)].first = temp_x-0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_v_index[wire_index_v.at(j)].second){
	    dead_v_index[wire_index_v.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_w_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      for (int j=0;j!=nwire_w_vec->at(i);j++){
	//	if (wire_index_w.at(j)>350&&wire_index_w.at(j)<500) 
	//  std::cout << wire_index_w.at(j) << std::endl;
	if (dead_w_index.find(wire_index_w.at(j))==dead_w_index.end()){
	  dead_w_index[wire_index_w.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_w_index[wire_index_w.at(j)].first){
	    dead_w_index[wire_index_w.at(j)].first = temp_x-0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_w_index[wire_index_w.at(j)].second){
	    dead_w_index[wire_index_w.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
      }
    }
    for (int j=0;j!=nwire_u_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(j));
      mcell->AddWire(wire,WirePlaneType_t(0),wire_charge_u.at(j),wire_charge_err_u.at(j));
    }
    for (int j=0;j!=nwire_v_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(j));
      mcell->AddWire(wire,WirePlaneType_t(1),wire_charge_v.at(j),wire_charge_err_v.at(j));
    }
    for (int j=0;j!=nwire_w_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(j));
      mcell->AddWire(wire,WirePlaneType_t(2),wire_charge_w.at(j),wire_charge_err_w.at(j));
    }
    mcells.push_back(mcell);
    
    //    map_mcell_cluster_id[mcell]=cluster_id;

    if (cluster_id != prev_cluster_id){
      cluster = new WireCellPID::PR3DCluster(cluster_id);
      map_cluster_parent_id[cluster] = parent_cluster_id->at(i);
      live_clusters.push_back(cluster);
      if (map_parentid_clusters.find(parent_cluster_id->at(i)) == map_parentid_clusters.end()){
	std::vector<WireCellPID::PR3DCluster*> temp_clusters;
	temp_clusters.push_back(cluster);
	map_parentid_clusters[parent_cluster_id->at(i)] = temp_clusters;
      }else{
	map_parentid_clusters[parent_cluster_id->at(i)].push_back(cluster);
      }
    }
    cluster->AddCell(mcell,time_slice);


    
    prev_cluster_id = cluster_id;
    ident++;
  }
  
  prev_cluster_id = -1;
  cluster_id_vec->clear();
  wire_index_u_vec->clear();
  wire_index_v_vec->clear();
  wire_index_w_vec->clear();
  nwire_u_vec->clear();
  nwire_v_vec->clear();
  nwire_w_vec->clear();
  flag_u_vec->clear();
  flag_v_vec->clear();
  flag_w_vec->clear();

  TDC->GetEntry(0);
  for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    std::vector<int> time_slices = time_slices_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);
    
    
    double temp_x1 = (time_slices.front()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    double temp_x2 = (time_slices.back()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    // std::cout << temp_x1/units::cm << " " << temp_x2/units::cm << std::endl;
    
    if (flag_u_vec->at(i)==0){
    
      for (int j=0;j!=nwire_u_vec->at(i);j++){
	if (dead_u_index.find(wire_index_u.at(j))==dead_u_index.end()){
	  dead_u_index[wire_index_u.at(j)] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_u_index[wire_index_u.at(j)].first){
	    dead_u_index[wire_index_u.at(j)].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_u_index[wire_index_u.at(j)].second){
	    dead_u_index[wire_index_u.at(j)].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_v_vec->at(i)==0){
    
      for (int j=0;j!=nwire_v_vec->at(i);j++){
	if (dead_v_index.find(wire_index_v.at(j))==dead_v_index.end()){
	  dead_v_index[wire_index_v.at(j)] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_v_index[wire_index_v.at(j)].first){
	    dead_v_index[wire_index_v.at(j)].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_v_index[wire_index_v.at(j)].second){
	    dead_v_index[wire_index_v.at(j)].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_w_vec->at(i)==0){
    
      for (int j=0;j!=nwire_w_vec->at(i);j++){
	if (dead_w_index.find(wire_index_w.at(j))==dead_w_index.end()){
	  dead_w_index[wire_index_w.at(j)] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_w_index[wire_index_w.at(j)].first){
	    dead_w_index[wire_index_w.at(j)].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_w_index[wire_index_w.at(j)].second){
	    dead_w_index[wire_index_w.at(j)].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
    
    ident++;
  }
 
  
  // veto 16 channels in U ...
  for (int i=2080; i!=2096;i++){
    if (dead_u_index.find(i)==dead_u_index.end()){
      dead_u_index[i] = std::make_pair((0*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm-0.1*units::cm, (2400*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm+0.1*units::cm);
    }
  }

  // Load T_ch_bad tree ...
  TTree *T_bad_ch = (TTree*)file->Get("T_bad_ch");
  if (T_bad_ch!=0){
    Int_t chid, plane;
    Int_t start_time,end_time;
    T_bad_ch->SetBranchAddress("chid",&chid);
    T_bad_ch->SetBranchAddress("plane",&plane);
    T_bad_ch->SetBranchAddress("start_time",&start_time);
    T_bad_ch->SetBranchAddress("end_time",&end_time);
    
    for (int i=0;i!=T_bad_ch->GetEntries();i++){
      T_bad_ch->GetEntry(i);
      double temp_x1 = (start_time/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
      double temp_x2 = (end_time/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
      if (plane==1){
	chid -= 2400;
      }else if (plane==2){
	chid -=4800;
      }
      if (plane==0){
	if (dead_u_index.find(chid)==dead_u_index.end()){
	  dead_u_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_u_index[chid].first){
	    dead_u_index[chid].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_u_index[chid].second){
	    dead_u_index[chid].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }else if (plane==1){
	if (dead_v_index.find(chid)==dead_v_index.end()){
	  dead_v_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	  //std::cout << plane << " " << chid << std::endl;
	}else{
	  if (temp_x1-0.1*units::cm < dead_v_index[chid].first){
	    dead_v_index[chid].first = temp_x1 - 0.1*units::cm;
	    //std::cout << plane << " a " << chid << std::endl;
	  }else if (temp_x2+0.1*units::cm > dead_v_index[chid].second){
	    dead_v_index[chid].second = temp_x2 + 0.1*units::cm;
	    //std::cout << plane << " b " << chid << std::endl;
	  }
	}
      }else if (plane==2){
	if (dead_w_index.find(chid)==dead_w_index.end()){
	  dead_w_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	  //std::cout << plane << " " << chid << std::endl;
	}else{
	  if (temp_x1-0.1*units::cm < dead_w_index[chid].first){
	    dead_w_index[chid].first = temp_x1 - 0.1*units::cm;
	    //std::cout << plane << " a " << chid << std::endl;
	  }else if (temp_x2+0.1*units::cm > dead_w_index[chid].second){
	    dead_w_index[chid].second = temp_x2 + 0.1*units::cm;
	    //std::cout << plane << " b " << chid << std::endl;
	  }
	}
      }
    }
    
  }
  cout << em("load clusters from file") << endl;

  
  // form a global map with the current map information
  std::map<int,std::map<const GeomWire*, SMGCSelection > > global_wc_map;
  for (size_t i=0; i!=live_clusters.size();i++){
    WireCellPID::PR3DCluster *cluster = live_clusters.at(i);
    SMGCSelection& mcells = cluster->get_mcells();
    for (auto it = mcells.begin(); it!= mcells.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int time_slice = mcell->GetTimeSlice();
      if (global_wc_map.find(time_slice)==global_wc_map.end()){
	std::map<const GeomWire*, SMGCSelection> temp_wc_map;
	global_wc_map[time_slice] = temp_wc_map;
      }
      std::map<const GeomWire*, SMGCSelection>& timeslice_wc_map = global_wc_map[time_slice];
      
      GeomWireSelection& uwires = mcell->get_uwires();
      GeomWireSelection& vwires = mcell->get_vwires();
      GeomWireSelection& wwires = mcell->get_wwires();
      std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
      if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(0))==bad_planes.end()){
	for (int j=0;j!=uwires.size();j++){
	  const GeomWire *wire = uwires.at(j);
	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
	    SMGCSelection temp_mcells;
	    temp_mcells.push_back(mcell);
	    timeslice_wc_map[wire] = temp_mcells;
	  }else{
	    timeslice_wc_map[wire].push_back(mcell);
	  }
	}
      }
      if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(1))==bad_planes.end()){
	for (int j=0;j!=vwires.size();j++){
	  const GeomWire *wire = vwires.at(j);
	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
	    SMGCSelection temp_mcells;
	    temp_mcells.push_back(mcell);
	     timeslice_wc_map[wire] = temp_mcells;
	  }else{
	    timeslice_wc_map[wire].push_back(mcell);
	  }
	}
      }
      if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(2))==bad_planes.end()){
	for (int j=0;j!=wwires.size();j++){
	  const GeomWire *wire = wwires.at(j);
	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
	    SMGCSelection temp_mcells;
	    temp_mcells.push_back(mcell);
	    timeslice_wc_map[wire] = temp_mcells;
	  }else{
	    timeslice_wc_map[wire].push_back(mcell);
	  }
	}
      }
    }
  }
  
  
  // replace by the new sampling points ...
  for (size_t i=0; i!=live_clusters.size();i++){
    WireCellPID::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
    live_clusters.at(i)->Create_point_cloud();
  }
  cout << em("Add X, Y, Z points") << std::endl;
  
  double first_t_dis = live_clusters.at(0)->get_mcells().front()->GetTimeSlice()*time_slice_width - live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().x;
  double offset_t = first_t_dis/time_slice_width;
  
  ToyCTPointCloud ct_point_cloud(0,2399,2400,4799,4800,8255, // channel range
				 offset_t, -first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w, // offset
				 1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
				 angle_u,angle_v,angle_w// angle
				 );
  ct_point_cloud.AddPoints(timesliceId,timesliceChannel,raw_charge,raw_charge_err);
  ct_point_cloud.AddDeadChs(dead_u_index, dead_v_index, dead_w_index);
  ct_point_cloud.build_kdtree_index();
  
  std::map<WireCellPID::PR3DCluster*, WireCellPID::PR3DCluster*> old_new_cluster_map;
  
  //std::cout << saved_parent_tpc_cluster_ids.size() << std::endl;
  for (size_t i=0; i!=live_clusters.size();i++){
    // if (live_clusters.at(i)->get_cluster_id() != 34
    	//     	&& live_clusters.at(i)->get_cluster_id() != 2
     	//&& live_clusters.at(i)->get_cluster_id() != 6
     	//&& live_clusters.at(i)->get_cluster_id() != 21 
    	// 	&& live_clusters.at(i)->get_cluster_id() != 80 
    // ) continue;
    
    if (live_clusters.at(i)->get_num_points() <= 2) continue;

    // no matched flash 
    if (map_tpc_flash_ids.find(map_cluster_parent_id[live_clusters.at(i)]) == map_tpc_flash_ids.end()) continue;
    double flash_time = map_flash_info[map_tpc_flash_ids[map_cluster_parent_id[live_clusters.at(i)]]].second;

    if (flag_in_time_only && (flash_time < lowerwindow || flash_time > upperwindow)) continue;

    std::cout << live_clusters.at(i)->get_cluster_id() << " " << map_cluster_parent_id[live_clusters.at(i)] << " " << flash_time << std::endl;

    if (flag_main_cluster_only){
      if (live_clusters.at(i)->get_cluster_id() == map_cluster_parent_id[live_clusters.at(i)]){
	live_clusters.at(i)->create_steiner_graph(ct_point_cloud, gds, nrebin, frame_length, unit_dis);
      live_clusters.at(i)->recover_steiner_graph();
      }
    }else{
      live_clusters.at(i)->create_steiner_graph(ct_point_cloud, gds, nrebin, frame_length, unit_dis);
      live_clusters.at(i)->recover_steiner_graph();
    }
  }
  cout << em("Build graph for all clusters") << std::endl;
  
  
  TFile *file1 = new TFile(Form("tracking_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  Trun->CloneTree(-1,"fast");
  if (T_bad_ch!=0){
     T_bad_ch->CloneTree(-1,"fast");
  }
  if (flag_bee_output){
    T_match->CloneTree(-1,"fast");
    T_flash->CloneTree(-1,"fast");
  }
  
  TTree *T_cluster ;
  Double_t x,y,z,q,nq;
  Int_t ncluster;
  Int_t temp_time_slice, ch_u, ch_v, ch_w;
  T_cluster = new TTree("T_cluster","T_cluster");
  T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
  T_cluster->Branch("x",&x,"x/D");
  T_cluster->Branch("y",&y,"y/D");
  T_cluster->Branch("z",&z,"z/D");
  T_cluster->Branch("q",&q,"q/D");
  T_cluster->Branch("nq",&nq,"nq/D");
  T_cluster->Branch("time_slice",&temp_time_slice,"time_slice/I");
  T_cluster->Branch("ch_u",&ch_u,"ch_u/I");
  T_cluster->Branch("ch_v",&ch_v,"ch_v/I");
  T_cluster->Branch("ch_w",&ch_w,"ch_w/I");
  T_cluster->SetDirectory(file1);
  
  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){
    // if (old_new_cluster_map.find(*it)==old_new_cluster_map.end()) continue;
    // WireCellPID::PR3DCluster* new_cluster = old_new_cluster_map[*it];

    WireCellPID::PR3DCluster* new_cluster = *it;  
    ncluster = map_cluster_parent_id[new_cluster]; //new_cluster->get_cluster_id();
    
    ToyPointCloud *pcloud = new_cluster->get_point_cloud();
    //ToyPointCloud *pcloud = new_cluster->get_point_cloud_steiner();
    //ToyPointCloud *pcloud = new_cluster->get_point_cloud_steiner_terminal();

    if (pcloud!=0){
      WireCell::WCPointCloud<double>& cloud = pcloud->get_cloud();
      
      for (size_t i=0;i!=cloud.pts.size();i++){
	x = cloud.pts[i].x/units::cm;
	y = cloud.pts[i].y/units::cm;
	z = cloud.pts[i].z/units::cm;
	SlimMergeGeomCell *mcell = cloud.pts[i].mcell;
	ch_u = cloud.pts[i].index_u;
	ch_v = cloud.pts[i].index_v;
	ch_w = cloud.pts[i].index_w;

	if (mcell==0){
	  temp_time_slice = -1;
	  q = 1;
	  nq = 1;
	}else{
	  temp_time_slice = mcell->GetTimeSlice();
	  nq = 1;
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),ch_u);
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),ch_v);
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),ch_w);
	  q = (mcell->Get_Wire_Charge(uwire) + mcell->Get_Wire_Charge(vwire) + mcell->Get_Wire_Charge(wwire))/3.;
	}
	T_cluster->Fill();
      }
    }
      
    // SMGCSelection& mcells = new_cluster->get_mcells();
    // for (size_t i=0;i!=mcells.size();i++){
    //   SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)mcells.at(i);
    //   PointVector ps = mcell->get_sampling_points();
    //   int time_slice = mcell->GetTimeSlice();
    //   if (ps.size()==0){
    // 	std::cout << "zero sampling points!" << std::endl;
    //   }else{
    // 	q = mcell->get_q() / ps.size();
    // 	nq = ps.size();
    // 	double offset_x = 0;
    // 	for (int k=0;k!=ps.size();k++){
    // 	  x = (ps.at(k).x- offset_x)/units::cm ;
    // 	  y = ps.at(k).y/units::cm;
    // 	  z = ps.at(k).z/units::cm;
    // 	  //	std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(ps.at(k));
    // 	  //temp_time_slice = time_chs.at(0);
    // 	  //ch_u = time_chs.at(1);
    // 	  //ch_v = time_chs.at(2);
    // 	  //ch_w = time_chs.at(3);
	  
    // 	  T_cluster->Fill();
    // 	}
    //   }
    // }
    

    
  }

  Double_t pu, pv, pw, pt;
  Double_t charge_save=1, ncharge_save=1, chi2_save=1, ndf_save=1;
  TTree *T_rec = new TTree("T_rec","T_rec");
  T_rec->Branch("x",&x,"x/D");
  T_rec->Branch("y",&y,"y/D");
  T_rec->Branch("z",&z,"z/D");
  T_rec->Branch("q",&charge_save,"q/D");
  T_rec->Branch("nq",&ncharge_save,"nq/D");
  T_rec->Branch("chi2",&chi2_save,"chi2/D");
  T_rec->Branch("ndf",&ndf_save,"ndf/D");
  T_rec->Branch("pu",&pu,"pu/D");
  T_rec->Branch("pv",&pv,"pv/D");
  T_rec->Branch("pw",&pw,"pw/D");
  T_rec->Branch("pt",&pt,"pt/D");
  T_rec->SetDirectory(file1);

  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  t_rec_charge->SetDirectory(file1);
  t_rec_charge->Branch("x",&x,"x/D");
  t_rec_charge->Branch("y",&y,"y/D");
  t_rec_charge->Branch("z",&z,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
  t_rec_charge->Branch("pu",&pu,"pu/D");
  t_rec_charge->Branch("pv",&pv,"pv/D");
  t_rec_charge->Branch("pw",&pw,"pw/D");
  t_rec_charge->Branch("pt",&pt,"pt/D");
  
  TTree *T_proj_data = new TTree("T_proj_data","T_proj_data");
  std::vector<int> *proj_data_cluster_id = new std::vector<int>;
  std::vector<std::vector<int>> *proj_data_cluster_channel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_data_cluster_timeslice= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_data_cluster_charge= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_data_cluster_charge_err= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_data_cluster_charge_pred= new std::vector<std::vector<int>>;
  
  T_proj_data->Branch("cluster_id",&proj_data_cluster_id);
  T_proj_data->Branch("channel",&proj_data_cluster_channel);
  T_proj_data->Branch("time_slice",&proj_data_cluster_timeslice);
  T_proj_data->Branch("charge",&proj_data_cluster_charge);
  T_proj_data->Branch("charge_err",&proj_data_cluster_charge_err);
  T_proj_data->Branch("charge_pred",&proj_data_cluster_charge_pred);
  T_proj_data->SetDirectory(file1);
  
  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){
    //     if (old_new_cluster_map.find(*it)==old_new_cluster_map.end()) continue;
   //  WireCellPID::PR3DCluster* new_cluster = old_new_cluster_map[*it];
//     //    new_cluster->establish_same_mcell_steiner_edges(gds, false);
// std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = new_cluster->get_two_boundary_wcps(1);
//     new_cluster->dijkstra_shortest_paths(wcps.first,1); 
//     new_cluster->cal_shortest_path(wcps.second,1);
//     //new_cluster->remove_same_mcell_steiner_edges();
    
    WireCellPID::PR3DCluster* new_cluster = *it;
    if (new_cluster->get_point_cloud_steiner()==0) continue;

    if (new_cluster->get_point_cloud_steiner()->get_num_points() >= 2){
      std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = new_cluster->get_two_boundary_wcps(2); 
      // std::cout << wcps.first.x/units::cm << " " << wcps.first.y/units::cm << " " << wcps.first.z/units::cm << std::endl;
      // std::cout << wcps.second.x/units::cm << " " << wcps.second.y/units::cm << " " << wcps.second.z/units::cm << std::endl;
      new_cluster->dijkstra_shortest_paths(wcps.first,2); 
      new_cluster->cal_shortest_path(wcps.second,2);
    }
    new_cluster->collect_charge_trajectory(ct_point_cloud);


    double flash_time = map_flash_info[map_tpc_flash_ids[map_cluster_parent_id[*it]]].second;
    
    std::cout << new_cluster->get_cluster_id() << std::endl;
    
    new_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
    
    
    ndf_save = new_cluster->get_cluster_id();
    charge_save = 0;
    ncharge_save = 0;
    chi2_save = 0;
    
    std::list<WCPointCloud<double>::WCPoint>& wcps_list = new_cluster->get_path_wcps();
    for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
       x = (*it).x/units::cm;
       y = (*it).y/units::cm;
       z = (*it).z/units::cm;

       Point temp_p((*it).x,(*it).y,(*it).z);
       std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
       pt = time_chs.at(0);
       pu = time_chs.at(1);
       pv = time_chs.at(2);
       pw = time_chs.at(3);
       
       T_rec->Fill();
     }
    
  }

   cout << em("shortest path ...") << std::endl;

  
  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){
    WireCellPID::PR3DCluster* cluster = *it;

    ndf_save = cluster->get_cluster_id();

    // original
    PointVector& pts = cluster->get_fine_tracking_path();
    //std::vector<double> dQ, dx;
    std::vector<double>& dQ = cluster->get_dQ();
    std::vector<double>& dx = cluster->get_dx();
    std::vector<double>& tpu = cluster->get_pu();
    std::vector<double>& tpv = cluster->get_pv();
    std::vector<double>& tpw = cluster->get_pw();
    std::vector<double>& tpt = cluster->get_pt();
    //hack for now 
    // for (auto it = pts.begin(); it!=pts.end(); it++){
    //   dQ.push_back(0);
    //   dx.push_back(1);
    //   // std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(*it);
    //   // tpt.push_back(time_chs.at(0));
    //   // tpu.push_back(time_chs.at(1));
    //   // tpv.push_back(time_chs.at(2));
    //   // tpw.push_back(time_chs.at(3));
    // }
    
    //hack for now ...
    //std::list<WCPointCloud<double>::WCPoint>& wcps_list = cluster->get_path_wcps();
    //PointVector pts;
    //std::vector<double> dQ, dx, tpu, tpv, tpw, tpt;
    
    // for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
    //   Point p((*it).x, (*it).y,  (*it).z);
    //   pts.push_back(p);
    //   dQ.push_back(0);
    //   dx.push_back(1);
    //   std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(p);
    //   tpt.push_back(time_chs.at(0));
    //   tpu.push_back(time_chs.at(1));
    //   tpv.push_back(time_chs.at(2));
    //   tpw.push_back(time_chs.at(3));
    // }
    
    if (pts.size()!=dQ.size() || pts.size()==0) continue;
    
    for (size_t i=0; i!=pts.size(); i++){
      x = pts.at(i).x/units::cm;
      y = pts.at(i).y/units::cm;
      z = pts.at(i).z/units::cm;
      charge_save = dQ.at(i);
      ncharge_save = dx.at(i)/units::cm;
      pu = tpu.at(i);
      pv = tpv.at(i);
      pw = tpw.at(i);
      pt = tpt.at(i);
      t_rec_charge->Fill();
    }
    
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & proj_data_u_map = cluster->get_proj_data_u_map();
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & proj_data_v_map = cluster->get_proj_data_v_map();
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & proj_data_w_map = cluster->get_proj_data_w_map();
    
    proj_data_cluster_id->push_back(ndf_save);
    std::vector<int> temp_channel;
    std::vector<int> temp_timeslice;
    std::vector<int> temp_charge;
    std::vector<int> temp_charge_err;
    std::vector<int> temp_charge_pred;
    for (auto it = proj_data_u_map.begin(); it!=proj_data_u_map.end(); it++){
      temp_channel.push_back(it->first.first);
      temp_timeslice.push_back(it->first.second);
      temp_charge.push_back(std::get<0>(it->second));
      temp_charge_err.push_back(std::get<1>(it->second));
      temp_charge_pred.push_back(std::get<2>(it->second));
    }
    for (auto it = proj_data_v_map.begin(); it!=proj_data_v_map.end(); it++){
      temp_channel.push_back(it->first.first);
      temp_timeslice.push_back(it->first.second);
      temp_charge.push_back(std::get<0>(it->second));
      temp_charge_err.push_back(std::get<1>(it->second));
      temp_charge_pred.push_back(std::get<2>(it->second));
    }
    for (auto it = proj_data_w_map.begin(); it!=proj_data_w_map.end(); it++){
      temp_channel.push_back(it->first.first);
      temp_timeslice.push_back(it->first.second);
      temp_charge.push_back(std::get<0>(it->second));
      temp_charge_err.push_back(std::get<1>(it->second));
      temp_charge_pred.push_back(std::get<2>(it->second));
    }
    proj_data_cluster_channel->push_back(temp_channel);
    proj_data_cluster_timeslice->push_back(temp_timeslice);
    proj_data_cluster_charge->push_back(temp_charge);
    proj_data_cluster_charge_err->push_back(temp_charge_err);
    proj_data_cluster_charge_pred->push_back(temp_charge_pred);

  }
  T_proj_data->Fill();

  
  // now save the original projected charge information
  // fill the bad channels ...
  TTree *T_proj = new TTree("T_proj","T_proj");
  std::vector<int> *proj_cluster_id = new std::vector<int>;
  std::vector<std::vector<int>> *proj_cluster_channel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_timeslice= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge_err= new std::vector<std::vector<int>>;
  T_proj->Branch("cluster_id",&proj_cluster_id);
  T_proj->Branch("channel",&proj_cluster_channel);
  T_proj->Branch("time_slice",&proj_cluster_timeslice);
  T_proj->Branch("charge",&proj_cluster_charge);
  T_proj->Branch("charge_err",&proj_cluster_charge_err);
  T_proj->SetDirectory(file1);

  for (auto it = map_parentid_clusters.begin(); it!=map_parentid_clusters.end(); it++){
    int cluster_id = it->first;
    std::vector<int> proj_channel;
    std::vector<int> proj_timeslice;
    std::vector<int> proj_charge;
    std::vector<int> proj_charge_err;
    std::vector<int> proj_flag;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      WireCellPID::PR3DCluster *cluster = (*it1);
      cluster->get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);
    }
    proj_cluster_id->push_back(cluster_id);
    proj_cluster_channel->push_back(proj_channel);
    proj_cluster_timeslice->push_back(proj_timeslice);
    proj_cluster_charge->push_back(proj_charge);
    proj_cluster_charge_err->push_back(proj_charge_err);
  }
  
  T_proj->Fill();
  
  file1->Write();
  file1->Close();
  
  
  return 0;
}
