#include "WCPSst/GeomDataSource.h"

#include "WCPData/SlimMergeGeomCell.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"
#include "WCPData/ToyCTPointCloud.h"

#include "WCPPID/ToyFiducial.h"
#include "WCPPID/CalcPoints.h"
#include "WCPPID/ProtectOverClustering.h"

#include "WCPPID/ExecMon.h"
#include "WCPPID/NeutrinoID.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-prod-nue /path/to/ChannelWireGeometry.txt /path/to/matching.root #entry -d[data=0, overlay=1, full mc=2]" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  

  bool flag_debug_output = true; // output
  int datatier = 0; // data=0, overlay=1, full mc=2
  int flag_calib_corr = 1;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'd':
      datatier = atoi(&argv[i][2]);
      break;
    case 'q':
      flag_calib_corr = atoi(&argv[i][2]);
      break;
    }
  }

  int flag_data = 1; // data
  if(datatier==1 || datatier==2) flag_data=0; // overlay, full mc
  
  WCPPID::ExecMon em("starting");
  cout << em("load geometry") << endl;
  
  WCPSst::GeomDataSource gds(argv[1]);
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
  std::cout << argv[2] << " " << argv[3] << std::endl;
  TString filename = argv[2];
  int entry_no = atoi(argv[3]);
  
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");

  if (entry_no >=Trun->GetEntries()) return 0;
  
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
    
  Trun->GetEntry(entry_no);
  double lowerwindow = 0;
  double upperwindow = 0;
  if((triggerbits>>11) & 1U) { lowerwindow=3.1875; upperwindow=4.96876;} // bnb  
  if((triggerbits>>9) & 1U) { lowerwindow=3.5625; upperwindow=5.34376; } // extbnb
  
  
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
  if (flag_calib_corr==1)
    mp.init_corr_files();

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
  Double_t low_time, high_time;
  Int_t type;
  Int_t flash_id;
  Int_t temp_run_no, temp_subrun_no, temp_event_no;
  Double_t total_PE;
  Double_t PE[32],PE_err[32];
  // std::vector<int> matched_tpc_ids;
  
  T_flash->SetBranchAddress("runNo",&temp_run_no);
  T_flash->SetBranchAddress("subRunNo",&temp_subrun_no);
  T_flash->SetBranchAddress("eventNo",&temp_event_no);
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("flash_id",&flash_id);
  T_flash->SetBranchAddress("low_time",&low_time);
  T_flash->SetBranchAddress("high_time",&high_time);
  T_flash->SetBranchAddress("total_PE",&total_PE);
  T_flash->SetBranchAddress("PE",PE);
  T_flash->SetBranchAddress("PE_err",PE_err);
  std::vector<int> *fired_channels = new std::vector<int>;
  std::vector<double> *l1_fired_time = new std::vector<double>;
  std::vector<double> *l1_fired_pe = new std::vector<double>;
  T_flash->SetBranchAddress("fired_channels",&fired_channels);
  T_flash->SetBranchAddress("l1_fired_time",&l1_fired_time);
  T_flash->SetBranchAddress("l1_fired_pe",&l1_fired_pe);
  
  TTree *T_match = (TTree*)file->Get("T_match");
  T_match->SetBranchAddress("runNo",&temp_run_no);
  T_match->SetBranchAddress("subRunNo",&temp_subrun_no);
  T_match->SetBranchAddress("eventNo",&temp_event_no);
  Int_t tpc_cluster_id;
  Int_t event_type;
  T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id); // parent cluster id
  T_match->SetBranchAddress("flash_id",&flash_id);  // flash id 
  T_match->SetBranchAddress("event_type",&event_type);

  Double_t strength;
  Double_t pe_pred[32];
  Double_t pe_meas[32];
  Double_t pe_meas_err[32];
  T_match->SetBranchAddress("strength",&strength);
  T_match->SetBranchAddress("pe_pred",pe_pred);
  T_match->SetBranchAddress("pe_meas",pe_meas);
  T_match->SetBranchAddress("pe_meas_err",pe_meas_err);
  T_match->SetBranchAddress("event_type",&event_type);
  bool flag_close_to_PMT;
  bool flag_at_x_boundary;
  double ks_dis;
  double chi2;
  int ndf;
  double cluster_length;
  
  T_match->SetBranchAddress("ks_dis",&ks_dis);
  T_match->SetBranchAddress("chi2",&chi2);
  T_match->SetBranchAddress("ndf",&ndf);
  T_match->SetBranchAddress("cluster_length",&cluster_length);
  
  
  std::map<int, std::pair<int, double> > map_flash_info;
  std::map<int, int> map_flash_tpc_ids;
  std::map<int, int> map_tpc_flash_ids;
  std::map<std::pair<int, int>, int> map_flash_tpc_pair_type;
  
  
  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    if (temp_run_no!=run_no || temp_subrun_no!=subrun_no || temp_event_no != event_no) continue;
    map_flash_info[flash_id] = std::make_pair(type, time);
  }
  for (int i=0;i!=T_match->GetEntries();i++){
    T_match->GetEntry(i);
    if (temp_run_no!=run_no || temp_subrun_no!=subrun_no || temp_event_no != event_no) continue;
    if (flash_id==-1) continue;
    map_flash_tpc_ids[flash_id] = tpc_cluster_id;
    map_tpc_flash_ids[tpc_cluster_id] = flash_id;
    map_flash_tpc_pair_type[std::make_pair(flash_id, tpc_cluster_id)] = event_type;
  }
  // no need to reconstruct flash ...
  

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

  WCPPID::ToyFiducial *fid = new WCPPID::ToyFiducial(3,800,-first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w,
								   1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
								   angle_u,angle_v,angle_w,// angle
								   3*units::cm, 117*units::cm, -116*units::cm, 0*units::cm, 1037*units::cm, 0*units::cm, 256*units::cm, flag_data);
  
  // load cells ... 
  GeomCellSelection mcells;
  //  CellIndexMap map_mcell_cluster_id;
  WCPPID::PR3DClusterSelection live_clusters;
  WCPPID::PR3DCluster *cluster;
  std::map<WCPPID::PR3DCluster*, int> map_cluster_parent_id; // cluster to main cluster
  std::map<int, std::vector<WCPPID::PR3DCluster*> > map_parentid_clusters; // main cluster to clusters

  int prev_cluster_id=-1;
  int ident = 0;
  TC->GetEntry(entry_no);
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
    
    if (cluster_id != prev_cluster_id){
      cluster = new WCPPID::PR3DCluster(cluster_id);
      map_cluster_parent_id[cluster] = parent_cluster_id->at(i);
      live_clusters.push_back(cluster);
      if (map_parentid_clusters.find(parent_cluster_id->at(i)) == map_parentid_clusters.end()){
	std::vector<WCPPID::PR3DCluster*> temp_clusters;
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

  
  TDC->GetEntry(entry_no);
  for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    std::vector<int> time_slices = time_slices_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);

    mcell->SetTimeSlice(time_slices.at(0));
    
    double temp_x1 = (time_slices.front()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    double temp_x2 = (time_slices.back()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    
    if (flag_u_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
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
      mcell->add_bad_planes(WirePlaneType_t(1));
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
      mcell->add_bad_planes(WirePlaneType_t(2));
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
    for (int j=0;j!=nwire_u_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(j));
      mcell->AddWire(wire,WirePlaneType_t(0));
    }
    for (int j=0;j!=nwire_v_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(j));
      mcell->AddWire(wire,WirePlaneType_t(1));
    }
    for (int j=0;j!=nwire_w_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(j));
      mcell->AddWire(wire,WirePlaneType_t(2));
    }

    // figure out the boundary here ...
    PointVector pcell;
    if (flag_u_vec->at(i)==0 && flag_v_vec->at(i)==0){
      const GeomWire *uwire_1 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.front());
      const GeomWire *uwire_2 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.back());
      float dis_u[3];
      float u_pitch = gds.pitch(kUwire);
      dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
      dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;

      const GeomWire *vwire_1 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.front());
      const GeomWire *vwire_2 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.back());
      float dis_v[3];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;

      std::vector<Vector> pcross(4);
      bool flag1 = gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, pcross[0]); // check the inner point
      
      if (flag1){
	// fill the outer point
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;

	  if (gds.crossing_point(dis_u[2],dis_v[0],kUwire,kVwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(k));
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_v[2],kUwire,kVwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[1],kUwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(wire_index_v.size()-1-k));
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_v[2],kUwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(wire_index_u.size()-1-k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[0],kUwire,kVwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 =  gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(k));
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_v[2],kUwire,kVwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(wire_index_u.size()-1-k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[1],kUwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(wire_index_v.size()-1-k));
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_v[2],kUwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
      
      mcell->AddBoundary(pcell);
      
    }else if (flag_u_vec->at(i)==0 && flag_w_vec->at(i)==0){
      const GeomWire *uwire_1 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.front());
      const GeomWire *uwire_2 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.back());
      float dis_u[2];
      float u_pitch = gds.pitch(kUwire);
      dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
      dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;

      const GeomWire *wwire_1 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.front());
      const GeomWire *wwire_2 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.back());
      float dis_w[2];
      float w_pitch = gds.pitch(kYwire);
      dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
      dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;

      std::vector<Vector> pcross(4);

      bool flag1 = gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[0],kUwire,kYwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(k));
	  dis_w[2] = gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_w[2],kUwire,kYwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[1],kUwire,kYwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(wire_index_w.size()-1-k));
	  dis_w[2] = gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_w[2],kUwire,kYwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 =gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(wire_index_u.size()-1-k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[0],kUwire,kYwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 =gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(k));
	  dis_w[2] = gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_w[2],kUwire,kYwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=wire_index_u.size();k++){
	  const GeomWire *uwire_3 =gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(wire_index_u.size()-1-k));
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[1],kUwire,kYwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(wire_index_w.size()-1-k));
	  dis_w[2] = gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_w[2],kUwire,kYwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
      mcell->AddBoundary(pcell);
      
    }else if (flag_v_vec->at(i)==0 && flag_w_vec->at(i)==0){
      const GeomWire *wwire_1 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.front());
      const GeomWire *wwire_2 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.back());
      float dis_w[2];
      float w_pitch = gds.pitch(kYwire);
      dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
      dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;
      
      const GeomWire *vwire_1 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.front());
      const GeomWire *vwire_2 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.back());
      float dis_v[2];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;
      
      std::vector<Vector> pcross(4);

       bool flag1 = gds.crossing_point(dis_w[0],dis_v[0],kYwire,kVwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(k));
	  dis_w[2] =  gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[0],kYwire,kVwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(k));
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_w[0],dis_v[2],kYwire,kVwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_w[0],dis_v[1],kYwire,kVwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(k));
	  dis_w[2] =  gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[1],kYwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(wire_index_v.size()-1-k));
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_w[0],dis_v[2],kYwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_w[1],dis_v[0],kYwire,kVwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(wire_index_w.size()-1-k));
	  dis_w[2] =  gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[0],kYwire,kVwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(k));
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_w[1],dis_v[2],kYwire,kVwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_w[1],dis_v[1],kYwire,kVwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=wire_index_w.size();k++){
	  const GeomWire *wwire_3 = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(wire_index_w.size()-1-k));
	  dis_w[2] =  gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[1],kYwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=wire_index_v.size();k++){
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(wire_index_v.size()-1-k));
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_w[1],dis_v[2],kYwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
      mcell->AddBoundary(pcell);
    }
    // figure out the boundary ...
    
    fid->AddDeadRegion(mcell,time_slices);
    ident++;
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
    Int_t temp_run_no, temp_subrun_no, temp_event_no;
    T_bad_ch->SetBranchAddress("runNo",&temp_run_no);
    T_bad_ch->SetBranchAddress("subRunNo",&temp_subrun_no);
    T_bad_ch->SetBranchAddress("eventNo",&temp_event_no);
    
    for (int i=0;i!=T_bad_ch->GetEntries();i++){
      T_bad_ch->GetEntry(i);
      if (temp_run_no != run_no || temp_subrun_no!=subrun_no || temp_event_no!=event_no) continue;
	  
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
	}else{
	  if (temp_x1-0.1*units::cm < dead_v_index[chid].first){
	    dead_v_index[chid].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_v_index[chid].second){
	    dead_v_index[chid].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }else if (plane==2){
	if (dead_w_index.find(chid)==dead_w_index.end()){
	  dead_w_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_w_index[chid].first){
	    dead_w_index[chid].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_w_index[chid].second){
	    dead_w_index[chid].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
  }
  cout << em("load clusters from file") << endl;
  
  // form a global map with the current map information
  std::map<int,std::map<const GeomWire*, SMGCSelection > > global_wc_map;
  for (size_t i=0; i!=live_clusters.size();i++){
    WCPPID::PR3DCluster *cluster = live_clusters.at(i);
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
    WCPPID::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
    live_clusters.at(i)->Create_point_cloud();
  }
  cout << em("Add X, Y, Z points") << std::endl;

   double first_t_dis = live_clusters.at(0)->get_mcells().front()->GetTimeSlice()*time_slice_width - live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().x;
  double offset_t = first_t_dis/time_slice_width;

  // test the fiducial volume cut 
  fid->set_offset_t(offset_t);
  
  ToyCTPointCloud ct_point_cloud(0,2399,2400,4799,4800,8255, // channel range
				 offset_t, -first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w, // offset
				 1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
				 angle_u,angle_v,angle_w// angle
				 );
  ct_point_cloud.AddPoints(timesliceId,timesliceChannel,raw_charge,raw_charge_err);
  ct_point_cloud.AddDeadChs(dead_u_index, dead_v_index, dead_w_index);
  ct_point_cloud.build_kdtree_index();

  // examine the clustering ... 
  std::vector<int> to_be_checked;
  for (auto it = map_flash_tpc_ids.begin(); it!=map_flash_tpc_ids.end(); it++){
    double flash_time = map_flash_info[it->first].second;
    if ( (flash_time < lowerwindow || flash_time > upperwindow)) continue;
    to_be_checked.push_back(it->second);
  }
  WCPPID::Protect_Over_Clustering(to_be_checked, live_clusters, map_cluster_parent_id, map_parentid_clusters, ct_point_cloud);
  
  for (size_t i=0; i!=live_clusters.size();i++){
    if (live_clusters.at(i)->get_point_cloud()==0){
      WCPPID::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
      live_clusters.at(i)->Create_point_cloud();
    }
  }
  
  
  std::vector<WCPPID::NeutrinoID *> neutrino_vec;
  // Code to select nue ...
  for (auto it = map_flash_tpc_ids.begin(); it!=map_flash_tpc_ids.end(); it++){
    double flash_time = map_flash_info[it->first].second;
    //std::cout << flash_time << " " << triggerbits << " " << lowerwindow << " " << upperwindow << std::endl;
    if ( (flash_time < lowerwindow || flash_time > upperwindow)) continue;
    if (map_parentid_clusters.find(it->second) == map_parentid_clusters.end()) continue;

    // hack for now ...
    //continue;
    
    std::vector<WCPPID::PR3DCluster*> temp_clusters = map_parentid_clusters[it->second];
    WCPPID::PR3DCluster* main_cluster = 0;
    for (auto it1 = temp_clusters.begin(); it1!=temp_clusters.end();it1++){
      if ((*it1)->get_cluster_id() == it->second){
	main_cluster = *it1;
	break;
      }
    }
    std::vector<WCPPID::PR3DCluster*> additional_clusters;
    for (auto it1 = temp_clusters.begin(); it1!=temp_clusters.end();it1++){
      if (*it1 != main_cluster)
	additional_clusters.push_back(*it1);
    }


    
    WCPPID::NeutrinoID *neutrino = new WCPPID::NeutrinoID(main_cluster, additional_clusters, live_clusters, gds, nrebin, frame_length, unit_dis, &ct_point_cloud, global_wc_map, flash_time);
    neutrino_vec.push_back(neutrino);
    
  }

  
  // start saving ...

  

  TFile *file1 = new TFile(Form("nue_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  TTree *Trun1 = new TTree("Trun","Trun");
  Trun1->SetDirectory(file1);
  Trun1->Branch("eventNo",&event_no,"eventNo/I");
  Trun1->Branch("runNo",&run_no,"runNo/I");
  Trun1->Branch("subRunNo",&subrun_no,"subRunNo/I");
  
  Trun1->Branch("triggerBits",&triggerbits,"triggerBits/i");
  Trun1->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun1->Branch("frame_length",&frame_length,"frame_length/I");
  Trun1->Branch("eve_num",&eve_num,"eve_num/I");
  Trun1->Branch("nrebin",&nrebin,"nrebin/I");
  Trun1->Branch("time_offset",&time_offset,"time_offset/I");
  int detector = 0; // MicroBooNE
  Trun1->Branch("detector",&detector,"detector/I");
  // Trun1->Branch("timesliceId",&timesliceId);
  // Trun1->Branch("timesliceChannel",&timesliceChannel);
  // Trun1->Branch("raw_charge",&raw_charge);
  // Trun1->Branch("raw_charge_err",&raw_charge_err);

  Double_t dQdx_scale = 0.1;
  Trun1->Branch("dQdx_scale",&dQdx_scale,"dQdx_scale/D");
  Double_t dQdx_offset = -1000;
  Trun1->Branch("dQdx_offset",&dQdx_offset,"dQdx_offset/D");
  
  Trun1->Fill();

  if (T_bad_ch!=0)
    T_bad_ch->CloneTree(-1,"fast");
  
  TTree *T_flash1 = new TTree("T_flash","T_flash");
  T_flash1->SetDirectory(file1);
  T_flash1->Branch("type",&type);
  T_flash1->Branch("flash_id",&flash_id);
  T_flash1->Branch("low_time",&low_time);
  T_flash1->Branch("high_time",&high_time);
  T_flash1->Branch("time",&time);
  T_flash1->Branch("total_PE",&total_PE);
  T_flash1->Branch("PE",PE,"PE[32]/D");
  T_flash1->Branch("PE_err",PE_err,"PE_err[32]/D");
  T_flash1->Branch("fired_channels",&fired_channels);
  T_flash1->Branch("l1_fired_time",&l1_fired_time);
  T_flash1->Branch("l1_fired_pe",&l1_fired_pe);
  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    if (temp_run_no!=run_no || temp_subrun_no!=subrun_no || temp_event_no != event_no) continue;
    T_flash1->Fill();
  }

  
  
  TTree *T_match1 = new TTree("T_match","T_match");
  T_match1->SetDirectory(file1);
  T_match1->Branch("tpc_cluster_id",&tpc_cluster_id,"tpc_cluster_id/I"); // parent cluster id
  T_match1->Branch("flash_id",&flash_id,"flash_id/I");  // flash id 
  
  T_match1->Branch("strength",&strength,"strength/D");
  T_match1->Branch("pe_pred",pe_pred,"pe_pred[32]/D");
  T_match1->Branch("pe_meas",pe_meas,"pe_meas[32]/D");
  T_match1->Branch("pe_meas_err",pe_meas_err,"pe_meas_err[32]/D");
  T_match1->Branch("event_type",&event_type,"event_type/I");
  T_match1->Branch("ks_dis",&ks_dis,"ks_dis/D");
  T_match1->Branch("chi2",&chi2,"chi2/D");
  T_match1->Branch("ndf",&ndf,"ndf/I");
  T_match1->Branch("cluster_length",&cluster_length,"cluster_length/D");

  for (int i=0;i!=T_match->GetEntries();i++){
    T_match->GetEntry(i);
    if (temp_run_no!=run_no || temp_subrun_no!=subrun_no || temp_event_no != event_no) continue;
    T_match1->Fill();
  }

  TTree *T_vtx = new TTree("T_vtx","T_vtx");
  T_vtx->SetDirectory(file1);
  Double_t x_vtx, y_vtx, z_vtx;
  Int_t type_vtx=0;
  // 1: Steiner-inspired graph
  // 2: Steiner terminals
  // 3: end point (extreme point, connecting to one object)
  // 4: kink (connecting to two objects
  // 5: vertex (connecting to three or more objects)
  Int_t flag_main_vtx=0; // main vertex or not
  Int_t cluster_id_vtx=-1; // -1 if not belonging to any cluster ...
  std::vector<int> *sub_cluster_ids_vtx = new std::vector<int>;
  
  T_vtx->Branch("x",&x_vtx,"x/D");
  T_vtx->Branch("y",&y_vtx,"y/D");
  T_vtx->Branch("z",&z_vtx,"z/D");
  T_vtx->Branch("type",&type_vtx,"type/I");
  T_vtx->Branch("flag_main",&flag_main_vtx,"flag_main/I");
  T_vtx->Branch("cluster_id",&cluster_id_vtx,"cluster_id/I");
  T_vtx->Branch("sub_cluster_ids",&sub_cluster_ids_vtx);
  
  for (size_t i=0;i!=neutrino_vec.size();i++){
    WCPPID::Map_Proto_Vertex_Segments& map_vertex_segments = neutrino_vec.at(i)->get_map_vertex_segments();
    //    WCPPID::Map_Proto_Segment_Vertices& map_segment_vertices = neutrino.get_map_segment_vertices();
    for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end();it++){
      WCPPID::ProtoVertex *vtx = it->first;
      x_vtx = vtx->get_fit_pt().x/units::cm;
      y_vtx = vtx->get_fit_pt().y/units::cm;
      z_vtx = vtx->get_fit_pt().z/units::cm;
      if (it->second.size()==1){
      	type_vtx = 3;
      }else if (it->second.size()==2){
      	type_vtx = 4;
      }else{
      	type_vtx = 5;
      }
      flag_main_vtx = vtx->get_flag_neutrino_vertex();
      cluster_id_vtx = vtx->get_cluster_id();
      sub_cluster_ids_vtx->clear();
      for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      	WCPPID::ProtoSegment *sg = *it1;
      	sub_cluster_ids_vtx->push_back(cluster_id_vtx*1000 + sg->get_id());
      }
      T_vtx->Fill();
    }
  }


  TTree *TMC = new TTree("TMC","TMC");
  
  WCPPID::WCRecoTree reco_tree;
  reco_tree.mc_daughters = new std::vector<std::vector<int> >;
  
  TMC->Branch("mc_Ntrack", &reco_tree.mc_Ntrack);  // number of tracks in MC
  TMC->Branch("mc_id", &reco_tree.mc_id, "mc_id[mc_Ntrack]/I");  // track id; size == mc_Ntrack
  TMC->Branch("mc_pdg", &reco_tree.mc_pdg, "mc_id[mc_Ntrack]/I");  // track particle pdg; size == mc_Ntrack
  TMC->Branch("mc_process", &reco_tree.mc_process, "mc_process[mc_Ntrack]/I");  // track generation process code; size == mc_Ntrack
  TMC->Branch("mc_mother", &reco_tree.mc_mother, "mc_mother[mc_Ntrack]/I");  // mother id of this track; size == mc_Ntrack
  TMC->Branch("mc_daughters", reco_tree.mc_daughters);  // daughters id of this track; vector
  TMC->Branch("mc_startXYZT", &reco_tree.mc_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
  TMC->Branch("mc_endXYZT", &reco_tree.mc_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
  TMC->Branch("mc_startMomentum", &reco_tree.mc_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
  TMC->Branch("mc_endMomentum", &reco_tree.mc_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
  TMC->SetDirectory(file1);

  for (size_t i=0; i!= neutrino_vec.size();i++){
    neutrino_vec.at(i)->fill_reco_simple_tree(reco_tree);
    TMC->Fill();
  }


  
  TTree *T_cluster ;
  Double_t x,y,z,q,nq;
  Int_t ncluster;
  Int_t temp_time_slice, ch_u, ch_v, ch_w;
  Int_t real_cluster_id;
  Int_t sub_cluster_id;
  T_cluster = new TTree("T_cluster","T_cluster");
  T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
  T_cluster->Branch("real_cluster_id",&real_cluster_id,"real_cluster_id/I");
  T_cluster->Branch("sub_cluster_id",&sub_cluster_id,"sub_cluster_id/I");
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
    WCPPID::PR3DCluster* new_cluster = *it;  
    ncluster = map_cluster_parent_id[new_cluster];
    real_cluster_id = new_cluster->get_cluster_id();
    ToyPointCloud *pcloud = new_cluster->get_point_cloud();
          
    if (pcloud!=0){
      WCP::WCPointCloud<double>& cloud = pcloud->get_cloud();
      std::vector<int>& point_sub_cluster_ids = new_cluster->get_point_sub_cluster_ids();
      for (size_t i=0;i!=cloud.pts.size();i++){
	x = cloud.pts[i].x/units::cm;
	y = cloud.pts[i].y/units::cm;
	z = cloud.pts[i].z/units::cm;
	
	if (point_sub_cluster_ids.size() == cloud.pts.size()){
	  if (point_sub_cluster_ids.at(i)==-1){
	    real_cluster_id = -1;
	    continue; // skip -1 points ... deghosting ...
	  }else{
	    real_cluster_id =new_cluster->get_cluster_id()*1000 + point_sub_cluster_ids.at(i);
	  }
	}
	//	real_cluster_id =new_cluster->get_cluster_id();
	
	sub_cluster_id =  new_cluster->get_cluster_id();//*1000;	  
	
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
  }

  Double_t pu, pv, pw, pt;
  Double_t charge_save=1, ncharge_save=1, chi2_save=1, ndf_save=1;
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
  Double_t reduced_chi2;
  t_rec_charge->Branch("reduced_chi2",&reduced_chi2,"reduced_chi2/D");
  Int_t flag_vertex_save;
  t_rec_charge->Branch("flag_vertex",&flag_vertex_save,"flag_vertex/I");
  t_rec_charge->Branch("cluster_id",&ncluster,"cluster_id/I");
  t_rec_charge->Branch("real_cluster_id",&real_cluster_id,"real_cluster_id/I");
  t_rec_charge->Branch("sub_cluster_id",&sub_cluster_id,"sub_cluster_id/I");

 
  // extra save ...
  TTree *t_rec_simple = new TTree("T_rec","T_rec");
  t_rec_simple->SetDirectory(file1);
  t_rec_simple->Branch("x",&x,"x/D");
  t_rec_simple->Branch("y",&y,"y/D");
  t_rec_simple->Branch("z",&z,"z/D");
  t_rec_simple->Branch("cluster_id",&ncluster,"cluster_id/I");
  t_rec_simple->Branch("real_cluster_id",&real_cluster_id,"real_cluster_id/I");
  t_rec_simple->Branch("sub_cluster_id",&sub_cluster_id,"sub_cluster_id/I");
  t_rec_simple->SetDirectory(file1);
  
  // extra save ...
  TTree *t_rec_deblob = new TTree("T_rec_charge_blob","T_rec_charge_blob");
  t_rec_deblob->SetDirectory(file1);
  t_rec_deblob->Branch("x",&x,"x/D");
  t_rec_deblob->Branch("y",&y,"y/D");
  t_rec_deblob->Branch("z",&z,"z/D");
  t_rec_deblob->Branch("q",&charge_save,"q/D");
  t_rec_deblob->Branch("nq",&ncharge_save,"nq/D");
  t_rec_deblob->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_deblob->Branch("ndf",&ndf_save,"ndf/D");
  t_rec_deblob->Branch("cluster_id",&ncluster,"cluster_id/I");
  t_rec_deblob->Branch("real_cluster_id",&real_cluster_id,"real_cluster_id/I");
  t_rec_deblob->Branch("sub_cluster_id",&sub_cluster_id,"sub_cluster_id/I");
  
  
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

  // // save the Steiner-inspired Graph ...
  // for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){
  //   WCPPID::PR3DCluster* cluster = *it;
  //   ndf_save = cluster->get_cluster_id();
  //   if (cluster->get_point_cloud_steiner()!=0){
  //     WCP::WCPointCloud<double>& cloud = cluster->get_point_cloud_steiner()->get_cloud();
  //     std::vector<bool>& flag_tagged = cluster->get_flag_tagged_steiner_graph();
  //     for (size_t i=0;i!=cloud.pts.size();i++){
  // 	x = cloud.pts[i].x/units::cm;
  // 	y = cloud.pts[i].y/units::cm;
  // 	z = cloud.pts[i].z/units::cm;
  // 	t_rec_simple->Fill();

  // 	charge_save = 1;
  // 	ncharge_save = 1;
  // 	chi2_save = 1;
  
  // 	if (flag_tagged.size()>0){
  // 	  if (!flag_tagged.at(i)) t_rec_deblob->Fill();
  // 	}
  //     }
      
  //   }
  // }

  

  // original save ...
  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){
    WCPPID::PR3DCluster* cluster = *it;
    ndf_save = map_cluster_parent_id[cluster];
    //    ndf_save = cluster->get_cluster_id();
    // original
    PointVector& pts = cluster->get_fine_tracking_path();
    //std::vector<double> dQ, dx;
    std::vector<double>& dQ = cluster->get_dQ();
    std::vector<double>& dx = cluster->get_dx();
    std::vector<double>& tpu = cluster->get_pu();
    std::vector<double>& tpv = cluster->get_pv();
    std::vector<double>& tpw = cluster->get_pw();
    std::vector<double>& tpt = cluster->get_pt();
    std::vector<double>& Vreduced_chi2 = cluster->get_reduced_chi2();
    std::vector<bool>& Vflag_vertex = cluster->get_flag_vertex();
    std::vector<int>& Vsub_cluster_id = cluster->get_sub_cluster_id();
    std::vector<bool>& Vflag_shower = cluster->get_flag_shower();
    
    if (pts.size()!=dQ.size() || pts.size()==0) continue;
    
    for (size_t i=0; i!=pts.size(); i++){
      x = pts.at(i).x/units::cm;
      y = pts.at(i).y/units::cm;
      z = pts.at(i).z/units::cm;
      charge_save = dQ.at(i) * dQdx_scale + dQdx_offset; // for display purpose ...
      ncharge_save = dx.at(i)/units::cm;
      pu = tpu.at(i);
      pv = tpv.at(i);
      pw = tpw.at(i);
      pt = tpt.at(i);
      reduced_chi2 = Vreduced_chi2.at(i);
      flag_vertex_save = Vflag_vertex.at(i);
      ncluster = ndf_save;
      real_cluster_id = Vsub_cluster_id.at(i);
      sub_cluster_id = Vsub_cluster_id.at(i);
      t_rec_charge->Fill();

      if (Vflag_shower.at(i))
	real_cluster_id = 1;
      else
	real_cluster_id = 4;
      
      
      t_rec_deblob->Fill();
    }

    ToyPointCloud *pcloud = cluster->get_point_cloud();          
    if (pcloud!=0){
      WCP::WCPointCloud<double>& cloud = pcloud->get_cloud();
      std::vector<int>& point_sub_cluster_ids = cluster->get_point_sub_cluster_ids();
      std::vector<bool>& point_flag_showers = cluster->get_point_flag_showers();
      for (size_t i=0;i!=cloud.pts.size();i++){
	x = cloud.pts[i].x/units::cm;
	y = cloud.pts[i].y/units::cm;
	z = cloud.pts[i].z/units::cm;
	if (point_sub_cluster_ids.size() == cloud.pts.size()){
	  if (point_sub_cluster_ids.at(i)==-1){
	    real_cluster_id = -1;
	    continue; // skip -1 points ... deghosting ...
	  }else{
	    
	    if (point_flag_showers.at(i))
	      real_cluster_id = 1;
	    else
	      real_cluster_id = 4;
	  }
	}
	sub_cluster_id =  cluster->get_cluster_id();
	t_rec_simple->Fill();
      }
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
      // if (it->first.first > 8256) std::cout << cluster->get_cluster_id() << " U " << it->first.first << std::endl;

      temp_channel.push_back(it->first.first);
      temp_timeslice.push_back(it->first.second);
      temp_charge.push_back(std::get<0>(it->second));
      temp_charge_err.push_back(std::get<1>(it->second));
      temp_charge_pred.push_back(std::get<2>(it->second));
    }
    for (auto it = proj_data_v_map.begin(); it!=proj_data_v_map.end(); it++){
      // if (it->first.first > 8256) std::cout << cluster->get_cluster_id() << " V " << it->first.first << std::endl;
      
      temp_channel.push_back(it->first.first);
      temp_timeslice.push_back(it->first.second);
      temp_charge.push_back(std::get<0>(it->second));
      temp_charge_err.push_back(std::get<1>(it->second));
      temp_charge_pred.push_back(std::get<2>(it->second));
    }
    for (auto it = proj_data_w_map.begin(); it!=proj_data_w_map.end(); it++){
      //  if (it->first.first > 8256) std::cout << cluster->get_cluster_id() << " W " << it->first.first << std::endl;
      
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
      WCPPID::PR3DCluster *cluster = (*it1);
      cluster->get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);
    }
    proj_cluster_id->push_back(cluster_id);
    proj_cluster_channel->push_back(proj_channel);
    proj_cluster_timeslice->push_back(proj_timeslice);
    proj_cluster_charge->push_back(proj_charge);
    proj_cluster_charge_err->push_back(proj_charge_err);
  }
  T_proj->Fill();
  
  
  TTree *t_bad = new TTree("T_bad","T_bad");
  t_bad->SetDirectory(file1);
  Int_t bad_npoints;
  Double_t bad_y[100],bad_z[100];
  t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
  t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
  t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");

  WCP::SMGCSelection& dead_mcells = fid->get_Dead_mcells();
  for (auto it = dead_mcells.begin(); it!=dead_mcells.end(); it++){
    const SlimMergeGeomCell *cell = (SlimMergeGeomCell*)(*it);
    PointVector ps = cell->boundary();
    bad_npoints = ps.size();
    //    std::cout << bad_npoints << std::endl;
    for (int j=0;j!=bad_npoints;j++){
      bad_y[j] = ps.at(j).y/units::cm;
      bad_z[j] = ps.at(j).z/units::cm;
    }
    t_bad->Fill();
  }
  
  file1->Write();
  file1->Close();
  
  
  return 0;
}
