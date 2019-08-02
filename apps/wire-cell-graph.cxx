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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/matching.root" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

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
      live_clusters.push_back(cluster);
    }
    cluster->AddCell(mcell,time_slice);


    
    prev_cluster_id = cluster_id;
    ident++;
  }
  
  prev_cluster_id = -1;

  cout << em("load clusters from file") << endl;

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
  
  
  for (size_t i=0; i!=live_clusters.size();i++){
     live_clusters.at(i)->create_steiner_graph(ct_point_cloud, gds, nrebin, frame_length, unit_dis);
     //ToyPointCloud *pcloud = live_clusters.at(i)->get_point_cloud_steiner();
     // std::cout << pcloud << std::endl;
    
     //old_new_cluster_map[live_clusters.at(i)] = new_cluster;
  }
  cout << em("Build graph for all clusters") << std::endl;
  
  
  TFile *file1 = new TFile(Form("graph_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  Trun->CloneTree()->Write();

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
    //WireCellPID::PR3DCluster* new_cluster = old_new_cluster_map[*it];
    WireCellPID::PR3DCluster* new_cluster = *it;  
    ncluster = new_cluster->get_cluster_id();
    
    ToyPointCloud *pcloud = new_cluster->get_point_cloud();
    //ToyPointCloud *pcloud = new_cluster->get_point_cloud_steiner();
    // ToyPointCloud *pcloud = new_cluster->get_point_cloud_steiner_terminal();

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
  
  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){

    WireCellPID::PR3DCluster* new_cluster = *it;
    //WireCellPID::PR3DCluster* new_cluster = old_new_cluster_map[*it];

    if (new_cluster->get_point_cloud_steiner()->get_num_points() >= 2){
      std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = new_cluster->get_two_boundary_wcps(2); 
      new_cluster->dijkstra_shortest_paths(wcps.first,2); 
      new_cluster->cal_shortest_path(wcps.second,2);
    }
    
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
    
    // // std::set<int> steiner_terminals = new_cluster->get_steiner_terminals();
    // std::set<int> steiner_terminals = new_cluster->get_selected_terminals();
    // WireCell::ToyPointCloud* point_cloud = new_cluster->get_point_cloud();
    // WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
    
    // for (auto it1 = steiner_terminals.begin();it1!=steiner_terminals.end();it1++){
    //   x = cloud.pts[*it1].x/units::cm;
    //   y = cloud.pts[*it1].y/units::cm;
    //   z = cloud.pts[*it1].z/units::cm;
    //   pu = cloud.pts[*it1].index_u;
    //   pv = cloud.pts[*it1].index_v;
    //   pw = cloud.pts[*it1].index_w;
    //   pt = cloud.pts[*it1].mcell->GetTimeSlice();
    //   ndf_save = (*it)->get_cluster_id();
    //   charge_save = 0;
    //   ncharge_save = 0;
    //   chi2_save = 0;
    //   T_rec->Fill();
    // }
  }

  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){

    WireCellPID::PR3DCluster* new_cluster = *it;
    //WireCellPID::PR3DCluster* new_cluster = old_new_cluster_map[*it];

    std::set<int> steiner_terminals = new_cluster->get_steiner_terminals();
    //std::set<int> steiner_terminals = new_cluster->get_selected_terminals();
    WireCell::ToyPointCloud* point_cloud = new_cluster->get_point_cloud();
    WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
    
    for (auto it1 = steiner_terminals.begin();it1!=steiner_terminals.end();it1++){
      x = cloud.pts[*it1].x/units::cm;
      y = cloud.pts[*it1].y/units::cm;
      z = cloud.pts[*it1].z/units::cm;
      pu = cloud.pts[*it1].index_u;
      pv = cloud.pts[*it1].index_v;
      pw = cloud.pts[*it1].index_w;
      pt = cloud.pts[*it1].mcell->GetTimeSlice();
      ndf_save = (*it)->get_cluster_id();
      charge_save = 0;
      ncharge_save = 0;
      chi2_save = 0;
      t_rec_charge->Fill();
    }
  }
  
  
  
  file1->Write();
  file1->Close();
  
  
  return 0;
}
