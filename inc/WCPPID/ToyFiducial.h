#ifndef WIRECELLPID_UBOONE_TOYFIDUCIAL_H
#define WIRECELLPID_UBOONE_TOYFIDUCIAL_H

#include "WCPData/Units.h"
#include "WCPData/Point.h"
#include "WCPData/SlimMergeGeomCell.h"
//#include "WCPData/FlashTPCBundle.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCPData/LMBDT.h"
#include "WCPData/Opflash.h"
#include "WCPData/PhotonLibrary.h"

#include "WCPPID/PR3DCluster.h"

#include "TVector3.h"
#include "TGraph.h"
#include "TFile.h"

#include <vector>
#include <map>

namespace WCPPID{
  class ToyFiducial{
  public:
    ToyFiducial(int dead_region_ch_ext = 3, double offset_t=800, double offset_u=0, double offset_v=0, double offset_w=0, double slope_t=1./2*units::mm, double slope_u=1./(3*units::mm), double slope_v=1./(3*units::mm), double slope_w=1./(3*units::mm), double angle_u=-1.047198, double angle=1.047198, double angle_w=0,
		double boundary_dis_cut=2*units::cm, double top=117*units::cm, double bottom=-116*units::cm, double upstream=0*units::cm, double downstream=1037*units::cm, double anode = 0*units::cm, double cathode=256*units::cm, int flag_data=1);
    ~ToyFiducial();

  
    
    
    void set_offset_t(double value){offset_t=value;};
    
    // helper functions ...
    bool inside_fiducial_volume(WCP::Point& p, double offset_x=0);
    bool inside_dead_region(WCP::Point& p);
    bool check_dead_volume(WCP::Point& p, TVector3& dir, double step = 1.0*units::cm, double offset_x=0);
    bool check_signal_processing(WCP::Point& p, TVector3& dir, WCP::ToyCTPointCloud& ct_point_cloud, double step = 1.0*units::cm, double offset_x=0);

    // TGM tagger ...
    bool check_neutrino_candidate(WCPPID::PR3DCluster *main_cluster, WCP::WCPointCloud<double>::WCPoint& wcp1, WCP::WCPointCloud<double>::WCPoint& wcp2, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud, bool flag_2view_check = true);
    bool check_tgm(WCPPID::PR3DCluster* main_cluster,WCP::Opflash* main_flash, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud);
    
    // fully contained tagger ...
    bool check_fully_contained(WCPPID::PR3DCluster* main_cluster, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud);
    
    // London's new tagger
    std::tuple<int, WCPPID::PR3DCluster*, WCP::Opflash*> cosmic_tagger(WCP::OpflashSelection& flashes, WCPPID::PR3DCluster* main_cluster, WCP::Opflash* main_flash, std::tuple<int, double, double, int>& bundle_info, WCP::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud, int run_no, int subrun_no, int event_no, bool flag_data, bool debug_tagger=false);

    // check STM code ...
    bool check_stm(WCPPID::PR3DCluster* cluster, std::vector<WCPPID::PR3DCluster*>& additional_clusters, double offset_x, double flash_time, WCP::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, int& event_type);

    bool eval_stm(WCPPID::PR3DCluster* main_cluster, int kink_num, double peak_range = 40*units::cm, double offset_x = 0*units::cm, double com_range = 35*units::cm);

    bool detect_proton(WCPPID::PR3DCluster* main_cluster, int kink_num);

    int find_first_kink(WCPPID::PR3DCluster* main_cluster);

    bool check_full_detector_dead();
    
    void AddDeadRegion(WCP::SlimMergeGeomCell* mcell, std::vector<int>& time_slices);

    WCP::SMGCSelection& get_Dead_mcells(){return mcells;};

    
    
  protected:
    // boundary
    double m_top; // top distance
    double m_bottom; // bottom distance
    double m_upstream;
    double m_downstream;
    double m_anode;
    double m_cathode;
    
    // space charge boundary
    double m_sc_bottom_1_x, m_sc_bottom_1_y;
    double m_sc_bottom_2_x, m_sc_bottom_2_y;

    double m_sc_top_1_x, m_sc_top_1_y;
    double m_sc_top_2_x, m_sc_top_2_y;

    double m_sc_upstream_1_x, m_sc_upstream_1_z;
    double m_sc_upstream_2_x, m_sc_upstream_2_z;

    double m_sc_downstream_1_x, m_sc_downstream_1_z;
    double m_sc_downstream_2_x, m_sc_downstream_2_z;

    std::vector<double> boundary_xy_x, boundary_xy_y;
    std::vector<double> boundary_xz_x, boundary_xz_z;
    
    // dead regions ... 
    WCP::SMGCSelection mcells;
    std::map<WCP::SlimMergeGeomCell*, std::pair<int,int>> mcell_time_map;
    std::map<int, std::set<WCP::SlimMergeGeomCell*>> ch_mcell_set_map;
    
    // conversion between positions to the channel and time ???

    // convert time into a position
    // (time_slice - offset_t) / slope_t = position_x 
    double offset_t, slope_t;
    // convert u wire number into a position
    // (u_index -offset_u) / slope_u = position_u
    double offset_u, slope_u;
    // convert v wire number into a position
    double offset_v, slope_v;
    // convert w wire number into a position 
    double offset_w, slope_w;
    double angle_u, angle_v, angle_w;

    int dead_region_ch_ext;

    TFile *file;
    TGraph *g_muon;
    TGraph *g_pion;
    TGraph *g_kaon;
    TGraph *g_proton;
    TGraph *g_electron;
    
  };
}

#endif
