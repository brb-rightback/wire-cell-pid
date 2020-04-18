#ifndef WCPPID_PROTOSEGMENT_H
#define WCPPID_PROTOSEGMENT_H

#include "WCPPID/ProtoVertex.h"
#include "WCPData/ToyPointCloud.h"

namespace WCPPID{
  class ProtoSegment{
  public:
    // initial creation
    ProtoSegment(int id, std::list<WCP::WCPointCloud<double>::WCPoint >& path_wcps , int cluster_id);
    ~ProtoSegment();
    
    
    std::vector<WCP::WCPointCloud<double>::WCPoint >& get_wcpt_vec(){return wcpt_vec;};
    std::vector<WCP::Point >& get_point_vec(){return fit_pt_vec;};
    void set_point_vec(std::vector<WCP::Point >& tmp_pt_vec);
    std::vector<double>& get_dQ_vec(){return dQ_vec;};
    std::vector<double>& get_dx_vec(){return dx_vec;};
    std::vector<double>& get_dQ_dx_vec(){return dQ_dx_vec;};
    std::vector<double>& get_pu_vec(){return pu_vec;};
    std::vector<double>& get_pv_vec(){return pv_vec;};
    std::vector<double>& get_pw_vec(){return pw_vec;};
    std::vector<double>& get_pt_vec(){return pt_vec;};
    std::vector<double>& get_reduced_chi2_vec(){return reduced_chi2_vec;};

    void build_pcloud_fit();
    
    void set_fit_associate_vec(std::vector<WCP::Point >& tmp_fit_pt_vec, std::vector<int>& tmp_fit_index, std::vector<bool>& tmp_fit_skip);
    void reset_fit_prop();
    std::vector<int>& get_fit_index_vec(){return fit_index_vec;};
    std::vector<bool>& get_fit_flag_skip(){return fit_flag_skip;};

    double get_length();
    double get_direct_length();
    double get_direct_length(int n1, int n2);
    double get_length(int n1, int n2);

    double get_direct_length(int n1, int n2, TVector3 dir_perp);
    double get_length(int n1, int n2, TVector3 dir_perp);

    double get_medium_dQ_dx();
    double get_medium_dQ_dx(int n1, int n2);
    double get_rms_dQ_dx();

    WCP::WCPointCloud<double>::WCPoint get_closest_wcpt(WCP::Point& test_p);

    int get_id(){return id;};
    int get_cluster_id(){return cluster_id;};
    // break a segment into two from the middle
    std::tuple<ProtoSegment*, ProtoVertex*, ProtoSegment*> break_segment_at_point(WCP::Point& p, int& acc_segment_id, int& acc_vertex_id);

    
    // after fit
    void clear_fit();
    // set the data ...
    void set_fit_vec(std::vector<WCP::Point >& tmp_fit_pt_vec, std::vector<double>& tmp_dQ_vec, std::vector<double>& tmp_dx_vec, std::vector<double>& tmp_pu_vec, std::vector<double>& tmp_pv_vec, std::vector<double>& tmp_pw_vec, std::vector<double>& tmp_pt_vec, std::vector<double>& tmp_reduced_chi2_vec);
    

    // if fit_flag is set, the fit_pts are useful ...
    bool get_fit_flag(){return flag_fit;};
    void set_fit_flag(bool flag){flag_fit = flag;};

    // get direction
    void print_dis();
    
    // get point
    std::pair<double, WCP::Point> get_closest_point(WCP::Point &p);
    std::tuple<double, double, double> get_closest_2d_dis(WCP::Point &p);
    double get_closest_2d_dis(double x, double y, int plane);
    
    // search for kinks ...  return  position, direction ...
    std::tuple<WCP::Point, TVector3, TVector3, bool> search_kink(WCP::Point& start_p);

    void reset_associate_points();
    WCP::ToyPointCloud* get_fit_pcloud(){return pcloud_fit;};
    WCP::ToyPointCloud* get_associated_pcloud(){return pcloud_associated;};
    WCP::ToyPointCloud* get_associated_pcloud_steiner(){return pcloud_associated_steiner;};
    
    void add_associate_point(WCP::WCPointCloud<double>::WCPoint& wcp, WCP::WC2DPointCloud<double>::WC2DPoint& wcp_u, WCP::WC2DPointCloud<double>::WC2DPoint& wcp_v, WCP::WC2DPointCloud<double>::WC2DPoint& wcp_w);
    void add_associate_point_steiner(WCP::WCPointCloud<double>::WCPoint& wcp);
    
    //
    bool is_shower_trajectory(double step_size =10.*units::cm);
    bool is_shower_topology();
    
    bool get_flag_shower();
    void set_flag_shower_trajectory(bool val){flag_shower_trajectory = val;};
    void set_flag_shower_topology(bool val){flag_shower_topology = val;};
    bool get_flag_shower_trajectory(){return flag_shower_trajectory;};
    bool get_flag_shower_topology(){return flag_shower_topology;};
    bool get_flag_shower_dQdx();

    bool is_dir_weak();
    void set_dir_weak(bool val){dir_weak = val;};
    int get_flag_dir(){return flag_dir;};
    void set_flag_dir(int val){flag_dir = val;};
    int get_particle_type();
    
    int set_particle_type(int type){particle_type = type;};
    void set_particle_mass(double val){particle_mass = val;};

    void set_particle_score(double val){particle_score = val;};
    double get_particle_score(){return particle_score;};
    double get_particle_mass(){return particle_mass;};
    double get_particle_4mom(int num){return particle_4mom[num];};

    
    void determine_dir_track(int start_n, int end_n, bool flag_print = false);
    bool do_track_pid(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range = 35*units::cm, double offset_length = 0*units::cm);
    std::vector<double> do_track_comp(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range = 35*units::cm, double offset_length = 0*units::cm);
    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2);
    void cal_4mom();
    TVector3 cal_dir_3vector();
    TVector3 cal_dir_3vector(WCP::Point& p, double dis_cut);
    
    double cal_kine_range();
    double cal_kine_dQdx();
    
    double cal_kine_range(double length);
    double cal_kine_dQdx(std::vector<double>& vec_dQ, std::vector<double>& vec_dx);
    
    double get_kine_charge(){return kenergy_charge;};
    void set_kine_charge(double val){kenergy_charge = val;};
    //    double cal_kine_charge(std::map<std::pair<int,int>, std::pair<double,double> >& charge_2d_u, std::map<std::pair<int,int>, std::pair<double,double> >& charge_2d_v, std::map<std::pair<int,int>, std::pair<double,double> >& charge_2d_w, WCP::ToyCTPointCloud* ct_point_cloud); // for shower ...
    
    
    void determine_dir_shower_trajectory(int start_n, int end_n, bool flag_print = false);
    void determine_dir_shower_topology(int start_n, int end_n, bool flag_print = false);
    
    
  protected:
    int id;
    int cluster_id;
    std::vector<WCP::WCPointCloud<double>::WCPoint > wcpt_vec;

    std::vector<WCP::Point > fit_pt_vec;
    std::vector<double> dQ_vec;
    std::vector<double> dx_vec;
    std::vector<double> dQ_dx_vec;
    std::vector<double> pu_vec;
    std::vector<double> pv_vec;
    std::vector<double> pw_vec;
    std::vector<double> pt_vec;
    std::vector<double> reduced_chi2_vec;

    std::vector<int> fit_index_vec;
    std::vector<bool> fit_flag_skip; // vertex???
    
    // point cloud ...
    WCP::ToyPointCloud* pcloud_fit;

    WCP::ToyPointCloud* pcloud_associated;
    //std::vector<bool> flag_good_associated;
    WCP::ToyPointCloud* pcloud_associated_steiner;
    //std::vector<bool> flag_good_associated_steiner;
    //    std::vector<WCP::Point > associated_points;

    bool flag_shower_trajectory;
    bool flag_shower_topology;
    bool flag_fit;
    
    
    int flag_dir;
    bool dir_weak;
    int particle_type; // -1 undetermined,
    double particle_score;
    // e- 11  e+ -11
    // muon- 13  muon+ -13
    // gamma 22
    // pi+ 211, pi- 111, pi- -211
    // kaon+ 321, K- -321
    // p  2212
    // n 2112
    double particle_mass;
    double particle_4mom[4];
    double kenergy_charge;
  };
  typedef std::vector<ProtoSegment*> ProtoSegmentSelection;
  
  struct ProtoSegmentCompare{
    bool operator() (ProtoSegment *a, ProtoSegment *b) const {
      
      if (a->get_id() < b->get_id()){
	return true;
      }else if (a->get_id() > b->get_id()){
	return false;
      }else if (a->get_id() == b->get_id()){
	return a < b;
      }
      return a < b;
    }
  };
  typedef std::set<ProtoSegment*, ProtoSegmentCompare> ProtoSegmentSet;
}

#endif 
