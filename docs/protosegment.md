# ProtoSegment Class Documentation

## Overview
The `ProtoSegment` class is part of the WCPPID namespace and represents a segment of a particle track or shower in a Wire-Cell detector. It handles track/shower properties, geometry, particle identification, and various analysis functions.

## Class Members

### Key Properties
- `id`: Unique identifier for the segment
- `cluster_id`: ID of the cluster this segment belongs to
- `particle_type`: Type of particle (-1 undetermined)
  - e-: 11, e+: -11
  - μ-: 13, μ+: -13
  - γ: 22
  - π+: 211, π⁰: 111, π-: -211
  - K+: 321, K-: -321
  - p: 2212
  - n: 2112
- `particle_mass`: Mass of the identified particle
- `particle_score`: Confidence score of particle identification (0-100)
- `particle_4mom`: 4-momentum vector components
- `kenergy_best`: Best estimate of kinetic energy
- `kenergy_charge`: Kinetic energy from charge measurement

### Geometry Data
- `wcpt_vec`: Vector of WCPoints representing the segment points
- `fit_pt_vec`: Vector of fitted 3D points
- `point_cloud`: Point cloud representation
- `pcloud_fit`: Fitted point cloud
- `pcloud_associated`: Associated points cloud
- `pcloud_associated_steiner`: Steiner tree associated points

### Track/Shower Properties
- `dQ_vec`: Charge deposition vector
- `dx_vec`: Distance intervals vector
- `dQ_dx_vec`: dQ/dx measurements vector
- `pu_vec`, `pv_vec`, `pw_vec`: Projections on U, V, W wire planes
- `pt_vec`: Time projections
- `reduced_chi2_vec`: χ² values for fits
- `flag_shower_trajectory`: Indicates shower-like trajectory
- `flag_shower_topology`: Indicates shower-like topology
- `flag_fit`: Indicates if segment has been fitted
- `flag_dir`: Direction flag (1: forward, -1: backward, 0: undetermined)
- `dir_weak`: Indicates if direction determination is weak

## Major Functions

### Constructor/Destructor
```cpp
ProtoSegment(int id, std::list<WCP::WCPointCloud<double>::WCPoint>& path_wcps, int cluster_id);
~ProtoSegment();
```

### Point Cloud Management
```cpp
void build_pcloud_fit();  // building a point cloud for fitted results
void add_associate_point(WCPointCloud<double>::WCPoint& wcp, ...); // add associated points, for EM shower for example
void reset_associate_points();
```

### Geometry Analysis
```cpp
double get_length();
double get_direct_length();
std::pair<double, WCP::Point> get_closest_point(WCP::Point& p);
std::tuple<double, double, double> get_closest_2d_dis(WCP::Point& p);
double WCPPID::ProtoSegment::get_max_deviation(int n1, int n2)
```
[details on get_length](./pattern_recognition/protosegment_get_length.md)

[details on get_direct_length](./pattern_recognition/protosegment_get_direct_length.md)

Get the maximum deviation for track segment to a straight line (n1 to n2).

### Track/Shower Classification
```cpp
bool is_shower_trajectory(double step_size = 10.*units::cm);
bool is_shower_topology(bool val = false);
bool determine_shower_direction();
bool get_flag_shower();
```

[see details about shower](./pattern_recognition/protosegment_is_shower.md)

[see details on determine shower direction](./pattern_recognition/protosegment_determine_shower_direction.md)

### Particle Identification
```cpp
void determine_dir_track(int start_n, int end_n, bool flag_print = false);
bool do_track_pid(std::vector<double>& L, std::vector<double>& dQ_dx, ...);
int get_particle_type();
```

[more details](./pattern_recognition/protosegment_track.md)

### Energy/Momentum Analysis
```cpp
double cal_kine_range();
double cal_kine_dQdx();
void cal_4mom();
TVector3 cal_dir_3vector();
```
[more details](./pattern_recognition/protosegment_kinematics.md)

### Track Breaking/Splitting
```cpp
std::tuple<ProtoSegment*, ProtoVertex*, ProtoSegment*> 
    break_segment_at_point(WCP::Point& p, int& acc_segment_id, int& acc_vertex_id);
std::tuple<WCP::Point, TVector3, TVector3, bool> WCPPID::ProtoSegment::search_kink(Point& start_p){
```

[more details on break_segment](./pattern_recognition/protosegment_break_segment_at_point.md)

[more details on search kink](./pattern_recognition/protosegment_search_kink.md)

### Data Access
```cpp
std::vector<WCP::WCPointCloud<double>::WCPoint>& get_wcpt_vec();
std::vector<WCP::Point>& get_point_vec();
std::vector<double>& get_dQ_vec();
// ... and many more getters for internal vectors
```

[details on get_dQ_dx](./pattern_recognition/protosegment_get_dQ_dx.md)



## Usage Notes

1. **Track/Shower Discrimination**:
   - Uses multiple criteria including trajectory shape, topology, and dQ/dx
   - Can identify shower-like segments using both geometric and calorimetric information

2. **Particle Identification**:
   - Employs dQ/dx patterns for particle type determination
   - Uses Kolmogorov-Smirnov tests to compare with expected patterns
   - Handles both tracks (μ, π, p) and showers (e, γ)

3. **Energy Reconstruction**:
   - Combines range-based and calorimetric methods
   - Handles both track-like and shower-like energy deposits
   - Accounts for detector effects and reconstruction uncertainties

4. **Point Cloud Operations**:
   - Maintains multiple point cloud representations
   - Supports efficient spatial searches and associations
   - Handles both raw and fitted point collections

5. **Direction Finding**:
   - Uses multiple algorithms for direction determination
   - Considers both geometric and calorimetric information
   - Handles special cases like short tracks and ambiguous directions