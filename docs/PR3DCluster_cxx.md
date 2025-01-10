# PR3DCluster Function Documentation

## 1. Cluster Management Functions

### AddCell
```cpp
PR3DCluster::AddCell(SlimMergeGeomCell* mcell, int time_slice)
```
- Adds a merge cell to the cluster
- Maintains bidirectional mappings:
  - Between cells and their time slices
  - Between time slices and their cells
- Updates both `cell_times_set_map` and `time_cells_set_map`

### Create_point_cloud
```cpp
PR3DCluster::Create_point_cloud(ToyPointCloud *global_point_cloud)
```
- Creates a point cloud representation of the cluster
- Samples points from all merge cells in the cluster
- Builds KD-tree index for efficient spatial queries
- Optionally adds points to a global point cloud

## 2. Position/Geometry Functions

### calc_ave_pos
```cpp
PR3DCluster::calc_ave_pos(Point& p, double dis)
```
- Calculates charge-weighted average position within radius `dis` of point `p`
- Uses point cloud for nearest neighbor search
- Returns weighted centroid based on cell charges

### Calc_PCA  [more details](./track_fitting/calc_PCA.md)
```cpp
PR3DCluster::Calc_PCA()
```
- Performs Principal Component Analysis on cluster points
- Calculates:
  - Cluster center
  - Principal axes (3 orthogonal directions)
  - Principal values (variance along each axis)
- Used for trajectory direction analysis and shape characterization

### Hough Transform [more details](./track_fitting/Hough.md)
The PR3DCluster class implements several Hough transform functions used for direction finding and pattern recognition in 3D space.

## 3. Tracking Functions

### do_tracking   [more details](./track_fitting/do_tracking.md)
```cpp
PR3DCluster::do_tracking(WCP::ToyCTPointCloud& ct_point_cloud, 
                        std::map<int,std::map<const GeomWire*, SMGCSelection>>& global_wc_map,
                        double time, bool flag_dQ_dx_fit_reg, bool flag_dQ_dx_fit)
```
- Main tracking function to reconstruct 3D trajectory
- Process:
  1. Organizes initial path points
  2. Performs trajectory fitting
  3. Calculates dQ/dx
- Handles both initial tracking and subsequent refinement passes
- Accounts for detector configuration and measurement uncertainties

### trajectory_fit [mode details](./track_fitting/trajectory_fit.md)
```cpp
PR3DCluster::trajectory_fit(WCP::PointVector& ps_vec, /*...multiple map parameters...*/)
```
- Fits continuous 3D trajectory through provided points
- Uses charge distribution in 2D wire plane projections
- Implements regularization for smooth trajectories
- Handles measurement uncertainties and dead channels

## 4. Charge Analysis Functions 

### dQ_dx_fit [more details](./track_fitting/dQ_dx_fit.md)
```cpp
PR3DCluster::dQ_dx_fit(std::map<int,std::map<const GeomWire*, SMGCSelection>>& global_wc_map,
                       /*...charge map parameters...*/,
                       double flash_time, double dis_end_point_ext, 
                       bool flag_dQ_dx_fit_reg)
```
- Calculates charge deposition (dQ/dx) along trajectory
- Accounts for detector effects:
  - Electron diffusion
  - Electron lifetime
  - Dead channels
- Uses regularization for stability
- Handles induction/collection plane differences

### cal_gaus_integral [more details](./track_fitting/cal_gaus_integral.md)
```cpp
PR3DCluster::cal_gaus_integral(int tbin, int wbin, double t_center, 
                              double t_sigma, double w_center, 
                              double w_sigma, int flag, double nsigma)
```
- Helper function for charge calculations
- Computes Gaussian integrals for charge spread
- Different modes for:
  - Collection plane (flag = 0)
  - Induction plane (flag = 1)
- Handles boundary effects

## 5. Path Organization Functions

### organize_wcps_path  [more details](./track_fitting/organize_wcps_path.md)
```cpp
PR3DCluster::organize_wcps_path(WCP::ToyCTPointCloud& ct_point_cloud,
                               std::list<WCPointCloud<double>::WCPoint>& path_wcps_list,
                               double low_dis_limit, double end_point_limit)
```
- Organizes wire cell points into continuous path
- Handles path endpoints specially
- Enforces minimum/maximum point spacing
- Returns organized point vector

### organize_ps_path [more details](./track_fitting/organize_ps_path.md)
```cpp
PR3DCluster::organize_ps_path(WCP::ToyCTPointCloud& ct_point_cloud,
                             WCP::PointVector& pts,
                             double low_dis_limit, double end_point_limit)
```
- Refines point spacing along path
- Handles endpoints with specified extension
- Ensures consistent point density
- Validates points against detector geometry

## 6. Data Preparation Functions

### prepare_data [more details](./track_fitting/prepare_data.md)
```cpp
PR3DCluster::prepare_data(WCP::ToyCTPointCloud& ct_point_cloud,
                         std::map<int,std::map<const GeomWire*, SMGCSelection>>& global_wc_map,
                         /*...charge map parameters...*/)
```
- Prepares 2D charge data for tracking
- Creates charge maps for each wire plane
- Handles dead channels and uncertainties
- Processes overlap regions between cells

### form_map  [more details](./track_fitting/form_map.md)
```cpp
PR3DCluster::form_map(WCP::ToyCTPointCloud& ct_point_cloud,
                      WCP::PointVector& pts,
                      /*...multiple map parameters...*/)
```
- Creates mappings between:
  - 3D points and 2D projections
  - Different wire plane views
- Handles charge sharing between adjacent wires
- Manages uncertainties in projections

## 7. Charge Collection and Projection Functions

### collect_charge_trajectory [more details](./track_fitting/collect_charge_trajectory.md)
```cpp
PR3DCluster::collect_charge_trajectory(ToyCTPointCloud& ct_point_cloud, 
                                     double dis_cut, double range_cut)
```
- Collects charge information along reconstructed trajectory
- Key operations:
  1. Creates set of existing time-channel pairs from cluster
  2. Forms trajectory from fine tracking or path points
  3. For each trajectory point:
     - Searches for nearby points within range_cut
     - Collects charge from U, V, W planes
     - Checks against existing time-channel pairs
  4. Stores collected charge in `collected_charge_map`
- Parameters:
  - `dis_cut`: Distance cut for trajectory point spacing
  - `range_cut`: Range for collecting nearby charges
- Handles:
  - Dead/noisy channels
  - Charge sharing between wires
  - Time slice boundaries

### get_projection [more details](./track_fitting/get_projection.md)
```cpp
PR3DCluster::get_projection(std::vector<int>& proj_channel,
                          std::vector<int>& proj_timeslice,
                          std::vector<int>& proj_charge,
                          std::vector<int>& proj_charge_err,
                          std::vector<int>& proj_flag,
                          std::map<int,std::map<const GeomWire*, SMGCSelection>>& global_wc_map)
```
- Creates 2D projections of the cluster onto wire planes
- Key features:
  1. Projects cluster onto U, V, W planes
  2. Handles wire overlaps and charge sharing
  3. Assigns flags to indicate charge quality:
     - flag 0: dead channels
     - flag 1: good channels or overlapping channels
     - flag 2: isolated channels
     - flag 3: additional channels
  4. Calculates charge uncertainties
- Operations:
  1. Identifies cluster merge cells
  2. Projects each cell onto wire planes
  3. Handles charge sharing for overlapping wires
  4. Calculates charge errors
  5. Fills projection vectors with results
- Handles special cases:
  - Dead channels
  - Overlapping cells
  - Channel boundaries
  - Noise thresholds

## Technical Notes

The code is part of a 3D track reconstruction system that:
1. Takes 2D wire plane measurements as input
2. Builds 3D space points from coincidences
3. Fits continuous trajectories through the points
4. Calculates charge deposition along tracks
5. Handles various detector effects and uncertainties

The architecture leverages:
- Point clouds for spatial organization
- Graph algorithms for connectivity
- Regularized fitting for robustness
- KD-trees for efficient nearest neighbor searches

The system is designed to handle common challenges in wire chamber reconstruction:
- Dead/noisy channels
- Charge sharing between wires
- Electron diffusion
- Signal induction effects
- Multiple scattering
- Delta rays