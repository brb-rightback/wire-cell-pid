# PR3DCluster Path Functions Analysis

This document details the functionality of various path-finding and point analysis methods in the PR3DCluster class.

## Point Analysis Functions

### Boundary Point Functions

1. **get_highest_lowest_wcps**
   - Purpose: Finds points with extreme Y coordinates
   - Parameters: `flag` (determines point cloud type)
   - Returns: Pair of points (highest, lowest)
   - Algorithm:
     - Iterates through cloud points
     - Tracks maximum and minimum Y coordinates
     - Considers excluded points

2. **get_front_back_wcps**
   - Purpose: Finds points with extreme Z coordinates
   - Similar to highest/lowest but for Z axis
   - Used for front/back boundary detection

3. **get_earliest_latest_wcps**
   - Purpose: Finds points with extreme X coordinates
   - Important for drift time analysis
   - Returns pair of points (earliest, latest)

### Complex Point Analysis

1. **get_extreme_wcps**
   - Purpose: Comprehensive extreme point analysis
   - Key Features:
     - Analyzes points along multiple axes
     - Considers time cells map
     - Groups nearby points
     - Returns vector of point vectors
   - Uses PCA for main axis determination

2. **get_main_axis_wcps**
   - Purpose: Finds extreme points along principal axis
   - Algorithm:
     - Calculates PCA
     - Projects points onto main axis
     - Returns extreme points pair

## Path Finding Methods

### Basic Path Functions

1. **get_local_extension**
   - Purpose: Extends path locally from given point
   - Uses VHoughTrans for direction finding
   - Considers drift direction constraints
   - Maximum extension of 10cm

[detailed algorithm](./track_fitting/get_local_extension.md)

2. **get_num_outside_range_points**
   - Purpose: Counts points outside main path range
   - Checks:
     - X coordinate range
     - Wire indices (U,V,W planes)
   - Used for path validation

[detailed algorithm](./track_fitting/get_num_outside_range.md)

### Advanced Path Analysis

1. **get_two_boundary_wcps**
   - Purpose: Finds optimal boundary points
   - Complex Algorithm:
     - Considers wire plane indices
     - Analyzes live channels
     - Uses cosmic ray flag for different metrics
   - Returns most separated valid point pair

2. **get_two_extreme_points**
   - Purpose: Finds maximally separated points
   - Algorithm:
     - Finds extremes in X,Y,Z
     - Calculates distances between candidates
     - Averages positions in 5cm radius

## Path Finding Core

### Dijkstra Implementation

1. **dijkstra_shortest_paths**
   - Purpose: Implements Dijkstra's algorithm
   - Features:
     - Uses Boost graph library
     - Supports different graph types
     - Stores parents and distances

2. **cal_shortest_path**
   - Purpose: Calculates path between points
   - Outputs:
     - Ordered list of points (path_wcps)
     - Associated cells (path_mcells)
   - Uses previously calculated Dijkstra results

## Usage Notes

1. **Point Cloud Selection**
   - flag=1: Regular point cloud
   - flag=2: Steiner point cloud

2. **Graph Initialization**
   - Requires Create_graph() call
   - Maintains separate graphs for regular/Steiner points

3. **Performance Considerations**
   - Caches results where possible
   - Uses excluded_points set for filtering
   - Maintains separate paths for different purposes
