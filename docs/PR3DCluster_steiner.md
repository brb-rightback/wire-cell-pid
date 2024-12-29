# PR3DCluster Steiner Tree Implementation Documentation

## Overview
The PR3DCluster Steiner tree implementation provides functionality for analyzing 3D clusters in wire chamber detectors using Steiner tree algorithms. This document describes the key components and functions of the implementation.

## Key Internal Variables

### Graph Related
- `MCUGraph *graph_steiner`: The main Steiner tree graph structure
- `WCP::ToyPointCloud *point_cloud_steiner`: Point cloud for Steiner tree vertices
- `WCP::ToyPointCloud *point_cloud_steiner_terminal`: Point cloud for terminal vertices
- `std::vector<bool> flag_steiner_terminal`: Flags indicating which points are terminals
- `std::vector<bool> flag_tagged_steiner_graph`: Tracking of tagged points in Steiner graph

### Terminal Related
- `std::set<int> steiner_terminal_indices`: Set of indices for terminal points
- `std::set<int> selected_terminal_indices`: Set of indices for selected terminals
- `std::set<int> steiner_graph_terminal_indices`: Terminal indices in final graph
- `std::set<int> steiner_graph_selected_terminal_indices`: Selected terminal indices in final graph

### Cell Mapping
- `std::map<WCP::SlimMergeGeomCell*, std::set<int>> cell_point_indices_map`: Maps cells to their point indices
- `std::vector<edge_descriptor> same_mcell_steiner_edges`: Edges connecting points in same merged cell
- `std::set<int> excluded_points`: Points excluded from analysis

## Key Functions

### create_steiner_graph()
Creates a Steiner tree graph structure from cluster data. For more details, refer to the [create_steiner_graph documentation](./steiner/create_steiner_graph.md).

**Inputs:**
- `WCP::ToyCTPointCloud& ct_point_cloud`: Point cloud data
- `WCPSst::GeomDataSource& gds`: Geometry data source
- `int nrebin`: Rebinning parameter
- `int frame_length`: Frame length
- `double unit_dis`: Unit distance

**Effects:**
- Creates and initializes graph_steiner
- Establishes point clouds and terminal flags
- Sets up edge connections

### recover_steiner_graph()
Reconstructs a Steiner graph from terminal indices.

**Effects:**
- Updates steiner_graph_terminal_indices
- Creates terminal graph structure
- Updates edge connections and weights
- Identifies minimal spanning tree

### Create_steiner_tree()
Creates a new Steiner tree from point cloud data. For more details, refer to the [create_steiner_tree documentation](./steiner/create_steiner_tree.md).

**Inputs:**
- `WCP::ToyPointCloud *point_cloud_steiner`: Target point cloud
- `std::vector<bool>& flag_steiner_terminal`: Terminal flags
- `WCP::GeomDataSource& gds`: Geometry data
- `WCP::SMGCSelection& old_mcells`: Merged cell selection
- `bool flag_path`: Path flag
- `bool disable_dead_mix_cell`: Dead cell handling flag

**Returns:**
- `MCUGraph*`: New Steiner tree graph

**Effects:**
- Creates terminal points
- Establishes graph structure
- Incorporates charge information

### find_steiner_terminals()
Identifies terminal points for Steiner tree construction. For more details, refer to the [find_steiner_terminals documentation](./steiner/find_steiner_terminals.md).

**Inputs:**
- `WCP::GeomDataSource& gds`: Geometry data
- `bool disable_dead_mix_cell`: Dead cell handling flag

**Effects:**
- Updates steiner_terminal_indices
- Creates cell-point mappings

### find_peak_point_indices()
Identifies peak points within merged geometry cells. For more details, refer to the [find_peak_point_indices documentation](./steiner/find_peak_point_indices.md).

**Inputs:**
- `SMGCSelection mcells`: Merged cell selection
- `WCP::GeomDataSource& gds`: Geometry data
- `bool disable_dead_mix_cell`: Dead cell handling flag
- `int nlevel`: Level parameter for search

**Returns:**
- `std::set<int>`: Set of peak point indices

### calc_charge_wcp()
Calculates charge for wire chamber points. For more details, refer to the [calc_charge_wcp documentation](./steiner/calc_charge_wcp.md).

**Inputs:**
- `WCP::WCPointCloud<double>::WCPoint& wcp`: Wire chamber point
- `WCP::GeomDataSource& gds`: Geometry data
- `bool disable_dead_mix_cell`: Dead cell handling flag
- `double charge_cut`: Charge threshold

**Returns:**
- `std::pair<bool,double>`: Validity flag and calculated charge

### form_cell_points_map()
Creates mappings between cells and their points. For more details, refer to the [form_cell_points_map documentation](./steiner/form_cell_points_map.md).


**Effects:**
- Populates cell_point_indices_map
- Organizes point indices by cell

## Algorithm Flow
1. **Initialization**
   - Create point clouds and graph structures
   - Initialize terminal flags and mappings

2. **Terminal Selection**
   - Identify significant points using charge information
   - Create terminal set based on peak points
   - Map terminals to graph vertices

3. **Graph Construction**
   - Create edges between terminals
   - Assign weights based on distances and charges
   - Build minimal spanning tree

4. **Optimization**
   - Refine terminal selection
   - Adjust edge weights
   - Create final Steiner tree structure

## Usage Notes
- The implementation requires proper initialization of geometry and point cloud data
- Charge calculations consider multiple wire planes (u, v, w)
- Dead/disabled channels are handled through specific flags
- Graph operations utilize the Boost Graph Library
- Terminal selection considers both spatial and charge information

## Dependencies
- Boost Graph Library
- WCP Data Structures
- Eigen Library (for matrix operations)
- PAAL Library (for Steiner tree algorithms)