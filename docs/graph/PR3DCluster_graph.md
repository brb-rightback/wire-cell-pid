# PR3DCluster Graph Documentation

## Overview
The PR3DCluster class implements graph-based functionality for 3D particle reconstruction in wire chamber detectors. It uses the Boost Graph Library to create and manipulate graphs representing connected points in 3D space.

## Key Internal Variables

### Graph-Related Variables
- `MCUGraph *graph`: Main graph structure using Boost's adjacency list
- `WCP::ToyPointCloud *point_cloud`: Collection of 3D points in the cluster 
- `std::vector<int> point_sub_cluster_ids`: IDs for point sub-clusters
- `std::vector<bool> point_flag_showers`: Flags indicating shower-like points

### Steiner Tree Variables
- `MCUGraph *graph_steiner`: Graph structure for Steiner tree calculations
- `WCP::ToyPointCloud *point_cloud_steiner`: Points used in Steiner tree
- `std::vector<bool> flag_steiner_terminal`: Flags for terminal points in Steiner tree
- `std::set<int> steiner_terminal_indices`: Indices of terminal points
- `std::vector<edge_descriptor> same_mcell_steiner_edges`: Edges connecting points in same merged cell

## Key Functions

### Del_graph()
**Purpose**: Deletes and cleans up the graph structure

**Implementation Details**:
- Checks if graph exists 
- Deletes graph pointer if it exists
- Sets pointer to null

### remove_same_mcell_steiner_edges(int flag)
**Purpose**: Removes edges between points in the same merged cell

**Input Parameters**:
- `flag`: Determines which graph to modify (1 for main graph, 2 for Steiner graph)

**Implementation Details**:
- Iterates through stored same-cell edges
- Removes edges from specified graph
- Clears the edge list

### establish_same_mcell_steiner_edges(WCP::GeomDataSource& gds, bool disable_dead_mix_cell, int flag)
**Purpose**: Creates edges between points within the same merged cell

**Input Parameters**:
- `gds`: Geometry data source
- `disable_dead_mix_cell`: Flag to disable mixed dead cells
- `flag`: Specifies which graph to modify

**Implementation Details**:
- Creates graph if it doesn't exist
- Maps points to their merged cells
- Creates edges between points in same cell
- Assigns edge weights based on distance

For more details, refer to the [establish_same_mcell_steiner_edges documentation](../steiner/establish_same_mcell_steiner_edges.md).

### Connect_graph(WCP::ToyCTPointCloud& ct_point_cloud, WCP::ToyPointCloud* ref_point_cloud)
**Purpose**: Connects separate components of the graph

**Input Parameters**:
- `ct_point_cloud`: Point cloud for connectivity testing
- `ref_point_cloud`: Reference point cloud (optional)

**Implementation Details**:
- Identifies disconnected components
- Creates temporary point clouds for each component
- Finds closest points between components
- Adds edges to connect components
- Uses MST (Minimum Spanning Tree) algorithm to optimize connections

For more details, refer to the [connect_graph documentation](./connect_graph.md).


### Create_graph(WCP::ToyCTPointCloud& ct_point_cloud, WCP::ToyPointCloud* ref_point_cloud)
**Purpose**: Creates initial graph structure

**Input Parameters**:
- `ct_point_cloud`: Point cloud for connectivity testing  
- `ref_point_cloud`: Reference point cloud (optional)

**Implementation Details**:
- Creates point cloud if needed
- Initializes graph with vertices for each point
- Establishes close connections
- Connects components

For more details, refer to the [create_graph documentation](./create_graph.md).

### Establish_close_connected_graph()
**Purpose**: Creates initial connections between nearby points

**Implementation Details**:
- Maps points to wire indices
- Creates edges between points sharing wire indices
- Optimizes edge creation for performance
- Handles special cases for different wire planes

For more details, refer to the [establish_close_connected_graph documentation](./establish_close_connected_graph.md).


### search_for_connection_isochronous()
**Purpose**: Searches for connections in isochronous tracks

**Input Parameters**:
- Multiple parameters for search criteria and geometry

**Implementation Details**:
- Handles special case of tracks parallel to wire planes
- Uses drift direction for calculations
- Finds connections based on geometric criteria

For more details, refer to the [search_for_connection_isochronous documentation](./search_for_connection_isochronous.md).


## Common Usage Patterns

1. Graph Creation:
```cpp
PR3DCluster* cluster = new PR3DCluster(id);
cluster->Create_graph(ct_point_cloud, ref_point_cloud);
```

2. Graph Modification:
```cpp
cluster->establish_same_mcell_steiner_edges(gds, true, 1);
cluster->Connect_graph(ct_point_cloud);
```

3. Cleanup:
```cpp
cluster->Del_graph();
```

## Important Notes

- The class heavily relies on Boost's graph library for core functionality
- Graph operations can be computationally intensive for large point clouds
- Care must be taken with memory management and graph cleanup
- Many functions have internal optimizations for handling large datasets
- Edge weights are typically based on 3D distances between points

## Performance Considerations

- Graph operations scale with number of points and edges
- MST calculations can be expensive for large disconnected components
- Memory usage increases with point cloud size
- Some operations use temporary data structures that should be cleaned up