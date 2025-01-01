# Understanding the shower_clustering_in_other_clusters Function

## Overview
The `shower_clustering_in_other_clusters` function is part of the neutrino event reconstruction in liquid argon time projection chambers. Its main purpose is to identify and cluster shower-like segments that are not part of the main cluster but should be associated with existing showers.

## Key Components and Data Structures

### Input Components:
- Main vertex
- Collection of vertices from main cluster
- Collection of vertices from existing showers
- Other clusters that haven't been processed yet
- Map of cluster to length relationships
- Map of vertex to shower relationships

### Important Data Structures:
```cpp
// Stores vertices from main cluster and existing showers
std::vector<WCPPID::ProtoVertex*> vertices;

// Maps clusters to their closest vertex information
std::map<WCPPID::PR3DCluster*, WCPPID::ProtoVertex*> map_cluster_associated_vertex;

// Tracks which clusters have been used in showers
std::set<int> used_shower_clusters;
```

## Algorithm Flow

1. **Initial Vertex Collection**
   ```cpp
   for (auto it = map_vertex_segments.begin(); it!=map_vertex_segments.end(); it++){
     WCPPID::ProtoVertex *vtx = it->first;
     if (vtx->get_cluster_id() == main_cluster->get_cluster_id() || 
         map_vertex_in_shower.find(vtx) != map_vertex_in_shower.end()){
       vertices.push_back(vtx);
     }
   }
   ```
   This step collects all vertices that are either:
   - Part of the main cluster
   - Already associated with existing showers

2. **Cluster Processing Loop**
   For each cluster that hasn't been processed yet:

   a. **Distance Calculation**
   ```cpp
   double min_dis = 1e9;
   WCPPID::ProtoVertex *min_vertex = 0;
   double main_dis;
   
   // Find closest vertex
   for (size_t i=0; i!=vertices.size(); i++){
     double dis = pcloud->get_closest_dis(vertices.at(i)->get_fit_pt());
     if (dis < min_dis){
       min_dis = dis;
       min_vertex = vertices.at(i);
     }
     if (vertices.at(i) == main_vertex)
       main_dis = dis;
   }
   ```

   b. **Main Vertex Preference**
   ```cpp
   if (min_dis > 0.8 * main_dis){
     min_dis = main_dis;
     min_vertex = main_vertex;
   }
   ```
   If the closest vertex is not significantly closer than the main vertex, prefer the main vertex.

3. **Shower Creation Process**
   For each qualifying segment:
   
   a. **Shower Initialization**
   ```cpp
   WCPPID::WCShower *shower = new WCPPID::WCShower();
   shower->set_start_vertex(min_vertex, connection_type);
   shower->set_start_segment(sg);
   ```

   b. **Direction Setting**
   - Determines the direction of the shower based on vertex relationships
   - Sets flag_dir based on distance comparisons

   c. **Particle Type Assignment**
   ```cpp
   if (sg->get_particle_type()==0 || 
       (fabs(sg->get_particle_type())==13 && sg->get_length() < 40*units::cm && sg->is_dir_weak())){
     sg->set_particle_type(11); // Set as electron
     sg->set_particle_mass(mp.get_mass_electron());
     sg->cal_4mom();
   }
   ```

4. **Shower Structure Completion**
   ```cpp
   std::set<WCPPID::ProtoSegment*> used_segments;
   shower->complete_structure_with_start_segment(map_vertex_segments, map_segment_vertices, used_segments);
   ```
   - Completes the shower structure by adding associated segments
   - Updates the relationships between vertices and segments

## Connection Types

The function uses different connection types:
- Type 3: Standard connection for nearby vertices
- Type 4: Long-distance connection (> 80 cm)

## Flow Diagram
```
[Start]
   │
   ▼
[Collect Vertices]
   │
   ▼
[For Each Unprocessed Cluster]
   │
   ├──> [Find Closest Vertex]
   │    │
   │    ├──> [Check Main Vertex Distance]
   │    │
   │    ▼
   ├──> [Create New Shower]
   │    │
   │    ├──> [Set Direction]
   │    │
   │    ├──> [Assign Particle Type]
   │    │
   │    ▼
   └──> [Complete Shower Structure]
   
[End]
```

## Key Considerations

1. **Distance Thresholds**
   - 80 cm threshold for standard connections
   - Different connection types based on distance

2. **Particle Type Conversion**
   - Converts certain muon-like segments to electron-like if they meet specific criteria
   - Length < 40 cm
   - Weak directionality

3. **Vertex Preference**
   - Prefers main vertex if distance difference is less than 20%
   - Helps maintain coherent event structure

## Function Calls

The function makes several important calls to other methods in the codebase:

1. **Point Cloud Operations**
   - `cluster->get_point_cloud()`: Gets the point cloud representation of a cluster
   - `pcloud->get_closest_dis(point)`: Calculates closest distance between points

2. **Vertex/Segment Operations**
   - `get_start_end_vertices(segment)`: Gets start and end vertices for a segment
   - `map_vertex_segments[vertex]`: Accesses segments connected to a vertex
   - `find_other_vertex(segment, vertex)`: Finds the opposite vertex of a segment 

3. **Shower Management** [WCShower](../wcshower.md)
   - `shower->set_start_vertex(vertex, connection_type)`: Sets initial vertex for shower
   - `shower->set_start_segment(segment)`: Sets initial segment for shower
   - `shower->set_start_point(point)`: Sets starting point for shower
   - `shower->complete_structure_with_start_segment(...)`: Builds complete shower structure
   - `shower->add_segment(segment, map_segment_vertices)`: Adds segment to shower

4. **Segment Properties** [ProtoSegment](../protosegment.md)
   - `segment->get_length()`: Gets segment length
   - `segment->get_cluster_id()`: Gets cluster ID
   - `segment->get_particle_type()`: Gets particle type
   - `segment->set_particle_type(type)`: Sets particle type
   - `segment->set_particle_mass(mass)`: Sets particle mass
   - `segment->cal_4mom()`: Calculates 4-momentum
   - `segment->get_point_vec()`: Gets vector of points
   - `segment->is_dir_weak()`: Checks if direction is weakly determined

5. **Map Updates**
   - `update_shower_maps()`: Updates global shower mapping relationships [more details](./update_shower_maps.md)
   - Used to maintain consistency after shower modifications

6. **TPCParams Access**
   - `Singleton<TPCParams>::Instance()`: Gets TPC parameters instance
   - Used for accessing physical constants like electron mass

## Usage Example

```cpp
// Example usage in neutrino event reconstruction
WCPPID::NeutrinoID neutrino_id(...);
neutrino_id.shower_clustering_in_other_clusters(true);  // true for saving results

// The function will:
// 1. Identify shower-like segments in other clusters
// 2. Create new showers and associate them with appropriate vertices
// 3. Update global shower maps and relationships
```