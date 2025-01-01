# Understanding calc_conflict_maps Function

The `calc_conflict_maps` function is designed to evaluate the consistency and conflicts in a vertex-segment network structure in neutrino interaction reconstruction. It calculates a conflict score based on various geometric and topological criteria.

## Overview

The function takes a `ProtoVertex` as input and returns a float value representing the total number of conflicts detected. A higher conflict score indicates more inconsistencies in the reconstructed topology.

## Key Data Structures

```cpp
std::map<WCPPID::ProtoSegment*, std::pair<ProtoVertex*, ProtoVertex*>> map_seg_dir;
std::set<WCPPID::ProtoVertex*> used_vertices;
```

- `map_seg_dir`: Maps segments to their start and end vertices
- `used_vertices`: Tracks visited vertices during traversal

## Main Processing Steps

### 1. Graph Traversal

The function starts from the input vertex and traverses the connected segments and vertices:

```cpp
std::vector<std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoSegment*>> segments_to_be_examined;
// Start with all segments connected to input vertex
for (auto it = map_vertex_segments[temp_vertex].begin(); it != map_vertex_segments[temp_vertex].end(); it++) {
    segments_to_be_examined.push_back(std::make_pair(temp_vertex, *it));
}
used_vertices.insert(temp_vertex);
```

### 2. Conflict Detection

The function checks for several types of conflicts:

#### a. Direction Conflicts (Weight: 0.5-1.0)
```cpp
if (sg->get_flag_dir()!=0 && (sg->get_flag_shower() && sg->get_length()>5*units::cm 
    || (!sg->get_flag_shower()))) {
    if (flag_start && sg->get_flag_dir()==-1 ||
        (!flag_start) && sg->get_flag_dir() == 1) {
        if (!sg->is_dir_weak()) num_conflicts++;
        else num_conflicts += 0.5;
    }
}
```

- Checks if segment directions are consistent with their connection to vertices
- Higher penalty for strong directional conflicts
- Lower penalty for weak directional conflicts

#### b. Vertex-Level Topology Conflicts

For each vertex, the function checks:

1. Multiple Incoming Particles (Weight: 1.0 or 0.5)
```cpp
if (n_in > 1) {
    if (n_in != n_in_shower) {
        num_conflicts += (n_in-1);
    } else {
        num_conflicts += (n_in-1)/2.;
    }
}
```

2. Shower-Track Mixing Conflicts (Weight: 1.0)
```cpp
if (n_in_shower > 0 && n_out_tracks > 0) {
    num_conflicts += std::min(n_in_shower, n_out_tracks);
}
```

#### c. Angular Conflicts

The function evaluates angular relationships between segments:

```cpp
if (max_angle < 35) {
    num_conflicts += 5;
} else if (max_angle < 70) {
    num_conflicts += 3;
} else if (max_angle < 85) {
    num_conflicts += 1;
} else if (max_angle < 110) {
    num_conflicts += 0.25;
}
```

Penalties based on angle ranges:
- < 35°: Very severe conflict (5 points)
- 35-70°: Severe conflict (3 points)
- 70-85°: Moderate conflict (1 point)
- 85-110°: Minor conflict (0.25 points)

### 3. Beam Direction Penalties

Additional penalties for segments going backwards relative to beam direction:

```cpp
if (angle_beam < 60 && max_angle < 110) {
    num_conflicts += 1;
} else if (angle_beam < 45 && max_angle < 70) {
    num_conflicts += 3;
}
```

## Example Scenarios

### 1. Clean Topology
```
Vertex
  ↑ (incoming track)
  |
  ↓ (outgoing track at 160°)
  |
```
Low conflict score: The angles and directions are consistent.

### 2. Problematic Topology
```
Vertex
  ↑ (incoming track)
  |
  ← (outgoing track at 90°)
  ↑ (another incoming track)
```
High conflict score due to:
- Multiple incoming tracks
- Sharp angle between segments
- Inconsistent directions

## Function Dependencies

### Called Functions [ProtoSegment](../protosegment.md)

1. `find_vertices(ProtoSegment* sg)` 
   - Returns pair of vertices (start and end) connected to a segment
   - Used for topology reconstruction
   - Returns: `std::pair<WCPPID::ProtoVertex*, WCPPID::ProtoVertex*>`

2. `cal_dir_3vector(Point& pt, double dis)`
   - Calculates direction vector at a point with given distance
   - Used for angular calculations
   - Returns: `TVector3`

3. `get_flag_shower()`
   - Checks if segment is shower-like
   - Returns: `bool`

4. `get_flag_shower_trajectory()`
   - More specific shower check based on trajectory
   - Returns: `bool`

5. `get_flag_dir()`
   - Gets direction flag of segment (-1, 0, or 1)
   - Returns: `int`

6. `is_dir_weak()`
   - Checks if direction determination is weak/uncertain
   - Returns: `bool`

7. `get_length()`
   - Gets length of segment
   - Returns: `double`

### Helper Functions Used In Maps

1. `map_vertex_segments`
   - Maps vertices to connected segments
   - Type: `std::map<ProtoVertex*, ProtoSegmentSet>`

2. `map_segment_vertices`
   - Maps segments to connected vertices
   - Type: `std::map<ProtoSegment*, ProtoVertexSet>`

### Core Math Functions

1. `TVector3::Angle()`
   - Calculates angle between two vectors
   - From ROOT framework
   - Returns: `double` (in radians)

## Usage in Neutrino Reconstruction

The function is primarily used to:
1. Evaluate candidate neutrino vertices
2. Identify problematic topological configurations
3. Help select the most likely neutrino interaction vertex

The conflict score helps discriminate between well-reconstructed vertices and those with inconsistent track/shower patterns.