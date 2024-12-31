# WCShower Class Documentation

## Overview
The WCShower class is designed to represent and analyze shower-like particle trajectories in a Wire-Cell detector. It handles both shower and track-like particles, managing their geometry, kinematics, and various properties through a connected network of vertices and segments.

## Class Members

### Key Data Members

#### Particle Properties
- `particle_type`: Integer identifying the type of particle
- `flag_shower`: Boolean indicating if the structure is shower-like
- `flag_kinematics`: Boolean indicating if kinematics have been calculated

#### Kinematic Properties
- `kenergy_range`: Kinetic energy calculated from range
- `kenergy_dQdx`: Kinetic energy calculated from dQ/dx
- `kenergy_charge`: Kinetic energy calculated from charge
- `kenergy_best`: Best estimate of kinetic energy
- `start_point`, `end_point`: WCP::Point objects defining trajectory endpoints
- `init_dir`: TVector3 representing initial direction

#### Structural Components
- `start_vertex`: Starting ProtoVertex pointer
- `start_connection_type`: Integer defining connection type
  - 1: Direct connection
  - 2: Indirect connection with gap
  - 3: Associated connection (ambiguous)
- `start_segment`: Starting ProtoSegment pointer
- `map_vtx_segs`: Map of vertices to connected segments
- `map_seg_vtxs`: Map of segments to connected vertices

#### Point Clouds
- `pcloud_fit`: ToyPointCloud for fitted points
- `pcloud_associated`: ToyPointCloud for associated points

## Major Functions

### Initialization and Setup
```cpp
void set_start_vertex(ProtoVertex* vertex, int type);
void set_start_segment(ProtoSegment* seg);
void set_start_point(WCP::Point p);
add_segment(...);
add_shower(...);
```
These functions initialize the shower structure with starting points and segments.

[add segment or shower](./pattern_recognition/wcshower_addition.md)

### Kinematics Calculation

The `calculate_kinematics()` function is a central piece that computes particle properties. The logic flow is:

1. Single Segment Case:
   ```mermaid
   graph TD
   A[Start] --> B{Single Segment?}
   B -->|Yes| C[Calculate Basic Properties]
   C --> D[Set Energy Based on Length]
   D --> E[Set Direction]
   B -->|No| F[Multiple Segment Processing]
   F --> G[Calculate Connected Components]
   G --> H[Determine Shower Properties]
   H --> I[Calculate Energies]
   ```

2. Multiple Segment Case:
   - Analyzes connectivity between segments
   - Determines if structure is track-like or shower-like
   - Computes cumulative properties

[see details](./pattern_recognition/wcshower_kinematics.md)

### Point Cloud Management [more details](./pattern_recognition/wcshower_point_clouds.md)
```cpp
void build_point_clouds();
void rebuild_point_clouds();
```
These functions manage point cloud representations of the shower structure:
- Creates separate clouds for fitted and associated points
- Builds KD-tree indices for efficient spatial queries
- Handles point addition and updates

### Structure Analysis

#### Connection Analysis [more details](./pattern_recognition/wcshower_get_info.md)
```cpp
std::pair<std::set<ProtoSegment*>, std::set<ProtoVertex*>> 
get_connected_pieces(ProtoSegment* seg);
```
Analyzes the connectivity of segments and vertices:
1. Starts from given segment
2. Recursively explores connected components
3. Returns sets of connected segments and vertices

#### Distance Calculations [more details](./pattern_recognition/wcshower_get_distance.md)
```cpp
double get_closest_dis(ProtoSegment* seg);
std::pair<double, WCP::Point> get_closest_point(WCP::Point& p);
```
Provides spatial analysis capabilities:
- Finds nearest points between segments
- Calculates distances between shower components
- Supports trajectory analysis

### Shower/Track Properties

#### Length Calculations [more details](./pattern_recognition/wcshower_get_length.md)
```cpp
double get_total_length();
double get_total_track_length();
```
Computes various length measurements:
- Total path length of all segments
- Length of track-like segments only
- Specific cluster lengths

#### dQ/dx Analysis [more details](./pattern_recognition/wcshower_get_info.md)
```cpp
std::vector<double> get_stem_dQ_dx(ProtoVertex* vertex, 
                                  ProtoSegment* sg, 
                                  int limit = 20);
```
Analyzes charge deposition:
- Calculates dQ/dx along segments
- Handles multiple segment transitions
- Supports particle identification

## get_num_main_segments Function [more details](./pattern_recognition/wcshower_get_info.md)

### Purpose
Counts the number of segments that belong to the same cluster as the start segment.

### Function Signature
```cpp
int get_num_main_segments()
```
## update_particle_type Function [more details](./pattern_recognition/wcshower_get_info.md)


### Purpose
Updates the particle type based on the characteristics of connected segments.

### Function Signature
```cpp
void update_particle_type()
```

## complete_structure_with_start_segment() [more details](./pattern_recognition/wcshower_pattern.md)

### Purpose
Builds a complete graph structure starting from an initial segment, connecting related vertices and segments.

## get_last_segment_vertex_long_muon() [more details](./pattern_recognition/wcshower_pattern.md)

### Purpose
Identifies the last segment and vertex in a muon track by traversing the connected segments.

## cal_dir_3vector() [more details](./pattern_recognition/wcshower_pattern.md)

### Purpose
Calculates the direction vector at a given point by averaging nearby points within a distance cut.

## Usage Example

```cpp
WCPPID::WCShower shower;

// Initialize
shower.set_start_vertex(vertex_ptr, 1);
shower.set_start_segment(segment_ptr);

// Calculate properties
shower.calculate_kinematics();

// Access results
if (shower.get_flag_shower()) {
    double energy = shower.get_kine_best();
    TVector3 direction = shower.get_init_dir();
    // Process shower properties
}
```

## Notes and Considerations

1. Memory Management
   - Class owns and manages point cloud objects
   - Proper cleanup in destructor
   - Careful handling of pointers to vertices and segments

2. Performance Considerations
   - Uses KD-trees for efficient spatial queries
   - Caches computed properties
   - Rebuilds point clouds only when necessary

3. Algorithm Selection
   - Adapts calculations based on structure type
   - Uses different energy calculations depending on particle type
   - Handles both track-like and shower-like patterns

4. Error Handling
   - Validates input parameters
   - Handles edge cases in connectivity
   - Provides fallback calculations when optimal methods fail