# WCShower Class Functions Documentation

## Table of Contents
1. [get_stem_dQ_dx Function](#get_stem_dqdx)
2. [get_num_main_segments Function](#get_num_main_segments)
3. [get_connected_pieces Function](#get_connected_pieces)
4. [update_particle_type Function](#update_particle_type)

## get_stem_dQ_dx Function <a name="get_stem_dqdx"></a>

### Purpose
Calculates the dQ/dx (charge deposition per unit length) values for a given vertex and segment combination, with an optional limit on the number of values to return.

### Function Signature
```cpp
std::vector<double> get_stem_dQ_dx(WCPPID::ProtoVertex *vertex, 
                                  WCPPID::ProtoSegment *sg, 
                                  int limit = 20)
```

### Logic Flow

```mermaid
graph TD
    A[Start] --> B{Check vertex position}
    B -->|Vertex at front| C[Process from start]
    B -->|Vertex at back| D[Process from end]
    C --> E[Calculate dQ/dx]
    D --> E
    E --> F{Is sg == start_segment?}
    F -->|No| G[Return results]
    F -->|Yes| H{Results < limit?}
    H -->|No| G
    H -->|Yes| I[Find next vertex & segment]
    I --> J{Found next?}
    J -->|No| G
    J -->|Yes| K[Add more dQ/dx values]
    K --> L{Results >= limit?}
    L -->|Yes| G
    L -->|No| I

```

### Key Steps
1. Initial dQ/dx calculation for given segment
   - Gets dQ and dx vectors from segment
   - Calculates dQ/dx based on vertex position (front/back)
   - Normalizes by detector constant (43e3/units::cm)

2. Extension for start segment (if needed)
   - Finds connected vertices and segments
   - Checks angle conditions (< 25 degrees)
   - Adds additional dQ/dx values up to limit
   - Maximum of 3 connected segments checked

## get_num_main_segments Function <a name="get_num_main_segments"></a>

### Purpose
Counts the number of segments that belong to the same cluster as the start segment.

### Function Signature
```cpp
int get_num_main_segments()
```

### Logic Flow

```mermaid
graph TD
    A[Start] --> B[Initialize count = 0]
    B --> C[Iterate through map_seg_vtxs]
    C --> D{Current segment cluster_id == start_segment cluster_id?}
    D -->|Yes| E[Increment count]
    D -->|No| F[Next segment]
    E --> F
    F --> G{More segments?}
    G -->|Yes| C
    G -->|No| H[Return count]
```

## get_connected_pieces Function <a name="get_connected_pieces"></a>

### Purpose
Returns sets of connected segments and vertices starting from a given segment.

### Function Signature
```cpp
std::pair<std::set<WCPPID::ProtoSegment*>, std::set<WCPPID::ProtoVertex*>> 
get_connected_pieces(WCPPID::ProtoSegment* tseg)
```

### Logic Flow

```mermaid
graph TD
    A[Start] --> B[Initialize sets]
    B --> C[Add initial segment]
    C --> D{New vertices or segments?}
    D -->|Yes| E{Process new vertex?}
    E -->|Yes| F[Get connected segments]
    E -->|No| G[Process new segment]
    F --> H[Add new segments]
    G --> I[Get connected vertices]
    I --> J[Add new vertices]
    H --> D
    J --> D
    D -->|No| K[Return sets]
```

### Implementation Details
1. Uses two vectors for processing:
   - new_segments: Segments to be processed
   - new_vertices: Vertices to be processed
   
2. Uses two sets for results:
   - used_segments: All connected segments
   - used_vertices: All connected vertices

3. Processing continues until no new connections found

## update_particle_type Function <a name="update_particle_type"></a>

### Purpose
Updates the particle type based on the characteristics of connected segments.

### Function Signature
```cpp
void update_particle_type()
```

### Logic Flow

```mermaid
graph TD
    A[Start] --> B{Multiple segments?}
    B -->|No| C[End]
    B -->|Yes| D[Calculate lengths]
    D --> E[Sum shower_length]
    D --> F[Sum track_length]
    E --> G{shower_length > track_length?}
    F --> G
    G -->|Yes| H[Set particle type = electron]
    G -->|No| C
    H --> I[Update particle mass]
    I --> J[Recalculate 4-momentum]
    J --> C
```

### Key Points
1. Only processes if multiple segments exist
2. Classifies segments into:
   - Shower-like: `seg->get_flag_shower() || fabs(seg->get_particle_type())!=2212`
   - Track-like: Otherwise

3. Updates properties if shower-like characteristics dominate:
   - Sets particle type to electron (11)
   - Updates mass to electron mass
   - Recalculates 4-momentum

### Related Constants
```cpp
// Electron particle type
const int ELECTRON_TYPE = 11;

// TPCParams electron mass
// Accessed via: mp.get_mass_electron()
```