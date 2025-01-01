# Understanding improve_maps Functions in NeutrinoID

The improve_maps functions in NeutrinoID are a set of algorithms that handle particle track and shower classification and direction determination in neutrino interactions. Here's a detailed breakdown of the key functions and their logic.

## 1. improve_maps_one_in

This function enforces the physics constraint that typically only one particle should enter a vertex in neutrino interactions.

### Main Algorithm Flow

#### Initialization
```cpp
bool flag_update = true;
std::set<WCPPID::ProtoVertex* > used_vertices;
std::set<WCPPID::ProtoSegment* > used_segments;
```

#### Main Processing Loop
```cpp
while(flag_update) {
    flag_update = false;
    for (auto it = map_vertex_segments.begin(); 
         it!=map_vertex_segments.end(); 
         it++) {
        // Process each vertex
        WCPPID::ProtoVertex *vtx = it->first;
        if (used_vertices.find(vtx)!=used_vertices.end()) 
            continue;
```

### Key Components

#### 1. Incoming Particle Counter
```cpp
int n_in = 0;
std::map<WCPPID::ProtoSegment*, bool> map_sg_dir;
for (auto it1 = it->second.begin(); 
     it1 != it->second.end(); 
     it1++) {
    WCPPID::ProtoSegment *sg = (*it1);
    if (used_segments.find(sg) != used_segments.end()) 
        continue;
    
    bool flag_start = determine_start_flag(sg, vtx);
    
    if ((flag_start && sg->get_flag_dir()==-1 || 
        (!flag_start) && sg->get_flag_dir()==1)) {
        if (flag_strong_check) {
            if (!sg->is_dir_weak()) n_in ++;
        } else {
            n_in ++;
        }
    }
}
```

#### 2. Direction Assignment Logic
```cpp
if (n_in > 0) {
    // Update directions for undetermined segments
    for (auto it1 = map_sg_dir.begin(); 
         it1!=map_sg_dir.end(); 
         it1++) {
        WCPPID::ProtoSegment *sg = it1->first;
        bool flag_start = it1->second;
        
        // Set outgoing direction
        if (flag_start) {
            sg->set_flag_dir(1);
        } else {
            sg->set_flag_dir(-1);
        }
        sg->cal_4mom();
        sg->set_dir_weak(true);
        used_segments.insert(sg);
        flag_update = true;
    }
}
```

### Special Cases

#### 1. Strong Direction Check
When flag_strong_check is true:
- Only considers segments with strong directional evidence
- Requires more confidence before making changes
```cpp
if (flag_strong_check && !sg->is_dir_weak()) {
    n_in++;
} else if (!flag_strong_check) {
    n_in++;
}
```

#### 2. Vertex Skipping
```cpp
if (map_sg_dir.size() == 0) {
    used_vertices.insert(vtx); // No segments to update
}
```

### Update Propagation
1. Changes are made one vertex at a time
2. Updates can trigger further updates at connected vertices
3. Process continues until no more changes are needed

### Termination Conditions
- No more direction updates possible
- All relevant vertices processed
- All undetermined segments assigned directions

## 2. improve_maps_shower_in_track_out

This function corrects unphysical configurations where electromagnetic showers appear to enter a vertex while tracks exit. Such configurations typically indicate incorrect direction assignments.

### Main Algorithm Flow

#### Initialization
```cpp
for (auto it = map_vertex_segments.begin(); 
     it!=map_vertex_segments.end();
     it++) {
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster_id) continue;
    if (it->second.size()==1) continue;
```

#### Segment Analysis Per Vertex
```cpp
std::vector<WCPPID::ProtoSegment*> in_showers;
bool flag_turn_shower_dir = false;

// Analyze each segment at vertex
for (auto it1 = it->second.begin(); 
     it1 != it->second.end(); 
     it1++) {
    WCPPID::ProtoSegment *sg = (*it1);
    bool flag_start = determine_segment_direction(sg, vtx);
    
    // Check incoming showers
    if (is_incoming(flag_start, sg)) {
        if (sg->get_flag_shower()) {
            in_showers.push_back(sg);
        }
    }
    // Check outgoing tracks
    else if (is_outgoing(flag_start, sg)) {
        if (!sg->get_flag_shower()) {
            if (!sg->is_dir_weak()) 
                flag_turn_shower_dir = true;
        }
    }
}
```

### Key Decision Points

#### 1. Configuration Detection
```cpp
bool should_reverse_showers(ProtoSegment* sg, 
                          ProtoVertex* vtx) {
    // Check if we have:
    // 1. Incoming shower(s)
    // 2. Outgoing track(s)
    // 3. Strong directional evidence on tracks
    return has_incoming_showers(sg, vtx) && 
           has_outgoing_tracks(sg, vtx) &&
           tracks_have_strong_direction(sg, vtx);
}
```

#### 2. Direction Reversal Logic
```cpp
if (flag_turn_shower_dir) {
    for (auto it1 = in_showers.begin(); 
         it1!= in_showers.end(); 
         it1++) {
        // Reverse shower direction
        (*it1)->set_flag_dir((*it1)->get_flag_dir() * (-1));
        (*it1)->set_dir_weak(true);
    }
}
```

### Special Cases

#### 1. Single Segment Vertices
- Skipped as no direction conflict possible
```cpp
if (it->second.size()==1) continue;
```

#### 2. Strong Direction Tracks
- Only reverse showers if tracks have strong directional evidence
```cpp
if (!sg->is_dir_weak()) {
    flag_turn_shower_dir = true;
}
```

### Physics Constraints Enforced

1. **Shower Development**
   - EM showers typically develop away from vertex
   - Direction should align with energy flow

2. **Track-Shower Relationships**
   - Tracks more reliable for direction
   - Showers adapt to track directions

3. **Energy Conservation**
   - Energy flow should be consistent
   - Maintains causal relationships

### Update Propagation
- Changes are local to vertex
- No cascading updates needed
- One-time direction correction

## 3. improve_maps_multiple_tracks_in

This function handles the unphysical case of multiple tracks entering a single vertex, typically indicating misidentified particles or incorrect direction assignments.

### Main Algorithm Flow

#### Initialization and Control Loop
```cpp
bool flag_update = true;
std::set<WCPPID::ProtoVertex* > used_vertices;
std::set<WCPPID::ProtoSegment* > used_segments;

while(flag_update) {
    flag_update = false;
```

#### Vertex Analysis
```cpp
for (auto it = map_vertex_segments.begin(); 
     it!=map_vertex_segments.end();
     it++) {
    WCPPID::ProtoVertex *vtx = it->first;
    if (vtx->get_cluster_id() != temp_cluster_id) continue;
    if (it->second.size()==1) continue;
    if (used_vertices.find(vtx)!=used_vertices.end()) continue;
    
    // Count incoming particles
    int n_in = 0;
    int n_in_shower = 0;
    std::vector<ProtoSegment* > in_tracks;
```

### Key Components

#### 1. Track Counter and Classifier
```cpp
for (auto it1 = it->second.begin(); 
     it1 != it->second.end(); 
     it1++) {
    WCPPID::ProtoSegment *sg = (*it1);
    bool flag_start = determine_start_flag(sg, vtx);
    
    if (is_incoming(flag_start, sg)) {
        n_in++;
        if (sg->get_flag_shower()) {
            n_in_shower++;
        } else {
            in_tracks.push_back(sg);
        }
    }
}
```

#### 2. Reclassification Logic
```cpp
if (n_in > 1 && n_in != n_in_shower) {
    for (auto it1 = in_tracks.begin(); 
         it1!=in_tracks.end(); 
         it1++) {
        WCPPID::ProtoSegment *sg1 = *it1;
        
        // Reclassify as electron
        sg1->set_particle_type(11);
        sg1->set_particle_mass(mp.get_mass_electron());
        if (sg1->get_particle_4mom(3)>0) 
            sg1->cal_4mom();
        
        flag_update = true;
    }
}
```

### Special Cases

#### 1. Multiple Shower Exemption
```cpp
if (n_in == n_in_shower) {
    // Multiple incoming showers allowed
    continue;
}
```

#### 2. Single Track Processing
```cpp
if (n_in == 1) {
    // Single incoming track is physically valid
    continue;
}
```

### Physics Rules Enforced

1. **Conservation Laws**
   - Enforce charge conservation
   - Maintain energy-momentum conservation

2. **Topology Constraints**
   - Single incoming track preferred
   - Multiple incoming tracks suggest misidentification

3. **Particle Identification**
   - Tracks reclassified as electrons when topology suggests shower
   - Maintains consistent event interpretation

### Update Management

#### Track Status Tracking
```cpp
used_vertices.insert(vtx);
// Update tracking for modified segments
for (auto sg : modified_segments) {
    used_segments.insert(sg);
}
```

#### Update Propagation
```cpp
if (flag_update) {
    // Continue iteration if changes made
    // Allows cascade effects through event
}
```

Example scenario showing transformation:
```
Before:                   After:
  μ-->                    e-->
  p-->  Vertex  π-->      e-->  Vertex  π-->
  μ-->                    e-->
```

### Termination Conditions
1. No more multiple incoming tracks
2. All vertices processed
3. All reclassifications complete

## 4. improve_maps_no_dir_tracks

This complex function handles tracks with no clear direction, implementing a comprehensive algorithm for particle identification and direction determination.

### Main Algorithm Flow

#### Initialization Phase:
1. Creates sets to track processed vertices and segments
2. Initializes flag_update for iterative processing
3. Begins main processing loop that continues until no more updates are needed

#### Per-Segment Processing:
For each segment in the cluster:
```cpp
// Check if segment needs processing
if (sg->get_flag_dir()==0 || sg->is_dir_weak() || 
    fabs(sg->get_particle_type())==2212) {
    // Process segment
}
```

### Classification Rules

#### 1. Multiple Shower Connection Rule
```cpp
if (nshowers[0] + nshowers[1] >2 && length < 5*units::cm ||
    nshowers[0]+1 == map_vertex_segments[two_vertices.first].size() && 
    (nshowers[1]+1 == map_vertex_segments[two_vertices.second].size()) && 
    nshowers[0]>0 && nshowers[1]>0 && length<5*units::cm) {
    // Reclassify as electron
    sg->set_particle_type(11);
    sg->set_particle_mass(mp.get_mass_electron());
}
```

#### 2. Shower End Point Rule
When a segment connects to vertex with multiple showers:
```cpp
// Check shower endpoint conditions
if (nshowers[0]+1 == map_vertex_segments[two_vertices.first].size() && 
    nshowers[0]>=2 && sg->get_particle_type() == 2212) {
    // Additional angle checks
    TVector3 v1 = sg->cal_dir_3vector(two_vertices.first->get_fit_pt(), 5*units::cm);
    double min_angle = calculate_minimum_angle(v1, two_vertices.first);
    
    if (meets_angle_criteria(min_angle, dQ_dx_rms)) {
        reclassify_as_electron(sg);
    }
}
```

#### 3. Long Track Analysis
For longer tracks (>10cm):
```cpp
if (current_sg->get_particle_type()==11 && 
    num_daughter_showers<=2 && 
    (!flag_shower_in) && 
    (!current_sg->get_flag_shower_topology()) && 
    (!current_sg->get_flag_shower_trajectory()) && 
    length > 10*units::cm && 
    (!flag_only_showers)) {
    
    // Check if should be reclassified as muon
    if (current_sg->get_particle_score()<=100) {
        double direct_length = current_sg->get_direct_length();
        if (direct_length >= 34*units::cm || 
            (direct_length < 34*units::cm && direct_length > 0.93 * length)) {
            current_sg->set_particle_type(13);
            current_sg->set_particle_mass(mp.get_mass_muon());
        }
    }
}
```

### Key Decision Points

#### 1. Energy Deposition Analysis
- Checks dQ/dx patterns:
  ```cpp
  if (sg->get_medium_dQ_dx()/(43e3/units::cm) > 1.4) {
      // High dQ/dx suggests proton
      sg->set_particle_type(2212);
  } else {
      // Lower dQ/dx suggests muon
      sg->set_particle_type(13);
  }
  ```

#### 2. Topology Analysis
Examines connected vertices and segments:
- Counts number of connected showers/tracks
- Analyzes angles between segments
- Checks for vertex activity
- Evaluates segment length ratios

#### 3. Angular Analysis
```cpp
bool check_angle_criteria(TVector3 dir1, TVector3 dir2, double length) {
    double angle = dir1.Angle(dir2)/3.1415926*180.;
    if (length > 50*units::cm) {
        // Different criteria for long segments
        TVector3 dir3 = calculate_long_direction(dir1, 50);
        TVector3 dir4 = calculate_long_direction(dir2, 50);
        angle = compare_angles(dir3, dir4);
    }
    return evaluate_angle_threshold(angle, length);
}
```

### Special Considerations

1. **Short Track Handling**
   - Tracks < 5cm get special treatment
   - More likely to be reclassified as showers
   - Higher sensitivity to shower connections

2. **Long Track Stability**
   - Tracks > 35cm maintain classification more readily
   - Requires stronger evidence for reclassification
   - Direction determination more stable

3. **Vertex Activity**
   - Considers activity at both endpoints
   - Weighs multiple particle connections
   - Special handling for vertex regions

### Update Propagation

The function implements iterative updates:
```cpp
while(flag_update) {
    flag_update = false;
    // Process all segments
    for (auto it = map_segment_vertices.begin(); 
         it != map_segment_vertices.end(); 
         it++) {
        // If any updates made, flag_update set true
        // Continue until no more updates needed
    }
}
```

This ensures all changes are properly propagated through the event structure.

## Common Patterns Across Functions

All improve_maps functions share some common patterns:

1. **Iterative Processing**
   - Continue until no more changes can be made
   - Each pass may enable more improvements

2. **Confidence Levels**
   - Strong vs weak directional information
   - Different thresholds for different scenarios

3. **Physics Constraints**
   - Conservation of charge
   - Particle behavior patterns
   - Energy deposition patterns

## Important Variables

- `flag_dir`: Direction flag (-1: incoming, 0: unknown, 1: outgoing)
- `dir_weak`: Boolean indicating confidence in direction
- `particle_type`: PDG code (11: electron, 13: muon, 2212: proton)
- `map_vertex_segments`: Maps vertices to connected segments
- `map_segment_vertices`: Maps segments to connected vertices

## Key Considerations

1. **Physics Consistency**
   - Maintain conservation laws
   - Respect particle behavior patterns
   - Consider detector effects

2. **Algorithmic Robustness**
   - Handle edge cases gracefully
   - Prevent infinite loops
   - Maintain consistency across updates

3. **Performance**
   - Efficient graph traversal
   - Minimal redundant calculations
   - Clear termination conditions

These functions work together to improve particle track and shower identification in complex neutrino interactions by enforcing physical constraints and topology rules.

## Function Dependencies

### improve_maps_one_in
Called functions:
- `find_other_vertex(ProtoSegment*, ProtoVertex*)`: Finds the other vertex connected to a segment
- `cal_4mom()`: Calculates 4-momentum of segment
- `set_flag_dir()`: Sets direction flag
- `set_dir_weak()`: Sets confidence in direction determination

### improve_maps_shower_in_track_out
Called functions:
- `get_cluster_id()`: Gets cluster identifier
- `get_wcpt_vec().front()/back()`: Gets endpoint wireCell points
- `get_wcpt().index`: Gets index of vertex point
- `get_flag_dir()`: Gets direction flag
- `is_dir_weak()`: Checks if direction determination is weak
- `get_flag_shower()`: Checks if segment is shower-like
- `set_flag_dir()`: Sets direction flag
- `set_dir_weak()`: Sets confidence in direction determination

### improve_maps_multiple_tracks_in  
Called functions:
- `get_cluster_id()`: Gets cluster identifier
- `get_wcpt_vec().front()/back()`: Gets endpoint wireCell points
- `get_wcpt().index`: Gets index of vertex point
- `get_flag_dir()`: Gets direction flag
- `get_flag_shower()`: Checks if segment is shower-like
- `set_particle_type()`: Sets particle PDG code
- `set_particle_mass()`: Sets particle mass
- `cal_4mom()`: Calculates 4-momentum

### improve_maps_no_dir_tracks
Called functions:
- `calculate_num_daughter_showers()`: Calculates number of daughter showers [see details](./calculate_num_daughter_showers.md)
- `calculate_num_daughter_tracks()`: Calculates number of daughter tracks
- `get_medium_dQ_dx()`: Gets average dQ/dx
- `cal_dir_3vector()`: Calculates direction vector
- `get_fit_pt()`: Gets fitted point
- `find_vertices()`: Finds vertices connected to segment
- `get_length()`: Gets segment length
- `get_direct_length()`: Gets direct length between endpoints
- `set_particle_type()`: Sets particle PDG code
- `set_particle_mass()`: Sets particle mass
- `set_flag_dir()`: Sets direction flag
- `set_dir_weak()`: Sets confidence in direction determination
- `cal_4mom()`: Calculates 4-momentum
- `get_rms_dQ_dx()`: Gets RMS of dQ/dx
- `get_particle_score()`: Gets particle identification score
- `get_closest_2d_dis()`: Gets closest 2D distance to points

### Helper Functions Called Across Multiple improve_maps Functions:
- `find_other_vertex()`: Used to navigate connected vertices [see_details](../NeutrinoID/find_other_vertex.md)
- `cal_4mom()`: Used to update particle kinematics
- `set_flag_dir()`: Used to update direction information
- `set_dir_weak()`: Used to update confidence in direction
- `get_cluster_id()`: Used for cluster identification
- `get_length()`: Used for segment length calculations