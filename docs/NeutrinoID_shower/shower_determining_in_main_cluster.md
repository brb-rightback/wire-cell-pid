# Analysis of shower_determining_in_main_cluster Function

## Overview
The `shower_determining_in_main_cluster` function is a key component in neutrino event reconstruction, specifically handling the identification and classification of particle showers within the main cluster of detector hits. Its primary purpose is to analyze track and shower patterns to determine particle types and their interactions.

## Function Flow

1. Initial Setup and Track Examination
   ```cpp
   examine_good_tracks(temp_cluster->get_cluster_id());
   ```
   - First evaluates all tracks in the cluster to identify high-quality tracks
   - Good tracks are typically long, straight tracks with consistent energy deposition

2. Multiple Track Handling
   ```cpp
   fix_maps_multiple_tracks_in(temp_cluster->get_cluster_id());
   ```
   - Handles cases where multiple tracks enter the same vertex
   - Sets tracks as "undetermined" if their direction or particle type is ambiguous

3. Shower-Track Direction Fix
   ```cpp
   fix_maps_shower_in_track_out(temp_cluster->get_cluster_id());
   ```
   - Analyzes cases where a shower enters and a track exits
   - May reverse shower direction based on topology

4. One-In Many-Out Processing
   ```cpp
   improve_maps_one_in(temp_cluster);
   ```
   - If one clear track/shower enters, helps determine direction of other tracks
   - Important for neutrino interaction vertex identification

## Key Steps in Track/Shower Determination

### 1. Shower-Track Conversion
```cpp
improve_maps_shower_in_track_out(temp_cluster->get_cluster_id());
```
- Converts tracks to showers if they're connected to a primary shower
- Uses topology and energy deposition patterns
- May iterate multiple times to refine classifications

### 2. Direction Determination
```cpp
improve_maps_no_dir_tracks(temp_cluster->get_cluster_id());
```
- Determines track directions based on:
  * Connection to vertex
  * Energy deposition patterns
  * Topology relative to other tracks/showers

### 3. Multiple Track Analysis
```cpp
improve_maps_multiple_tracks_in(temp_cluster->get_cluster_id());
```
- Handles complex vertex topologies
- May convert tracks to showers based on:
  * Angular distribution
  * Energy patterns
  * Interaction likelihood

## Logical Flow Diagram

```
┌────────────────────┐
│ Initial Track      │
│ Examination        │
└─────────┬──────────┘
          ▼
┌────────────────────┐
│ Fix Multiple       │
│ Tracks In          │
└─────────┬──────────┘
          ▼
┌────────────────────┐
│ Fix Shower/Track   │
│ Directions         │
└─────────┬──────────┘
          ▼
┌────────────────────┐
│ Improve One-In     │
│ Many-Out           │
└─────────┬──────────┘
          ▼
┌────────────────────┐
│ Convert Ambiguous  │
│ Tracks to Showers  │
└─────────┬──────────┘
          ▼
┌────────────────────┐
│ Final Direction    │
│ Determination      │
└────────────────────┘
```

## Example Scenarios

### Scenario 1: Simple Track to Shower Conversion
```cpp
if (sg->get_particle_type()==0 || fabs(sg->get_particle_type())==13 
    && sg->get_length() < 40*units::cm && sg->is_dir_weak()) {
    sg->set_particle_type(11);  // Convert to electron
    sg->set_particle_mass(mp.get_mass_electron());
    sg->cal_4mom();
}
```
This shows how an ambiguous or muon-like track might be converted to an electron shower based on length and direction strength.

### Scenario 2: Multiple Track Handling
```cpp
if (map_vertex_segments[vtx].size() > 1) {
    bool flag_non_ele = false;
    for (auto it2 = map_vertex_segments[vtx].begin(); 
         it2 != map_vertex_segments[vtx].end(); it2++) {
        if (sg2 == sg1) continue;
        if (!sg2->get_flag_shower()) flag_non_ele = true;
    }
    if ((!flag_non_ele) && map_vertex_segments[vtx].size()<=3) 
        flag_good_track = true;
}
```
This demonstrates how the function handles vertices with multiple segments, checking if they should be classified as tracks or showers.

## Key Features

1. **Iterative Refinement**
   - Makes multiple passes to improve classification
   - Each pass can affect subsequent decisions

2. **Topological Analysis**
   - Uses vertex connections
   - Considers track/shower patterns
   - Analyzes energy deposition

3. **Direction Determination**
   - Based on vertex connections
   - Uses energy gradients
   - Considers overall event topology

4. **Particle Identification**
   - Distinguishes electrons from muons
   - Identifies electromagnetic showers
   - Handles ambiguous cases

## Function Call Hierarchy

The function makes several key function calls in sequence. Here's the complete list of functions it invokes:

1. `examine_good_tracks(temp_cluster->get_cluster_id())`  [more details](./examine_good_tracks.md)
   - Evaluates tracks in the cluster to identify high-quality tracks
   - First function called in the sequence

2. `fix_maps_multiple_tracks_in(temp_cluster->get_cluster_id())` [more details](./fix_maps.md)
   - Makes tracks undetermined if multiple tracks enter vertex
   - Called after initial track examination

3. `fix_maps_shower_in_track_out(temp_cluster->get_cluster_id())` [more details](./fix_maps.md)
   - Called twice in the sequence
   - Reverses shower direction if needed based on track patterns

4. `improve_maps_one_in(temp_cluster)` [more details](./improve_maps.md)
   - Handles cases where one track/shower enters vertex
   - Helps determine direction of other particles

5. `improve_maps_shower_in_track_out(temp_cluster->get_cluster_id())` [more details](./improve_maps.md)
   - Called twice with different parameters
   - Uses shower information to determine remaining particle classifications

6. `improve_maps_no_dir_tracks(temp_cluster->get_cluster_id())` [more details](./improve_maps.md)
   - Helps change tracks around showers to showers
   - Assists with directionality determination

7. `improve_maps_multiple_tracks_in(temp_cluster->get_cluster_id())` [more details](./improve_maps.md)
   - Changes track classification if multiple tracks enter
   - Called after shower/track relationship analysis

8. `judge_no_dir_tracks_close_to_showers(temp_cluster->get_cluster_id())` [more details](./judge_no_dir_tracks_close_to_showers.md)
   - Makes final judgments about tracks near showers
   - Called near end of sequence

9. `examine_maps(temp_cluster)` [more details](./examine_maps.md)
   - Final examination of all mappings
   - Verifies consistency of classifications

10. `examine_all_showers(temp_cluster)` [more details](./examine_all_showers.md)
    - Final check of all identified showers
    - Last major function called

Call sequence diagram:
```
examine_good_tracks()
    ↓
fix_maps_multiple_tracks_in()
    ↓
fix_maps_shower_in_track_out()
    ↓
improve_maps_one_in()
    ↓
improve_maps_shower_in_track_out() [First call]
    ↓
improve_maps_no_dir_tracks()
    ↓
improve_maps_shower_in_track_out() [Second call]
    ↓
improve_maps_multiple_tracks_in()
    ↓
fix_maps_shower_in_track_out()
    ↓
judge_no_dir_tracks_close_to_showers()
    ↓
examine_maps()
    ↓
examine_all_showers()
```

## Conclusion

The function employs a sophisticated multi-step approach to determine particle types and directions in neutrino interactions. It's designed to handle complex topologies and make decisions based on multiple criteria, iteratively improving its classifications through several passes of analysis. The sequence of function calls shows how it progressively refines its analysis, starting with basic track identification and moving through increasingly sophisticated shower and track determinations.