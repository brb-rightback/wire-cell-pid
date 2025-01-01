# NeutrinoID Function Documentation

## examine_merge_showers()

This function handles the merging of showers that originate from the main neutrino vertex. 

### Key Steps:

1. **Identify Main Vertex Showers**
   - Iterates through showers connected to the main vertex
   - Only processes showers with start_vertex.second == 1 (direct connection)

2. **Direction Analysis**
   - For each shower, calculates directional vectors
   - Compares angles between potential shower pairs

3. **Merging Logic**
```cpp
// Angle threshold check for merging
if (dir1.Angle(dir2)/3.1415926*180. < 10) {
    // Add second shower to merge candidates
    map_shower_merge_showers[shower1].insert(shower2);
    used_showers.insert(shower2);
}
```

4. **Shower Combination**
   - Merges showers that share similar directions
   - Updates particle types and kinematics
   - Recalculates charge

### Flow Diagram:
```
Start
  ↓
Check Main Vertex Showers
  ↓
Calculate Directions
  ↓
Compare Angles (<10°)
  ↓
Group Similar Showers
  ↓
Merge & Update Properties
  ↓
End
```

