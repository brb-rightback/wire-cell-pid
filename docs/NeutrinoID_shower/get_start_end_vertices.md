## get_start_end_vertices()

A utility function that returns the start and end vertices of a given segment.

### Function Operation:

1. **Input**: Takes a ProtoSegment pointer
2. **Processing**: 
   - Compares segment endpoint indices with vertex indices
   - Uses map_segment_vertices for lookup
3. **Output**: Returns pair of ProtoVertex pointers (start, end)

### Example Usage:
```cpp
auto vertices = get_start_end_vertices(segment);
ProtoVertex* start = vertices.first;
ProtoVertex* end = vertices.second;

// Useful for direction determination
if (start->get_wcpt().index == segment->get_wcpt_vec().front().index) {
    // Forward direction
} else {
    // Reverse direction
}
```

### Key Points:
- Essential for determining segment orientation
- Used in shower reconstruction and topology analysis
- Critical for maintaining connectivity information

## Integration of Functions

These functions work together in the following typical sequence:

1. Initial shower identification/creation
2. `calculate_shower_kinematics()` computes basic properties
3. `examine_merge_showers()` looks for merger opportunities
4. `update_shower_maps()` maintains consistency
5. Process repeats as needed

This integrated approach ensures proper shower reconstruction and characterization while maintaining data consistency throughout the process.