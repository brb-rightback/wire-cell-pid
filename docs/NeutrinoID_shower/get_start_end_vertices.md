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

