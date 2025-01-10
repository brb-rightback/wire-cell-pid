## update_shower_maps()

This function maintains the consistency of various mapping relationships between vertices, segments, and showers.

### Key Maps Updated:

1. **map_vertex_to_shower**
   - Maps vertices to their associated showers
   - Critical for vertex-based shower lookup

2. **map_vertex_in_shower** & **map_segment_in_shower**
   - Track which vertices/segments belong to which showers
   - Prevents double-counting/overlap

3. **used_shower_clusters**
   - Tracks cluster IDs that have been incorporated into showers

### Example Usage:
```cpp
// After modifying shower structure
shower->add_segment(new_segment, map_segment_vertices);
update_shower_maps();  // Ensures consistency

// Now safe to query
auto associated_shower = map_segment_in_shower[segment];
```

