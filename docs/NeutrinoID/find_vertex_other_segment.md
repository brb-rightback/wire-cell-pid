# Analysis of Vertex Finding Functions in NeutrinoID

## Function: `find_vertex_other_segment`

### Purpose
This function attempts to find or create a vertex for a given endpoint of a track/segment by examining geometrical relationships with existing vertices and segments.

### Key Parameters
- `temp_cluster`: The PR3DCluster being analyzed
- `flag_forward`: Boolean indicating whether to check front or back of tracking path
- `wcp`: WCPoint reference point to start searching from

### Logical Flow
1. First attempts to find vertex through `check_end_point` with default parameters
2. If unsuccessful, makes two more attempts with progressively relaxed parameters:
   - Second try: vtx_cut1=1.2cm, vtx_cut2=2.5cm
   - Third try: vtx_cut1=1.5cm, vtx_cut2=3.0cm

3. After finding candidate vertices/segments, follows one of three paths:

[check_end_point alg](./check_end_point.md)

```cpp
if (existing_vertex) {
    return existing_vertex;
} 
else if (found_segment) {
    // Create new vertex from break point
    // Update segment connectivity
    // Return new vertex
}
else {
    // Create brand new vertex at wcp
    return new_vertex;
}
```

### Function Call Tree
```
find_vertex_other_segment
├── check_end_point
│   ├── get_closest_point
│   └── proto_extend_point
├── get_closest_wcpoint
└── add_proto_connection
```

