flowchart TD
    A[Start] --> B[Clear existing cell_point_indices_map]
    B --> C[Get reference to point cloud]
    C --> D[Get next merged cell from mcells]
    D --> E{Is there a cell?}
    E -->|No| K[End]
    E -->|Yes| F[Get vector of point indices for this cell]
    F --> G{Does cell exist in map?}
    G -->|No| H[Create empty set for cell]
    G -->|Yes| I[Get existing set]
    H --> I
    I --> J[Add point indices to set]
    J --> D