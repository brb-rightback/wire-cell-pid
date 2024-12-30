# Energy and Momentum Analysis in ProtoSegment

## Overview
The ProtoSegment class implements several methods for calculating kinetic energy and momentum of particle trajectories. The main functions can be grouped into three categories:

1. Kinetic Energy Calculation
   - cal_kine_range()
   - cal_kine_dQdx()
2. Direction Vector Calculation
   - cal_dir_3vector()
   - cal_dir_3vector(Point& p, double dis_cut)
3. Four-Momentum Assembly
   - cal_4mom()

## Kinetic Energy Calculation

### Strategy Selection Logic

```mermaid
graph TD
    A[Calculate Kinetic Energy] --> B{Track Length}
    B -->|< 4cm| C[Use dQ/dx Method]
    B -->|≥ 4cm| D{Is Shower Trajectory?}
    D -->|Yes| C
    D -->|No| E[Use Range Method]
```

### 1. Range-Based Method (cal_kine_range)

The range-based method calculates kinetic energy based on the particle's total path length using particle-specific range-to-energy conversion graphs.

Key points:
- Uses predefined range-to-energy curves from TPCParams
- Different curves for different particle types (e-, μ±, π±, K±, p)
- More reliable for longer tracks (≥ 4cm)
- Formula: KE = graph->Eval(length/units::cm) * units::MeV

### 2. dQ/dx-Based Method (cal_kine_dQdx)

This method calculates kinetic energy by integrating the energy loss along the track using charge deposition measurements.

```mermaid
graph LR
    A[Input dQ/dx] --> B[Convert to dE/dx]
    B --> C[Integrate over track length]
    C --> D[Total Kinetic Energy]
    
    subgraph Conversion
    B --"Modified Box Model"--> E["dE/dx = (exp(dQ*β) - α)/(β/1.38/0.273)"]
    end
```

Key features:
- Uses modified box model for dQ/dx to dE/dx conversion
- Parameters: α = 1.0, β = 0.255
- Handles edge effects at track endpoints
- More reliable for short tracks (< 4cm)

## Direction Vector Calculation

### 1. Basic Direction (cal_dir_3vector)

```mermaid
graph TD
    A[Determine Direction] --> B{Check flag_dir}
    B -->|+1| C[Forward Direction]
    B -->|-1| D[Backward Direction]
    C --> E[Average first 5 points]
    D --> F[Average last 5 points]
    E --> G[Normalize Vector]
    F --> G
```

Features:
- Uses flag_dir to determine forward/backward direction
- Averages multiple points for stability
- Returns unit vector

### 2. Local Direction (cal_dir_3vector with point)

Calculates direction vector around a specific point within a distance cut:
- Finds all points within dis_cut radius
- Computes average position relative to reference point
- Returns normalized direction vector

## Four-Momentum Assembly (cal_4mom)

The cal_4mom function combines kinetic energy and direction calculations to construct the full 4-momentum vector.

```mermaid
graph TD
    A[Calculate 4-Momentum] --> B[Determine KE Method]
    B --> C[Calculate Total Energy]
    C --> D[Calculate Momentum Magnitude]
    D --> E[Get Direction Vector]
    E --> F[Construct 4-Vector]
    
    subgraph Energy Choice
    B --> G{Length < 4cm?}
    G -->|Yes| H[Use dQ/dx]
    G -->|No| I{Is Shower?}
    I -->|Yes| H
    I -->|No| J[Use Range]
    end
```

Process:
1. Calculates kinetic energy (KE) using appropriate method
2. Total energy = KE + rest mass
3. Momentum magnitude = sqrt(E² - m²)
4. Direction from cal_dir_3vector()
5. Constructs components:
   - E = KE + mass
   - p_x = p * dir.X()
   - p_y = p * dir.Y()
   - p_z = p * dir.Z()

## Usage Considerations

1. Energy Method Selection:
   - Short tracks (< 4cm): Always use dQ/dx method
   - Shower-like segments: Use dQ/dx method
   - Long tracks: Use range method unless shower-like

2. Direction Determination:
   - Relies on flag_dir from pattern recognition
   - Forward (+1): Start → End
   - Backward (-1): End → Start
   - Weak directions stored in dir_weak flag

3. Particle Type Dependencies:
   - Different range curves for different particles
   - Mass values affect momentum calculation
   - Supported types: e±, μ±, π±, K±, p

4. Calibration Constants:
   - Box model parameters (α, β)
   - Energy units conversion
   - Charge collection calibration