# ProtoSegment `get_direct_length()` Functions Documentation

## Overview

The `get_direct_length()` functions in the ProtoSegment class calculate straight-line distances between points in different contexts. There are three overloads of this function, each serving a specific purpose in analyzing the geometry of particle tracks.

## Function Variants

### 1. Basic Direct Length
```cpp
double get_direct_length()
```

#### Purpose
Calculates the straight-line distance between the first and last points of the segment.

#### Implementation
```cpp
double length = sqrt(pow(fit_pt_vec.front().x - fit_pt_vec.back().x, 2) + 
                    pow(fit_pt_vec.front().y - fit_pt_vec.back().y, 2) + 
                    pow(fit_pt_vec.front().z - fit_pt_vec.back().z, 2)); 
```

#### Usage
- Used to get a quick measure of the overall segment length
- Useful for comparing with the integrated path length to assess track straightness
- Common in particle identification and track characterization

### 2. Direct Length Between Points
```cpp
double get_direct_length(int n1, int n2)
```

#### Purpose
Calculates the straight-line distance between two arbitrary points in the segment specified by their indices.

#### Parameters
- `n1`: Starting point index
- `n2`: Ending point index

#### Implementation Details
- Performs bounds checking on indices
- Uses the fitted points vector (`fit_pt_vec`)
- Returns Euclidean distance between the specified points

#### Code
```cpp
if (n1 < 0) n1 = 0;  
if (n1+1 > fit_pt_vec.size()) n1 = int(fit_pt_vec.size())-1;
if (n2 < 0) n2 = 0;  
if (n2+1 > fit_pt_vec.size()) n2 = int(fit_pt_vec.size())-1;

double length = sqrt(pow(fit_pt_vec.at(n1).x - fit_pt_vec.at(n2).x, 2) +
                    pow(fit_pt_vec.at(n1).y - fit_pt_vec.at(n2).y, 2) +
                    pow(fit_pt_vec.at(n1).z - fit_pt_vec.at(n2).z, 2));
```

### 3. Direct Length With Direction Constraint
```cpp
double get_direct_length(int n1, int n2, TVector3 dir_perp)
```

#### Purpose
Calculates the direct length between two points, considering only the components perpendicular to a specified direction.

#### Parameters
- `n1`: Starting point index
- `n2`: Ending point index
- `dir_perp`: Direction vector to project against

#### Implementation Details
- Performs bounds checking on indices
- Projects the displacement vector onto the plane perpendicular to `dir_perp`
- Uses vector algebra to remove the parallel component

#### Code
```cpp
if (n1 < 0) n1 = 0;  
if (n1+1 > fit_pt_vec.size()) n1 = int(fit_pt_vec.size())-1;
if (n2 < 0) n2 = 0;  
if (n2+1 > fit_pt_vec.size()) n2 = int(fit_pt_vec.size())-1;

TVector3 temp_dir(fit_pt_vec.at(n1).x - fit_pt_vec.at(n2).x,
                  fit_pt_vec.at(n1).y - fit_pt_vec.at(n2).y,
                  fit_pt_vec.at(n1).z - fit_pt_vec.at(n2).z);
double length = sqrt(pow(temp_dir.Mag(), 2) - 
                    pow(temp_dir.Dot(dir_perp.Unit()), 2));
```

## Common Applications

### Track Analysis
- Used to calculate track straightness by comparing direct length with integrated path length
- Helps identify kinks or curves in particle tracks
- Important for particle identification (PID) algorithms

### Shower Analysis
- Used with direction constraints to analyze shower spread
- Helps characterize shower development perpendicular to the main axis
- Useful in electromagnetic shower identification

### Track/Shower Discrimination
The various `get_direct_length()` functions are often used together for:
- Calculating track straightness metrics
- Analyzing shower development patterns
- Determining particle interaction types
- Supporting particle identification algorithms

## Usage Examples

1. Basic Track Straightness:
```cpp
double direct_length = segment.get_direct_length();
double path_length = segment.get_length();
double straightness = direct_length / path_length;
```

2. Sectional Analysis:
```cpp
// Analyze middle third of track
int n_points = segment.get_point_vec().size();
double mid_straightness = segment.get_direct_length(n_points/3, 2*n_points/3);
```

3. Shower Analysis:
```cpp
TVector3 shower_axis = calculate_shower_axis();
double transverse_spread = segment.get_direct_length(0, -1, shower_axis);
```

## Notes and Caveats

1. **Index Handling**
   - All functions include bounds checking
   - Invalid indices are automatically clamped to valid ranges
   - Negative indices are set to 0
   - Indices beyond vector size are set to last valid index

2. **Numerical Precision**
   - Functions use standard floating-point arithmetic
   - Results are in the same units as input coordinates
   - No specific handling of very small distances

3. **Performance Considerations**
   - Calculations are relatively lightweight
   - No dynamic memory allocation
   - Suitable for frequent calls in analysis loops