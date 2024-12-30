# ProtoSegment Length Functions Documentation

## Overview

The ProtoSegment class implements several length calculation functions that serve different purposes in trajectory analysis. These functions operate on the `fit_pt_vec` member variable, which contains a vector of 3D points representing the fitted trajectory.

## Core Length Functions

### `get_length()`

```cpp
double get_length()
```

The base length calculation function that measures the total path length of the segment by:
- Iterating through adjacent points in `fit_pt_vec`
- Calculating the Euclidean distance between each pair of consecutive points
- Summing these distances to get the total length
- Returns the accumulated length in the native distance units

This represents the true path length following all points in the trajectory.

### `get_direct_length()`

```cpp
double get_direct_length()
```

Calculates the straight-line distance between the start and end points of the segment by:
- Taking only the first and last points in `fit_pt_vec`
- Computing the 3D Euclidean distance between them
- Returns this direct distance in native units

This is useful for comparing to the full path length to measure how straight or curved a segment is.

## Parameterized Length Functions 

### `get_length(int n1, int n2)`

```cpp
double get_length(int n1, int n2)
```

Similar to base `get_length()` but operates on a subsection of the trajectory:
- Takes start index `n1` and end index `n2` as parameters
- Performs bounds checking to ensure valid indices
- Calculates cumulative path length between specified points
- Returns the partial path length in native units

Useful for analyzing specific portions of a segment.

### `get_direct_length(int n1, int n2)`

```cpp
double get_direct_length(int n1, int n2)
```

Similar to base `get_direct_length()` but for a subsection:
- Takes start index `n1` and end index `n2` as parameters 
- Performs bounds checking
- Calculates straight-line distance between specified points
- Returns the direct distance in native units

### `get_length(int n1, int n2, TVector3 dir_perp)`

```cpp
double get_length(int n1, int n2, TVector3 dir_perp)
```

Specialized length calculation that:
- Takes start/end indices and a direction vector
- Projects the path onto plane perpendicular to direction vector
- Calculates cumulative length of projected path
- Returns projected length in native units

Used for measuring length perpendicular to a specific direction.

### `get_direct_length(int n1, int n2, TVector3 dir_perp)`

```cpp
double get_direct_length(int n1, int n2, TVector3 dir_perp)
```

Similar to above but:
- Projects start and end points onto perpendicular plane
- Returns direct distance between projected points
- Used for analyzing straightness perpendicular to direction

## Common Uses

These length functions serve several key purposes in trajectory analysis:

1. Measuring total path length vs direct distance to quantify curvature
2. Analyzing subsections of trajectories independently
3. Calculating projected lengths for specific analysis needs
4. Supporting PID and tracking algorithms that need length information

## Implementation Notes

- All functions operate on the fitted points in `fit_pt_vec`
- Bounds checking is performed to handle invalid indices
- Native units (typically cm) are preserved
- 3D Euclidean distances are used throughout
- Results can be used for track/shower discrimination
- Projected lengths help with specific physics analyses

## Example Usage

```cpp
ProtoSegment* seg = /* ... */;

// Get full path length
double total_length = seg->get_length();

// Get direct length between points 5 and 10
double partial_direct = seg->get_direct_length(5, 10);

// Get projected length perpendicular to beam direction
TVector3 beam_dir(0,0,1);
double proj_length = seg->get_length(0, -1, beam_dir);
```