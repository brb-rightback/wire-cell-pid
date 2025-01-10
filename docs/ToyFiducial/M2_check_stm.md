# Comparison of `check_stm` and `M2_check_stm`

This document compares two similar functions in the ToyFiducial namespace:
- `WCPPID::ToyFiducial::check_stm`
- `WCPPID::ToyFiducial::M2_check_stm`

## Overview

Both functions serve to identify Single-Track Michel (STM) events but have slight differences in their classification approach. The main workflow is similar, with key differences in the final classification steps.

## Common Algorithm Steps

1. **Initial Validation**
   - Check for empty Steiner point cloud
   - Compute and validate boundary points
   - Identify candidate exit points

2. **Path Analysis**
   - Construct rough path between endpoints
   - Collect charge along trajectory
   - Perform fine tracking
   - Identify kinks in dQ/dx distribution
   - Calculate path segments

3. **Cluster Analysis**
   - Check for additional clusters
   - Search for other tracks
   - Validate against fiducial volume

## Key Differences

1. **Event Classification**
   - `check_stm`: Uses event_type bitmask for classification
   - `M2_check_stm`: No event_type parameter, simpler return logic

2. **Return Behavior**
   - `check_stm`:
     ```cpp
     event_type |= 1UL << 3;
     return true;
     ```
   - `M2_check_stm`:
     ```cpp
     return false;  // For similar conditions
     ```

## Common Return Conditions

Both functions return `false` when:
- Steiner point cloud is empty
- Track is fully contained in fiducial volume
- Other specific validation checks fail

## Conclusion

The primary difference lies in the classification mechanism:
- `check_stm` uses event type flags for detailed classification
- `M2_check_stm` employs simplified boolean returns

The core algorithm logic remains largely identical between both versions.