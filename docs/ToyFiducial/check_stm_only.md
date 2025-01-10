# Comparison of `check_stm` and `check_stm_only`

## Overview

Both functions are designed to identify and classify tracks in the detector, but `check_stm_only` is a modified version of `check_stm` that specifically ignores Tagged Gamma (TGM) classifications.

## Key Differences

### 1. TGM Handling
- `check_stm`: Sets event type and returns true for TGM cases
  ```cpp
  if ((exit_L < 3*units::cm || left_L < 3*units::cm || flag_TGM_anode)) {
      event_type |= 1UL << 3;
      return true;
  }
  ```
- `check_stm_only`: Ignores TGM cases by commenting out the event type setting
  ```cpp
  if ((exit_L < 3*units::cm || left_L < 3*units::cm || flag_TGM_anode)) {
      // event_type |= 1UL << 3;
      // return true;
  }
  ```

### 2. Event Classification
- `check_stm`: Uses event_type bitmask for classification
- `check_stm_only`: Removes event type classification, focusing only on STM identification

### 3. Return Behavior
- `check_stm`: Returns true for both STM and TGM cases
- `check_stm_only`: Only returns true for STM cases, ignoring TGM scenarios

## Common Features

Both functions share:
1. Initial cluster validation
2. Boundary point checking
3. Path analysis and tracking
4. Charge collection and analysis
5. Track evaluation logic

## Implementation Details

1. **Initial Checks**
   - Both perform identical Steiner point cloud validation
   - Both use the same boundary point detection

2. **Track Analysis**
   - Same path construction logic
   - Identical charge collection mechanisms
   - Same kink detection algorithms

3. **Classification Logic**
   - Both use similar track quality metrics
   - Share the same proton detection logic
   - Use identical fiducial volume checks

## Conclusion

The main distinction between these functions is their handling of TGM events:
- `check_stm` provides comprehensive classification including TGMs
- `check_stm_only` focuses solely on STM identification, ignoring TGM cases

This modification makes `check_stm_only` more specialized but potentially less complete in terms of overall event classification.
