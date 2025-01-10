# ToyFiducial Class Overview

## Purpose
The ToyFiducial class is part of the WCPPID namespace and provides functionality for fiducial volume checks and cosmic ray tagging in a wire chamber particle detector. It helps identify and classify different types of particle tracks, particularly through-going muons (TGM) and stopping muons (STM).

## Key Features

### Initialization & Configuration
- Configurable detector boundaries and geometry
- Support for both Monte Carlo and real data configurations
- Dead region handling and channel mapping
- Configurable space charge boundary definitions

### Core Functionality

#### Boundary Checking
- `inside_fiducial_volume()`: Checks if a point is inside the detector's fiducial volume [detailed documentation](./ToyFiducial/inside_fiducial_volume.md)
- `inside_dead_region()`: Determines if a point is in a detector dead region [detailed documentation](./ToyFiducial/inside_dead_region.md)
- `check_dead_volume()`: Examines dead volume along a track direction [detailed documentation](./ToyFiducial/check_dead_region.md)
- `check_signal_processing()`: Verifies signal processing quality along tracks [detailed documentation](./ToyFiducial/check_signal_processing.md)
- `inside1_outside0_SCB()`: inside the space charge boundary [detailed documentation](./ToyFiducial/inside1_outside0_SCB.md)
- `check_boundary()`: check boundary [detailed documentation](./ToyFiducial/check_boundary.md)
- `inside_x_region()`: check boundary [detailed documentation](./ToyFiducial/inside_x_region.md)

#### Track Classification

##### Through-Going Muon (TGM) Detection
- `check_tgm()`: Main TGM identification function [detailed documentation](./ToyFiducial/check_tgm.md)
- `M2_check_tgm()`: Enhanced TGM checking algorithm [detailed documentation](./ToyFiducial/M2_check_tgm.md)
- Handles both flash-matched and flash-independent TGM identification

##### Stopping Muon (STM) Detection
- `check_stm()`: Primary STM identification function [detailed documentation](./ToyFiducial/check_stm.md)
- `check_stm_only()`: STM identification function [detailed documentation](./ToyFiducial/check_stm_only.md)
- `M2_check_stm()`: Enhanced STM checking algorithm [detailed documentation](./ToyFiducial/M2_check_stm.md)
- `eval_stm()`: Evaluates track characteristics for STM classification [detailed documentation](./ToyFiducial/eval_stm.md)

#### Flash Matching
- `calculate_pred_pe()`: Calculates predicted photoelectron response [detailed documentation](./ToyFiducial/calculate_pred_pe.md)
- Provides flash-to-track matching capabilities
- Handles PMT response predictions

### Advanced Features

### Add Dead Region
- `AddDeadRegion`: add dead region [detailed documentation](./ToyFiducial/AddDeadRegion.md)

#### Neutrino candidate
- `check_neutrino_candidate()`: identify neutrino candidate [detailed documentation](./ToyFiducial/check_neutrino_candidate.md)
- `check_fully_contained():` check if the activity is fully contained [detailed documentation](./ToyFiducial/check_fully_contained.md)

#### Cosmic Ray Tagging
- `glm_tagger()`: Global cosmic ray tagging algorithm [detailed documentation](./ToyFiducial/glm_tagger.md)
- `M2_cosmic_tagger()`: Enhanced cosmic ray identification [detailed documentation](./ToyFiducial/M2_cosmic_tagger.md)
- Combines multiple detection strategies:
    - Flash matching
    - Track topology analysis
    - Boundary intersection checks

### Full Detector dead
- `check_full_detector_dead()`: check if the full detector is dead [detailed documentation](./ToyFiducial/check_full_detector_dead.md)

#### Track Analysis
- `find_first_kink()`: Identifies track bending points [detailed documentation](./ToyFiducial/find_first_kink.md)
- `detect_proton()`: Proton track identification [detailed documentation](./ToyFiducial/detect_proton.md)
- Track trajectory and dE/dx analysis

## Key Data Structures

### Geometry Definition
```cpp
double m_top;        // Top boundary
double m_bottom;     // Bottom boundary
double m_upstream;   // Upstream boundary
double m_downstream; // Downstream boundary
double m_anode;      // Anode boundary
double m_cathode;    // Cathode boundary
```

### Space Charge Boundaries
```cpp
double m_sc_bottom_1_x, m_sc_bottom_1_y;
double m_sc_bottom_2_x, m_sc_bottom_2_y;
double m_sc_top_1_x, m_sc_top_1_y;
double m_sc_top_2_x, m_sc_top_2_y;
```

## Usage Examples

### Basic Fiducial Volume Check
```cpp
Point p(x, y, z);
bool is_inside = toy_fiducial.inside_fiducial_volume(p, offset_x);
```

### Cosmic Ray Tagging
```cpp
auto [tag_type, cluster, flash] = toy_fiducial.M2_cosmic_tagger(
    eventTime, flashes, main_cluster, additional_clusters, 
    main_flash, bundle_info, pl, time_offset, nrebin, 
    unit_dis, ct_point_cloud, run_no, subrun_no, event_no, 
    flag_data, global_wc_map
);
```

## Return Values

### Tag Types
- 0: Normal event (not cosmic)
- 1: TGM with flash
- 2: TGM without flash
- 3: STM candidate

## Dependencies
- Relies on WCP (Wire Cell Processing) framework
- Uses Boost Graph Library for track analysis
- Requires ROOT for histogram operations
- Uses Eigen for matrix operations

## Performance Considerations
- Computationally intensive track analysis
- Multiple boundary checking operations
- Complex flash matching calculations
- Memory-intensive point cloud operations

## Best Practices
- Initialize with appropriate boundary values
- Consider both MC and data configurations
- Handle flash matching carefully
- Account for dead regions and signal quality
- Consider multiple track topologies
- Validate results with multiple methods