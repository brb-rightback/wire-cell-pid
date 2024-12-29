# GLM Tagger Function Explanation

The `glm_tagger()` function is a cosmic ray tagging algorithm that analyzes particle tracks in a liquid argon time projection chamber (TPC) detector. Its main purpose is to identify and classify different types of cosmic ray tracks.

## Overview

The function takes inputs like:
- Flash/light information from photon detectors
- 3D reconstructed particle tracks/clusters
- Detector geometry information
- Timing information

And outputs:
- Tag type (0-3 indicating track classification)
- Main cluster pointer
- Associated flash pointer

## Key Track Classifications

1. **Normal Track (0)**: Tracks that don't meet special criteria
2. **TGM with Flash (1)**: Through-going muon with associated light flash
3. **TGM without Flash (2)**: Through-going muon without associated light flash  
4. **STM (3)**: Stopping muon candidate

## Main Algorithm Steps

1. Initial Checks:
```cpp
// Check if track is long enough to analyze
Point p0(extreme_points[0][0].x, extreme_points[0][0].y, extreme_points[0][0].z);
Point p1(extreme_points[1][0].x, extreme_points[1][0].y, extreme_points[1][0].z);
double distance = sqrt(pow(p0.x-p1.x,2) + pow(p0.y-p1.y,2) + pow(p0.z-p1.z,2));

if(distance <= 10*units::cm) {
    return std::make_tuple(0, main_cluster, main_flash); // Too short to analyze
}
```

2. Flash Matching Analysis:
```cpp
// Analyze predicted vs measured light patterns
std::vector<double> pred_pmt_light = calculate_pred_pe(
    eventTime, run_no, offset_x, pl, 
    main_cluster, additional_clusters, 
    flash, flag_match_data, flag_timestamp
);

// Compare predicted and measured light patterns
double chi2_user = 0;
int ndf_user = 0;
for(int i=0; i<pred_pmt_light.size(); i++) {
    double flash_pe = flash->get_PE(i);
    double pred_pe = pred_pmt_light[i];
    double pe_err = flash->get_PE_err(i);
    
    chi2_user += pow(pred_pe - flash_pe,2)/pow(pe_err,2);
    if(pred_pe!=0 || flash_pe!=0) {
        ndf_user++;
    }
}
```

3. Boundary Checks:
```cpp
// Check if track endpoints are near detector boundaries
bool inside_boundary = inside_fiducial_volume(point, offset_x, tolerance_vec);

// Check specific regions like anode
bool near_anode = false;
if(point.x - offset_x < dist2anode_cut) {
    near_anode = true;
}
```

4. Track Classification:
```cpp
// Example of TGM classification
if(flag_small_SCB_p0==false && flag_small_SCB_p1==false && 
   flag_large_SCB_p0==true && flag_large_SCB_p1==true) {
    
    if(flag_NotNear_anode && chi2_user < bundle_chi2) {
        bool flag_tgm = check_tgm(main_cluster, flash, offset_x, ct_point_cloud);
        if(flag_tgm) {
            flag_M2 = 1; // TGM with flash
            new_flash = flash;
        }
    }
}
```

5. Final Output:
```cpp
return std::make_tuple(flag_M2, main_cluster, new_flash);
```

## Key Tolerances and Parameters

The algorithm uses several key parameters to tune its sensitivity:

```cpp
// Example tolerance parameters
double chi2_stm_tol = 1;        // Chi2 tolerance for STM
double chi2_tgm_tol = 3;        // Chi2 tolerance for TGM
double bundle_ks_stm_tol = 0.04;  // KS test tolerance for STM
double bundle_ks_tgm_tol = 0.03;  // KS test tolerance for TGM
double stm_flash_tol = 1.2*units::m;  // Flash position tolerance for STM
double tgm_flash_tol = 1.6*units::m;  // Flash position tolerance for TGM
```

## Mermaid Flow Diagram

```mermaid
graph TD
    A[Start] --> B{Track > 10cm?}
    B -->|No| C[Return Normal 0]
    B -->|Yes| D[Calculate Light Pattern]
    D --> E{Check Original KS}
    E -->|Bad Match| F[Check TGM]
    E -->|Good Match| G[Check STM]
    F --> H{Near Boundary?}
    H -->|Yes| I{Check Light Match}
    I -->|Good| J[Tag as TGM w/Flash 1]
    I -->|Bad| K[Tag as TGM wo/Flash 2]
    H -->|No| L[Continue Checks]
    G --> M{At Cathode?}
    M -->|Yes| N{Check Light Match}
    N -->|Good| O[Tag as STM 3]
    N -->|Bad| P[Continue Normal]
    M -->|No| P
    P --> Q[Return Normal 0]
```

## Important Notes

1. The function uses multiple criteria to make classifications:
   - Track geometry
   - Light pattern matching
   - Detector boundary proximity
   - Track behavior (stopping vs through-going)

2. The tolerances and parameters can be tuned based on detector conditions

3. The algorithm prioritizes classifications in this order:
   - TGM with flash
   - TGM without flash 
   - STM
   - Normal tracks

4. Flash matching quality is assessed using both:
   - Chi-square comparison
   - Kolmogorov-Smirnov test

5. Multiple boundary checks are performed to ensure track classifications are robust