## calculate_shower_kinematics()

This function handles the calculation of shower kinematic properties, treating different particle types appropriately.

### Main Operations:

1. **Regular Showers (non-muon)**
   - Calculates basic kinematics
   - Computes charge deposition
   - Sets kinematics flag

2. **Long Muons**
   - Special handling for muon-like tracks
   - Uses segments_in_long_muon for calculations
   - Different charge calculation approach

### Implementation Pattern:
```cpp
for (auto shower : showers) {
    if (!shower->get_flag_kinematics()) {
        if (abs(shower->get_particle_type()) != 13) {
            // Regular shower calculations
            shower->calculate_kinematics();
            double charge = cal_kine_charge(shower);
        } else {
            // Long muon calculations
            shower->calculate_kinematics_long_muon(segments_in_long_muon);
            double charge = cal_kine_charge(shower);
        }
        shower->set_flag_kinematics(true);
    }
}
```

## Integration of Functions

These functions work together in the following typical sequence:

1. Initial shower identification/creation
2. `calculate_shower_kinematics()` computes basic properties
3. `examine_merge_showers()` looks for merger opportunities [see details](./examine_merge_showers.md)
4. `update_shower_maps()` maintains consistency [see details](./update_shower_maps.md)
5. Process repeats as needed

This integrated approach ensures proper shower reconstruction and characterization while maintaining data consistency throughout the process.