# NeutrinoID BDT Functions Documentation

## Main BDT Functions

### cal_bdts_xgboost()
- Primary function for calculating XGBoost-based BDT scores
- Uses a large number of variables (>200) for classification
- Uses XGB_nue_seed2_0923.xml model
- Returns a log-transformed score: log10((1+val)/(1-val))

### cal_bdts()
- Alternative BDT calculation function
- Calculates multiple component BDT scores
- Uses BDTcombine_BDT800_3.weights.xml model
- Returns final combined score

## Component BDT Functions

### MIP (Minimum Ionizing Particle) Related
- `cal_mipid_bdt()`: Calculates MIP identification score
- `cal_mipquality_bdt()`: Evaluates quality of MIP signals
- Uses variables like:
  - Energy
  - dQ/dx measurements
  - Track lengths
  - Angular measurements

### Topology Related

#### Track-Shower Separation
- `cal_tro_1_bdt()`: Track/shower characteristics evaluation
- `cal_tro_2_bdt()`: Secondary track/shower analysis
- `cal_tro_4_bdt()`: Direction and energy-based track/shower analysis
- `cal_tro_5_bdt()`: Vertex-based track/shower classification

#### Shower Analysis
- `cal_stw_spt_bdt()`: Shower pattern and trajectory analysis
- `cal_stw_2_bdt()`: Secondary shower characteristics
- `cal_stw_3_bdt()`: Tertiary shower analysis
- `cal_stw_4_bdt()`: Additional shower pattern recognition

### Pattern Recognition

#### Gap Analysis
- `cal_gap_bdt()`: Analyzes gaps in particle tracks
- Considers:
  - Track prolongation
  - Parallel tracks
  - Energy measurements

#### Vertex Analysis
- `cal_vis_1_bdt()`: Primary vertex analysis
- `cal_vis_2_bdt()`: Secondary vertex characteristics
- Evaluates:
  - Number of segments at vertex
  - Angular distributions
  - Energy measurements

### Physics Process Identification

#### Pi0 Related
- `cal_pio_1_bdt()`: Primary π⁰ identification
- `cal_pio_2_bdt()`: Secondary π⁰ characteristics
- Analyzes:
  - Mass reconstruction
  - Energy measurements
  - Angular distributions

#### Signal Related
- `cal_sig_1_bdt()`: Primary signal characteristics
- `cal_sig_2_bdt()`: Secondary signal patterns
- Focuses on:
  - Energy distributions
  - Shower angles
  - Pattern recognition

### Low Energy Effects
- `cal_lol_1_bdt()`: Primary low energy analysis
- `cal_lol_2_bdt()`: Secondary low energy patterns
- Analyzes:
  - Energy thresholds
  - Track segments
  - Angular distributions

### Bad Reconstruction Identification
- `cal_br1_bdt()`: Primary bad reconstruction patterns
- `cal_br3_bdt()`: Tertiary reconstruction analysis
- `cal_br3_3_bdt()`: Specific reconstruction pattern analysis
- `cal_br3_5_bdt()`: Complex reconstruction pattern analysis
- `cal_br3_6_bdt()`: Additional reconstruction quality checks

### Miscellaneous Analysis
- `cal_hol_lol_bdt()`: High/low overlap analysis
- `cal_cme_anc_bdt()`: Cosmic muon energy analysis
- `cal_mgo_mgt_bdt()`: Multiple gamma analysis
- `cal_stemdir_br2_bdt()`: Stem direction analysis
- `cal_trimuon_bdt()`: Three-muon event analysis
- `cal_br4_tro_bdt()`: Combined reconstruction and track analysis

## Implementation Details

All BDT functions:
1. Initialize a TMVA::Reader
2. Add relevant variables
3. Load corresponding XML weight files
4. Calculate and return scores
5. Include default values for error cases

Most functions include checks for filled variables and return default values if required variables are not available.

## Common Features

### Input Processing
- Most functions handle vector inputs
- Process multiple measurements per event
- Return worst-case (minimum) scores for vector inputs

### Error Handling
- Default values provided for missing data
- Boundary checks for vector operations
- Protection against undefined states

### Performance Considerations
- Optimized for speed with early returns
- Efficient handling of vector operations
- Minimal memory allocation