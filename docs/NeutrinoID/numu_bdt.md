# NeutrinoID NuMu BDT Functions Documentation

## Overview
This file contains functions for calculating various Boosted Decision Tree (BDT) scores used in neutrino muon (NuMu) identification.

## Main Function

### `cal_numu_bdts_xgboost()`
Primary function that calculates the overall NuMu BDT score using XGBoost model. It:
- Calculates individual BDT scores for various sub-classifiers
- Uses TMVA reader with multiple variables related to cosmic tracks, particle types, and geometrical properties
- Returns a transformed score using log10((1+val)/(1-val))

## Sub-Classifier BDT Functions

### `cal_numu_bdts()`
Legacy BDT calculation function that:
- Calculates scores for multiple sub-classifiers
- Uses TMVA reader with different set of variables
- Returns raw BDT score without transformation

### `cal_cosmict_2_4_bdt(float default_val)`
Calculates cosmic tagger BDT score for type 2 and 4 events:
- Uses variables like particle type, muon tracks, shower length
- Only calculates if cosmict_2 is filled
- Returns default value if calculation cannot be performed

### `cal_cosmict_3_5_bdt(float default_val)` 
Calculates cosmic tagger BDT score for type 3 and 5 events:
- Uses variables related to track properties inside detector
- Only calculates if cosmict_3 is filled
- Returns default value if calculation cannot be performed

### `cal_cosmict_6_bdt(float default_val)`
Calculates cosmic tagger BDT score for type 6 events:
- Uses minimal set of variables including direction and angle
- Only calculates if cosmict_6 is filled
- Returns default value if calculation cannot be performed

### `cal_cosmict_7_bdt(float default_val)`
Calculates cosmic tagger BDT score for type 7 events:
- Uses extensive set of variables including shower properties
- Only calculates if cosmict_7 is filled
- Returns default value if calculation cannot be performed

### `cal_cosmict_8_bdt(float default_val)`
Calculates cosmic tagger BDT score for type 8 events:
- Uses basic variables like muon length and accumulative length
- Only calculates if cosmict_8 is filled
- Returns default value if calculation cannot be performed

### `cal_cosmict_10_bdt(float default_val)`
Calculates cosmic tagger BDT score for type 10 events:
- Processes vector of events and returns minimum score
- Uses variables related to vertex position and shower properties
- Returns default value if no events to process

### `cal_numu_1_bdt(float default_val)`
Calculates NuMu classifier type 1 BDT score:
- Processes vector of events and returns maximum score
- Uses variables related to particle properties and track characteristics
- Returns default value if no events to process

### `cal_numu_2_bdt(float default_val)`
Calculates NuMu classifier type 2 BDT score:
- Processes vector of events and returns maximum score
- Uses track length and daughter particle information
- Returns default value if no events to process

### `cal_numu_3_bdt(float default_val)`
Calculates NuMu classifier type 3 BDT score:
- Uses single event variables including particle type and track lengths
- Returns raw BDT score
- Returns default value if calculation cannot be performed

## Input Files
All BDT functions rely on pre-trained models stored in XML files:
- Main XGBoost model: "numu_scalars_scores_0923.xml"
- Cosmic taggers: "cos_tagger_X.weights.xml" (where X is tagger type)
- NuMu classifiers: "numu_taggerX.weights.xml" (where X is classifier type)

## Important Notes
1. Functions handle invalid values (NaN) by setting them to 0
2. Many functions operate on vectors of events, choosing maximum or minimum scores
3. All functions accept a default value parameter returned when calculation cannot be performed
4. Score transformations are applied in some cases to normalize outputs