# Climate variability thresholds govern transitions between vegetation growth patterns

## Description of the data and file structure

This repository contains the source code and processed datasets underlying the analysis presented in the manuscript titled above. The study investigates the global spatiotemporal dynamics of vegetation responses to climate variability. By integrating Vector Autoregression (VAR) models, sliding-window comparison, segmented regression, and Machine Learning (Random Forest + SHAP), the analysis identifies the climate variability thresholds associated with transitions between vegetation growth patterns.

## Files and variables

**File:** `Result1_Parallel_SlidingWindow_VAR_GPP_History.R`
**Description:** Fits Vector Autoregression (VAR) models in sliding windows for global grid cells and calculates Impulse Response Functions (IRF) to capture the dynamic legacy effects of climate on vegetation growth (GPP). (Language: R)

**File:** `Result1_VegetationResponseOutcomes.R`
**Description:** Classifies each grid cell into vegetation response outcomes from the between-period change in the IRF, and produces the threshold sensitivity panels and classification maps. (Language: R)

**File:** `Result2_Thershold_sensitivity.R`
**Description:** Fits segmented regression models of vegetation against climate per grid cell and climate factor to detect the breakpoints where climate variability alters vegetation growth variability. (Language: R)

**File:** `Result3_Prepare_Threshold_SHAP_data.py`
**Description:** Assembles the per-grid climate thresholds and response-type labels into a single table for the SHAP analysis. (Language: Python)

**File:** `Result3_RandomForest_SHAP_Thresholds.py`
**Description:** Trains a Random Forest classifier to predict the vegetation response type from the climate thresholds, and uses SHAP values to interpret feature importance and extract the physiological thresholds of specific climate variables. (Language: Python)

**File:** `Result4_Validation_of_legacy_transitions.R`
**Description:** Validates the response-type transitions against climate threshold-crossing events, and produces the global maps, frequency charts, and significance tests of state transitions. (Language: R)

## Code/software

The analysis and data processing were performed using free and open-source software:

* R (version 4.4.3 or higher recommended)
* Python (version 3.12 or higher recommended)

## Access information
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
All observational data that support the findings of this study are available as follows. The temperature, precipitation, and vapor pressure deficit from the CRU v4.0.7 data are available at https://crudata.uea.ac.uk/cru/data/hrg/. The short-wave radiation flux data from ECMWF Reanalysis 5 are available at https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview. The historical and SSP245 scenario atmospheric CO2 concentration data, derived from CMIP6, are available at https://aims2.llnl.gov/search/cmip6/. The root-zone and surface soil moisture data from GLEAM version 3.7 are available at https://www.gleam.eu/. The GLASS AVHRR GPP data are available at http://www.glass.umd.edu/GPP/AVHRR/. The GPP data simulated by the TL-LUE model are available at https://datadryad.org/stash/dataset/doi:10.5061/dryad.dfn2z352k. The AVHRR GIMMS NDVI3g V1.2 data are available at https://zenodo.org/records/8253971. The AVHRR GIMMS LAI4g V1.2 data are available at https://zenodo.org/records/8281930. The LCSIF product are available at https://zenodo.org/records/7916851 and https://zenodo.org/records/7916879. The VIIRS GLSP product VNP22C2 data are available at https://viirsland.gsfc.nasa.gov/Products/NASA/PhenologyESDR.html. The global vegetation map GLC-FCS30D is available at https://zenodo.org/records/8239305. The elevation data from WorldClim 2.1 is available at https://www.worldclim.org/data/worldclim21.html. The world map of Köppen-Geiger climate classification is available at http://koeppen-geiger.vu-wien.ac.at/present.htm. The SPEI global drought monitor is available at https://spei.csic.es/.  
