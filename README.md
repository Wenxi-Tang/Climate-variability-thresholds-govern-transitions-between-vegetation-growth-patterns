# Climate variability thresholds govern transitions between vegetation growth patterns

## Description of the data and file structure

This repository contains the source code and processed datasets underlying the analysis presented in the manuscript titled above. The study investigates the global spatiotemporal dynamics of vegetation responses to climate variability. By integrating Vector Autoregression (VAR) models, Sliding Window approaches, and Machine Learning (Random Forest + SHAP), the analysis identifies specific climate variability thresholds that trigger transitions between vegetation growth patterns.

### Files and variables

#### File: 01_Parallel_SlidingWindow_VAR_GPP.R

**Description:** Fits Vector Autoregression (VAR) models in sliding windows for global grid cells. Calculates Impulse Response Functions (IRF) to capture the dynamic legacy effects of climate on vegetation growth (GPP). (**Language:** R)

#### File: 02_Parallel_Response_Pattern_Classification.R

**Description:** Classifies vegetation response patterns based on the VAR and IRF results. (**Language:** R)

#### File: 03_Threshold_Analysis_Segmented_Regression.R

**Description:** Uses segmented regression models to detect breakpoints where the standard deviation of climate drivers significantly alters the standard deviation of vegetation growth. (**Language:** R)

#### File: 04_RandomForest_SHAP_Classification_Thresholds.py

**Description:** Performs the driver analysis using Machine Learning. (1) Trains a Random Forest classifier to predict the vegetation response type based on climate factors. (2) Uses SHAP values to interpret feature importance and extract physiological thresholds for specific climate variables. (**Language:** Python)

#### File: 05_Comprehensive_Visualization_Threshold_Crossing.R

**Description:** Visualizes the global spatial distribution of threshold crossing events, creates frequency bar charts of state transitions, and produces heatmaps of significance tests. (**Language:** R)

## Code/software

The analysis and data processing were performed using free and open-source software:

* **R** (version 4.4.3 or higher recommended)
* **Python** (version 3.12 or higher recommended)

## Access information
All observational data that support the findings of this study are available as follows. The temperature, precipitation, and vapor pressure deficit from the CRU v4.0.7 data are available at https://crudata.uea.ac.uk/cru/data/hrg/. The short-wave radiation flux data from ECMWF Reanalysis 5 are available at https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview. The historical and SSP245 scenario atmospheric CO₂ concentration data, derived from CMIP6, are available at https://aims2.llnl.gov/search/cmip6/. The root-zone and surface soil moisture data from GLEAM version 3.7 are available at https://www.gleam.eu/. The GLASS AVHRR GPP data are available at http://www.glass.umd.edu/GPP/AVHRR/. The GPP data simulated by the TL-LUE model are available at https://datadryad.org/stash/dataset/doi:10.5061/dryad.dfn2z352k. The AVHRR GIMMS NDVI3g V1.2 data are available at https://zenodo.org/records/8253971. The AVHRR GIMMS LAI4g V1.2 data are available at https://zenodo.org/records/8281930. The LCSIF product are available at https://zenodo.org/records/7916851 and https://zenodo.org/records/7916879. The VIIRS GLSP product VNP22Q2 data are available at https://viirsland.gsfc.nasa.gov/Products/NASA/PhenologyESDR.html. The global vegetation map GLC 2000 is available at https://forobs.jrc.ec.europa.eu/products/glc2000/glc2000.php. The elevation data from WorldClim 2.1 is available at https://www.worldclim.org/data/worldclim21.html. The world map of Köppen-Geiger climate classification was freely obtained from http://koeppen-geiger.vu-wien.ac.at/present.htm.
