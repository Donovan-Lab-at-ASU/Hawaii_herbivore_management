[![DOI](https://zenodo.org/badge/573900863.svg)](https://zenodo.org/badge/latestdoi/573900863)

This repository includes data and analysis scripts to accompany:

# Evidence for managing herbivores for reef resilience
*under consideration for pubilication*

### Author of analysis and code: Mary K Donovan
-----

### Description:
This work analyzes data from the Hawaii Monitoring and Reporting Collaborative (HIMARC) to investigate multiple lines of evidence underpinning herbivore management for reef resilience.

### Contents:
#### Scripts:
* **1_herbivore_driver_model.R:** R script that defines a JAGS model for understanding the relationship between drivers and herbivore biomass, and executes the model.
* **1b_herbivore_driver_model_plot.R:** R script that reads in output of 1_herbivore_driver_model.R and creates figures and tables
* **2_map_and_summarize_herbivores.R:** R script that reads in output of 1_herbivore_driver_model.R and creates predictive map and spatial summaries of posterior estimates
* **3_herbivore-benthic.R:** R script for two analyses of herbivore-benthic relationships
* **3b_fishing_herbivores_benthos.R:** R script for predictions of fishing-herbivore and heribvore-benthic patterns
* **4_case_studies.R:** R script with analyses of parrotfish biomass across four fisheries management actions

#### /data:
Data files with herbivore and benthic surveys from the Hawaii Monitoring and Reporting Collaborative (HIMARC)
