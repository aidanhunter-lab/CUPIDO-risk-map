# Mapping microplastic's ecological impact in the Southern Ocean

## Overview
<img src="misc/CUPIDO-logo.png" align="right" width="200">
This work is output from the UKRI Future Leaders Fellowship project <a href="https://www.bas.ac.uk/project/cupido/">CUPIDO</a> (MR/T020962/1), that is an investigation of how zooplankton interact with marine microplastic and the potential consequences for the carbon cycle.
The repository contains the code required to reproduce the study entitled <i>Microplastic hotspots mapped across the Southern Ocean reveal areas of potential ecological impact</i>, published in Nature Scientific Reports by Hunter et al. (2024).
The study objective was to assess and visualise the spatial distribution of ecological risk posed by microplastics in the Southern Ocean.
The approach taken was a spatial overlap analysis of multiple datasets comprising plankton, physical ocean variables, ship traffic, research stations, and <i>in situ</i> microplastic samples.
These data may be explored, filtered, and downloaded using an interactive web application: <a href="https://southernoceanplastic.data.bas.ac.uk/">southernoceanplastic.data.bas.ac.uk</a>.

## Code summary
#### Data
All compiled data required to generate output maps are stored in CUPIDO-risk-map/data/. The original data sets are not included, but may be downloaded via some of the MatLab scripts or directly from source. The compiled data may then be reproduced from MatLab scripts included here. The original data may be included at a later date (except for the proprietary ship traffic data which are conditionally available upon direct request).

#### Software
All scripts for downloading data and preparing them for mapping were written in MatLab 2023b. At some point, for consistency and to run entirely on open source software, these may be converted into R.

Code required to generate all maps, including the interactive map, was written in R 4.4.1.
The simplest way to view the interactive map is from the weblink (<a href="https://southernoceanplastic.data.bas.ac.uk/">southernoceanplastic.data.bas.ac.uk</a>), but it may also be launched directly from this code base by opening an R session from within the CUPIDO-risk-map/R/ directory then inputting the command `source('run_app.R')`.
<!-- The simplest way to view the interactive map is from the weblink ([southernoceanplastic.data.bas.ac.uk](southernoceanplastic.data.bas.ac.uk)), but it may also be launched directly from this code base by opening an R session from within the CUPIDO-risk-map/R/ directory then inputting the command `source('run_app.R')`. -->
Running the interactive map locally requires installation of several software dependencies (see below), and it may take a little longer to load the first time it is opened.

#### Packages
* R packages and the version numbers used: cowplot 1.1.3, DT 0.33, flextable 0.9.6, ggiraph 0.8.10, ggnewscale 0.4.10, ggplot2 3.5.1, grid 4.4.1, gridExtra 2.3, gtable 0.3.5, mapdata 2.3.1, maps 3.4.2, patchwork 1.2.0, RColorBrewer 1.1-3, remotes 2.5.0, reshape2 1.4.4, scales 1.3.0, sf 1.0-16, shiny 1.8.1.1, sp 2.1-4, this.path 2.4.0, viridisLite 0.4.2
* Other software dependencies and the version numbers used: GEOS 3.10.2, GDAL 3.4.1, PROJ 8.2.1

