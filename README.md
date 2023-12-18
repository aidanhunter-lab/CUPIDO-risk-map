# CUPIDO-risk-map
**Mapping the distribution of microplastic's ecological impact throughout the Southern Ocean.**
![project logo](misc/CUPIDO-logo.png)

**License to download, share, fork and modify this software is temporarily revoked. Express permission must be granted by the code developer. The open source license to share will be reinstated as soon as possible.**

Code required to generate all maps, including the interactive map, is written in R. To run the interactive map, open an R session from within the CUPIDO-risk-map/R/ directory. Then, if using RStudio, open app.R and run the code; if R is running in a terminal then type 'library(shiny)' to load the rshiny package then 'runApp()' to launch the interactive map. Opening the interactive map may take a little longer the first time it's opened.

All scripts for downloading data and preparing it for mapping are written in MatLab. At some point, for consistency and to run entirely on open source software, these may be converted into R.

All compiled data required to generate maps are stored in CUPIDO-risk-map/data/. The original data sets are not included, but may be downloaded via some of the MatLab scripts or directly from source. The compiled data may then be reproduced from MatLab scripts included here. All files necessary to reproduce results are included in this repository.
The original data may be included at a later date (except for the propriatory ship traffic data which may only be available upon direct request).






