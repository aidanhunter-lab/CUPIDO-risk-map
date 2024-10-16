# CUPIDO-risk-map
**Mapping the distribution of microplastic's ecological impact throughout the Southern Ocean.**
![project logo](misc/CUPIDO-logo.png)

**License to download, share, fork and modify this software is temporarily revoked. Express permission must be granted by the code developer. The open source license to share will be reinstated as soon as possible.**

### Code summary

All compiled data required to generate maps are stored in CUPIDO-risk-map/data/. The original data sets are not included, but may be downloaded via some of the MatLab scripts or directly from source. The compiled data may then be reproduced from MatLab scripts included here. All files necessary to reproduce results are included in this repository.
The original data may be included at a later date (except for the propriatory ship traffic data which may only be available upon direct request).

All scripts for downloading data and preparing it for mapping were written in MatLab 2023b. At some point, for consistency and to run entirely on open source software, these may be converted into R.

Code required to generate all maps, including the interactive map, was written in R 4.4.1. The R Shiny interactive map will soon be publicly available from weblink: [southernoceanplastic.data.bas.ac.uk](southernoceanplastic.data.bas.ac.uk). To launch the interactive map directly from this code base, open an R session from within the CUPIDO-risk-map/R/ directory, then input the command `source('run_app.R')`. The interactive map may take a little longer to load the first time it is opened.

### Required code packages

* R packages and the version numbers used: cowplot 1.1.3, DT 0.33, flextable 0.9.6, ggiraph 0.8.10, ggnewscale 0.4.10, ggplot2 3.5.1, grid 4.4.1, gridExtra 2.3, gtable 0.3.5, mapdata 2.3.1, maps 3.4.2, patchwork 1.2.0, RColorBrewer 1.1-3, remotes 2.5.0, reshape2 1.4.4, scales 1.3.0, sf 1.0-16, shiny 1.8.1.1, sp 2.1-4, this.path 2.4.0, viridisLite 0.4.2
* Other software dependencies and the version numbers used: GEOS 3.10.2, GDAL 3.4.1, PROJ 8.2.1

