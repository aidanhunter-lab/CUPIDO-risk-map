library(shiny)
library(sp)
library(sf)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)
library(ggnewscale)
library(RColorBrewer)
# library(plotly) # doesn't work with layering grobs...
library(ggiraph) # this is good for interactivity, but cannot be converted to grobs for layout...
# library(ggpubr) # maybe this package's 'ggarrange' function can help.

library(reshape2)
library(flextable) # for displaying tables when mouse hovers over mapped data points
# library(tidyverse)
library(DT)

# library(shinybrowser) # get browser dimension & other info

# CAN THE CODE EFFICIENCY (MAP PLOTTING SPEED) BE INCREASED BY DEFINING THE
# TOOLTIP DATA TABLE IN THE PREAMBLE? THEN R-SHINY FILTERS WILL SUBSET THIS
# PRE-CALCULATED TABLE, RATHER THAN RECREATING THE TOOLTIP TABLE WITH EACH NEW
# SELECTION OF FILTERING VARIABLES...
# THIS WILL REQUIRE CASTING TO WIDE-FORM THEN MELTING BACK TO LONG-FORM...

# getwd()
setwd('~/Documents/Git Repos/CUPIDO-risk-map')

# Load map shape files -----------------------------------------------------
verbose <- FALSE

# nc1 <- st_read('data/map/Antarctic coastline polygons/high res/add_coastline_high_res_polygon_v7_6.shp', quiet = !verbose)
nc1 <- st_read('data/map/Antarctic coastline polygons/medium res/add_coastline_medium_res_polygon_v7_6.shp', quiet = !verbose)
nc2 <- st_read('data/map/sub Antarctic coastline polygons/sub_antarctic_coastline_high_res_polygon_v1.0.shp', quiet = !verbose)

# Merge coastline data sets
nc1$location = 'Antarctica'
nc2$source <- NULL
nc2$surface <- 'land'
nc1 <- nc1[,c('location','surface','geometry')]
nc2 <- nc2[,c('location','surface','geometry')]
nc <- rbind(nc1, nc2)

# Coordinate reference system
crs_world <- 4326
crs_use <- 3031
# attr(crs_world, 'class') <- 'crs'

# Reset coastline bounding box to map limits -- defined by most northern latitude
latmax <- -50
mapLimits <- data.frame(Longitude = -180:180, 
                        Latitude = rep(latmax, 361))
mapLimits <- st_as_sf(mapLimits, coords = 1:2, crs = crs_world)
mapLimits <- st_transform(mapLimits, crs_use)
attr(st_geometry(nc), 'bbox') <- st_bbox(mapLimits)

# bb_nc1 <- st_bbox(nc1)
# bb_nc2 <- st_bbox(nc2)
# bb_nc <- setNames(
#   c(min(bb_nc1['xmin'], bb_nc2['xmin']), min(bb_nc1['ymin'], bb_nc2['ymin']),
#     max(bb_nc1['xmax'], bb_nc2['xmax']), max(bb_nc1['ymax'], bb_nc2['ymax'])),
#   c('xmin','ymin','xmax','ymax'))
# attr(bb_nc, 'class') <- 'bbox'
# attr(st_geometry(nc), 'bbox') <- bb_nc

  
# nc <- st_read('data/map/add_coastline_medium_res_polygon_v7.3.shp', quiet = !verbose)
nc$surface = factor(nc$surface, levels = unique(nc$surface))


# Load ecoregion shape files -----------------------------------------------------
# verbose = TRUE
eco <- st_read('data/marine ecoregions/data/commondata/data0/meow_ecos_expl_clipped_expl.shp', quiet = !verbose)

# Filter ecoregions by province
keepProvinces <- c(60, 61, 48, 59, 62) # 60=Scotia Sea, 61=Continental High Antarctic, 48=Magellenic, 59=Subantarctic Islands, 62=Subantarctic New Zealand
eco <- eco[eco$PROV_CODE %in% keepProvinces,]

# Change coordinate reference system
eco <- st_transform(eco, crs_use)

# Match bounding box to coastline data
attr(st_geometry(eco), 'bbox') <- st_bbox(nc)

# Crop ecoregions to map limits
eco <- st_crop(eco, st_bbox(eco))




# Load shipping data  -----------------------------------------------------

ship <- read.table('data/shipping/McCarthy_2022/ship activity in ecoregions_figure 3', sep = ';', header = TRUE)

ship$time_total <- ship$Number.of.visits * ship$Mean.time.in.ecoregion..days.

# Match to spatial data on ecoregions
ship_poly <- eco[eco$ECO_CODE_X %in% ship$Ecoregion.number,]
ship_poly$Number.of.visits <- sapply(ship_poly$ECO_CODE_X, function(z) ship$Number.of.visits[ship$Ecoregion.number == z])
ship_poly$Mean.time.in.ecoregion <- sapply(ship_poly$ECO_CODE_X, function(z) ship$Mean.time.in.ecoregion..days.[ship$Ecoregion.number == z])
ship_poly$time_total <- sapply(ship_poly$ECO_CODE_X, function(z) ship$time_total[ship$Ecoregion.number == z])

# Load plastic data -----------------------------------------------------

# filename = 'plastic_quantity.csv'
filename = 'plastic_quantity_new2.csv'
filepath = 'data/plastic_quantity/'
DATA = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing ----

# Include spaces in names of data for better display
DATA_names_orig <- names(DATA)
DATA_names_spaces <- DATA_names_orig
DATA_names_spaces[c(2:9, 11:12)] <- c('Sample type','Sample gear','Litter ID method','Litter category','Litter scale','Plastic form','Plastic size','Sample at station','Site category','Sample ID')
DATA_name_swap <- as.data.frame(cbind(original = DATA_names_orig, spaces = DATA_names_spaces))
# names(DATA) <- DATA_name_swap$spaces

# Date format
DATA$Date_1 = as.character(as.Date(DATA$Date_1, '%d-%b-%Y'))
DATA$Date_2 = as.character(as.Date(DATA$Date_2, '%d-%b-%Y'))
DATA$Date_1[is.na(DATA$Date_1)] = ''
DATA$Date_2[is.na(DATA$Date_2)] = ''
# Date column shows specific sample date or a date range
DATA$`Sample date` <- sapply(1:nrow(DATA), function(z){
  d1 <- DATA$Date_1[z]
  d2 <- DATA$Date_2[z]
  if(d1 != '' & d2 != ''){
    paste0(d1, '\n', '-> ', d2)} else{
      if(d1 != '' & d2 == ''){
        d1} else{
          'not specified'}}
})

# Include Year column
DATA$Year = as.integer(substr(DATA$Date_1, 1, 4))
# DATA <- subset(DATA, select = -c(Date_1, Date_2))

# Sample type -- sub/surface water, sediment, etc
s <- as.character(DATA$SampleType)
marine <- s == 'seawater'
DATA$Depth <- as.character(DATA$Depth)
dep <- DATA$Depth
surf <- dep %in% c('5m','<1m', '5 m', '<1 m', 'surface') | suppressWarnings(as.numeric(dep) <= 5)
subsurf <- dep %in% c('subsurface')  | suppressWarnings(as.numeric(dep) > 5)
surf <- marine & surf
subsurf <- marine & subsurf
surf[is.na(surf)] <- FALSE
subsurf[is.na(subsurf)] <- FALSE
s[surf] <- 'marine surface'
s[subsurf] <- 'marine subsurface'
DATA$SampleType <- factor(s, levels = unique(s))


# Reorder data sources by publication date
sources <- levels(DATA$Source)
nsources <- length(sources)
ncs <- nchar(sources)
sourceYear <- suppressWarnings(
  as.integer(sapply(1:nsources, FUN = function(z) substr(sources[z], ncs[z]-4, ncs[z]-1)))
)
DATA$Source <- factor(DATA$Source, levels = sources[order(sourceYear)])

# Column for interactive (mouse scrolling) labelling of position -- samples at a
# fixed station have single lat-lon coordinates; longer tows have start/end coords.
# Fixed stations...
west_ind <- DATA$Longitude < 0
east_ind <- DATA$Longitude >= 0
lon_lab <- abs(round(DATA$Longitude, 2))
lon_lab[west_ind] <- paste(lon_lab[west_ind], 'W')
lon_lab[east_ind] <- paste(lon_lab[east_ind], 'E')
lat_lab <- paste(abs(round(DATA$Latitude, 2)), 'S')
coord_lab <- paste0('(', lon_lab, ', ', lat_lab, ')')
DATA$Coordinates[DATA$SampleAtStation] <- coord_lab[DATA$SampleAtStation]
# Longer tows...
notAtStation <- !DATA$SampleAtStation
west_ind1 <- DATA$Longitude_start < 0 & notAtStation
west_ind2 <- DATA$Longitude_end < 0 & notAtStation
east_ind1 <- DATA$Longitude_start >= 0 & notAtStation
east_ind2 <- DATA$Longitude_end >= 0 & notAtStation
lon_lab1 <- abs(round(DATA$Longitude_start, 2))
lon_lab2 <- abs(round(DATA$Longitude_end, 2))
lon_lab1[west_ind1] <- paste(lon_lab1[west_ind1], 'W')
lon_lab2[west_ind2] <- paste(lon_lab2[west_ind2], 'W')
lon_lab1[east_ind1] <- paste(lon_lab1[east_ind1], 'E')
lon_lab2[east_ind2] <- paste(lon_lab2[east_ind2], 'E')
lat_lab1 <- paste(abs(round(DATA$Latitude_start, 2)), 'S')
lat_lab2 <- paste(abs(round(DATA$Latitude_end, 2)), 'S')
coord_lab1 <- paste0('(', lon_lab1, ', ', lat_lab1, ')')
coord_lab2 <- paste0('(', lon_lab2, ', ', lat_lab2, ')')
DATA$`Coordinates (start)`[notAtStation] <- coord_lab1[notAtStation]
DATA$`Coordinates (end)`[notAtStation] <- coord_lab2[notAtStation]


DATA$Variable <- as.character(DATA$Variable)
DATA$Variable[DATA$Variable == 'massDensity'] <- 'mass density'
DATA$Variable[DATA$Variable == 'massConcentration'] <- 'mass concentration'
DATA$Variable <- factor(DATA$Variable, levels = c('concentration', 'density', 'mass concentration', 'mass density'))


# plastic_units <- unique(DATA[,c('Variable','Unit')])
plastic_units <- unique(DATA[,c('SampleType', 'Variable','Unit')])
plastic_units$SampleType <- as.character(plastic_units$SampleType)
plastic_units$SampleType[grepl('water', plastic_units$SampleType) | 
                           grepl('marine', plastic_units$SampleType)] <- 'water'
plastic_units <- unique(plastic_units)
plastic_units$SampleType <- factor(plastic_units$SampleType, levels = unique(plastic_units$SampleType))

# DATA <- subset(DATA, select = - Unit)


# More variables that may need filtering here: SampleGear, LitterIDMethod,
# LitterCategory, LitterScale, PlasticForm, PlasticSize, Site, SiteCategory

# PlasticForm is certainly important -- R shiny filtering should be among a few
# groups (fragment, fibre, film, other/unspecified), but the data display should
# show the original categories...
DATA$PlasticForm <- as.character(DATA$PlasticForm)
DATA$PlasticForm[DATA$PlasticForm == ''] = 'unspecified'
allPlasticForms <- unique(DATA$PlasticForm)
fragments <- c('fragment', 'sphere', 'particle', 'flake', 'granule', 'pellet')
fibres <- c('fibre', 'line', 'line/fibre', 'filament')
film <- c('film')
other <- allPlasticForms[!allPlasticForms %in% c(fragments, fibres, film, 'all')]
DATA$PlasticForm_grouped <- DATA$PlasticForm
DATA$PlasticForm_grouped[DATA$PlasticForm %in% fragments] = 'fragment'
DATA$PlasticForm_grouped[DATA$PlasticForm %in% fibres] = 'fibre'
DATA$PlasticForm_grouped[DATA$PlasticForm %in% film] = 'film'
DATA$PlasticForm_grouped[DATA$PlasticForm %in% other] = 'other/unspecified'

DATA$PlasticForm <- factor(DATA$PlasticForm, levels = c(fragments, fibres, film, other, 'all'))
DATA$PlasticForm_grouped <- factor(DATA$PlasticForm_grouped, levels = c('fragment','fibre','film','other/unspecified','all'))
# names(DATA)[names(DATA) == 'PlasticForm'] <- 'Plastic form' # rename column for better presentation in tooltip tables


# Plastic scale -- SCALES NEED SPECIFIED BETTER IN DATA SET SO THE BELOW WILL NEED ADJUSTED...
DATA$LitterScale <- as.character(DATA$LitterScale)
DATA$LitterScale[grepl('micro and meso', DATA$LitterScale)] <- 'micro and meso'
DATA$LitterScale[DATA$LitterScale %in% c('see size data','')] <- 'unspecified'
DATA$LitterScale <- factor(DATA$LitterScale, levels = c('micro', 'meso', 'micro and meso', 'macro', 'unspecified'))
# names(DATA)[names(DATA) == 'LitterScale'] <- 'Litter scale' # rename column for better presentation in tooltip tables


# Plastic size
DATA$PlasticSize <- as.character(DATA$PlasticSize)
DATA$PlasticSize[DATA$PlasticSize == ''] = 'unspecified'
DATA$PlasticSize <- factor(DATA$PlasticSize, levels = unique(DATA$PlasticSize))
# names(DATA)[names(DATA) == 'PlasticSize'] <- 'Plastic size' # rename column for better presentation in tooltip tables


# Litter category -- only include plastic litter
# unique(DATA$LitterCategory)
DATA$LitterCategory <- as.character(DATA$LitterCategory)
DATA <- subset(DATA, LitterCategory == 'plastic')

# # How many LitterScales are in each sample?
# ns = rep(0, length(unique(DATA$data_id)))
# for(i in 1:length(ns)){
#   d = DATA[DATA$data_id == paste('sample', i),]
#   ns[i] = length(unique(d$LitterScale))
# }
# x = cbind(ns, unique(DATA$data_id))
# DATA[DATA$data_id %in% x[x[,1] == '3',2],]

# Measurement statistic
DATA$Statistic <- as.character(DATA$Statistic)
DATA$Statistic[DATA$Statistic == 'stddev'] = 's.d.'
DATA$Statistic <- factor(DATA$Statistic, levels = c('raw','mean','s.d.','min','max'))


# Order the data
DATA <- DATA[order(DATA$Source, DATA$SampleID, DATA$LitterScale, DATA$PlasticSize, 
                   DATA$PlasticForm, DATA$Variable, DATA$Statistic, DATA$Replicate),]
# DATA <- DATA[order(DATA$Source, DATA$SampleID, DATA$`Litter scale`, DATA$`Plastic size`, 
#                    DATA$`Plastic form`, DATA$Variable, DATA$Statistic, DATA$Replicate),]



# Include a data_id variable unique to each sample
#head(DATA)
DATA$order <- 1:nrow(DATA)
d <- unique(DATA[c('Source', 'SampleAtStation',  'SampleID')])
d$data_id <- paste('sample', 1:nrow(d))

DATA <- merge(DATA, d, sort = FALSE)
DATA <- DATA[order(DATA$order),]
DATA <- DATA[,names(DATA) != 'order']

# for(i in 1:nrow(DATA_name_swap)){
#   names(DATA)[names(DATA) == DATA_name_swap$original[i]] = DATA_name_swap$spaces[i]
# }


# Mapping coordinates
DATA_sf = st_as_sf(DATA, coords = c('Longitude', 'Latitude'), crs = crs_world)
DATA_sf = st_transform(DATA_sf, crs_use)

# Load research station data ----------------------------------------------
filename = 'COMNAP_Antarctic_Facilities_Master.csv'
filepath = 'data/research stations/COMNAP/'
STATIONS = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing ----
STATIONS$Type = reorder(STATIONS$Type, rep(1, nrow(STATIONS)), sum, decreasing = TRUE)

# Define tooltip variable for interactive info...

# Attribute values cannot contain single quotes "'", so remove them... I want to
# find a way to include these...
STATIONS$Official_Name <- gsub("'", "", STATIONS$Official_Name)
stationTypes <- levels(STATIONS$Type)
STATIONS$order <- 1:nrow(STATIONS)

sl <- lapply(stationTypes, FUN = function(z){
  n <- gsub(z, '', STATIONS$Official_Name[STATIONS$Type == z])
  o <- STATIONS$order[STATIONS$Type == z]
  matrix(c(o,n), ncol = 2)
})
d <- as.data.frame(do.call('rbind', sl))
d <- d[order(as.numeric(d[,1])),]
STATIONS$Official_Name <- d[,2]
# remove spaces at end of names
n <- nchar(STATIONS$Official_Name)
i <- substr(STATIONS$Official_Name, n, n) == ' '  
w <- any(i)
while(w){
  STATIONS$Official_Name[i] <- substr(STATIONS$Official_Name[i], 1, n[i] - 1)
  n <- nchar(STATIONS$Official_Name)
  i <- substr(STATIONS$Official_Name, n, n) == ' '  
  w <- any(i)
}
# remove spaces at start of names
n <- nchar(STATIONS$Official_Name)
i <- substr(STATIONS$Official_Name, 1, 1) == ' '  
w <- any(i)
while(w){
  STATIONS$Official_Name[i] <- substr(STATIONS$Official_Name[i], 2, n[i])
  n <- nchar(STATIONS$Official_Name)
  i <- substr(STATIONS$Official_Name, 1, 1) == ' '  
  w <- any(i)
}
STATIONS$Official_Name <- gsub('- ', ' ', STATIONS$Official_Name)

STATIONS$tooltip <- paste0(STATIONS$Official_Name, ' ', STATIONS$Type, '\n ',
                           'Primary operator: ', STATIONS$Operator_primary, '\n ',
                           'Peak population size = ', STATIONS$Peak_Population)

STATIONS$Record_ID <- paste('facility', STATIONS$Record_ID)

# Mapping coordinates
STATIONS_sf = st_as_sf(STATIONS, coords = c('Longitude_DD', 'Latitude_DD'), crs = crs_world)
STATIONS_sf = st_transform(STATIONS_sf, crs_use)
xy = matrix(unlist(STATIONS_sf$geometry), 2, nrow(STATIONS))
STATIONS$x = xy[1,]
STATIONS$y = xy[2,]





# Load krill data ---------------------------------------------------------

filename = 'krill_data_mapped.csv'
filepath = 'MatLab/temp/'
KRILL = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing -----------------------------------------------------
KRILL <- KRILL[!is.nan(KRILL$value),] # omit missing data
# Mean abundance across months for each grid cell
KRILL_allMonths <- aggregate(value ~ lonmin + lonmax + latmin + latmax, KRILL, 'mean')
KRILL_allMonths$month <- 'all'
KRILL_allMonths <- KRILL_allMonths[,c(1:4,6,5)]
KRILL <- rbind(KRILL, KRILL_allMonths)
# Round to 2 significant figures
sigDigits <- 2
KRILL$value <- signif(KRILL$value, sigDigits)

# I think this function is too complicated -- maybe only need call to Polygon function, and not the Polygons function...
createPolygons <- function(dat){
  # Returns list of polygons describing lon/lat bounding box for each row of dat
  n = nrow(dat)
  lapply(1:n, FUN = function(z){
    x = dat[z,]
    Polygons(
      list(
        Polygon(
          cbind(
            c(x$lonmin, x$lonmax, x$lonmax, x$lonmin, x$lonmin),
            c(x$latmin, x$latmin, x$latmax, x$latmax, x$latmin)
          ))), paste0('row', z)
    )}
  )
}

krill_poly <- createPolygons(KRILL) # generate polygons
krill_poly <- SpatialPolygons(krill_poly, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
krill_poly <- st_as_sf(krill_poly)
krill_poly <- st_transform(krill_poly, crs = crs_use) # convert coordinates
krill_poly$month <- KRILL$month
krill_poly$value <- KRILL$value

# create a discrete colour scheme based on krill (log scale) abundance (value column)
ncols <- 10
clims <- 2 ^ {1:ncols-1}
clims[1] <- 0
colgroup <- paste(clims[1:{ncols-1}], clims[2:ncols], sep = '-')
colgroup[ncols] <- paste0('>', clims[ncols])
# colgroup[ncols+1] <- 'no data'
clims[ncols+1] <- Inf
for(i in 1:ncols){
  j <- clims[i] <=  krill_poly$value & krill_poly$value < clims[i+1]
  krill_poly$colourgroup[j] = colgroup[i]
}
# krill_poly$colourgroup[is.na(krill_poly$value)] <- colgroup[ncols+1]
krill_poly$colourgroup <- factor(krill_poly$colourgroup, levels = colgroup)



# Load chlorophyll data ---------------------------------------------------------

# The high-res data are too big for fast processing... It's quite slow to load
# and later functions are VERY slow -- I don't know where the bottle neck is.
# Perhaps dplyr could help... but it's probably best to simply reduce the
# resolution.

# filename = 'chl_data_mapped.csv'
filename = 'chl_data_mapped_highRes_3x1.csv'
# filename = 'chl_data_mapped_highRes_1x0.33.csv'
filepath = 'MatLab/temp/'
CHL = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing -----------------------------------------------------
CHL <- CHL[!is.nan(CHL$value),] # omit missing data
# Mean abundance across months for each grid cell
CHL_allMonths <- aggregate(value ~ lonmin + lonmax + latmin + latmax, CHL, 'mean')
CHL_allMonths$month <- 'all'
CHL_allMonths <- CHL_allMonths[,c(1:4,6,5)]
CHL <- rbind(CHL, CHL_allMonths)

# Round to 2 significant figures
CHL$value <- signif(CHL$value, sigDigits)

chl_poly <- createPolygons(CHL) # generate polygons
chl_poly <- SpatialPolygons(chl_poly, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
chl_poly <- st_as_sf(chl_poly)
chl_poly <- st_transform(chl_poly, crs = crs_use) # convert coordinates
chl_poly$month <- CHL$month
chl_poly$value <- CHL$value


# Plotting parameters -----------------------------------------------------
# Plot symbols - store in separate data frame
# Research stations use (mostly) unfilled symbols
pltShapes <- c(1, 0, 5, 2, 6, 19) 
stationTypes <- levels(STATIONS$Type)
nstationTypes <- length(stationTypes)
pltSymbols <- data.frame(Class = rep('ResearchStation', nstationTypes), Type = stationTypes, symbol = pltShapes[1:nstationTypes])
# Plastic samples use filled plot symbols (21:25) that differ according to sample type.
pltShapes <- 21:25
sampleTypes <- levels(DATA$SampleType)
# sampleTypes <- levels(DATA$'Sample type')
nsampleTypes <- length(sampleTypes)
pltSymbols <- rbind(pltSymbols,
                    data.frame(Class = rep('PlasticSample', nsampleTypes), Type = sampleTypes, symbol = pltShapes[1:nsampleTypes]))
pltSymbols$Type <- factor(pltSymbols$Type, levels = pltSymbols$Type)

# Plot colours
# Use a qualitative palette for plastic sources
pltColours <- brewer.pal(12, 'Paired') # the Set3 and Paired palettes has a maximum of 12 colours
pltColours <- pltColours[c(seq(2, 12, 2), seq(1, 11, 2))]
sources <- levels(DATA$Source)
nsources <- length(sources)
pltColours <- data.frame(Source = sources, colour = pltColours[1:nsources])

# Plot bounding box
bbox_map <- st_bbox(nc)
# st_bbox(eco)
bbox_dat <- st_bbox(DATA_sf)
bbox_krill <- st_bbox(krill_poly)
bbox_chl <- st_bbox(chl_poly)
# bbox <- setNames(c(min(bbox_map$xmin, bbox_dat$xmin, bbox_krill$xmin, bbox_chl$xmin),
#                    min(bbox_map$ymin, bbox_dat$ymin, bbox_krill$ymin, bbox_chl$ymin),
#                    max(bbox_map$xmax, bbox_dat$xmax, bbox_krill$xmax, bbox_chl$xmax),
#                    max(bbox_map$ymax, bbox_dat$ymax, bbox_krill$ymax, bbox_chl$ymax)), names(bbox_map))

bbox <- bbox_map

aspectRatio = unname(diff(bbox[c(1,3)]) / diff(bbox[c(2,4)]))

linebreaks <- function(n){HTML(strrep(br(), n))} # convenience function


# Plotting function -------------------------------------------------------

# Some nice examples, worth a look, presented here https://davidgohel.github.io/ggiraph/articles/offcran/shiny.html

# Options for flextables
set_flextable_defaults(
  font.size = 10,
  font.family = 'Roboto',
  theme_fun = 'theme_vanilla'
)


fun_flextable <- function(x, longVars, singleRowVars, sampleType){
  # Function defining tables to interactively display when mouse hovers over mapped (plastic) data points.
  # Inputs: x = data, pu = plastic_units, singleRowVars = variables only taking one row of table
  anyMeasures <- any(!is.na(x[,!names(x) %in% singleRowVars]))
  if(anyMeasures){
    # Arrange data frame
    x <-suppressMessages(melt(x, id.vars = singleRowVars)) # change to long-form
    x <- x[!is.na(x$value),]
    # Disentangle variables...
    l <- strsplit(as.character(x$variable), '_')
    y <- setNames(as.data.frame(do.call('rbind', l)), longVars)
    x <- cbind(x, y)
    # Create empty Value column -- this stores numeric measurement and unit
    x$Value <- x$value
    # Remove/compress unnecessary/redundant variables...
    # Measurement statistics
    if(all(c('min','max') %in% x$Statistic)){ # convert min/max into range
      x_ <- x[x$Statistic %in% c('min', 'max'),]
      n <- nrow(x_)
      v <- sapply(seq(1, n, 2), function(z) paste0('[', x_$value[z], ', ', x_$value[z+1], ']'))
      x$Value[x$Statistic == 'min'] <- v
      x$Statistic[x$Statistic == 'min'] <- 'range'
      x <- x[x$Statistic != 'max',]
    }
    if(all(c('mean','s.d.') %in% x$Statistic)){ # convert mean/s.d. into mean +- s.d.
      x_ <- x[x$Statistic %in% c('mean', 's.d.'),]
      n <- nrow(x_)
      v <- sapply(seq(1, n, 2), function(z) paste0(x_$value[z], ' +- ', x_$value[z+1]))
      # v <- sapply(seq(1, n, 2), function(z) bquote(.(x_$value[z]) %+-% .(x_$value[z+1])), simplify = 'array')
      x$Value[x$Statistic == 's.d.'] <- v
      x$Statistic[x$Statistic == 's.d.'] <- 'mean +- s.d.'
      # x$Statistic[x$Statistic == 's.d.'] <- sapply(1:{n/2}, function(...) bquote('mean' %+-% 's.d.'))
      x <- x[x$Statistic != 'mean',]
    }
    # Sample replicates
    if(length(unique(x$Replicate)) == 1){
      x <- x[,-which(names(x) == 'Replicate')]
      longVars <- longVars[longVars != 'Replicate']}
    # Include unit if Depth is numeric
    depth <- suppressWarnings(as.numeric(x$Depth))
    if(!any(is.na(depth)))
      x$Depth <- paste(depth, 'm')
    x$Measure <- x$Variable # MAYBE RENAME THE Variable COLUMN AS Measure IN THE MAIN DATA SET...
    # Identify correct unit and create Value column
    # pu <- pu[pu$SampleType == sampleType,]
    # w <- outer(pu$Variable, x$Measure, '==')
    n <- nrow(x)
    # x$Value <- sapply(1:n, function(z) paste(x$Value[z], pu$Unit[w[,z]]))
    x$Value <- paste(x$Value, x$Unit)
    if(x$SampleAtStation[1]){
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates (start)', 'Coordinates (end)', 'SampleAtStation')]
    }else{
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates', 'SampleAtStation')]}
    if(n > 1) x[singleRowVars][2:n,] <- ''
    
    # Include spaces in column names
    newnames <- sapply(1:ncol(x), function(z){
      n <- names(x)[z]
      if(n %in% DATA_name_swap$original) m <- DATA_name_swap$spaces[n == DATA_name_swap$original] else m <- n
      c(n,m)
    })
    names(x) <- newnames[2,]
    singleRowVars <- newnames[2,][newnames[1,] %in% singleRowVars]
    longVars <- newnames[2,][newnames[1,] %in% longVars]
    # # Remove quotes from column names
    # names(x) <- gsub('`', '', names(x))
    # singleRowVars <- gsub('`', '', singleRowVars)
    # longVars <- gsub('`', '', longVars)
    
    # Create flextable
    ft <- flextable(x, col_keys = c(singleRowVars, longVars, 'Value'))#, 
    # ft <- flextable(x, col_keys = c(singleRowVars, 'Plastic type', 'Measure', 'Value'))#, 
    ft <- bold(ft, part = 'header', bold = TRUE)
    ft <- set_table_properties(ft, layout = 'autofit')
    as.character(htmltools_value(ft, ft.shadow = FALSE))
  }else{
    if(x$SampleAtStation[1]){
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates (start)', 'Coordinates (end)', 'SampleAtStation')]
    }else{
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates', 'SampleAtStation')]}
    x <- cbind(x[singleRowVars], Measure = 'none')
    # Include spaces in column names
    newnames <- sapply(1:ncol(x), function(z){
      n <- names(x)[z]
      if(n %in% DATA_name_swap$original) m <- DATA_name_swap$spaces[n == DATA_name_swap$original] else m <- n
      c(n,m)
    })
    names(x) <- newnames[2,]
    # # Remove quotes from column names
    # names(x) <- gsub('`', '', names(x))
    ft <- flextable(x)#, 
    ft <- bold(ft, part = 'header', bold = TRUE)
    ft <- set_table_properties(ft, layout = 'autofit')
    as.character(htmltools_value(ft, ft.shadow = FALSE))
  }
}

# change function structure -- better efficiency, avoids errors...
make_plot <- function(dat, background = 'none', displayEcoregions = FALSE, components = 'all', ptSize = 6, legPtSize = 4, alpha = 0.6, polyLineWidth = 0.75){#, pu = plastic_units){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to generate map plot, called from inside the server function after
  # data have been filtered by user selected inputs.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Inputs:
  # dat               = list of all required (reactively filtered) data sets
  # background        = character specifying which background layer to plot
  # displayEcoregions = TRUE/FALSE: plot distinct ecoregions -- required for shipping data from McCarthy
  # components        = 'main', 'legend' or 'all' to output the map, the legend, or a list of both
  # alpha             = point transparency
  # Output:
  # A ggplot object containing interactive geoms compatible with girafe
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dat_background = dat$background
  dat_plastic = dat$plastic
  dat_stations = dat$stations
  symbols = dat$symbols
  
  anyPlastic <- nrow(dat_plastic) > 0
  anyStations <- nrow(dat_stations) > 0
  anyBackground <- background != 'none'
  
  # Background layer
  switch(background,
         none = {
           plt_background <- ggplot()
         },
         krill = {
           plt_background <-
             ggplot() +
             geom_sf(data = dat_background, aes(fill = colourgroup)) +
             scale_fill_viridis_d(option = 'plasma',
                                  name = bquote(atop(Krill, individuals / m^2)))
           },
         chl = {
           plt_background <-
             ggplot() +
             geom_sf(data = dat_background, aes(fill = value)) +
             scale_fill_viridis_c(option = 'viridis', trans = 'log10',
                                  name = bquote(atop(Chlorophyll, mg / m^3)))
         },
         ship = {
           plt_background <- 
             ggplot() +
             geom_sf(data = dat_background, aes(fill = time_total)) + 
             scale_fill_viridis_c(option = 'mako', trans = 'log10',
                                  # scale_fill_viridis_c(option = 'cividis', trans = 'log10',
                                  name = bquote(atop(Ship ~ time, days)))
         }
  )
  
  plt_background <- plt_background +
    theme_void() +
    theme(panel.background = element_rect(fill = 'white', colour = 'white') #,
#          legend.direction = 'horizontal',
          # plot.margin = unit(rep(-0.2,4), 'cm'),
#          legend.position = "bottom" #, legend.key.size = unit(0.4, 'cm')
    )
  
  
  switch(components,
         
         main = {
           # Produce main map -- no legend
           
           # Map layers
           plt_map <- plt_background + guides(fill = 'none') + theme(legend.position = 'right')
           
           # Ecoregions
           if(displayEcoregions){
             plt_map <- plt_map +
               geom_sf(data = eco, 
                       linewidth = polyLineWidth, colour = 'black', fill = alpha('white', 0))
           }
           
           # Coastline
           plt_map <- plt_map +
             new_scale('fill') +
             geom_sf(data = nc,
                     aes(fill = surface),
                     show.legend = FALSE) +
             scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
           
           if(anyStations | anyPlastic){
             plt_map <- plt_map +
               scale_shape_manual(values = setNames(symbols$symbol, symbols$Type))
           }
           
           # Research stations
           if(anyStations){
             plt_map <- plt_map +
               geom_sf_interactive(data = dat_stations,
                                   aes(shape = Type, colour = Seasonality, data_id = Record_ID, tooltip = tooltip),
                                   alpha = alpha, size = ptSize, stroke = 1, show.legend = FALSE) +
               scale_colour_manual(values = c('forestgreen','firebrick'))
           }
           
           # Plastic samples
           if(anyPlastic){
             plt_map <- plt_map + 
               new_scale('fill') +
               geom_sf_interactive(data = dat_plastic,
                                   aes(fill = Source, shape = SampleType, data_id = data_id, tooltip = tooltip),
                                   alpha = alpha, size = ptSize, show.legend = FALSE) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
           }
           
           plt_map <- plt_map +
             coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax']))# +

           output_plot <- ggdraw(plt_map)
           return(output_plot)
           
         },
  
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
         legend = {
           
           # Research stations
           if(anyStations){
             symbols_ <- symbols[symbols$Class == 'ResearchStation',]
             plt_stations <-
               ggplot() +
               geom_sf(data = dat_stations,
                       aes(shape = Type, colour = Seasonality),
                       alpha = 1, size = ptSize, stroke = 1) +
               scale_colour_manual(values = c('forestgreen','firebrick')) +
               scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Facility', override.aes = list(size = legPtSize)),
                      colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
             leg_stations <- get_legend(plt_stations)
           }else{
             leg_stations <- NULL
           }
           
           # Plastic samples
           if(anyPlastic){
             symbols_ = symbols[symbols$Class == 'PlasticSample',]
             plt_plastic_samples <-
               ggplot() +
               geom_sf(data = dat_plastic,
                       aes(fill = Source, shape = SampleType),
                       alpha = 1, size = ptSize) +
               # geom_sf_interactive(data = dat_plastic,
               #                     aes(fill = Source, shape = SampleType, data_id = SampleType, tooltip = SampleType),
               #                     alpha = 1, size = 4) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source)) +
               scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Sample type', override.aes = list(size = legPtSize)),
                      fill = guide_legend(override.aes = list(shape = c(21), size = legPtSize)))
             leg_plastic <- get_legend(plt_plastic_samples)
           }else{
             leg_plastic <- NULL
           }
           
           # Background
           if(anyBackground){
             leg_background <- get_legend(plt_background)
           }else{
             leg_background <- NULL
           }
           
           # output_legend <- ggdraw(
           #   plot_grid(
           #     plot_grid(
           #       leg_plastic, leg_stations,
           #       ncol = 2, align = 'h', axis = 't'),
           #     leg_background,
           #     ncol = 1, rel_heights = c(9, 1), axis = 'l', align = 'v'),
           #   clip = 'on')

           output_legend <- ggdraw(
             plot_grid(
               leg_background, leg_plastic, leg_stations,
               ncol = 3), #, align = 'h', axis = 't'),
           ) #clip = 'on')
           
           
           # library(patchwork)
           # output_legend <- (ggdraw(leg_plastic) + ggdraw(leg_stations)) / ggdraw(leg_background)
           # 
           # wrap_plots(ggdraw(leg_plastic), ggdraw(leg_stations), ggdraw(leg_background))

           return(output_legend)
           
         },
         
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         # components = 'all' is not necessary and included only for completeness
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         all = {
           # Main plot...
           # Map layers
           plt_map <- plt_background + guides(fill = 'none') + theme(legend.position = 'right')
           
           # Ecoregions
           if(displayEcoregions){
             plt_map <- plt_map +
               geom_sf(data = eco, 
                       linewidth = polyLineWidth, colour = 'black', fill = alpha('white', 0))
           }
           
           # Coastline
           plt_map <- plt_map +
             new_scale('fill') +
             geom_sf(data = nc,
                     aes(fill = surface),
                     show.legend = FALSE) +
             scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
           
           if(anyStations | anyPlastic){
             plt_map <- plt_map +
               scale_shape_manual(values = setNames(symbols$symbol, symbols$Type))
           }
           
           # Research stations
           if(anyStations){
             plt_map <- plt_map +
               geom_sf_interactive(data = dat_stations,
                                   aes(shape = Type, colour = Seasonality, data_id = Record_ID, tooltip = tooltip),
                                   alpha = alpha, size = ptSize, stroke = 1, show.legend = FALSE) +
               scale_colour_manual(values = c('forestgreen','firebrick'))
           }
           
           # Plastic samples
           if(anyPlastic){
             plt_map <- plt_map + 
               new_scale('fill') +
               geom_sf_interactive(data = dat_plastic,
                                   aes(fill = Source, shape = SampleType, data_id = data_id, tooltip = tooltip),
                                   alpha = alpha, size = ptSize, show.legend = FALSE) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
           }
           
           plt_map <- plt_map +
             coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax']))# +
           
           output_plot <- ggdraw(plt_map)
           
           # Legends...
           # Research stations
           if(anyStations){
             symbols_ <- symbols[symbols$Class == 'ResearchStation',]
             plt_stations <-
               ggplot() +
               geom_sf(data = dat_stations,
                       aes(shape = Type, colour = Seasonality),
                       alpha = 1, size = ptSize, stroke = 1) +
               scale_colour_manual(values = c('forestgreen','firebrick')) +
               scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Facility', override.aes = list(size = legPtSize)),
                      colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
             leg_stations <- get_legend(plt_stations)
           }else{
             leg_stations <- NULL
           }
           
           # Plastic samples
           if(anyPlastic){
             symbols_ = symbols[symbols$Class == 'PlasticSample',]
             plt_plastic_samples <-
               ggplot() +
               geom_sf(data = dat_plastic,
                       aes(fill = Source, shape = SampleType),
                       alpha = 1, size = ptSize) +
               # geom_sf_interactive(data = dat_plastic,
               #                     aes(fill = Source, shape = SampleType, data_id = SampleType, tooltip = SampleType),
               #                     alpha = 1, size = 4) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source)) +
               scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Sample type', override.aes = list(size = legPtSize)),
                      fill = guide_legend(override.aes = list(shape = c(21), size = legPtSize)))
             leg_plastic <- get_legend(plt_plastic_samples)
           }else{
             leg_plastic <- NULL
           }
           
           # Background
           if(anyBackground){
             leg_background <- get_legend(plt_background)
           }else{
             leg_background <- NULL
           }
           
           output_legend <- ggdraw(
             plot_grid(
               leg_background, leg_plastic, leg_stations,
               ncol = 3)
           )
           
           return(
             list(plot = output_plot, legend = output_legend)
           )
           
         }
  )
  
}

# Define UI for plastic data map app ----
ui <- fluidPage(
  
  # App title
  titlePanel('Mapping Southern Ocean Plastic Data'),
  
  fluidRow(
    column(width = 3,
           wellPanel(
             # Input: year range
             sliderInput('YearRange', 'Years:',
                         min = min(DATA$Year, na.rm = TRUE) , max = max(DATA$Year, na.rm = TRUE),
                         value = range(DATA$Year, na.rm = TRUE), step = 1, sep = ''),
             
             # Input: measurement variable
             selectInput('Variable', 'Measurement:',
                         c('Concentration (pieces/m3)' = 'concentration',
                           'Mass concentration (g/m3)' = 'mass concentration',
                           'Density (pieces/km2)' = 'density',
                           'Mass density (g/km2)' = 'mass density'),
                         multiple = TRUE, selected = c('concentration', 'mass concentration', 'density', 'mass density')),
             
             # Input: plastic type
             selectInput('PlasticForm_grouped', 'Plastic form:',
                         c('Fragment' = 'fragment',
                           'Fibre' = 'fibre',
                           'Film' = 'film',
                           'Other/unspecified' = 'other/unspecified'),
                         multiple = TRUE, selected = c('fragment', 'fibre', 'film', 'other/unspecified')),
             
             # selectInput('Type', 'Plastic type:',
             #             c('Fragment' = 'fragment',
             #               'Fibre' = 'fibre',
             #               'Film' = 'film'),
             #             multiple = TRUE, selected = c('fragment', 'fibre', 'film')),
             
             # Input: plastic scale
             selectInput('LitterScale', 'Plastic scale:',
                         c('Micro' = 'micro',
                           'Micro & meso' = 'micro and meso',
                           'Meso' = 'meso',
                           'Macro' = 'macro',
                           'Unspecified' = 'unspecified'),
                         multiple = TRUE, selected = c('micro', 'micro and meso', 'meso', 'macro', 'unspecified')),
             
             # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment)
             selectInput('SampleType', 'Sample type:',
                         c('Sea surface' = 'marine surface',
                           'Sea subsurface' = 'marine subsurface',
                           'Freshwater' = 'freshwater',
                           'Wastewater' = 'wastewater',
                           'Sediment' = 'sediment'),
                         multiple = TRUE, selected = c('marine surface','marine subsurface','freshwater','wastewater','sediment')),
             
             # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment)
             selectInput('StationType', 'Facility:',
                         c('Station' = 'Station',
                           'Camp' = 'Camp',
                           'Refuge' = 'Refuge',
                           'Airfield' = 'Airfield Camp',
                           'Laboratory' = 'Laboratory',
                           'Depot' = 'Depot'),
                         multiple = TRUE, selected = c('Station', 'Camp', 'Refuge', 'Airfield Camp', 'Laboratory', 'Depot')),

             # Input: display ecoregions
             checkboxInput('DisplayEcoregions','Display ecoregions',
                           value = FALSE),

             # Input: background layers
             radioButtons('background', 'Background data:',
                          c('None' = 'none',
                            'Krill (Jan--Mar)' = 'krill_all',
                            'Krill (Jan)' = 'krill_1',
                            'Krill (Feb)' = 'krill_2',
                            'Krill (Mar)' = 'krill_3',
                            'Chlorophyll (Jan--Mar)' = 'chl_all',
                            'Chlorophyll (Jan)' = 'chl_1',
                            'Chlorophyll (Feb)' = 'chl_2',
                            'Chlorophyll (Mar)' = 'chl_3',
                            'Shipping activity' = 'ship'
                          ),
                          selected = 'none')
             
             
             # # Input: plot value transformation
             # radioButtons('tran', 'Scale:',
             #              c('Natural' = 'norm',
             #                'Log' = 'log')),
             
             # # Input: Checkbox for whether outliers should be included
             # checkboxInput('outliers', 'Show outliers', TRUE),
             
             # width = 3
             
           )
    ),
    
    # Main panel for displaying outputs
    column(width = 7, #offset = 0, style='padding:0px;',
           tags$body(tags$div(id='ppitest', style="width:1in;visible:hidden;padding:0px")),
           tags$script('$(document).on("shiny:connected", function(e) {
                                    var d =  document.getElementById("ppitest").offsetWidth;
                                    Shiny.onInputChange("dpi", d);
                                });
                                $(window).resize(function(e) {
                                    var d =  document.getElementById("ppitest").offsetWidth;
                                    Shiny.onInputChange("dpi", d);
                                });
                            '),
           # h3("Data plot"),
           plotOutput('blank', width = '100%', height = '1px'),
           ggiraphOutput('plt')
    ),
    
    column(width = 2, offset = 0, style='padding:0px;',
           linebreaks(8), # add whitespace above legends
           plotOutput('legend')
    )
  ),
  
  fluidRow(
    column(width = 12, #offset = 3,
           
           wellPanel(
             h4('Selected plastic data'),
             dataTableOutput('datatab_plastic')
           )
           
    )
  ),
  
  fluidRow(
    column(width = 12,
           
           wellPanel(
             h4('Selected facilities'),
             dataTableOutput('datatab_stations')
           )
           
    )
  )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output, session) {
  
  # Updating some Shiny inputs takes a long time because plotting the plastic
  # data is costly. Use 'debounce' function to delay updates, preventing each
  # individual change being rendered provided input changes are made quickly.
  listInputs <- debounce({
    reactive({
      list(YearRange = input$YearRange,
           Variable = input$Variable,
           PlasticForm = input$PlasticForm_grouped,
           # PlasticForm_grouped = input$PlasticForm_grouped,
           LitterScale = input$LitterScale,
           SampleType = input$SampleType
      )
    })
  }, 2000)

  
  # Filter plastic data according to Shiny inputs
  filtered_plastic_data <-reactive({
    x <- listInputs()
    return(
      subset(DATA,
             x$YearRange[1] <= Year & Year <= x$YearRange[2] &
               Variable %in% x$Variable &
               PlasticForm_grouped %in% x$PlasticForm &
               LitterScale %in% x$LitterScale &
               SampleType %in% x$SampleType
      )
    )    
  })
  

  transformed_plastic_data <- reactive({
    d <- filtered_plastic_data()
    anyPlastic <- nrow(d) > 0
    if(anyPlastic){
      dn <- names(d)
      # Change to wide-form -- choose variables to produce one row per data_id
      dummyGroupingVariables <- c('PlasticForm_grouped')
      yvars <- c('LitterScale', 'PlasticSize', 'PlasticForm', 'Variable', 'Statistic', 'Replicate', 'Unit') # variables to change to wide format
      xvars <- dn[!dn %in% c(yvars, 'Value', dummyGroupingVariables)] # everything else, apart from Value and dummy grouping variables, stays in long format
      # Handle spaces in data column names
      xvars_ <- grepl(' ', xvars)
      xvars[xvars_] <- paste0('`', xvars[xvars_], '`')
      yvars_ <- grepl(' ', yvars)
      yvars[yvars_] <- paste0('`', yvars[yvars_], '`')
      dat_plastic_recast <- suppressMessages(
        dcast(d, list(xvars, yvars), value.var = 'Value'))
      wellOrganised <- nrow(dat_plastic_recast) == length(unique(d$data_id))
      if(!wellOrganised) warning('Data are not well organised! Some samples are associated with multiple rows of data: we require a single row per sample in the wide-form data set.')
      # Transform to correct CRS
      dat_plastic_recast <- st_as_sf(dat_plastic_recast, coords = c('Longitude', 'Latitude'), crs = crs_world)
      dat_plastic_recast <- st_transform(dat_plastic_recast, crs_use)
      xvars <- c(xvars, 'geometry')
      # Define interactive tooltip using flextable...
      # The yvars (now stretched into wide form) are displayed (long form) in
      # tooltip. Now also choose other variables to display in a single row.
      keepVars <- c('SampleGear', 'Depth', 'Sample date', 'Coordinates', 'Coordinates (start)', 'Coordinates (end)', 'SampleAtStation') # SampleAtStation is needed to choose appropriate coordinates display
      omitVars <- xvars[!xvars %in% keepVars]
      temp_dat <- as.data.frame(dat_plastic_recast)[!names(dat_plastic_recast) %in% omitVars]
      temp_dat <- split(temp_dat, seq_len(nrow(temp_dat)))
      sType <- as.character(dat_plastic_recast$SampleType) # units depend on SampleType
      sType[grepl('water', sType) | grepl('marine', sType)] <- 'water'
      dat_plastic_recast$tooltip <- sapply(1:length(temp_dat), function(z) fun_flextable(x = temp_dat[[z]], longVars = yvars, singleRowVars = keepVars, sampleType = sType[z])) # this is time-consuming and could maybe be moved to pre-amble/saved data set...
      return(
        dat_plastic_recast
      )
    }else{
      return(d)
    }
  })
  
  # Filter the research station data
  filtered_station_data = reactive({
      return(
        subset(STATIONS_sf,
               STATIONS_sf$Type %in% input$StationType
        )
      )
  })
  
  # Get background type
  which_background <- reactive({
      backgroundOptions <- c('none', 'krill', 'chl', 'ship')
      t <- sapply(backgroundOptions, FUN = function(z) grepl(z, input$background))
      background <- names(which(t)) # get the chosen background
      return(background)
  })
  
  # Filter background data
  filtered_background_dat <- reactive({
    background <- which_background()
    background_dat <- switch(background,
                             krill = krill_poly,
                             chl = chl_poly,
                             ship = ship_poly
    ) # get background data
    if(background != 'none'){
      x <- strsplit(input$background, '_')[[1]]
      if(length(x) == 2){
        m <- x[2] # get the chosen month
        background_dat <- subset(background_dat, month == m) # filter by month
      }
    }
    return(background_dat)
  })
  
  # Display ecoregions
  display_ecoregions <-reactive(
    return(input$DisplayEcoregions)
  )
  
  # Filter plot symbols based on sample & station types
  Symbols <- reactive({
    x <- listInputs()
      return(
        subset(pltSymbols, Type %in% c(x$SampleType, input$StationType))
        # subset(pltSymbols, Type %in% c(input$SampleType, input$StationType))
      )
  })
  
  # pltsymbols = subset(pltSymbols, Type %in% c(SampleTypes, as.character(pltSymbols$Type[1:6])))
  
  # Group all data in a named list
  listData <- reactive({
    # plastic <- filtered_plastic_data()
    plastic <- transformed_plastic_data()
    stations <- filtered_station_data()
    background <- filtered_background_dat()
    symbols <- Symbols()
    return(
      list(background = background, stations = stations, plastic = plastic, symbols = symbols)
    )
  })
  
    # Create blank plot to specify sizing/aspect ratio
    output$blank = renderPlot({
     ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point() + theme_void()
    })

    blankwidth = reactive({
      # this is the magic that makes it work
      bw = session$clientData$output_blank_width
      return(5/4 * bw) # the sizing on girafe functions is weird, this scaling by 5/4 is a bit of a hack method...
    })

    blankheight <- reactive({
      blankwidth() / aspectRatio
    })
  
  plot_main <- reactive({
    return(
      make_plot(dat = listData(), background = which_background(), displayEcoregions = display_ecoregions(), components = 'main')
    )
  })

  plot_legend <- reactive({
    return(
      make_plot(dat = listData(), background = which_background(), components = 'legend')
    )
  })

  
  output$plt <- renderGirafe({
    p <- plot_main()
    w <- blankwidth()
    h <- blankheight()
    x <- girafe(code = print(p),
                options = list(opts_sizing(rescale = FALSE),
                               # opts_tooltip(css = 'padding:5px;background:white;border-radius:2px 2px 2px 2px;font-size:12pt;'),
                               opts_tooltip(css = 'background:white;'),
                               opts_zoom(min = 0.5, max = 10),
                               opts_hover(css = 'opacity:1.0;stroke-width:4;cursor:pointer;', reactive = TRUE),
                               # opts_hover(css = 'opacity:1.0;stroke-width:2;r:6pt;width:12pt;height:12pt;cursor:pointer;', reactive = TRUE),
                               opts_hover_inv(css = 'opacity:0.2;cursor:pointer;'),
                               # opts_hover(css = 'fill:#FF3333;stroke:black;cursor:pointer;', reactive = TRUE)),
                               opts_selection(type = 'multiple', css = 'opacity:1.0;stroke-width:4;')),
                             # opts_selection(type = 'multiple', css = 'opacity:1.0;stroke-width:2;r:6pt;width:12pt;height:12pt;')),
                width_svg = (w / {input$dpi}),
                height_svg = (h / {input$dpi})
    )
    x
  }
  )
  
  output$legend <- renderPlot({
    p <- plot_legend()
    p
  }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Display table of selected plastic data
  
  # output$hover_data <- renderPrint({
  #   input$plt_hovered
  # })
  
  selected_data <- reactive({
    input$plt_selected
  })

  output$datatab_plastic <- renderDataTable({
    # d <- DATA
    d <- filtered_plastic_data()
    d <- d[d$data_id %in% selected_data(),]
    if(nrow(d) < 1) return(NULL)
    
    d <- d[DATA_name_swap$original] # display only columns appearing in original data set
    names(d) <- sapply(1:ncol(d), function(z){
      n <- names(d)[z]
      i <- DATA_name_swap$original == n
      if(any(i)) n <- DATA_name_swap$spaces[i]
      n
    })
    # d <- subset(d, select = -c(Year, Coordinates, data_id))
    row.names(d) <- NULL
    d <- datatable(d,
                   options = list(
                     paging = FALSE, searching = FALSE # ,fillContainer = TRUE,
                    # bPaginate = FALSE,  searching = FALSE, fillContainer = TRUE
                     # autoWidth = TRUE
                   ))
    d
  })

  output$datatab_stations <- renderDataTable({
    d <- filtered_station_data()
    d <- as.data.frame(d)
    d <- d[d$Record_ID %in% selected_data(),]
    if(nrow(d) < 1) return(NULL)
    d <- subset(d, select = -c(Record_ID, Elevation_Datum, tooltip, geometry, order, Webcam_URL))
    # d <- subset(d, select = -c(Elevation_Datum, tooltip, geometry, order, Webcam_URL))
    row.names(d) <- NULL
    d <- datatable(d,
                   options = list(
                     paging = FALSE, searching = FALSE, autoWidth = TRUE, scrollX = TRUE
                   ))
    d
  })
  
  
      
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Potential extra plots of data distributions...
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # output$summary <- renderPrint({
  #   summary(dat())
  # })
  # 
  # output$view <- renderTable({
  #   dat()
  #   # head(d())
  # })
  
  
  
  # ylabel <- reactive({
  #   return(
  #     switch(input$Variable,
  #            concentration = expression(concentration ~ (pieces / m^3)),
  #            density = expression(density ~ (pieces / km^2)),
  #            massDensity = expression(concentration ~ (g / km^2)))
  #   )
  # })
  
}

# Create Shiny app ----
shinyApp(ui, server)
