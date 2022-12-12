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


# getwd()
setwd('~/Documents/Git Repos/CUPIDO-risk-map')

# Load map shape file -----------------------------------------------------
verbose <- FALSE

# nc1 <- st_read('data/map/Antarctic coastline polygons/high res/add_coastline_high_res_polygon_v7_6.shp', quiet = !verbose)
nc1 <- st_read('data/map/Antarctic coastline polygons/medium res/add_coastline_medium_res_polygon_v7_6.shp', quiet = !verbose)
nc2 <- st_read('data/map/sub Antarctic coastline polygons/sub_antarctic_coastline_high_res_polygon_v1.0.shp', quiet = !verbose)

# merge coastline data sets
nc1$location = 'Antarctica'
nc2$source <- NULL
nc2$surface <- 'land'
nc1 <- nc1[,c('location','surface','geometry')]
nc2 <- nc2[,c('location','surface','geometry')]
nc <- rbind(nc1, nc2)

# nc <- st_read('data/map/add_coastline_medium_res_polygon_v7.3.shp', quiet = !verbose)
nc$surface = factor(nc$surface, levels = unique(nc$surface))
# Coordinate reference system
crs_world = 4326
crs_use = 3031



# Load plastic data -----------------------------------------------------

# filename = 'plastic_quantity.csv'
filename = 'plastic_quantity_new2.csv'
filepath = 'data/plastic_quantity/'
DATA = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing ----
# Date format
DATA$Date_1 = as.character(as.Date(DATA$Date_1, '%d-%b-%Y'))
DATA$Date_2 = as.character(as.Date(DATA$Date_2, '%d-%b-%Y'))
DATA$Date_1[is.na(DATA$Date_1)] = ''
DATA$Date_2[is.na(DATA$Date_2)] = ''
# Date column shows specific sample date or a date range
DATA$`Sample Date` <- sapply(1:nrow(DATA), function(z){
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
DATA <- subset(DATA, select = -c(Date_1, Date_2))

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

# Column for interactive (mouse scrolling) labelling of position
west_ind <- DATA$Longitude < 0
east_ind <- !west_ind
lon_lab <- abs(round(DATA$Longitude, 2))
lon_lab[west_ind] <- paste(lon_lab[west_ind], 'W')
lon_lab[east_ind] <- paste(lon_lab[east_ind], 'E')
lat_lab <- paste(abs(round(DATA$Latitude, 2)), 'S')
coord_lab <- paste0('(', lon_lab, ', ', lat_lab, ')')
DATA$Coordinates <- coord_lab

DATA$Variable <- as.character(DATA$Variable)
DATA$Variable[DATA$Variable == 'massDensity'] <- 'mass density'
DATA$Variable[DATA$Variable == 'massConcentration'] <- 'mass concentration'

# plastic_units <- unique(DATA[,c('Variable','Unit')])
plastic_units <- unique(DATA[,c('SampleType', 'Variable','Unit')])
plastic_units$SampleType <- as.character(plastic_units$SampleType)
plastic_units$SampleType[grepl('water', plastic_units$SampleType) | 
                           grepl('marine', plastic_units$SampleType)] <- 'water'
plastic_units <- unique(plastic_units)
plastic_units$SampleType <- factor(plastic_units$SampleType, levels = unique(plastic_units$SampleType))

DATA <- subset(DATA, select = - Unit)

# Include a data_id variable
#head(DATA)
DATA$order <- 1:nrow(DATA)
d <- unique(DATA[c('Source', 'SampleAtStation',  'SampleID')])
# d <- unique(DATA[c('Source','SampleID')])
d$data_id <- paste('sample', 1:nrow(d))

DATA <- merge(DATA, d, sort = FALSE)
DATA <- DATA[order(DATA$order),]
DATA <- DATA[,names(DATA) != 'order']


# More variables that may need filtering here: SampleGear, LitterIDMethod,
# LitterCategory, LitterScale, PlasticForm, PlasticSize, Site, SiteCategory

# PlasticForm is certainly important -- R shiny filtering should be among a few
# groups (fragment, fibre, film, other/unspecified), but the data display should
# show the original categories...
DATA$PlasticForm_grouped <- as.character(DATA$PlasticForm)
# unique(DATA$PlasticForm_grouped)
fragments <- DATA$PlasticForm_grouped %in% c('fragment', 'sphere', 'particle', 'flake', 'granule', 'pellet')
fibres <- DATA$PlasticForm_grouped %in% c('fibre', 'line', 'line/fibre', 'filament')
film <- DATA$PlasticForm_grouped %in% c('film')
other <- DATA$PlasticForm_grouped %in% c('other', 'foam', '')
DATA$PlasticForm_grouped[fragments] = 'fragment'
DATA$PlasticForm_grouped[fibres] = 'fibre'
DATA$PlasticForm_grouped[film] = 'film'
DATA$PlasticForm_grouped[other] = 'other/unspecified'


# Mapping coordinates
DATA_sf = st_as_sf(DATA, coords = c("Longitude", "Latitude"), crs = crs_world)
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
STATIONS_sf = st_as_sf(STATIONS, coords = c("Longitude_DD", "Latitude_DD"), crs = crs_world)
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
# Research stations use line-based symbols (crosses/hatches)
# pltShapes <- c(3,4,8,7,9,10,12,13,14,11)

pltShapes <- c(1, 0, 5, 2, 6, 19) 

stationTypes <- levels(STATIONS$Type)
nstationTypes <- length(stationTypes)
pltSymbols <- data.frame(Class = rep('ResearchStation', nstationTypes), Type = stationTypes, symbol = pltShapes[1:nstationTypes])
# Plastic samples use filled plot symbols (21:25) that differ according to sample type.
pltShapes <- 21:25
sampleTypes <- levels(DATA$SampleType)
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

# Set (minimum) size of plot bounding box
bbox_map <- st_bbox(nc)
bbox_dat <- st_bbox(DATA_sf)
bbox_krill <- st_bbox(krill_poly)
bbox_chl <- st_bbox(chl_poly)
bbox <- setNames(c(min(bbox_map$xmin, bbox_dat$xmin, bbox_krill$xmin, bbox_chl$xmin),
                   min(bbox_map$ymin, bbox_dat$ymin, bbox_krill$ymin, bbox_chl$ymin),
                   max(bbox_map$xmax, bbox_dat$xmax, bbox_krill$xmax, bbox_chl$xmax),
                   max(bbox_map$ymax, bbox_dat$ymax, bbox_krill$ymax, bbox_chl$ymax)), names(bbox_map))
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


#fun_flextable(z, plastic_units, singleRowVars)

fun_flextable <- function(x, pu, singleRowVars, sampleType){
  # Function defining tables to interactively display when mouse hovers over mapped (plastic) data points.
  # Inputs: x = data, pu = plastic_units, singleRowVars = variables only taking one row of table
  anyMeasures <- any(!is.na(x[,!names(x) %in% singleRowVars]))
  if(anyMeasures){
    # Arrange data frame
    x <-suppressMessages(melt(x, id.vars = singleRowVars))
    x <- x[!is.na(x$value),]
    l <- strsplit(as.character(x$variable), '_')
    x$`Plastic type` <- sapply(l, function(z) z[1])
    x$Measure <- sapply(l, function(z) z[2])
    pu <- pu[pu$SampleType == sampleType,]
    w <- outer(pu$Variable, x$Measure, '==')
    n <- nrow(x)
    x$Value <- sapply(1:n, function(z) paste(x$value[z], pu$Unit[w[,z]]))
    if(n > 1) x[singleRowVars][2:n,] <- ''
    # Create flextable
    ft <- flextable(x, col_keys = c(singleRowVars, "Plastic type", "Measure", "Value"))#, 
    ft <- bold(ft, part = "header", bold = TRUE)
    ft <- set_table_properties(ft, layout = "autofit")
    as.character(htmltools_value(ft, ft.shadow = FALSE))
  }else{
    x <- cbind(x[singleRowVars], Measure = 'none')
    ft <- flextable(x)#, 
    ft <- bold(ft, part = "header", bold = TRUE)
    ft <- set_table_properties(ft, layout = "autofit")
    as.character(htmltools_value(ft, ft.shadow = FALSE))
  }
}


make_plot <- function(dat, background = 'none', components = 'all', ptSize = 6, legPtSize = 4, alpha = 0.6){#, pu = plastic_units){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to generate map plot, called from inside the server function after
  # data have been filtered by user selected inputs.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Inputs:
  # dat        = list of all required (reactively filtered) data sets
  # background = character specifying which background layer to plot
  # components = 'main', 'legend' or 'all' to output the map, the legend, or a list of both
  # alpha      = point transparency
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
  
  switch(components,
         
         main = {
           
           # Background layer
           switch(background,
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
                  }
           )
           
           # Coastline
           if(background == 'none'){
             plt_map <-
               ggplot() +
               geom_sf(data = nc,
                       aes(fill = surface),
                       show.legend = FALSE) +
               scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
           }else{
             plt_map <-
               plt_background +
               new_scale('fill') +
               geom_sf(data = nc,
                       aes(fill = surface),
                       show.legend = FALSE) +
               scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
           }
           
           plt_map <- plt_map +
             theme_void()
           
           # Combine plot layers
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
               # geom_sf(data = dat_stations,
               #         aes(shape = Type, colour = Seasonality),
               #         alpha = 1, size = 4, stroke = 1, show.legend = FALSE) +
               scale_colour_manual(values = c('forestgreen','firebrick'))
           }
           
           # Plastic samples
           if(anyPlastic){
             plt_map <- plt_map +
               new_scale('fill') +
               # geom_sf(data = dat_plastic,
               #         aes(fill = Source, shape = SampleType),
               #         alpha = 1, size = 4, show.legend = FALSE) +
               geom_sf_interactive(data = dat_plastic,
                                   aes(fill = Source, shape = SampleType, data_id = data_id, tooltip = tooltip),
                                   alpha = alpha, size = ptSize, show.legend = FALSE) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
           }
           
           plt_map <- plt_map +
             coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax']))# +
             # theme_void()
           
           output_plot <- ggdraw(plt_map)
           return(output_plot)
           
         },
         
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
         legend = {
           
           # Research stations
           if(anyStations){
             plt_stations <-
               ggplot() +
               geom_sf(data = dat_stations,
                       aes(shape = Type, colour = Seasonality),
                       alpha = 1, size = ptSize, stroke = 1) +
               scale_colour_manual(values = c('forestgreen','firebrick')) +
               scale_shape_manual(values = symbols$symbol[symbols$Class == 'ResearchStation']) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Facility', override.aes = list(size = legPtSize)),
                      colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
           }
           
           # Plastic samples
           if(anyPlastic){
             plt_plastic_samples <-
               ggplot() +
               geom_sf(data = dat_plastic,
                       aes(fill = Source, shape = SampleType),
                       alpha = 1, size = ptSize) +
               # geom_sf_interactive(data = dat_plastic,
               #                     aes(fill = Source, shape = SampleType, data_id = SampleType, tooltip = SampleType),
               #                     alpha = 1, size = 4) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source)) +
               scale_shape_manual(values = symbols$symbol[symbols$Class == 'PlasticSample']) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Sample type', override.aes = list(size = legPtSize)),
                      fill = guide_legend(override.aes = list(shape = c(21), size = legPtSize)))
           }
           
           # Extract legends and make subplot
           if(anyStations & anyPlastic){
             leg_plastic <- get_legend(plt_plastic_samples)
             leg_stations <- get_legend(plt_stations)
             output_legend <- ggdraw(plot_grid(leg_plastic, leg_stations,
                                               ncol = 2, align = 'h'))
           }else{
             if(anyPlastic){
               leg_plastic <- get_legend(plt_plastic_samples)
               output_legend <- ggdraw(plot_grid(leg_plastic, NULL,
                                                 ncol = 2, align = 'h'))
             }else{
               if(anyStations){
                 leg_stations <- get_legend(plt_stations)
                 output_legend <- ggdraw(plot_grid(leg_stations, NULL,
                                                   ncol = 2, align = 'h'))
               }else{
                 output_legend <-ggdraw(plot_grid(NULL, NULL,
                                                  ncol = 2))
               }}}
           
           # # Extract legends and make subplot
           # if(anyStations & anyPlastic){
           #   leg_plastic <- get_legend(plt_plastic_samples)
           #   leg_stations <- get_legend(plt_stations)
           #   output_legend <- ggdraw(plot_grid(leg_plastic, leg_stations,
           #                                     ncol = 1, align = 'v'))
           # }else{
           #   if(anyPlastic){
           #     leg_plastic <- get_legend(plt_plastic_samples)
           #     output_legend <- ggdraw(plot_grid(leg_plastic, NULL,
           #                                       ncol = 1, align = 'v'))
           #   }else{
           #     if(anyStations){
           #       leg_stations <- get_legend(plt_stations)
           #       output_legend <- ggdraw(plot_grid(leg_stations, NULL,
           #                                         ncol = 1, align = 'v'))
           #     }else{
           #       output_legend <-ggdraw(plot_grid(NULL, NULL,
           #                                        ncol = 1))
           #     }}}
           
           return(output_legend)
           
         },
         
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         # components = 'all' is not necessary and included only for completeness
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         all = {
           
           # Background layer
           switch(background,
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
                  }
           )
           
           # Coastline
           if(background == 'none'){
             plt_map <-
               ggplot() +
               geom_sf(data = nc,
                       aes(fill = surface),
                       show.legend = FALSE) +
               scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
           }else{
             plt_map <-
               plt_background +
               new_scale('fill') +
               geom_sf(data = nc,
                       aes(fill = surface),
                       show.legend = FALSE) +
               scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
           }
           
           # Research stations
           if(anyStations){
             plt_stations <-
               ggplot() +
               geom_sf(data = dat_stations,
                       aes(shape = Type, colour = Seasonality),
                       alpha = 1, size = ptSize, stroke = 1) +
               scale_colour_manual(values = c('forestgreen','firebrick')) +
               scale_shape_manual(values = symbols$symbol[symbols$Class == 'ResearchStation']) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Facility', override.aes = list(size = legPtSize)),
                      colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
           }
           
           # Plastic samples
           if(anyPlastic){
             plt_plastic_samples <-
               ggplot() +
               geom_sf(data = dat_plastic,
                       aes(fill = Source, shape = SampleType),
                       alpha = 1, size = ptSize) +
               # geom_sf_interactive(data = dat_plastic,
               #                     aes(fill = Source, shape = SampleType, data_id = SampleType, tooltip = SampleType),
               #                     alpha = 1, size = 4) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source)) +
               scale_shape_manual(values = symbols$symbol[symbols$Class == 'PlasticSample']) +
               theme(legend.key = element_blank()) +
               guides(shape = guide_legend(title = 'Sample type', override.aes = list(size = legPtSize)),
                      fill = guide_legend(override.aes = list(shape = c(21), size = legPtSize)))
           }
           
           plt_map <- plt_map +
             theme_void()
           
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
               # geom_sf(data = dat_stations,
               #         aes(shape = Type, colour = Seasonality),
               #         alpha = 1, size = 4, stroke = 1, show.legend = FALSE) +
               scale_colour_manual(values = c('forestgreen','firebrick'))
           }
           
           # Plastic samples
           if(anyPlastic){
             plt_map <- plt_map +
               new_scale('fill') +
               # geom_sf(data = dat_plastic,
               #         aes(fill = Source, shape = SampleType),
               #         alpha = 1, size = 4, show.legend = FALSE) +
               geom_sf_interactive(data = dat_plastic,
                                   aes(fill = Source, shape = SampleType, data_id = data_id, tooltip = tooltip),
                                   alpha = alpha, size = ptSize, show.legend = FALSE) +
               scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
           }
           
           # Make main plot
           plt_map <- plt_map +
             coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax']))# +
             # theme_void()
           
           output_plot <- ggdraw(plt_map)
           
           
           # Extract legends and make subplot
           if(anyStations & anyPlastic){
             leg_plastic <- get_legend(plt_plastic_samples)
             leg_stations <- get_legend(plt_stations)
             output_legend <- ggdraw(plot_grid(leg_plastic, leg_stations,
                                               ncol = 1, align = 'v'))
           }else{
             if(anyPlastic){
               leg_plastic <- get_legend(plt_plastic_samples)
               output_legend <- ggdraw(plot_grid(leg_plastic, NULL,
                                                 ncol = 1, align = 'v'))
             }else{
               if(anyStations){
                 leg_stations <- get_legend(plt_stations)
                 output_legend <- ggdraw(plot_grid(leg_stations, NULL,
                                                   ncol = 1, align = 'v'))
               }else{
                 output_legend <-ggdraw(plot_grid(NULL, NULL,
                                                  ncol = 1))
               }}}
           
           return(
             list(plot = output_plot, legend = output_legend)
           )
           
         }
  )
  
}



# Define UI for plastic data map app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Mapping Southern Ocean Plastic Data"),
  
  fluidRow(
    column(width = 3,
           wellPanel(
             # Input: year range -------------------------------------------------------
             sliderInput("YearRange", "Years:",
                         min = min(DATA$Year, na.rm = TRUE) , max = max(DATA$Year, na.rm = TRUE),
                         value = range(DATA$Year, na.rm = TRUE), step = 1, sep = ''),
             
             # Input: measurement variable ----
             selectInput("Variable", "Measurement:",
                         c("Concentration (pieces/m3)" = "concentration",
                           "Mass concentration (g/m3)" = "mass concentration",
                           "Density (pieces/km2)" = "density",
                           "Mass density (g/km2)" = "mass density"),
                         multiple = TRUE, selected = c("concentration", "mass concentration", "density", "mass density")),
             
             # Input: plastic type ----
             # I NEED TO THINK ABOUT HOW TO HANDLE THE TOTAL COLUMN... OMIT IT FOR NOW
             # BUT MAYBE INCLUDE IT LATER WHEN GRAPHING DATA...
             selectInput("PlasticForm_grouped", "Plastic form:",
                         c("Fragment" = "fragment",
                           "Fibre" = "fibre",
                           "Film" = "film",
                           "Other/unspecified" = 'other/unspecified'),
                         multiple = TRUE, selected = c("fragment", "fibre", "film", "other/unspecified")),
             
             # selectInput("Type", "Plastic type:",
             #             c("Fragment" = "fragment",
             #               "Fibre" = "fibre",
             #               "Film" = "film"),
             #             multiple = TRUE, selected = c("fragment", "fibre", "film")),
             
             # # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment) ----
             # selectInput("SampleType", "Sample type:",
             #             c("Near-surface" = "surface",
             #               "Subsurface" = "subsurface"),
             #             multiple = TRUE, selected = c("surface", "subsurface")),
             
             # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment) ----
             selectInput("SampleType", "Sample type:",
                         c("Sea surface" = "marine surface",
                           "Sea subsurface" = "marine subsurface",
                           "Freshwater" = "freshwater",
                           "Wastewater" = "wastewater",
                           "Sediment" = "sediment"),
                         multiple = TRUE, selected = c("marine surface","marine subsurface","freshwater","wastewater","sediment")),
             
             # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment) ----
             selectInput("StationType", "Facility:",
                         c("Station" = "Station",
                           "Camp" = "Camp",
                           "Refuge" = "Refuge",
                           "Airfield" = "Airfield Camp",
                           "Laboratory" = "Laboratory",
                           "Depot" = "Depot"),
                         multiple = TRUE, selected = c("Station", "Camp", "Refuge", "Airfield Camp", "Laboratory", "Depot")),
             
             # Input: background layers ----
             radioButtons("background", "Background data:",
                          c("None" = "none",
                            "Krill (Jan--Mar)" = "krill_all",
                            "Krill (Jan)" = "krill_1",
                            "Krill (Feb)" = "krill_2",
                            "Krill (Mar)" = "krill_3",
                            "Chlorophyll (Jan--Mar)" = "chl_all",
                            "Chlorophyll (Jan)" = "chl_1",
                            "Chlorophyll (Feb)" = "chl_2",
                            "Chlorophyll (Mar)" = "chl_3"
                          ),
                          selected = 'none')
             
             
             # # Input: plot value transformation ----
             # radioButtons("tran", "Scale:",
             #              c("Natural" = "norm",
             #                "Log" = "log")),
             
             # # Input: Checkbox for whether outliers should be included ----
             # checkboxInput("outliers", "Show outliers", TRUE),
             
             # width = 3
             
           )
    ),
    
    # Main panel for displaying outputs ----
    column(width = 7, #offset = 0, style='padding:0px;',
           tags$body(tags$div(id="ppitest", style="width:1in;visible:hidden;padding:0px")),
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
             h4("Selected plastic data"),
             dataTableOutput("datatab_plastic")
           )
           
    )
  ),
  
  fluidRow(
    column(width = 12,
           
           wellPanel(
             h4("Selected facilities"),
             dataTableOutput("datatab_stations")
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
               SampleType %in% x$SampleType
      )
    )    
  })
  
  # # Filter plastic data according to Shiny inputs
  # filtered_plastic_data <-reactive({
  #   x <- listInputs()
  #   return(
  #     subset(DATA,
  #            x$YearRange[1] <= Year & Year <= x$YearRange[2] &
  #              DATA$Variable %in% x$Variable &
  #              DATA$Type %in% x$Type &
  #              DATA$SampleType %in% x$SampleType
  #     )
  #   )    
  # })
  
  # filtered_plastic_data <-reactive({
  #   return(
  #     subset(DATA,
  #            input$YearRange[1] <= Year & Year <= input$YearRange[2] &
  #              DATA$Variable %in% input$Variable &
  #              DATA$Type %in% input$Type &
  #              DATA$SampleType %in% SampleType_d()
  #            # DATA$SampleType %in% input$SampleType
  #     )
  #   )
  # })
  
  transformed_plastic_data <- reactive({
    d <- filtered_plastic_data()
    anyPlastic <- nrow(d) > 0
    if(anyPlastic){
      dat_plastic_recast <- dcast(d,
                                  data_id + Source + SampleID + Site + SampleType + Depth + Longitude + Latitude + Coordinates +
                                    `Sample Date` + Year ~ PlasticForm + Variable, fun.aggregate = mean, value.var = 'Value')
      # Transform to correct CRS
      dat_plastic_recast <- st_as_sf(dat_plastic_recast, coords = c("Longitude", "Latitude"), crs = crs_world)
      dat_plastic_recast <- st_transform(dat_plastic_recast, crs_use)
      # Define interactive tooltip using flextable
      omitVars <- c('Source', 'SampleID', 'Site', 'SampleType',
                    'Year', 'geometry', 'data_id')
      temp_dat <- as.data.frame(dat_plastic_recast)[!names(dat_plastic_recast) %in% omitVars]
      temp_dat <- split(temp_dat, seq_len(nrow(temp_dat)))
      singleRowVars <- c('Sample Date', 'Coordinates', 'Depth')
      sType <- as.character(dat_plastic_recast$SampleType)
      sType[grepl('water', sType) | grepl('marine', sType)] <- 'water'
      dat_plastic_recast$tooltip <-  sapply(1:length(temp_dat), function(z) fun_flextable(temp_dat[[z]], plastic_units, singleRowVars, sType[z]))
      # dat_plastic_recast$tooltip <-  sapply(temp_dat, function(z) fun_flextable(z, plastic_units, singleRowVars, sType[z]))
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
      backgroundOptions <- c('none', 'krill', 'chl')
      t <- sapply(backgroundOptions, FUN = function(z) grepl(z, input$background))
      background <- names(which(t)) # get the chosen background
      return(background)
  })
  
  # Filter background data
  filtered_background_dat <- reactive({
    background <- which_background()
    background_dat <- switch(background,
                             krill = krill_poly,
                             chl = chl_poly
    ) # get background data
    if(background != 'none'){
      m <- strsplit(input$background, '_')[[1]][2] # get the chosen month
      background_dat <- subset(background_dat, month == m) # filter by month
    }
    return(background_dat)
  })
  
  # Filter plot symbols based on sample & station types
  Symbols <- reactive({
    x <- listInputs()
      return(
        subset(pltSymbols, Type %in% c(x$SampleType, input$StationType))
        # subset(pltSymbols, Type %in% c(input$SampleType, input$StationType))
      )
  })
  
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
      make_plot(listData(), which_background(), 'main')
    )
  })

  plot_legend <- reactive({
    return(
      make_plot(listData(), which_background(), 'legend')
    )
  })

  
  output$plt <- renderGirafe({
    p <- plot_main()
    w <- blankwidth()
    h <- blankheight()
    x <- girafe(code = print(p),
                options = list(opts_sizing(rescale = FALSE),
                               # opts_tooltip(css = "padding:5px;background:white;border-radius:2px 2px 2px 2px;font-size:12pt;"),
                               opts_tooltip(css = "background:white;"),
                               opts_zoom(min = 0.5, max = 10),
                               opts_hover(css = "opacity:1.0;stroke-width:4;cursor:pointer;", reactive = TRUE),
                               # opts_hover(css = "opacity:1.0;stroke-width:2;r:6pt;width:12pt;height:12pt;cursor:pointer;", reactive = TRUE),
                               opts_hover_inv(css = "opacity:0.2;cursor:pointer;"),
                               # opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE)),
                               opts_selection(type = "multiple", css = "opacity:1.0;stroke-width:4;")),
                             # opts_selection(type = "multiple", css = "opacity:1.0;stroke-width:2;r:6pt;width:12pt;height:12pt;")),
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
    # d <- subset(d, select = -c(Stat, Year, Coordinates, data_id))
    d <- subset(d, select = -c(Year, Coordinates, data_id))
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
