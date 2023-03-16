# Individual map plots for paper

# Packages ----------------------------------------------------------------

# library(shiny)
library(sp)
library(sf)
# library(mapdata)
library(maptools) # for the ContourLines2SLDF() function used to map contour lines
library(ggplot2)
# library(gtable)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(scales)
# library(plotly) # doesn't work with layering grobs...
library(ggiraph) # this is good for interactivity, but cannot be converted to grobs for layout...
# library(ggpubr) # maybe this package's 'ggarrange' function can help.
library(reshape2)
# library(flextable) # for displaying tables when mouse hovers over mapped data points
# library(tidyverse)
# library(DT)
# library(shinybrowser) # get browser dimension & other info

# library(ggnewscale)
# Use devtools to install ggnewscale version 0.4.3 as there seems to be bug in the
# latest version -- this is annoying, but ggnewscale is vital for this map.
# See this solution from the package developer: https://github.com/eliocamp/ggnewscale/issues/45
# (I'm not convinced this helped! The issue seemed to solved by including '_new'
# as a suffix in arguments to the 'guides' function)
library(devtools)
ggnewscale_version <- '0.4.3'
ggnewscale_download_path <- paste0('eliocamp/ggnewscale@v', ggnewscale_version)
ggnewscale_available <- require(ggnewscale, quietly = TRUE)
if(!ggnewscale_available){
  install_github(ggnewscale_download_path)
  library(ggnewscale)
}else{
  if(packageVersion('ggnewscale') != ggnewscale_version){
    detach("package:ggnewscale", unload = TRUE)
    remove.packages('ggnewscale')
    install_github(ggnewscale_download_path)
    library(ggnewscale)
  }
}

# Directory info ----------------------------------------------------------

wd_orig <- getwd()
wd_base <- '~/Documents/Git Repos/CUPIDO-risk-map'
# setwd('~/Documents/Git Repos/CUPIDO-risk-map')

# Load data ---------------------------------------------------------------

source('functions.R', local = TRUE)

significantTrendsOnly <- TRUE

res <- '3x1'

# source('import_data.R') # function (get_data) to load/organise data
get_data(res = res, baseDirectory = wd_base, shinyDirectory = wd_orig, 
         sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
         sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly,
         roundShipTime = TRUE, indexGridCells = FALSE)

# Plotting parameters -----------------------------------------------------
set_plot_params(nc, DATA_sf, STATIONS_sf)

# Plots -------------------------------------------------------------------

# 1. background data

make_background_map <- function(background, group, backgroundOnly = TRUE, stationsOnly = FALSE, stationPopSize = FALSE, eco = FALSE, trendType = 'unspecified', ship_data = 'unspecified', pval_contours = FALSE, contourLineWidth = 0.25, ptStroke = 1){
  if(backgroundOnly & stationsOnly){
    warning('Cannot select both backgroundOnly and stationsOnly')
    return(NULL)}
  if(!backgroundOnly & !stationsOnly){
    warning('Either backgroundOnly or stationsOnly must be TRUE')
    return(NULL)}
  if(backgroundOnly){
    d <- get(paste(background, 'poly', sep = '_'))
    d <- switch(background,
                chl = subset(d, month == group),
                krill = subset(d, month == group),
                sst = subset(d, month == group),
                pH = subset(d, month == group),
                ship = switch(ship_data,
                              summary = subset(d, ship_class == group),
                              raw = subset(d, activity == group))
    )
    dat <- list(nc = nc, background = d, plastic = DATA_sf[0,],  stations = STATIONS_sf[0,], symbols = pltSymbols)
  }
  if(stationsOnly){
    d <- STATIONS_sf
    if(group %in% d$Type) d <- subset(d, Type == group)
    dat <- list(nc = nc, background = NULL, plastic = DATA_sf[0,],  stations = d, symbols = pltSymbols)
  }
  # dat <- list(nc = nc, background = d, plastic = DATA_sf[0,],  stations = STATIONS_sf[0,], pltSymbols = pltSymbols)
  mp <- make_plot(dat = dat, background = background, backgroundOnly = backgroundOnly, stationsOnly = stationsOnly, 
                  stationPopSize = stationPopSize, displayEcoregions = eco, plotSignificanceContours = pval_contours, contourLineWidth = contourLineWidth, ptStroke = ptStroke)
  # mp$plot
  mp
}

# Chlorophyll
background <- 'chl'
month <- 'all'
p_chl <- make_background_map(background, month)
p_chl$plot_complete

# Krill
background <- 'krill'
month <- 'all'
p_krill <- make_background_map(background, month)# + theme_nothing()
p_krill$plot_complete

# Shipping
background <- 'ship'
ship_data = 'raw'
ship_class <- 'all'
p_ship <- make_background_map(background, ship_class, eco = FALSE, ship_data = ship_data)# + theme_nothing()
p_ship$plot_complete

ship_classes <- levels(ship_poly$activity)
for(i in 1:length(ship_classes)){
  ship_class <- ship_classes[i]
  plt_name <- paste0('p_ship_', ship_class)
  p <- make_background_map(background, ship_class, eco = FALSE, ship_data = ship_data)# + theme_nothing()
  assign(plt_name, p)
}

npanels <- 6
rw <- c(0.9, 0.1)
ra <- 2
p_ship1 <- plot_grid(p_ship_fishing$plot, p_ship_fishing$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')), rel_widths = rw)
p_ship2 <- plot_grid(p_ship_tourism$plot, p_ship_tourism$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')), rel_widths = rw)
p_ship3 <- plot_grid(p_ship_supply$plot, p_ship_supply$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')), rel_widths = rw)
p_ship4 <- plot_grid(p_ship_research$plot, p_ship_research$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')), rel_widths = rw)
p_ship5 <- plot_grid(p_ship_other$plot, p_ship_other$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')), rel_widths = rw)
p_ship6 <- plot_grid(p_ship_all$plot, p_ship_all$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')), rel_widths = rw)

p_ship_combine <- plot_grid(p_ship1, p_ship2,
                   p_ship3, p_ship4,
                   p_ship5, p_ship6,
                   nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
                   labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9)  + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'), plot.background = element_rect(fill = 'white', colour = 'white'))

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.45 # scaling factors
sc_h <- 0.9
pw <- sc_w * A4_w
ph <- sc_h * A4_h

fileName <- paste0('ship traffic data_', res, '.png')
ggsave(fileName, plot = p_ship_combine, device = 'png', width = pw, height = ph, units = 'in')




# Sea surface temperature
background <- 'sst'
month <- 'all'
type <- 'trend'
pval_contours <- FALSE
contourLineWidth <- 0.25
p_sst <- make_background_map(background, month, trendType = type, pval_contours = pval_contours, contourLineWidth = contourLineWidth)# + theme_nothing()
p_sst$plot_complete

# ggsave('sst.png', plot = p_sst, device = 'png', width = 6, height = 4.15, units = 'in')

# pH
background <- 'pH'
month <- 'all'
type <- 'trend'
pval_contours <- FALSE
contourLineWidth <- 0.25
p_pH <- make_background_map(background, month, trendType = type, pval_contours = pval_contours, contourLineWidth = contourLineWidth)# + theme_nothing()
p_pH$plot_complete

# ggsave('pH.png', plot = p_pH, device = 'png', width = 6, height = 4.15, units = 'in')

# stations
background <- 'none'
group <- 'all'
p_stations <- make_background_map(background, group, backgroundOnly = FALSE, stationsOnly = TRUE, stationPopSize = TRUE)
p_stations$plot_complete


# combine all plots
npanels <- 6
rw <- c(0.8, 0.2)
p_chl_ <- plot_grid(p_chl$plot, p_chl$legend, rel_widths = rw) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_krill_ <- plot_grid(p_krill$plot, p_krill$legend, rel_widths = rw) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_sst_ <- plot_grid(p_sst$plot, p_sst$legend, rel_widths = rw) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_pH_ <- plot_grid(p_pH$plot, p_pH$legend, rel_widths = rw) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_ship_ <- plot_grid(p_ship$plot, p_ship$legend, rel_widths = rw) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_stations_ <- plot_grid(p_stations$plot, p_stations$legend, rel_widths = rw) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))

p_all <- plot_grid(p_chl_, p_krill_,
                   p_sst_, p_pH_,
                   p_ship_, p_stations_,
                   nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
                   labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9)  + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'), plot.background = element_rect(fill = 'white', colour = 'white'))

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.5 # scaling factors
sc_h <- 1.2
pw <- sc_w * A4_w
ph <- sc_h * A4_h

fileName <- paste0('environ background data_', res, '.png')
ggsave(fileName, plot = p_all, device = 'png', width = pw, height = ph, units = 'in')


# 2. plastic data

make_plastic_map <- function(plastic, stations, ptSize, Symbols, background = 'none', group = NA, eco = FALSE){
  if(background == 'none'){
    dat <- list(nc = nc, background = NULL, plastic = plastic, stations = stations, symbols = Symbols)
  } else return(NULL)
  mp <- make_plot(dat = dat, background = background, plasticOnly = TRUE, displayEcoregions = eco, ptSize = ptSize)
  mp
}


background <- 'none'
plastic <- DATA_sf
# Get line breaks in legend
sg <- levels(plastic$SampleType_grouped)
ll <- strsplit(sg, '/')
sgnew <- sapply(ll, function(z) if(length(z) == 1) return(z) else paste(z, collapse = '/\n'))
plastic$SampleType_grouped <- factor(plastic$SampleType_grouped, levels = sg, labels = sgnew)
Symbols <- pltSymbols
Symbols$Type <- as.character(Symbols$Type)
for(i in 1:nrow(Symbols)){
  j <- Symbols$Type[i]
  if(j %in% sg) Symbols$Type[i] <- sgnew[which(sg == j)]}
Symbols$Type <- as.factor(Symbols$Type)
stations <- NULL
ptSize <- 4
# Trim off some points at the map boundaries -- this is tricky, and should maybe be done later...

# i <- plastic$Source == 'Suaria (2020)'
# lat <- as.numeric(sapply(strsplit(plastic$Coordinates, ' '), function(z) z[3]))
# lat_ <- lat[!is.na(lat)]
# lat_ <- sort(lat_)
# lm <- lat_[1:2]
# plastic <- plastic[!{lat %in% lm},]

# lat_start <- as.numeric(sapply(strsplit(plastic$`Coordinates (start)`, ' '), function(z) z[3]))
# lat_end <- as.numeric(sapply(strsplit(plastic$`Coordinates (end)`, ' '), function(z) z[3]))
# lat_start_ <- lat_start[!is.na(lat_start)]
# lat_end_ <- lat_end[!is.na(lat_end)]
# lat_ <- cbind(lat_start_, lat_end_)
# lat_ <- apply(lat_, 1, min)
# lat_ <- sort(lat_)
# lm <- lat_[1:2]
# plastic <- plastic[!{{lat_start %in% lm} | {lat_end %in% lm}},]


p_pla <- make_plastic_map(plastic, stations, ptSize, Symbols)

p_plastic <- p_pla$plot_complete + 
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.1 # scaling factors
sc_h <- 0.8
pw <- sc_w * A4_w
ph <- sc_h * pw
# ph <- sc_h * A4_h

ggsave('plastic samples_trim.png', plot = p_plastic, device = 'png', width = pw, height = ph, units = 'in')




# background <- 'none'
# plastic <- DATA_sf
# stations <- STATIONS_sf
# ptSize <- 4
# p_pla <- make_plastic_map(plastic, stations, ptSize) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
# p_pla
# 
# ggsave('plastic samples and stations.png', plot = p_pla, device = 'png', width = 10, height = 8.3, units = 'in')










