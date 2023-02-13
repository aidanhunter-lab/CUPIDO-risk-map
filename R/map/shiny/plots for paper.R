# Individual map plots for paper

# Packages ----------------------------------------------------------------

library(shiny)
library(sp)
library(sf)
library(mapdata)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(scales)
# library(plotly) # doesn't work with layering grobs...
library(ggiraph) # this is good for interactivity, but cannot be converted to grobs for layout...
# library(ggpubr) # maybe this package's 'ggarrange' function can help.
library(reshape2)
library(flextable) # for displaying tables when mouse hovers over mapped data points
# library(tidyverse)
library(DT)
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

significantTrendsOnly <- FALSE

# source('import_data.R') # function (get_data) to load/organise data
get_data(res = '3x1', baseDirectory = wd_base, shinyDirectory = wd_orig, 
         sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
         sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly)


# Plotting parameters -----------------------------------------------------
# Plot symbols - store in separate data frame
# Research stations use (mostly) unfilled symbols
pltShapes <- c(1, 0, 5, 2, 6, 19) 
# stationTypes <- unique(STATIONS$Type)
stationTypes <- levels(STATIONS_sf$Type)
# stationTypes <- levels(STATIONS$Type)
nstationTypes <- length(stationTypes)
pltSymbols <- data.frame(Class = rep('ResearchStation', nstationTypes), Type = stationTypes, symbol = pltShapes[1:nstationTypes])
# Plastic samples use filled plot symbols (21:25) that differ according to sample type.
pltShapes <- 21:25
sampleTypes <- levels(DATA_sf$SampleType_grouped)
# sampleTypes <- levels(DATA$SampleType_grouped)
# sampleTypes <- levels(DATA$SampleType)
# sampleTypes <- levels(DATA$'Sample type')
nsampleTypes <- length(sampleTypes)
pltSymbols <- rbind(pltSymbols,
                    data.frame(Class = rep('PlasticSample', nsampleTypes), Type = sampleTypes, symbol = pltShapes[1:nsampleTypes]))
pltSymbols$Type <- factor(pltSymbols$Type, levels = pltSymbols$Type)

# Plot colours
# Use a qualitative palette for plastic sources -- the default is limited to 12 colours
ncol <- length(unique(DATA_sf$Source))
pltColours <- brewer.pal(12, 'Paired') # the Set3 and Paired palettes has a maximum of 12 colours
# plot(1:12,1:12, col=pltColours, pch = 19)
if(ncol < 13){
  pltColours <- pltColours[c(seq(2, 12, 2), seq(1, 11, 2))]
}else{
  # make new colours as gradient between existing Paired palette colours
  ne <- ceiling({ncol - 12} / 6)
  mc <- matrix(pltColours, 2)
  mc <- sapply(1:6, function(z){
    cf <- colorRampPalette(c(mc[1,z], mc[2,z]))
    cf <- cf(ne+2)
    cf[2:{ne+1}]
  })
  mc <- matrix(mc, ne, 6)
  pltColours <- c(
    pltColours[seq(2, 12, 2)],
    pltColours[seq(1, 11, 2)],
    as.vector(t(mc[seq(ne, 1, -1),]))
  )
}

# pltColours <- brewer.pal(12, 'Paired') # the Set3 and Paired palettes has a maximum of 12 colours
# pltColours <- pltColours[c(seq(2, 12, 2), seq(1, 11, 2))]

sources <- levels(DATA_sf$Source)
# sources <- levels(DATA$Source)
nsources <- length(sources)
pltColours <- data.frame(Source = sources, colour = pltColours[1:nsources])

du <- unique(data.frame(Source = DATA_sf$Source, URL = DATA_sf$URL))
# du <- unique(DATA[c('Source', 'URL')])
pltColours <- merge(pltColours, du, by = 'Source')

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

BBox <- bbox_map

aspectRatio = unname(diff(BBox[c(1,3)]) / diff(BBox[c(2,4)]))

linebreaks <- function(n){HTML(strrep(br(), n))} # convenience function



# Plots -------------------------------------------------------------------

# 1. background data

make_background_map <- function(background, group, eco = FALSE, trendType = 'unspecified', ship_data = 'unspecified', pval_contours = FALSE, contourLineWidth = 0.25){
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
  dat <- list(nc = nc, background = d, plastic = DATA_sf[0,],  stations = STATIONS_sf[0,], pltSymbols = pltSymbols)
  mp <- make_plot(dat = dat, background = background, backgroundOnly = TRUE, 
                  displayEcoregions = eco, plotSignificanceContours = pval_contours, contourLineWidth = contourLineWidth)
  mp$plot
}

# Chlorophyll
background <- 'chl'
month <- 'all'
p_chl <- make_background_map(background, month)
p_chl


# Krill
background <- 'krill'
month <- 'all'
p_krill <- make_background_map(background, month)# + theme_nothing()
p_krill

# Shipping
background <- 'ship'
ship_data = 'raw'
ship_class <- 'all'
p_ship <- make_background_map(background, ship_class, eco = FALSE, ship_data = ship_data)# + theme_nothing()
p_ship

ship_classes <- levels(ship_poly$activity)
for(i in 1:length(ship_classes)){
  ship_class <- ship_classes[i]
  plt_name <- paste0('p_ship_', ship_class)
  p <- make_background_map(background, ship_class, eco = FALSE, ship_data = ship_data)# + theme_nothing()
  assign(plt_name, p)
}

ncl <- 2
nrw <- 3
p_ship_combine <- grid.arrange(p_ship_all, p_ship_fishing, p_ship_tourism,
                               p_ship_supply, p_ship_research, p_ship_other,
                               nrow = nrw, ncol = ncl)
pw <- 6 * ncl
ph <- 4.15 * nrw
ggsave('ship traffic data.png', plot = p_ship_combine, device = 'png', width = pw, height = ph, units = 'in')



# Sea surface temperature
background <- 'sst'
month <- 'all'
type <- 'trend'
pval_contours <- TRUE
contourLineWidth <- 0.25
p_sst <- make_background_map(background, month, trendType = type, pval_contours = pval_contours, contourLineWidth = contourLineWidth)# + theme_nothing()
p_sst

ggsave('sst.png', plot = p_sst, device = 'png', width = 6, height = 4.15, units = 'in')

# pH
background <- 'pH'
month <- 'all'
type <- 'trend'
pval_contours <- TRUE
contourLineWidth <- 0.25
p_pH <- make_background_map(background, month, trendType = type, pval_contours = pval_contours, contourLineWidth = contourLineWidth)# + theme_nothing()
p_pH

ggsave('pH.png', plot = p_pH, device = 'png', width = 6, height = 4.15, units = 'in')



npanels <- 5
nrw <- ceiling(npanels ^ 0.5)
ncl <- ceiling(npanels / nrw)
p_all <- plot_grid(p_chl, p_krill, p_ship, p_sst, p_pH,
                   nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
                   labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))

p_all


pw <- 5 * ncl
ph <- 4.15 * nrw

ggsave('environ background data.png', plot = p_all, device = 'png', width = pw, height = ph, units = 'in')


# 2. plastic data

make_plastic_map <- function(plastic, stations, ptSize, background = 'none', group = NA, eco = FALSE){
  if(background == 'none'){
    dat <- list(background = NULL, plastic = plastic, stations = stations, symbols = pltSymbols)
  } else return(NULL)
  mp <- make_plot(dat = dat, background = background, plasticOnly = TRUE, displayEcoregions = eco, ptSize = ptSize)
  mp$plot
}

background <- 'none'
plastic <- DATA_sf
stations <- STATIONS_sf
ptSize <- 4
p_pla <- make_plastic_map(plastic, stations, ptSize) + theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_pla

ggsave('plastic samples and stations.png', plot = p_pla, device = 'png', width = 10, height = 8.3, units = 'in')










