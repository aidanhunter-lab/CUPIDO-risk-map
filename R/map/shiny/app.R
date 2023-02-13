#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Shiny app for interactive map of Southern Ocean displaying plastic & environmental data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Packages ----------------------------------------------------------------

library(shiny)
library(sp)
library(sf)
library(mapdata) # for Antarctic coastline data used to estimate distance of facilities from the coast (I would prefer a higher-res shape file...)
library(maptools) # for the ContourLines2SLDF() function used to map contour lines
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

# It will be nice to able to select between trend data and anomalies from within
# the Shiny app. This will involve returning both data sets from the get_data
# function, then including a selection option in the Shiny dashboard that chooses
# between them. This is not critical to include, but would be a nice extension
# to what I've already produced

if(significantTrendsOnly){
  # If all insignificant trends are omitted (plotted in white) then we don't need
  # to keep the p-values
  sst_poly <- sst_poly[sst_poly$metric != 'p-value',]
  pH_poly <- pH_poly[pH_poly$metric != 'p-value',]
}


# Flextable options ------------------------------------------------------

set_flextable_defaults(
  font.size = 10,
  font.family = 'Roboto',
  theme_fun = 'theme_vanilla'
)

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


# Plotting function -------------------------------------------------------

# make_plot <- function(dat, background = 'none', displayEcoregions = FALSE, backgroundOnly = FALSE, plasticOnly = FALSE, ptSize = 6, legPtSize = 4, alpha = 0.6, polyLineWidth = 0.75, contourLineWidth = 0.75, legWidth = 0.3, mapAspectRatio = aspectRatio, pHType = 'unspecified', sstType = 'unspecified', shipSummaryDataOrRaw = 'unspecified', plotSignificanceContours = FALSE, significanceContours = NULL, sst_sigContours = significanceContours, pH_sigContours = significanceContours){
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Function to generate map plot, called from inside the server function after
#   # data have been filtered by user selected inputs.
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Inputs:
#   # dat               = list of all required (reactively filtered) data sets
#   # background        = character specifying which background layer to plot
#   # displayEcoregions = TRUE/FALSE: plot distinct ecoregions -- required for shipping data from McCarthy
#   # alpha             = point transparency
#   # Output:
#   # A ggplot object containing interactive geoms compatible with girafe
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   dat_background = dat$background
#   dat_plastic = dat$plastic
#   dat_stations = dat$stations
#   symbols = dat$symbols
#   pltColours = dat$pltColours
# 
#   anyPlastic <- nrow(dat_plastic) > 0
#   anyStations <- nrow(dat_stations) > 0
#   anyBackground <- background != 'none'
# 
#   #~~~~~~~~~~~~~~~~~
#   # Background layer
#   #~~~~~~~~~~~~~~~~~
#   switch(background,
# 
#          none = {
#            plt_background <- ggplot()
#          },
# 
#          krill = {
#            plt_background <-
#              ggplot() +
#              geom_sf(data = dat_background, aes(fill = colourgroup)) +
#              scale_fill_viridis_d(option = 'plasma',
#                                   name = bquote(atop(Krill, (number / m^2))))
#          },
# 
#          chl = {
#            plt_background <-
#              ggplot() +
#              geom_sf(data = dat_background, aes(fill = value)) +
#              scale_fill_viridis_c(option = 'viridis', trans = 'log10',
#                                   name = bquote(atop(Chlorophyll, (mg / m^3))))
#          },
# 
#          sst = {
# 
#            legLabel <- switch(sstType,
#                               unspecified = bquote(SST),
#                               trend = bquote(atop(SST ~ trend, (degree * C / year))),
#                               anomaly = bquote(atop(SST ~ anomaly, (degree * C))))
#            nv <- dat_background$value[dat_background$value < 0]
#            pv <- dat_background$value[dat_background$value > 0]
#            any_positive <- length(pv) > 0
#            any_negative <- length(nv) > 0
#            ncol_sst <- 11
#            ncol_sst_ <- 0.5 * {ncol_sst - 1}
#            colfunc <- colorRampPalette(c('dodgerblue4', 'white', 'red3'))
#            Cols <- colfunc(ncol_sst)
#            cval <- c(quantile(nv, seq(0, 1, length = ncol_sst_ + 1))[1:ncol_sst_], 0, quantile(pv, seq(0, 1, length = ncol_sst_ + 1))[2:{ncol_sst_ + 1}])
#              if(any_positive & any_negative){
#              cval <- {cval - min(cval)} / diff(range(cval))
#              plt_background <-
#                ggplot() +
#                geom_sf(data = dat_background, aes(fill = value)) +
#                scale_fill_gradientn(colours = Cols, values = cval,
#                                     name = legLabel)
#            }else{if(any_negative){
#              cval <- {cval[1:{ncol_sst_+1}] - min(cval[1:{ncol_sst_+1}])} / diff(range(cval[1:{ncol_sst_+1}]))
#              cval <- {cval - min(cval)} / diff(range(cval))
#              plt_background <-
#                ggplot() +
#                geom_sf(data = dat_background, aes(fill = value)) +
#                scale_fill_gradientn(colours = Cols[1:{ncol_sst_+1}], values = cval,
#                                     name = legLabel)
#            }else{if(any_positive){
#              cval <- {cval[{ncol_sst_+1}:ncol_sst] - min(cval[{ncol_sst_+1}:ncol_sst])} / diff(range(cval[{ncol_sst_+1}:ncol_sst]))
#              cval <- {cval - min(cval)} / diff(range(cval))
#              plt_background <-
#                ggplot() +
#                geom_sf(data = dat_background, aes(fill = value)) +
#                scale_fill_gradientn(colours = Cols[{ncol_sst_+1}:ncol_sst], values = cval,
#                                     name = legLabel)
#            }}}
#          },
# 
#          pH = {
# 
#            legLabel <- switch(pHType,
#                               unspecified = bquote(pH),
#                               trend = bquote(atop(pH ~ trend, (1/year))),
#                               anomaly = bquote(pH ~ anomaly))
#            nv <- dat_background$value[dat_background$value < 0]
#            pv <- dat_background$value[dat_background$value > 0]
#            any_positive <- length(pv) > 0
#            any_negative <- length(nv) > 0
#            ncol_pH <- 11
#            ncol_pH_ <- 0.5 * {ncol_pH - 1}
#            colfunc <- colorRampPalette(c('darkorange2', 'white', 'darkorchid4'))
#            Cols <- colfunc(ncol_pH)
#            cval <- c(quantile(nv, seq(0, 1, length = ncol_pH_ + 1))[1:ncol_pH_], 0, quantile(pv, seq(0, 1, length = ncol_pH_ + 1))[2:{ncol_pH_ + 1}])
#            # Colour scale selection depends on range of values
#            if(any_positive & any_negative){
#              cval <- {cval - min(cval)} / diff(range(cval))
#              plt_background <-
#                ggplot() +
#                geom_sf(data = dat_background, aes(fill = value)) +
#                scale_fill_gradientn(colours = Cols, values = cval,
#                                     name = legLabel)
#            }else{if(any_negative){
#              cval <- {cval[1:{ncol_pH_+1}] - min(cval[1:{ncol_pH_+1}])} / diff(range(cval[1:{ncol_pH_+1}]))
#              cval <- {cval - min(cval)} / diff(range(cval))
#              plt_background <-
#                ggplot() +
#                geom_sf(data = dat_background, aes(fill = value)) +
#                scale_fill_gradientn(colours = Cols[1:{ncol_pH_+1}], values = cval,
#                                    name = legLabel)
#            }else{if(any_positive){
#              cval <- {cval[{ncol_pH_+1}:ncol_pH] - min(cval[{ncol_pH_+1}:ncol_pH])} / diff(range(cval[{ncol_pH_+1}:ncol_pH]))
#              cval <- {cval - min(cval)} / diff(range(cval))
#              plt_background <-
#                ggplot() +
#                geom_sf(data = dat_background, aes(fill = value)) +
#                scale_fill_gradientn(colours = Cols[{ncol_pH_+1}:ncol_pH], values = cval,
#                                     name = legLabel)
#            }}}
#          },
# 
#          ship = {
# 
#            plt_background <-switch(shipSummaryDataOrRaw,
# 
#                   summary = {
#                     shipType <- unique(dat_background$ship_class)
#                     ShipType <- paste0(toupper(substr(shipType,1,1)), substr(shipType,2,nchar(shipType)))
#                     leg_lab <- bquote(.(paste0('Ship time:', '\n', shipType, ' vessels', '\n', '(days)')))
#                     vesselsPresent <- !is.na(dat_background$total_time)
#                     plt <-
#                       ggplot() +
#                       geom_sf(data = dat_background[vesselsPresent,], aes(fill = total_time), colour = 'black', linewidth = polyLineWidth) +
#                       scale_fill_viridis_c(option = 'mako', trans = 'log10', direction = -1,
#                                            name = leg_lab) +
#                       geom_sf(data = dat_background[!vesselsPresent,], fill = 'white', colour = 'black', linewidth = polyLineWidth)
#                     plt
#                   },
# 
#                   raw = {
# 
#                     shipActivity <- unique(dat_background$activity)
#                     leg_lab <- bquote(.(paste0('Ship time:', '\n', shipActivity, ' vessels', '\n', '(days/yr)')))
#                     plt <-
#                       ggplot() +
#                       geom_sf(data = dat_background, aes(fill = value), colour = 'black', linewidth = polyLineWidth) +
#                       scale_fill_viridis_c(option = 'mako', trans = 'log10', direction = -1,
#                                            name = leg_lab) # +
#                     # geom_sf(data = dat_background[!vesselsPresent,], fill = 'white', colour = 'black', linewidth = polyLineWidth)
#                     plt
# 
#                   },
# 
#                   unspecified = {
#                     plt <- ggplot()
#                     plt
#                   }
#                   )
# 
# 
# 
#            # shipType <- unique(dat_background$ship_class)
#            # ShipType <- paste0(toupper(substr(shipType,1,1)), substr(shipType,2,nchar(shipType)))
#            # leg_lab <- bquote(.(paste0('Ship time:', '\n', shipType, ' vessels', '\n', '(days)')))
#            # vesselsPresent <- !is.na(dat_background$total_time)
#            #
#            # plt_background <-
#            #   ggplot() +
#            #   geom_sf(data = dat_background[vesselsPresent,], aes(fill = total_time), colour = 'black', linewidth = polyLineWidth) +
#            #   scale_fill_viridis_c(option = 'mako', trans = 'log10', direction = -1,
#            #                        name = leg_lab) +
#            #   geom_sf(data = dat_background[!vesselsPresent,], fill = 'white', colour = 'black', linewidth = polyLineWidth)
#            #
# 
# 
#          }
#   )
# 
#   if(plotSignificanceContours){
#     if(background %in% c('sst', 'pH')){
#       # Create contours in standard lat-lon coordinates
#       dat_ <- st_transform(dat_background, crs_world)
#       centroids <- as.data.frame(st_coordinates(st_centroid(dat_$geometry))) # get coordinates of polygon centroids
#       # Account for rounding error after coordinate transform
#       centroids <- round(centroids, 4)
#       centroids$pval <- dat_$pval
#       centroids$plevel <- dat_$plevel
#       xgrid <- sort(unique(centroids$X))
#       ygrid <- sort(unique(centroids$Y))
#       nx <- length(xgrid)
#       ny <- length(ygrid)
#       xygrid <- expand.grid(X = xgrid, Y = ygrid)
#       centroids <- merge(xygrid, centroids, all.x = TRUE)
#       centroids <- centroids[order(centroids$Y, centroids$X),]
#       zgrid <- matrix(centroids$pval, nx, ny)
# 
#       if(background == 'pH'){
#         if(is.null(pH_sigContours)){
#           log_pval <- log10(dat_$pval)
#           seq_pval <- seq(min(log_pval), max(log_pval), length = 5)
#           seq_pval <- pretty(seq_pval)
#           significanceContours <- 10 ^ seq_pval[2:4]
#         }
#       }
#       if(background == 'sst'){
#         if(is.null(sst_sigContours)){
#           significanceContours <- c(0.05, 0.25, 0.5)
#         }
#       }
#       cL <- contourLines(xgrid, ygrid, zgrid, levels = significanceContours)
#       # Transform contour lines to spatial object
#       library(maptools)
#       crs_world_ <- paste0('+init=epsg:', crs_world)
#       cLdf <- ContourLines2SLDF(cL, proj4string = CRS(crs_world_))
#       # Transform to sf
#       cLdf_sf <- st_as_sf(cLdf)
#       # Tranform to stereogrpahic coordinates
#       cLdf_sf <- st_transform(cLdf_sf, crs_use)
#       names(cLdf_sf)[names(cLdf_sf) == 'level'] <- 'p-value'
#       plt_background <- plt_background +
#         geom_sf(data = cLdf_sf, aes(linetype = `p-value`), linewidth = contourLineWidth)#, inherit.aes = FALSE)# +
#     }
#   }
# 
#   plt_background <- plt_background +
#     theme(
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       panel.background = element_rect(fill = 'white', colour = 'white'),
#       plot.background = element_rect(fill = 'white', colour = 'white')
#     )
# 
#   leg_background <- get_legend(plt_background)
# 
#   #~~~~~~~~~~~~~~~~~~~~~~
#   # Main map -- no legend
#   #~~~~~~~~~~~~~~~~~~~~~~
# 
#   plt_map <- plt_background + guides(fill = 'none') + guides(linetype = 'none')# + theme(legend.position = 'right')
# 
#   # Ecoregions
#   if(displayEcoregions){
#     plt_map <- plt_map +
#       geom_sf(data = eco,
#               linewidth = polyLineWidth, colour = 'black', fill = alpha('white', 0))
#   }
# 
#   # Coastline
#   plt_map <- plt_map +
#     new_scale('fill') +
#     geom_sf(data = nc,
#             aes(fill = surface),
#             show.legend = FALSE) +
#     scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
# 
#   plt_map <- plt_map + guides(fill_new = 'none')
# 
# 
#   if(backgroundOnly){
#     return(
#       list(plot = ggdraw(plot_grid(plt_map, leg_background, ncol = 2, rel_widths = c(0.8, 0.2))))
#       )
#   }
# 
#   if(anyStations | anyPlastic){
#     plt_map <- plt_map +
#       scale_shape_manual(values = setNames(symbols$symbol, symbols$Type))
#   }
# 
#   # Research stations
#   if(anyStations){
#     plt_map <- plt_map +
#       new_scale_colour() +
#       geom_sf_interactive(data = dat_stations,
#                           aes(shape = Type, colour = Seasonality, data_id = Record_ID, tooltip = tooltip),
#                           alpha = alpha, size = ptSize, stroke = 1, show.legend = FALSE) +
#       scale_colour_manual(values = c('forestgreen','firebrick'))
#   }
# 
#   # Plastic samples
#   if(anyPlastic){
#     plt_map <- plt_map +
#       new_scale_fill() +
#       geom_sf_interactive(data = dat_plastic,
#                           aes(fill = Source, shape = SampleType_grouped, data_id = data_id, tooltip = tooltip),
#                           alpha = alpha, size = ptSize, show.legend = FALSE) +
#       scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
#   }
# 
#   plt_map <- plt_map + guides(fill_new_new = 'none')
# 
#   plt_map <- plt_map +
#     coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax']))# +
# 
#   #~~~~~~~
#   # Legend
#   #~~~~~~~
# 
#   # Research stations
#   if(anyStations){
#     symbols_ <- symbols[symbols$Class == 'ResearchStation',]
#     plt_stations <-
#       ggplot() +
#       geom_sf(data = dat_stations,
#               aes(shape = Type, colour = Seasonality),
#               alpha = 1, size = ptSize, stroke = 1) +
#       scale_colour_manual(values = c('forestgreen','firebrick')) +
#       scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
#       theme(legend.key = element_blank()) +
#       guides(shape = guide_legend(title = 'Facility', override.aes = list(size = legPtSize)),
#              colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
#     leg_stations <- get_legend(plt_stations)
#     # get legend dimensions
#     lw = (leg_stations$widths)
#     lh = (leg_stations$heights)
#     lw = lw[grepl('cm', lw)]
#     lh = lh[grepl('cm', lh)]
#     leg_stations_width = sum(as.numeric(gsub('cm', '', lw)))
#     leg_stations_height = sum(as.numeric(gsub('cm', '', lh)))
#   }else{
#     leg_stations <- NULL
#     leg_stations_width = 0.5
#     leg_stations_height = 0.5
#   }
# 
#   if(anyPlastic){
#     symbols_ = symbols[symbols$Class == 'PlasticSample',]
#     plt_plastic_samples <-
#       ggplot() +
#       geom_sf(data = dat_plastic,
#               aes(fill = Source, shape = SampleType_grouped),
#               alpha = 1, size = ptSize) +
#       guides(
#         shape = guide_legend(
#           title = 'Sample type',
#           override.aes = list(size = legPtSize)),
#         fill = guide_legend_interactive(
#           override.aes = list(
#             shape = 21,
#             size = legPtSize))
#       ) +
#       scale_fill_manual_interactive(
#         values = setNames(pltColours$colour, pltColours$Source),
#         data_id = setNames(pltColours$Source, pltColours$Source),
#         tooltip = setNames(pltColours$URL, pltColours$Source),
#         labels = function(Source) {
#           lapply(Source, function(s) {
#             u <- unique(dat_plastic$URL[dat_plastic$Source == s])
#             label_interactive(
#               as.character(s),
#               data_id = as.character(s),
#               onclick = paste0("window.open(`", as.character(u), "`);"),
#               tooltip = as.character(u)
#             )
#           })
#         }
#       ) +
#       scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
#       theme(legend.key = element_blank())
#     leg_plastic <- get_legend(plt_plastic_samples)
#     # get legend dimensions
#     lw = (leg_plastic$widths)
#     lh = (leg_plastic$heights)
#     lw = lw[grepl('cm', lw)]
#     lh = lh[grepl('cm', lh)]
#     leg_plastic_width = sum(as.numeric(gsub('cm', '', lw)))
#     leg_plastic_height = sum(as.numeric(gsub('cm', '', lh)))
#   }else{
#     leg_plastic <- NULL
#     leg_plastic_width = 0.5
#     leg_plastic_height = 0.5
#   }
# 
#   # Background
#   if(anyBackground){
#     # leg_background <- get_legend(plt_background)
#     # get legend dimensions
#     lw = (leg_background$widths)
#     lh = (leg_background$heights)
#     lw = lw[grepl('cm', lw)]
#     lh = lh[grepl('cm', lh)]
#     leg_background_width = sum(as.numeric(gsub('cm', '', lw)))
#     leg_background_height = sum(as.numeric(gsub('cm', '', lh)))
#   }else{
#     leg_background <- NULL
#     leg_background_width = 0.5
#     leg_background_height = 0.5
#   }
# 
#   if(plasticOnly){
#     return(
#       list(plot = ggdraw(plot_grid(plt_map, leg_plastic, leg_stations, ncol = 3, rel_widths = c(0.7, 0.15, 0.15))))
#     )
#   }
# 
#   # Find size (cm) of combined legend -- include spacings
#   leg_width <- {2 * 1} + 3 * max(c(leg_plastic_width, leg_stations_width, leg_background_width))
#   leg_height <- {3 * 2} + max(c(leg_plastic_height, leg_stations_height, leg_background_height))
#   # Size of map & complete plot
#   tot_width = leg_width / legWidth
#   map_width = tot_width - leg_width
#   map_height <- map_width / mapAspectRatio
#   tot_height <- max(map_height, leg_height)
# 
#   # Combine all map components
#   leg_complete <- plot_grid(leg_background, leg_plastic, leg_stations,
#                             ncol = 3)
# 
#   plt_complete <- ggdraw(
#     plot_grid(
#       plt_map,
#       leg_complete,
#       ncol = 2, rel_widths = c(1-legWidth, legWidth), align = 'h', axis = 't'))
# 
#   return(
#     list(
#       plot = plt_complete, width = tot_width, height = tot_height
#       )
#     )
# }

# Define UI for plastic data map app ----

# Define lists of names for background data
backgroundData_choiceNames <- list(
  'None',
  'Krill (Jan--Mar)',
  'Krill (Jan)',
  'Krill (Feb)',
  'Krill (Mar)',
  'Chlorophyll (Jan--Mar)',
  'Chlorophyll (Jan)',
  'Chlorophyll (Feb)',
  'Chlorophyll (Mar)',
  paste('SST', sstType, '(Jan--Mar)', sep = ' '),
  paste('SST', sstType, '(Jan)', sep = ' '),
  paste('SST', sstType, '(Feb)', sep = ' '),
  paste('SST', sstType, '(Mar)', sep = ' '),
  paste('pH', pHType, '(Jan--Mar)', sep = ' '),
  paste('pH', pHType, '(Jan)', sep = ' '),
  paste('pH', pHType, '(Feb)', sep = ' '),
  paste('pH', pHType, '(Mar)', sep = ' ')
)

backgroundData_choiceValues = list(
  'none',
  'krill_all',
  'krill_1',
  'krill_2',
  'krill_3',
  'chl_all',
  'chl_1',
  'chl_2',
  'chl_3',
  'sst_all',
  'sst_1',
  'sst_2',
  'sst_3',
  'pH_all',
  'pH_1',
  'pH_2',
  'pH_3'
)

# These shipSummaryDataOrRaw switches have become redundant, but keep them in case
# the code is modified in future -- it may be useful again...
backgroundData_choiceNames<- switch(shipSummaryDataOrRaw,
                                    summary = c(backgroundData_choiceNames, c('Shipping (all)',
                                                                              'Shipping (fishing)',
                                                                              'Shipping (tourism)',
                                                                              'Shipping (supply)',
                                                                              'Shipping (research)',
                                                                              'Shipping (other)')),
                                    raw = c(backgroundData_choiceNames, c('Shipping (all)',
                                                                          'Shipping (fishing)',
                                                                          'Shipping (tourism)',
                                                                          'Shipping (supply)',
                                                                          'Shipping (research)',
                                                                          'Shipping (other)'))
)

backgroundData_choiceValues <- switch(shipSummaryDataOrRaw,
                                      summary = c(backgroundData_choiceValues, c('ship_all',
                                                                                 'ship_fishing',
                                                                                 'ship_tourism',
                                                                                 'ship_supply',
                                                                                 'ship_research',
                                                                                 'ship_other')),
                                      raw = c(backgroundData_choiceValues, c('ship_all',
                                                                             'ship_fishing',
                                                                             'ship_tourism',
                                                                             'ship_supply',
                                                                             'ship_research',
                                                                             'ship_other'))
)

ui <- fluidPage(
  
  # App title
  titlePanel('Mapping Southern Ocean Plastic Data'),
  
  fluidRow(
    column(width = 3,
           wellPanel(
             # Input: year range
             sliderInput('YearRange', 'Years:',
                         min = min(DATA_sf$Year, na.rm = TRUE) , max = max(DATA_sf$Year, na.rm = TRUE),
                         value = range(DATA_sf$Year, na.rm = TRUE), step = 1, sep = ''),
             # sliderInput('YearRange', 'Years:',
             #             min = min(DATA_sf$Year, na.rm = TRUE) , max = max(DATA_sf$Year, na.rm = TRUE),
             #             value = range(DATA_long$Year, na.rm = TRUE), step = 1, sep = ''),
             
             # Input: measurement variable
             selectInput('Variable', 'Measurement:',
                         c('Concentration (pieces/m3)' = 'concentration',
                           'Mass concentration (g/m3)' = 'mass concentration',
                           'Density (pieces/km2)' = 'density',
                           'Mass density (g/km2)' = 'mass density',
                           'Presence/absence' = 'presence/absence'),
                         multiple = TRUE, selected = c('concentration', 'mass concentration', 'density', 'mass density', 'presence/absence')),
             
             # Input: plastic type
             selectInput('PlasticForm_grouped', 'Plastic form:',
                         c('Fragment' = 'fragment',
                           'Fibre' = 'fibre',
                           'Film' = 'film',
                           'Other/unspecified' = 'other/unspecified'),
                         multiple = TRUE, selected = c('fragment', 'fibre', 'film', 'other/unspecified')),
             
             # Input: plastic scale
             selectInput('LitterScale', 'Plastic scale:',
                         c('Micro' = 'micro',
                           'Micro & meso' = 'micro and meso',
                           'Meso' = 'meso',
                           'Macro' = 'macro',
                           'Unspecified' = 'unspecified'),
                         multiple = TRUE, selected = c('micro', 'micro and meso', 'meso', 'macro', 'unspecified')),
             
             # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment)
             selectInput('SampleType_grouped', 'Sample type:',
                         c('Seawater' = 'seawater',
                           'Fresh/waste water' = 'freshwater/wastewater',
                           'Sediment' = 'sediment',
                           'Ice' = 'ice',
                           'Air' = 'air'),
                         multiple = TRUE, selected = c('seawater', 'freshwater/wastewater', 'sediment', 'ice', 'air')),

             # # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment)
             # selectInput('SampleType', 'Sample type:',
             #             c('Sea surface' = 'marine surface',
             #               'Sea subsurface' = 'marine subsurface',
             #               'Freshwater' = 'freshwater',
             #               'Wastewater' = 'wastewater',
             #               'Sediment' = 'sediment'),
             #             multiple = TRUE, selected = c('marine surface','marine subsurface','freshwater','wastewater','sediment')),
             
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

             # Input: display contours of linear model p-values
             checkboxInput('DisplaySignificanceContours','Display linear trend p-values',
                           value = FALSE),
             
             # Input: background layers
             radioButtons('background', 'Background data:',
                          choiceNames = backgroundData_choiceNames,
                          choiceValues = backgroundData_choiceValues,
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
    column(width = 9, #offset = 0, style='padding:0px;',
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
           plotOutput('blank', width = '100%', height = '10px'),
           ggiraphOutput('plt')
    ),
    
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

# Define server logic ----
server <- function(input, output, session) {
  
  # Updating some Shiny inputs takes a long time because plotting the plastic
  # data is costly. Use 'debounce' function to delay updates, preventing each
  # individual change being rendered provided input changes are made quickly.
  listInputs <- debounce({
    reactive({
      list(
        YearRange = input$YearRange,
        Variable = input$Variable,
        PlasticForm = input$PlasticForm_grouped,
        LitterScale = input$LitterScale,
        SampleType = input$SampleType_grouped,
        StationType = input$StationType
      )
    })
  }, 2000)
  # Filter data according to Shiny inputs...
  # Research station data
  filtered_station_data <- reactive({
    x <- listInputs()
    return(
      subset(STATIONS_sf, Type %in% x$StationType)
    )
  })
  # Plastic data
  filtered_plastic_data <-reactive({
    x <- listInputs()
    y <- subset(DATA_long,
                x$YearRange[1] <= Year & Year <= x$YearRange[2] &
                  Variable %in% x$Variable &
                  PlasticForm_grouped %in% x$PlasticForm &
                  LitterScale %in% x$LitterScale &
                  SampleType_grouped %in% x$SampleType
    )
    return(
      list(data = y, samples = unique(y$data_id))
    )
  })
  transformed_plastic_data <- reactive({
    fp <- filtered_plastic_data()
    return(
      subset(DATA_sf, data_id %in% fp$samples)
    )
  })
  
  
  
  # Get background type
  which_background <- reactive({
      backgroundOptions <- c('none', 'krill', 'chl', 'sst', 'pH', 'ship')
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
                             sst = sst_poly,
                             pH = pH_poly,
                             ship = ship_poly
    ) # get background data
    if(background != 'none'){
      # filter by month (chl/krill/sst) or ship type
      if(background %in% c('chl', 'krill', 'sst', 'pH')){
        x <- strsplit(input$background, '_')[[1]]
        if(length(x) == 2){
          m <- x[2] # get the chosen month(s)
          background_dat <- subset(background_dat, month == m) # filter by month
        } else warning('Background data cannot be filtered by month: check that R Shiny ui input options match the data.')
      }else{
        
        if(background == 'ship'){
          background_dat <- switch(shipSummaryDataOrRaw,
                                   summary = {
                                     x <- strsplit(input$background, '_')[[1]]
                                     if(x[2] %in% background_dat$ship_class){
                                       background_dat <- subset(background_dat, ship_class == x[2])
                                     } else{warning('Background data cannot be filtered by ship type: check that R Shiny input options match the data.')}
                                   },
                                   raw = {
                                     x <- strsplit(input$background, '_')[[1]]
                                     if(x[2] %in% background_dat$activity){
                                       background_dat <- subset(background_dat, activity == x[2])
                                     } else{warning('Background ship data cannot be filtered by activity: check that R Shiny input options match the data.')}
                                   })
        }
      }
    }
    return(background_dat)
  })
  
  # Display ecoregions
  display_ecoregions <-reactive(
    return(input$DisplayEcoregions)
  )
  
  # Display p-values of linear trends
  display_trend_pvals <-reactive(
    return(input$DisplaySignificanceContours)
  )
  
  
  
  # Filter plot symbols based on sample & station types
  Symbols <- reactive({
    x <- listInputs()
      return(
        subset(pltSymbols, Type %in% c(x$SampleType, x$StationType))
      )
  })
  
  # Group all data in a named list
  listData <- reactive({
    plastic <- transformed_plastic_data()
    stations <- filtered_station_data()
    background <- filtered_background_dat()
    symbols <- Symbols()
    return(
      list(nc = nc, background = background, stations = stations, plastic = plastic, symbols = symbols) #, pltColours = pltColours)
    )
  })
  
  # Create blank plot to specify sizing/aspect ratio of main plot
  output$blank = renderPlot({
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point() + theme_void()
  })
  blankwidth = reactive({
    # this is the magic that makes it work
    bw = session$clientData$output_blank_width
    return(5/4 * bw) # the sizing on girafe functions is weird, this scaling by 5/4 is a bit of a hack method...
  })
  # blankheight <- reactive({
  #   blankwidth() / aspectRatio
  # })
  
  plot_main <- reactive({
    return(
      make_plot(dat = listData(), background = which_background(), displayEcoregions = display_ecoregions(),
                     plotSignificanceContours = display_trend_pvals())
    )
  })

  output$plt <- renderGirafe({
    pm <- plot_main()
    p <- pm$plot
    w <- pm$width
    h <- pm$height
    bw <- blankwidth()
    bh <- bw / {w / h}
    # w <- blankwidth()
    # h <- blankheight()
    x <- girafe(code = print(p),
                width_svg = (bw / {input$dpi}),
                height_svg = (bh / {input$dpi})
    )
    x <- girafe_options(x,
                        opts_sizing(rescale = FALSE),
                        opts_tooltip(css = 'background:white;'),
                        opts_zoom(min = 0.5, max = 10),
                        opts_hover(css = 'opacity:1.0;stroke-width:4;cursor:pointer;', reactive = TRUE),
                        opts_hover_inv(css = 'opacity:0.2;cursor:pointer;'),
                        opts_selection(type = 'multiple', css = 'opacity:1.0;stroke-width:4;'),
                        opts_selection_key(css = girafe_css("stroke:red; stroke-width:2px",
                                                            text = "stroke:none;fill:red;font-size:12px")),
                        opts_hover(css = 'opacity:1.0;stroke-width:4;cursor:pointer;', reactive = TRUE)
    )
    x
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
    # fd <- filtered_data()
    # fp <- fd$plastic_wide
    fp <- filtered_plastic_data()
    d <- fp$data
    d <- d[d$data_id %in% selected_data(),]
    if(nrow(d) < 1) return(NULL)
    
    d <- d[DATA_name_swap$original] # display only columns appearing in original data set
    names(d) <- sapply(1:ncol(d), function(z){
      n <- names(d)[z]
      i <- DATA_name_swap$original == n
      if(any(i)) n <- DATA_name_swap$spaces[i]
      n
    })
    
    d$URL <- paste0("<a href='", d$URL, "'>", d$URL,"</a>")
    # d <- subset(d, select = -c(Year, Coordinates, data_id))
    row.names(d) <- NULL
    d <- datatable(d,
                   options = list(
                     paging = FALSE, searching = FALSE # ,fillContainer = TRUE,
                    # bPaginate = FALSE,  searching = FALSE, fillContainer = TRUE
                     # autoWidth = TRUE
                   ), escape = FALSE)
    d
  })

  output$datatab_stations <- renderDataTable({
    # fd <- filtered_data()
    # d <- fd$stations
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
