# Individual map plots for paper

# Packages ----------------------------------------------------------------

# library(sp)
library(sf)
# library(maptools) # for the ContourLines2SLDF() function used to map contour lines DEPRECATED PACKAGE
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(scales)
library(ggiraph) # this is good for interactivity, but cannot be converted to grobs for layout...
library(reshape2)
library(this.path)

# # Install ggnewscale version 0.4.3 from GitHub as there seems to be bug in the
# # latest version -- this is annoying, but ggnewscale is vital for this map.
# # See this solution from the package developer: https://github.com/eliocamp/ggnewscale/issues/45
# # (I'm not convinced this helped! The issue seemed to be solved by including '_new'
# # as a suffix in arguments to the 'guides' function)
# library(remotes)
# ggnewscale_version <- '0.4.10'
# ggnewscale_download_path <- paste0('eliocamp/ggnewscale@v', ggnewscale_version)
# ggnewscale_available <- require(ggnewscale, quietly = TRUE)
# if(!ggnewscale_available){
#   install_github(ggnewscale_download_path)
#   library(ggnewscale)
# }else{
#   if(packageVersion('ggnewscale') != ggnewscale_version){
#     detach("package:ggnewscale", unload = TRUE)
#     remove.packages('ggnewscale')
#     install_github(ggnewscale_download_path)
#     library(ggnewscale)
#   }
# }

# Directory info ----------------------------------------------------------

dir_base <- dirname(this.dir()) # project base directory
dir_data <- paste(dir_base, 'data', sep = '/')
dir_map <- paste(dir_data, 'map', sep = '/')
dir_plots <- paste(dir_base, 'plots', sep = '/')
if(!dir.exists(dir_plots)) dir.create(dir_plots)

# Load data ---------------------------------------------------------------

source('functions.R', local = TRUE)

significantTrendsOnly <- TRUE
# loadTooltipFromFile <- TRUE # this is only needed for the interactive map
displayAllLitterTypes <- TRUE # If FALSE then only plastics are displayed. Set to TRUE to show samples of non-plastics -- this 

roundShipTime <- FALSE

res_options <- c('9x3', '3x1')
res <- res_options[1]

get_data(
  baseDirectory = dir_base, dataDirectory = dir_data, mapDirectory = dir_map,
  res = res, allLitterTypes = displayAllLitterTypes, sstType = 'trend',
  pHType = 'trend', shipSummaryDataOrRaw = 'raw',
  sstTrend_significantOnly = significantTrendsOnly,
  pHTrend_significantOnly = significantTrendsOnly, roundShipTime = roundShipTime,
  indexGridCells = FALSE, loadTooltipFromFile = TRUE)


# Plots -------------------------------------------------------------------

# 1. background data

nc <- nc_cells
eco <- eco_cells
latlim <- lat_lim_cells
set_plot_params(nc, DATA_sf, STATIONS_sf)

make_background_map <- function(
    background, group, backgroundOnly = TRUE, stationsOnly = FALSE,
    stationPopSize = FALSE, eco = FALSE, trendType = 'unspecified',
    ship_data = 'unspecified', shipOrPersonTime = 'ship time',
    pval_contours = FALSE, contourLineWidth = 0.25,
    ptStroke = 1, axisTextSize = 4, discreteColourScheme = FALSE,
    overlay_labels = NULL, latlim = NULL, nColours = NULL,
    background.na.remove = FALSE, na.in.legend = TRUE,
    manual_legend_breaks = FALSE){
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
                              summary = subset(d, ship_class == group & variable == shipOrPersonTime),
                              raw = subset(d, activity == group & variable == shipOrPersonTime))
    )
    dat <- list(nc = nc, background = d, plastic = DATA_sf[0,],
                stations = STATIONS_sf[0,], symbols = pltSymbols)
  }
  if(stationsOnly){
    d <- STATIONS_sf
    if(group %in% d$Type) d <- subset(d, Type == group)
    dat <- list(nc = nc, background = NULL, plastic = DATA_sf[0,],
                stations = d, symbols = pltSymbols)
  }
  mp <- make_plot(
    dat = dat, background = background, backgroundOnly = backgroundOnly,
    stationsOnly = stationsOnly, stationPopSize = stationPopSize,
    displayEcoregions = eco, plotSignificanceContours = pval_contours,
    contourLineWidth = contourLineWidth, ptStroke = ptStroke,
    axisTextSize = axisTextSize, singleLatAxis = FALSE,
    discreteColourScheme = discreteColourScheme,
    overlay_labels = overlay_labels, latlim = latlim, nColours = nColours,
    background.na.remove = background.na.remove, na.in.legend = na.in.legend,
    manual_legend_breaks = manual_legend_breaks)

  mp
}

# Include NA values for empty grid cells in all background data sets. Use the
# sea surface temperature data to define the grid cells because this data set
# has values for every cell.
all_cells <- unique(sst_poly[,'geometry'])
# Chlorophyll
for(l in unique(chl_poly$month)){
  d <- chl_poly[chl_poly$month == l,]
  i <- !all_cells$geometry %in% d$geometry
  if(sum(i) == 0) next
  j <- all_cells[i,]
  j$month <- l
  j$value <- NA
  chl_poly <- rbind(chl_poly, j)
  chl_poly <- chl_poly[order(chl_poly$month),]
}
# Krill
for(l in unique(krill_poly$month)){
  d <- krill_poly[krill_poly$month == l,]
  i <- !all_cells$geometry %in% d$geometry
  if(sum(i) == 0) next
  j <- all_cells[i,]
  j$month <- l
  j$value <- NA
  j$colourgroup <- NA
  krill_poly <- rbind(krill_poly, j)
  krill_poly <- krill_poly[order(krill_poly$month),]
}
# Sea surface temperature
for(l in unique(sst_poly$metric)){
  d <- sst_poly[sst_poly$metric == l,]
  i <- !all_cells$geometry %in% d$geometry
  if(sum(i) == 0) next
  j <- all_cells[i,]
  j$metric <- l
  j$value <- NA
  j$pval <- NA
  j$plevel <- NA
  sst_poly <- rbind(sst_poly, j)
  sst_poly <- sst_poly[order(sst_poly$metric),]
}
# pH
for(l in unique(pH_poly$metric)){
  d <- pH_poly[pH_poly$metric == l,]
  i <- !all_cells$geometry %in% d$geometry
  if(sum(i) == 0) next
  j <- all_cells[i,]
  j$metric <- l
  j$value <- NA
  j$pval <- NA
  j$plevel <- NA
  pH_poly <- rbind(pH_poly, j)
  pH_poly <- pH_poly[order(pH_poly$metric),]
}
# Shipping
for(l in unique(ship_poly$activity)){
  for(m in unique(ship_poly$variable)){
    d <- ship_poly[ship_poly$activity == l & ship_poly$variable == m,]
    i <- !all_cells$geometry %in% d$geometry
    if(sum(i) == 0) next
    j <- all_cells[i,]
    j$year <- 'all'
    j$activity <- l
    j$variable <- m
    j$value <- NA
    ship_poly <- rbind(ship_poly, j)
    ship_poly <- ship_poly[order(ship_poly$activity, ship_poly$variable),]
  }
}

# Chlorophyll
background <- 'chl'
month <- 'all'
p_chl <- make_background_map(background, month, discreteColourScheme = TRUE, latlim = latlim, background.na.remove = TRUE)
p_chl$plot_complete

# Krill
background <- 'krill'
month <- 'all'
p_krill <- make_background_map(background, month, discreteColourScheme = TRUE, latlim = latlim)
p_krill$plot_complete

# Sea surface temperature
background <- 'sst'
month <- 'all'
type <- 'trend'
pval_contours <- FALSE
contourLineWidth <- 0.25
p_sst <- make_background_map(background, month, trendType = type,
                             pval_contours = pval_contours,
                             contourLineWidth = contourLineWidth, 
                             discreteColourScheme = TRUE, latlim = latlim)
p_sst$plot_complete

# pH
background <- 'pH'
month <- 'all'
type <- 'trend'
pval_contours <- FALSE
contourLineWidth <- 0.25
p_pH <- make_background_map(background, month, trendType = type,
                            pval_contours = pval_contours,
                            contourLineWidth = contourLineWidth,
                            discreteColourScheme = TRUE, latlim = latlim)
p_pH$plot_complete

# Shipping
background <- 'ship'
ship_data <- 'raw'
ship_class <- 'all'
metrics <- c('ship time', 'person time')
p_ship_s <- make_background_map(background, ship_class, ship_data = ship_data,
                              discreteColourScheme = TRUE, latlim = latlim,
                              shipOrPersonTime = metrics[1])
p_ship_p <- make_background_map(background, ship_class, ship_data = ship_data,
                                discreteColourScheme = TRUE, latlim = latlim,
                                shipOrPersonTime = metrics[2])
p_ship_s$plot_complete
p_ship_p$plot_complete

# Shipping by vessel type
axisTextSize <- 3
ship_classes <- levels(ship_poly$activity)
manual_legend_breaks_p <- c(0,5,50,500,2500,9.5e4) # choose a single break vector to identical scale for all vessel types
manual_legend_breaks_p_all <- c(0,5,50,500,2500,1.5e5)
manual_legend_breaks_s <- c(0,5,50,100,200,400)
manual_legend_breaks_s_all <- c(0,5,50,100,200,850)

for(i in 1:length(ship_classes)){
  ship_class <- ship_classes[i]
  plt_name0 <- paste0('p_ship_', ship_class)
  for(j in metrics){
    manual_legend_breaks <- switch(j,
                                   `ship time` = if(ship_classes[i] != 'all') manual_legend_breaks_s else manual_legend_breaks_s_all,
                                   `person time` = if(ship_classes[i] != 'all') manual_legend_breaks_p else manual_legend_breaks_p_all)
    plt_name <- paste0(plt_name0, '_', substr(j,1,1))
    p <- make_background_map(
      background, ship_class, ship_data = ship_data, axisTextSize = axisTextSize,
      discreteColourScheme = TRUE, latlim = latlim, shipOrPersonTime = j, nColours = 6,
      na.in.legend = FALSE, manual_legend_breaks = manual_legend_breaks)
    assign(plt_name, p)
  }
}

plotAllShipTypes <- function(metric){
  shipPlots <- switch(metric,
                      `ship time` = list(p_ship_fishing = p_ship_fishing_s,
                                         p_ship_tourism = p_ship_tourism_s,
                                         p_ship_supply = p_ship_supply_s,
                                         p_ship_research = p_ship_research_s,
                                         p_ship_other = p_ship_other_s,
                                         p_ship_all = p_ship_all_s),
                      `person time` = list(p_ship_fishing = p_ship_fishing_p,
                                           p_ship_tourism = p_ship_tourism_p,
                                           p_ship_supply = p_ship_supply_p,
                                           p_ship_research = p_ship_research_p,
                                           p_ship_other = p_ship_other_p,
                                           p_ship_all = p_ship_all_p))
  rw <- c(0.9, 0.1)
  ra <- 2
  return(list(
    fishing = plot_grid(
      shipPlots$p_ship_fishing$plot,
      shipPlots$p_ship_fishing$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')),
      rel_widths = rw),
    tourism = plot_grid(
      shipPlots$p_ship_tourism$plot,
      shipPlots$p_ship_tourism$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')),
      rel_widths = rw),
   supply = plot_grid(
      shipPlots$p_ship_supply$plot,
      shipPlots$p_ship_supply$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')),
      rel_widths = rw),
    research = plot_grid(
      shipPlots$p_ship_research$plot,
      shipPlots$p_ship_research$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')),
      rel_widths = rw),
    other = plot_grid(
      shipPlots$p_ship_other$plot,
      shipPlots$p_ship_other$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')),
      rel_widths = rw),
    all = plot_grid(
      shipPlots$p_ship_all$plot,
      shipPlots$p_ship_all$legend + theme(plot.margin = unit(c(0,0,0,-ra), 'cm')),
      rel_widths = rw)
  ))
}

allShipPlots_s <- plotAllShipTypes(metrics[1])
allShipPlots_p <- plotAllShipTypes(metrics[2])

npanels <- 6
nrw <- ceiling(npanels ^ 0.5)
ncl <- ceiling(npanels / nrw)

p_ship_combine_s <- plot_grid(
  allShipPlots_s$fishing, allShipPlots_s$tourism, allShipPlots_s$supply,
  allShipPlots_s$research, allShipPlots_s$other, allShipPlots_s$all,
  nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
  labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9,
  label_fontface = 2, label_fontfamily = 'serif')  + 
  theme(plot.margin = unit(c(0,0.5,0,0), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'))

p_ship_combine_p <- plot_grid(
  allShipPlots_p$fishing, allShipPlots_p$tourism, allShipPlots_p$supply,
  allShipPlots_p$research, allShipPlots_p$other, allShipPlots_p$all,
  nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
  labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9,
  label_fontface = 2, label_fontfamily = 'serif')  + 
  theme(plot.margin = unit(c(0,0.5,0,0), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'))

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.45 # scaling factors
sc_h <- 0.9
pw <- sc_w * A4_w
ph <- sc_h * A4_h

fileName <- paste0('ship traffic data_ship days_', res, '.png')
ggsave(fileName, plot = p_ship_combine_s, device = 'png', path = dir_plots, width = pw, height = ph, units = 'in')

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.45 # scaling factors
sc_h <- 0.9
pw <- sc_w * A4_w
ph <- sc_h * A4_h

fileName <- paste0('ship traffic data_person days_', res, '.png')
ggsave(fileName, plot = p_ship_combine_p, device = 'png', path = dir_plots, width = pw, height = ph, units = 'in')



# stations
# labelSeas <- data.frame(
#   # Coordinates taken from Wikipedia
#   lon = c(-45, -175, 22, 7, -40, 160, -65.9, 92, -40, 50, -160),
#   lat = c(-73, -75, -68, -68, -57.5, -67, -58.58, -66, -50, -50, -50),
#   Name = c('Weddell\nSea', 'Ross\nSea', 'Riiser-Larsen\nSea', 'Lazarev\nSea',
#             'Scotia\nSea', 'Somov\nSea', 'Drake\n      Passage', 'Davis\nSea',
#            'Atlantic\nSector', 'Indian Ocean\nSector', 'Pacific\nSector'
#            ),
#   hjust = c(0.5, 0.6, 0.5, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0.5),
#   vjust = c(0.5, 0.7, 0, 0.1, 0.5, 0.75, 0.4, 1, 0.5, 0.5, 0.5),
#   xnudge = c(2, 0, 2, -12, 2.5, 0, -16, 0, 0, 0, 0),
#   ynudge = c(0, 0, 1, 0, 0, 0, 0, 1.5, 0, 0, 0))
# labelSeas$lon <- labelSeas$lon + labelSeas$xnudge
# labelSeas$lat <- labelSeas$lat + labelSeas$ynudge
# labelSeas$family <- 'serif'
# labelSeas$fontface <- 4
# labelSeas$fontface[{nrow(labelSeas)-2}:nrow(labelSeas)] <- 3
# labelSeas <- st_as_sf(labelSeas, coords = c('lon', 'lat'),
#                       crs = crs_world, remove = FALSE)
# labelSeas <- st_transform(labelSeas, crs_use)
# xy <- matrix(unlist(labelSeas$geometry), 2, nrow(labelSeas))
# labelSeas$x <- xy[1,]
# labelSeas$y <- xy[2,]
# labelSeas$colour <- 'skyblue3'
# labelSeas$size <- 3
# labelSeas$size[{nrow(labelSeas)-2}:nrow(labelSeas)] <- 3.5

background <- 'none'
group <- 'all'
p_stations <- make_background_map(background, group, backgroundOnly = FALSE,
                                  stationsOnly = TRUE, stationPopSize = TRUE,
                                  overlay_labels = NULL, latlim = latlim)
p_stations$plot_complete


# combine all plots
npanels <- 6
rw <- c(0.75, 0.25)
p_chl_ <- plot_grid(p_chl$plot, p_chl$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_krill_ <- plot_grid(p_krill$plot, p_krill$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_sst_ <- plot_grid(p_sst$plot, p_sst$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_pH_ <- plot_grid(p_pH$plot, p_pH$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_ship_s_ <- plot_grid(p_ship_s$plot, p_ship_s$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_ship_p_ <- plot_grid(p_ship_p$plot, p_ship_p$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))
p_stations_ <- plot_grid(p_stations$plot, p_stations$legend, rel_widths = rw) + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'))

p_all_s <- plot_grid(
  p_chl_, p_krill_, p_sst_, p_pH_, p_ship_s_, p_stations_,
  nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
  labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9,
  label_fontface = 2, label_fontfamily = 'serif') + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'))

p_all_p <- plot_grid(
  p_chl_, p_krill_, p_sst_, p_pH_, p_ship_p_, p_stations_,
  nrow = nrw, ncol = ncl, align = 'hv', axis = 'tblr',
  labels = LETTERS[1:npanels], label_x = 0.05, label_y = 0.9,
  label_fontface = 2, label_fontfamily = 'serif') + 
  theme(plot.margin = unit(c(0,0,0,0), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'))

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.6 # scaling factors
sc_h <- 1.2
pw <- sc_w * A4_w
ph <- sc_h * A4_h

fileName <- paste0('environ background data_ship days_', res, '.png')
ggsave(fileName, plot = p_all_s, device = 'png', path = dir_plots, width = pw,
       height = ph, units = 'in')

fileName <- paste0('environ background data_person days_', res, '.png')
ggsave(fileName, plot = p_all_p, device = 'png', path = dir_plots, width = pw,
       height = ph, units = 'in')


# 2. plastic data

nc <- nc_plastic
eco <- eco_plastic
set_plot_params(nc, DATA_sf, STATIONS_sf)
latlim <- lat_lim_plastic

make_plastic_map <- function(
    plastic, stations, ptSize, Symbols, background = 'none', group = NA,
    eco = FALSE, alpha = 0.85, latlim = NULL, overlay_labels = NULL){
  if(background == 'none'){
    dat <- list(nc = nc, background = NULL, plastic = plastic, stations = stations, symbols = Symbols)
  } else return(NULL)
  mp <- make_plot(dat = dat, background = background, plasticOnly = TRUE, displayEcoregions = eco, ptSize = ptSize, alpha = alpha, latlim = latlim, overlay_labels = overlay_labels)
  mp
}

whiteSpace <- function(n) paste0(rep(' ', n), collapse = '')
labelSeas <- data.frame(
  lon = c(-43, -160, 10, -14, -44, 175, -69, 90.5, -30, 65, -135, 70, -10, -120, 10, 109, 110, 140),
  lat = c(-74, -76.5, -70.5, -66.5, -53.5, -59, -61, -57.75, latlim[2], latlim[2], latlim[2], -73, -55, -79.5, -72.2, -67.25, -58.5, -60),
  Name = c(
    paste0(whiteSpace(1), 'Weddell\n', 'Sea', whiteSpace(1)),
    paste0(whiteSpace(0), 'Ross\n', whiteSpace(0), 'Sea'),
    paste0('Riiser-\n', whiteSpace(2), 'Larsen\n', whiteSpace(7), 'Sea'),
    paste0(whiteSpace(16), 'Lazarev\n', whiteSpace(14), 'Sea'),
    'Scotia\n\nSea',
    paste0('Somov\n', whiteSpace(2), 'Sea'),
    paste0('Drake\n', whiteSpace(2), 'Passage'),
    paste0(whiteSpace(1), 'Davis\n', whiteSpace(1), 'Sea'),
    'Atlantic\nSector', 'Indian Ocean\nSector', 'Pacific\nSector',
    'Amery Ice\nShelf',
    paste0('South', whiteSpace(12), '\nGeorgia', whiteSpace(10)),
    paste0('Antarctic', whiteSpace(6), '\nPeninsula', whiteSpace(10)),
    paste0(whiteSpace(30), 'Lützow-\n', whiteSpace(30), 'Holm Bay'),
    'Vincennes\nBay',
    paste0('Mawson', '\nSea'),
    paste0(whiteSpace(4), "D'Urville\n", 'Sea')
  ),
  hjust = c(0.5, 0.5, 0, 0.5, 0.5, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.5, 0.5, 0.5, 0.5),
  vjust = c(0.5, 1, 0, 0.1, 0.5, 0, 0.4, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 1, 0, 0.5, 0.5),
  xnudge = c(2, 0, 2, -12, 2.5, 0, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  ynudge = c(0, 0, 1, 0, 0, 0, 0, 0, -10, -10, -10, 0, 0, 0, 0, 0, 0, 0))
labelSeas$lon <- labelSeas$lon + labelSeas$xnudge
labelSeas$lat <- labelSeas$lat + labelSeas$ynudge
labelSeas$family <- 'serif'
labelSeas$fontface <- 4
labelSeas$fontface[grepl('Sector', labelSeas$Name)] <- 3
labelSeas <- st_as_sf(labelSeas, coords = c('lon', 'lat'),
                      crs = crs_world, remove = FALSE)
labelSeas <- st_transform(labelSeas, crs_use)
xy <- matrix(unlist(labelSeas$geometry), 2, nrow(labelSeas))
labelSeas$x <- xy[1,]
labelSeas$y <- xy[2,]
labelSeas$colour <- 'navy' # skyblue3
labelSeas$size <- 3
labelSeas$size[grepl('Sector', labelSeas$Name)] <- 5
labelSeas$size[apply(
  Vectorize(grepl, 'pattern')(c('Ross','Weddell', 'Scotia'),
                              labelSeas$Name), 1, any)] <- 4
labelSeas$size[apply(
  Vectorize(grepl, 'pattern')(c('Amery','Georgia','Peninsula','Holm','Vincennes'),
                              labelSeas$Name), 1, any)] <- 2.5

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

p_pla <- make_plastic_map(plastic, stations, ptSize, Symbols, alpha = 0.9,
                          latlim = latlim, overlay_labels = labelSeas)

p_plastic <- p_pla$plot_complete + 
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))

A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.2 # scaling factors
sc_h <- 0.8
pw <- sc_w * A4_w
ph <- sc_h * pw

ggsave('plastic samples.png', plot = p_plastic, device = 'png', path = dir_plots,
       width = pw, height = ph, units = 'in')


