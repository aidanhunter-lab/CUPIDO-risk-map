#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Evaluate risk of potential plastic emissions upon planktonic biota by mapping
# overlap of biota, existing stresses, and human activity (shipping and
# terrestrial facilities) as a proxy for plastic emissions.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Summary -----------------------------------------------------------------
# The novelty of the work is assessing the potential of plastic to contribute to
# environmental stress.
# Consider two kinds of 'risk':
# (1) Risk of interaction between plastic and biota.
# (2) Plastic overlap with existing stresses -- multi-stressor risk.
# For (1) risk of interaction we use data on spatial distributions of primary
# productivity (chlorophyll) and of a keystone zooplankton species (krill) to
# assess where the distributions of these biota are likely to overlap with areas
# of relatively high plastic emissions.
# For (2) multi-stressor risk we use spatial distribution data on long-term
# trends of sea surface temperature and pH to investigate where existing areas
# of concern overlap with areas of relatively high plastics emissions. These
# variables were selected as many regard their long term variability to pose a
# risk to the existing ecosystem structure and function.

# Note that although temperature and pH may decrease or increase over time and
# that both increases and decreases were observed in different regions, we 
# concentrate only on the risk of increasing temperature and decreasing pH
# because these trends are currently the focus of much research effort, whereas
# the impacts of decreasing temperature and rising pH in Antarctic waters are
# lesser known due to much less research effort/funding. This is important to
# note with respect to temperature because there is a greater sea surface area
# of temperature decline than there is of temperature rise, whereas almost all
# Antarctic waters exhibit a decline in pH.

# Load packages -----------------------------------------------------------

# CHECK WHICH PACKAGES I CAN REMOVE

library(MASS) # for Box-Cox transform
library(car) # for Yeo-Johnson
library(ggplot2)
library(gridExtra)
library(cowplot)
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

# mmtable2 may be useful for displaying ranking values...
mmtable2_available <- require(mmtable2, quietly = TRUE)
if(!mmtable2_available){
  install_github("ianmoran11/mmtable2")
  library(mmtable2)
}

library(this.path)


# Directories -------------------------------------------------------------

dir_base <- dirname(this.dir())
dir_data <- paste(dir_base, 'data', sep = '/')
dir_map <- paste(dir_data, 'map', sep = '/')
dir_plots <- paste(dir_base, 'plots', sep = '/')

# Load functions ----------------------------------------------------------

source('functions.R') 

# Load and filter data ----------------------------------------------------

significantTrendsOnly <- TRUE
loadTooltipFromFile <- TRUE
displayAllLitterTypes <- TRUE

res <- '9x3'

get_data(
  baseDirectory = dir_base, dataDirectory = dir_data, mapDirectory = dir_map,
  res = res, allLitterTypes = displayAllLitterTypes, sstType = 'trend',
  pHType = 'trend', shipSummaryDataOrRaw = 'raw',
  sstTrend_significantOnly = significantTrendsOnly,
  pHTrend_significantOnly = significantTrendsOnly,
  loadTooltipFromFile = loadTooltipFromFile)

nc <- nc_cells

# Filter by month(s) -- summer months only.
mon <- 'all' # select month: 1, 2, 3, and/or 'all'

chl <- subset(chl_poly, month == mon)
krill <- subset(krill_poly, month == mon)

if(!is.null(sst_poly$month))
  sst <- subset(sst_poly, month == mon) else
    sst <- sst_poly

if(!is.null(pH_poly$month))
  pH <- subset(pH_poly, month == mon) else
    pH <- pH_poly
      
ship_type <- 'all'
ship <- subset(ship_poly, activity == ship_type)

nranks <- 5 # number of risk categories


# Chlorophyll -------------------------------------------------------------

# High productivity regions are associated with greater potential for environmental
# risk, so risk increases with chlorophyll concentration.

# Specify values for data rankings by observation of data distribution.

x <- {1:nrow(chl)} / nrow(chl)
y <- sort(chl$value)
plot(x, y)
plot(x, y, log = 'y')

yr <- y[-c({length(y)-1}:length(y))]
yl <- log10(yr)
r <- 10 ^ seq(min(yl), max(yl), length = nranks + 1) # even spread on log scale
r <- round(r / {orderOfMagnitude(r) / 10}) * {orderOfMagnitude(r) / 10}
r[c(1,nranks+1)] <- range(y)
rank_limits_chl <- r
abline(h = r)

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= chl$value & chl$value < r[i+1]} else{
      j <- r[i] <= chl$value & chl$value <= r[i+1]
    }
  v <- chl$value[j]
  p <- r[i:{i+1}]
  chl$score[j] <- {i-1} + {v - r[i]} / diff(p)
}

chl$rank <- score2rank(chl$score)

plt_chl <- plot_fun(chl, Title = 'Chlorophyll', legend_title = 'prevalence')
plt_chl$plot

# Nice plot of raw data
chl$order <- 1:nrow(chl)
chl <- chl[order(chl$value),]
chl$index <- 1:nrow(chl)
chl <- chl[order(chl$order),]

plt_chl_raw <- 
  ggplot(chl) + 
  geom_point(aes(x = index, y = value, fill = rank, colour = rank), shape = 21, size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = 'prevalence', guide = 'legend') +
  scale_fill_viridis_c(name = 'prevalence', guide = 'legend') +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0, hjust = 0) +
  xlab(expression(grid ~ cell ~ index)) +
  ylab(expression(chlorophyll ~ abundance ~ (mg / m^3))) +
  scale_y_continuous(trans = 'log10', labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())

# Krill -------------------------------------------------------------------

# Krill is a Southern Ocean key-stone species. Regions with high krill density 
# have greater potential for environmental risk.
# Take a similar approach to chlorophyll but note that there are grid cells with
# no data and grid cells with measurement zeros. Grid cells without data are left
# blank (white) while grid cells with measurement zeros are assigned the lowest rank.

# Specify values for data rankings by observation of data distribution.

x <- {1:nrow(krill)} / nrow(krill)
y <- sort(krill$value)
par(mfrow = c(1,1))
plot(x, y)
plot(x, y, log = 'y')

zero_krill <- krill$value == 0
observed_krill <- !zero_krill

krill_p <- krill[observed_krill,] # positive values

x_p <- {1:nrow(krill_p)} / nrow(krill_p)
y_p <- sort(krill_p$value)
plot(x_p, y_p)
plot(x_p, y_p, log = 'y')

# Rank as an even spread on log scale after removing outliers.
y_r <- y_p[-c(1,{length(y_p)-2}:length(y_p))]
x_r <- 1:length(y_r)
plot(x_r, y_r, log = 'y')
yl <- log10(y_r)

r <- 10 ^ seq(min(yl), max(yl), length = nranks + 1)
r <- round(r / {orderOfMagnitude(r) / 10}) * {orderOfMagnitude(r) / 10} # round to nice values
r[c(1, nranks+1)] <- c(0, max(y)) # adjust for outliers and measurement zeros
abline(h = r)

rank_limits_krill <- r

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= krill$value & krill$value < r[i+1]} else{
      j <- r[i] <= krill$value & krill$value <= r[i+1]
    }
  v <- krill$value[j]
  p <- r[i:{i+1}]
  krill$score[j] <- {i-1} + {v - r[i]} / diff(p)
}

krill$rank <- score2rank(krill$score)

plt_krill <- plot_fun(krill, Title = 'Krill', legend_title = 'prevalence')
plt_krill$plot


# Nice plot of raw data
krill$order <- 1:nrow(krill)
krill <- krill[order(krill$value),]
krill$index <- 1:nrow(krill)
krill <- krill[order(krill$order),]

plt_krill_raw <- 
  ggplot(krill) + 
  geom_point(aes(x = index, y = value, fill = rank, colour = rank), shape = 21, size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = 'prevalence', guide = 'legend') +
  scale_fill_viridis_c(name = 'prevalence', guide = 'legend') +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0, hjust = 0) +
  xlab(expression(grid ~ cell ~ index)) +
  ylab(expression(krill ~ abundance ~ (individuals / m^2))) +
  # scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) + 
  scale_y_continuous(trans = 'log10', labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())



# Temperature -------------------------------------------------------------

# Potential environmental risk due to temperature change is typically associated
# rising temperature, but the SST data record regions of both rising and falling
# temperatures.
# Set cells with declining or non-significant sst trends to minimum risk.

sst_pos <- sst$value > 0 # index grid cells with trend of increasing temperature (all p-values)
sst_neg <- !sst_pos # cells with negative values

insig <- sst$pval > 0.05 # insignificant trends
min_risk <- sst_neg | insig # cells to assign minimum risk from sst trends

sst_selected_cells <- !min_risk # other cells need ranked based on size of positive trend
sst_selected <- sst[sst_selected_cells,]
o <- order(sst_selected$value)

y <- sst_selected$value[o]
pv <- sst_selected$plevel[o]
x <- 1:length(y)
plot(x, y)

# trim outliers then group
ol <- c(1:4, {length(y)-3}:length(y))
# normalise with square root transform
yl <- y ^ 0.5
y_r <- yl[-ol]
yrange <- range(y_r)
# Set all positive trends in sst as rank 2-5, leaving rank 1 for the non-significant or negative trends.
r <- c(seq(yrange[1], yrange[2], length = nranks)) ^ 2
r <- c(min(sst$value), r) # account for negative trends
r[2:nranks] <- round(r[2:nranks] / {orderOfMagnitude(r[2:nranks]) / 10}) * {orderOfMagnitude(r[2:nranks]) / 10}
r[c(1,nranks+1)] <- range(sst$value)
rank_limits_sst <- r

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= sst$value & sst$value < r[i+1]} else{
      j <- r[i] <= sst$value & sst$value <= r[i+1]
    }
  v <- sst$value[j]
  p <- r[i:{i+1}]
  sst$score[j] <- {i-1} + {v - r[i]} / diff(p)
}

sst$rank <- score2rank(sst$score)

plt_sst <- plot_fun(sst, Title = 'SST rising trends', legend_title = 'stress')
plt_sst$plot

# Nice plot of raw data
sst$order <- 1:nrow(sst)
sst <- sst[order(sst$value),]
sst$index <- 1:nrow(sst)
sst <- sst[order(sst$order),]

plt_sst_raw <- 
  ggplot(sst) + 
  geom_point(aes(x = index, y = value, fill = rank, colour = rank), shape = 21, size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = 'stress', guide = 'legend') +
  scale_fill_viridis_c(name = 'stress', guide = 'legend') +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[3:nranks], label = as.character(r[3:nranks]), vjust = 0, hjust = 0) +
  annotate('text', x = 1, y = r[2], label = as.character(r[2]), vjust = 1, hjust = 0) +
  xlab(expression(grid ~ cell ~ index)) +
  ylab(expression(SST ~ trend ~ (degree * C / year))) +
  # scale_y_continuous(trans = 'sqrt', labels = function(x) format(x, scientific = FALSE)) +
  # scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())


# pH -------------------------------------------------------------

pH_neg <- pH$value < 0 # index grid cells with trend of decreasing pH (all p-values)
pH_pos <- !pH_neg # cells with positive values

pH_selected_cells <- pH_neg
pH_selected <- pH[pH_selected_cells,]

y <- sort(pH_selected$value)
x <- 1:length(y)
plot(x, y)

# trim outliers then group
ol <- c(1, length(y))
y_r <- y[-ol]
x_r <- 1:length(y_r)

yrange <- range(y_r)
rlo <- floor(yrange[1] / orderOfMagnitude(yrange[1]) * 10) * orderOfMagnitude(yrange[1]) / 10
rhi <- ceiling(yrange[2] / orderOfMagnitude(yrange[2]) *10) * orderOfMagnitude(yrange[2]) / 10
r <- seq(rlo, rhi, length = nranks+1)
r <- round(r / {orderOfMagnitude(r) / 10}) * {orderOfMagnitude(r) / 10}
r[1] <- min(y)
r[nranks + 1] <- max(y)
abline(h = r[2:nranks], lty = 3)

rank_limits_pH <- r

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= pH$value & pH$value < r[i+1]} else{
      j <- r[i] <= pH$value & pH$value <= r[i+1]
    }
  v <- pH$value[j]
  p <- r[i:{i+1}]
  pH$score[j] <- {nranks-i} + {r[i+1] - v} / diff(p)
}
pH <- pH[pH_selected_cells,]
pH$rank <- score2rank(pH$score)

plt_pH <- plot_fun(pH, Title = 'pH declining trends', legend_title = 'stress')
plt_pH$plot

# Nice plot of raw data
pH$order <- 1:nrow(pH)
pH <- pH[order(pH$value),]
pH$index <- 1:nrow(pH)
pH <- pH[order(pH$order),]

plt_pH_raw <- 
  ggplot(pH) + 
  geom_point(aes(x = index, y = value, fill = rank, colour = rank), shape = 21, size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = 'stress', guide = 'legend') +
  scale_fill_viridis_c(name = 'stress', guide = 'legend') +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0, hjust = 0) +
  xlab(expression(grid ~ cell ~ index)) +
  ylab(expression(pH ~ trend ~ (1 / year))) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())


# Shipping ----------------------------------------------------------------

x <- {1:nrow(ship)} / nrow(ship)
y <- sort(ship$value)
par(mfrow = c(1,1))
plot(x, y)
plot(x, y, log = 'y')

y_ <- y[y >= 1] # data less than 1 is lowest rank
x_ <- 1:length(y_)
plot(x_, y_, log = 'y')
y_r <- y_[-c({length(y_)-3}:length(y_))] # omit outliers
x_r <- 1:length(y_r)
plot(x_r, y_r, log = 'y')

yl <- log10(y_r)
r <- 10 ^ seq(min(yl), max(yl), length = nranks)
r[1] <- 1
r <- c(0, r)
r <- round(r / {orderOfMagnitude(r) / 10}) * {orderOfMagnitude(r) / 10} # round to nice values
r[nranks + 1] <- max(y) # include outliers
plot(x, y, log = 'y')
abline(h = r)

rank_limits_ship <- r

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= ship$value & ship$value < r[i+1]} else{
      j <- r[i] <= ship$value & ship$value <= r[i+1]
    }
  v <- ship$value[j]
  p <- r[i:{i+1}]
  ship$score[j] <- {i-1} + {v - r[i]} / diff(p)
}

ship$rank <- score2rank(ship$score)

plt_ship <- plot_fun(ship, Title = 'Ship traffic', legend_title = 'activity')
plt_ship$plot

# Nice plot of raw data
ship$order <- 1:nrow(ship)
ship <- ship[order(ship$value),]
ship$index <- 1:nrow(ship)
ship <- ship[order(ship$order),]

plt_ship_raw <- 
  ggplot(ship) + 
  geom_point(aes(x = index, y = value, fill = rank, colour = rank), shape = 21, size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = 'activity', guide = 'legend') +
  scale_fill_viridis_c(name = 'activity', guide = 'legend') +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0, hjust = 0) +
  xlab(expression(grid ~ cell ~ index)) +
  ylab(expression(ship ~ traffic ~ (ship ~ days / year))) +
  scale_y_continuous(trans = 'log10', labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())



# Facilities --------------------------------------------------------------

# Each facility is assigned a rank determined by its population size, then
# mapped grid cells are assigned a disturbance rank determined by their
# population density (people/km^2) calculated as a function of distance to
# facilities.

stations <- STATIONS_sf

stations$Seasonality <- factor(stations$Seasonality,
                               levels = c('Year-Round','Seasonal'),
                               labels = c('year-round', 'seasonal'))

# Population size determines the disturbance rank. 17 (15%) of the facilities do
# not have data on population size. These can either be omitted from the analyses,
# or infilled using mean values.
infillMissingPopSize <- TRUE
if(!infillMissingPopSize){
  # omit facilities with missing data
  stations <- stations[!is.na(stations$Peak_Population),]
}else{
  # or infill using mean values
  i <- is.na(stations$Peak_Population)
  typesMissing <- unique(stations$Type[i])
  for(j in 1:length(typesMissing)){
    m <- stations$Type == typesMissing[j]
    p <- stations$Peak_Population[m & !i]
    infill <- round(10 ^ mean(log10(p)))
    stations$Peak_Population[m & i] <- infill
  }
}

x <- {1:nrow(stations)} / nrow(stations)
y <- sort(stations$Peak_Population)
par(mfrow = c(1,1))
plot(x, y)
plot(x, y, log = 'y')

# Manually choose the cut-offs for risk ranks...
r <- c(min(y),10,25,50,100,max(y))
abline(h = r)

# Account for seasonality by reducing population size for seasonal facilities.
# Leave population size of year-round facilities unchanged.
stations$Peak_Population_adjusted <-
  stations$Peak_Population * 0.5 * {{stations$Seasonality == 'year-round'} + 1}
y <- sort(stations$Peak_Population_adjusted)
plot(x, y, log = 'y')
abline(h = r)

rank_limits_stationPopSize <- r

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= stations$Peak_Population_adjusted & stations$Peak_Population_adjusted < r[i+1]} else{
      j <- r[i] <= stations$Peak_Population_adjusted & stations$Peak_Population_adjusted <= r[i+1]
    }
  v <- stations$Peak_Population_adjusted[j]
  p <- r[i:{i+1}]
  stations$score_pop[j] <- {i-1} + {v - r[i]} / diff(p)
}

stations$rank_pop <- score2rank(stations$score_pop)

# Facility distance from the coast is a useful metric, but it's not needed for 
# calculating disturbance ranks.

# Distance from coast -- these distances were crudely estimated from lo-res coastline data
y <- sort(stations$Distance_From_Coast)
plot(x, y, log = 'y')

# manually select risk rank cut-offs for distance from coast
r <- c(max(y), 100, 50, 10, 5, min(y))
abline(h = r)

rank_limits_stationCoastDist <- r

for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] >= stations$Distance_From_Coast & stations$Distance_From_Coast > r[i+1]} else{
      j <- r[i] >= stations$Distance_From_Coast & stations$Distance_From_Coast >= r[i+1]
    }
  v <- stations$Distance_From_Coast[j]
  p <- r[i:{i+1}]
  stations$score_dist[j] <- {i-1} + {v - r[i]} / diff(p)
}

stations$rank_dist <- score2rank(stations$score_dist)


# For station's disturbance ranking use only population size
stations$rank <- stations$rank_pop

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now calculate disturbance rankings for mapped grid cells
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid_cells <- st_transform(grid_cells_all, crs = crs_world)
grid_centroids <- st_centroid(grid_cells)

# Get lat-lon coords of facilities
stations <- st_transform(stations, crs = crs_world)
xy <- stations$geometry
xy <- as.data.frame(st_coordinates(xy)) # get the lon-lats of each facility
stations$lon <- xy$X
stations$lat <- xy$Y

# Data dimensions
ngrid_cells <- nrow(grid_cells)
nfacilities <- nrow(stations)

# Find distances between facilities and grid cell centroids
dist2cells <- st_distance(grid_centroids, stations)
if(attr(dist2cells, 'units')$numerator == 'm') dist2cells <- dist2cells * 1e-3 # convert m to km
dist2cells <- matrix(as.numeric(dist2cells), nrow = ngrid_cells, ncol = nfacilities) # convert class to standard matrix
minDist2cell <- apply(dist2cells, 2, min) # distance (km) of facility to nearest grid cell
nearest_cell <- apply(dist2cells, 2, which.min)
stationsInCells <- st_intersection(stations, grid_cells) # all facilities within a mapped grid cell
ind <- stations$order %in% stationsInCells$order
# Indicate the grid cell containing each facility - NA indicates inland facilities not contained in a grid cell
stations$cell_index[ind] <- grid_centroids$cell_index[nearest_cell[ind]]

# List the facilities contained within each grid cell
cellsWithFacilities <- sort(unique(stations$cell_index))
stationsInCell <- vector('list', length = ngrid_cells)
for(i in 1:length(cellsWithFacilities)){
  cell <- cellsWithFacilities[i]
  facilities <- which(stations$cell_index == cell)
  stationsInCell[[cell]] <- facilities
}
grid_cells$stations <- stationsInCell # NULL values indicate grid cells not containing facilities

# Get grid cell areas
grid_cell_area <- st_area(grid_cells)
if(all(attr(grid_cell_area, 'units')$numerator == 'm')) grid_cell_area <-
  grid_cell_area * 1e-6 # convert m^2 to km^2
grid_cells$area <- as.numeric(grid_cell_area)

# For each grid cell we know the great circle distance from each facility. From
# this distance, and the Earth radius, we can calculate the area of the (globe surface)
# circle surrounding a facility and intersecting the grid cell centroid. Within
# this area, A, the population density associated with the facility is Pd = P/A,
# where P is facility population. To find the associated 'population density', Pg,
# within the intersecting grid cell, we find the ratio of grid cell area, Ag, to
# circle area, A, then multiply: Pg = (Ag/A)(P/A). Something to bear in mind is
# that, if the circle intersects the grid cell centroid, only a fraction of the
# grid cell is enclosed within the circle, so dividing the grid cell area by 2 will
# produce a better estimate => Pg = (1/2)(Ag/A)(P/A) people per km^2.
# Then, for each grid cell, these population densities are summed over all facilities.
# It may also be useful to include some limit on the radius used when calculating
# areas, say, 1000km. Also, circles that intersect large areas of land may require
# special attention.

circleAreaOnGlobe <- function(s, R = 6371, distanceLimit = NULL){
  # Returns surface area of circle defined by 'radius' s. The 'radius' is the arc
  # length between two points on the globe. R is average Earth radius (km); s
  # should have units of km; returned area has units of km^2.
  # If distanceLimit (km) is not NULL then cells further than distanceLimit from
  # facilities are not impacted by those facilities, simply by returning Inf.
  Area <- 2 * pi * R^2 * {1 - cos(s / R)}
  if(is.numeric(distanceLimit)) Area[s > distanceLimit] <- Inf
  return(Area)
}

# I'm not sure if it's sensible to set a distance limit -- it seems a bit ad-hoc...
distanceLimit <- NULL
# distanceLimit <- 1000 # set a 1000 km distance limit, beyond which facility populations do not have impact

# Area of the circles centred on facilities and intersecting grid cell centroids
station2cellCircleArea <- circleAreaOnGlobe(dist2cells, distanceLimit = distanceLimit) # km^2
# Population density (people/km^2) within circle centered at facility and intersecting grid cells
popDensity_circle <- sweep(1 / station2cellCircleArea, 2, stations$Peak_Population_adjusted, '*')
# Scaling factor -- ratio of grid cell area to circle area
cell2circleAreaRatio <- sweep(1 / station2cellCircleArea, 1, grid_cells$area, '*')
# Given that circles intersect centroids of grid cells, we need another scaling
# factor representing the approximate(?) fraction of grid cell contained in circle
modAreaRatio <- 0.393
# Each facility's contribution to grid cell population density (people/km^2)
popDensity_cell <- modAreaRatio * cell2circleAreaRatio * popDensity_circle

# The above method is appropriate only for grid cells that are distant from facilities.
# The method breaks down for the grid cells containing facilities, but for these
# cells we can directly calculate population density.
# For each cell containing facilities, recalculate the population density
for(i in 1:length(cellsWithFacilities)){
  cell <- cellsWithFacilities[i]
  ind <- grid_cells$cell_index == cell
  a <- grid_cells$area[ind]
  facilities <- grid_cells$stations[[cell]]
  for(j in 1:length(facilities)){
    facility <- facilities[j]
    pop <- stations$Peak_Population_adjusted[facility]
    popDensity <- pop / a
    popDensity_cell[ind, facility] <- popDensity
  }
}

# Sum over facilities to get population density per grid cell
popDensity_cell <- apply(popDensity_cell, 1, sum)
grid_cells$popDensity <- popDensity_cell
# Scale units to get people / 10,000 km^2
grid_cells$popDensity <- 1e4 * grid_cells$popDensity

# As with previous data, specify some ranking for disturbance in grid cells associated
# with population density.
y <- sort(grid_cells$popDensity) # people / 10000 km^2
x <- 1:length(y)

# Cells with zero population density are not assigned a rank and will not appear
# on the map. This only applies when distanceLimit is set.
y_ <- y[y>0]
x_ <- 1:length(y_)
plot(x_, y_, log = 'y')

y_r <- y_[1:{length(y_)-3}] # omit three large outliers
x_r <- 1:length(y_r)

rlo <- floor(y_r[1] / orderOfMagnitude(y_r[1])) * orderOfMagnitude(y_r[1])
rhi <- floor(y_r[length(y_r)] / orderOfMagnitude(y_r[length(y_r)])) * orderOfMagnitude(y_r[length(y_r)])

rank_limits_stations <- 10 ^ seq(log10(rlo), log10(rhi), length = nranks + 1)
rank_limits_stations <- round(rank_limits_stations / {orderOfMagnitude(rank_limits_stations) / 10}) * {orderOfMagnitude(rank_limits_stations) / 10}
rank_limits_stations[c(1,nranks+1)] <- range(y) # account for outliers

abline(h = rank_limits_stations)

# Manually specify ranking limits
if(is.null(distanceLimit))
  rank_limits_stations <- c(0.001, 0.01, 0.1, 1, 10, 200)

if(!is.null(distanceLimit))
  if(distanceLimit == 1000)
    rank_limits_stations <- c(0.0001, 0.01, 0.1, 1, 10, 200)

r <- rank_limits_stations

# Assign disturbance ranks to each grid cell
for(i in 1:nranks){
  if(i < nranks){
    j <- r[i] <= grid_cells$popDensity & grid_cells$popDensity < r[i+1]} else{
      j <- r[i] <= grid_cells$popDensity & grid_cells$popDensity <= r[i+1]
    }
  v <- grid_cells$popDensity[j]
  p <- r[i:{i+1}]
  grid_cells$score[j] <- {i-1} + {v - r[i]} / diff(p)
}
grid_cells$rank <- score2rank(grid_cells$score)

# Convert to coordinate reference system used in plots
stations <- st_transform(stations, crs = crs_use)
grid_cells <- st_transform(grid_cells, crs = crs_use)

# Omit cells with unspecified rank -- only applies when distanceLimit is set.
grid_cells_reduced <- grid_cells[grid_cells$rank %in% 1:nranks,]

# Make plots
pltStations <- plot_fun_stations(data_stations = stations, data_grid_cells = grid_cells_reduced, rank_limits_PopSize = rank_limits_stationPopSize, Title = 'Facilities')
pltStations$plot

grid.arrange(plt_ship$plot + theme(legend.position = 'none'), pltStations$plot + theme(legend.position = 'none'), pltStations$legend,
             ncol = 3, widths = c(0.425,0.425,0.15))


# Nice plot of raw data
grid_cells_reduced$order <- 1:nrow(grid_cells_reduced)
grid_cells_reduced <- grid_cells_reduced[order(grid_cells_reduced$popDensity),]
grid_cells_reduced$index <- 1:nrow(grid_cells_reduced)
grid_cells_reduced <- grid_cells_reduced[order(grid_cells_reduced$order),]

plt_stations_raw <- 
  ggplot(grid_cells_reduced) + 
  geom_point(aes(x = index, y = popDensity, fill = rank, colour = rank), shape = 21, size = 3, alpha = 0.8) +
  scale_colour_viridis_c(name = 'activity', guide = 'legend') +
  scale_fill_viridis_c(name = 'activity', guide = 'legend') +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0, hjust = 0) +
  xlab(expression(grid ~ cell ~ index)) +
  ylab(expression(pop. ~ density ~ (people / 10000 ~ km^2))) +
  scale_y_continuous(trans = 'log10', labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())


# Combine -----------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~
# Combine all plots
#~~~~~~~~~~~~~~~~~~

plt_chl_combined <- ggdraw(arrangeGrob(
  plt_chl_raw,
  plt_chl$plot + labs(title = element_blank()),
  nrow = 2, heights = c(0.35, 0.65), top = plt_chl$plot$labels$title))
plt_krill_combined <- ggdraw(arrangeGrob(
  plt_krill_raw,
  plt_krill$plot + labs(title = element_blank()),
  nrow = 2, heights = c(0.35, 0.65), top = plt_krill$plot$labels$title))
plt_sst_combined <- ggdraw(arrangeGrob(
  plt_sst_raw,
  plt_sst$plot + labs(title = element_blank()),
  nrow = 2, heights = c(0.35, 0.65), top = plt_sst$plot$labels$title))
plt_pH_combined <- ggdraw(arrangeGrob(
  plt_pH_raw,
  plt_pH$plot + labs(title = element_blank()),
  nrow = 2, heights = c(0.35, 0.65), top = plt_pH$plot$labels$title))
plt_ship_combined <- ggdraw(arrangeGrob(
  plt_ship_raw,
  plt_ship$plot + labs(title = element_blank()),
  nrow = 2, heights = c(0.35, 0.65), top = plt_ship$plot$labels$title))
plt_stations_combined <- ggdraw(arrangeGrob(
  plt_stations_raw,
  pltStations$plot + labs(title = element_blank()),
  nrow = 2, heights = c(0.35, 0.65), top = pltStations$plot$labels$title))


plt_biota <- ggdraw(arrangeGrob(
  plt_chl$plot + theme(legend.position = 'none'),
  plt_krill$plot + theme(legend.position = 'none'),
  plt_chl$legend,
  nrow = 1, widths = c(0.425, 0.425, 0.15), top = 'Biota'))
plt_stress <- ggdraw(arrangeGrob(
  plt_sst$plot + theme(legend.position = 'none'),
  plt_pH$plot + theme(legend.position = 'none'),
  plt_pH$legend,
  nrow = 1, widths = c(0.425, 0.425, 0.15), top = 'Stressors'))
plt_activity <- ggdraw(arrangeGrob(
  plt_ship$plot + theme(legend.position = 'none'),
  pltStations$plot + theme(legend.position = 'none'),
  pltStations$legend,
  nrow = 1, widths = c(0.425, 0.425, 0.15), top = 'Activity'))

# Save these plots for supplementary material
# Chlorophyll and krill
plt_raw1 <- plt_chl_raw + theme(legend.position = 'none')
plt_raw2 <- plt_krill_raw + theme(legend.position = 'none')
leg_raw <- ggdraw(get_legend(plt_chl_raw))
plt_map1 <- plt_chl$plot + theme(legend.position = 'none') + 
  labs(title = element_blank())
plt_map2 <- plt_krill$plot + theme(legend.position = 'none') + 
  labs(title = element_blank())
leg_map <- plt_chl$legend

rh <- c(0.4, 0.6) # relative heights of panels
lx <- c(0.05, 0.15) # label adjustment factors
ly <- c(1, 0.85)

label_fontface <- 2
label_fontfamily <- 'serif'

title_left <- ggdraw() +
  draw_label('Chlorophyll', fontface = label_fontface,
             fontfamily = label_fontfamily, x = 0.55)
title_right <- ggdraw() +
  draw_label('Krill', fontface = label_fontface,
             fontfamily = label_fontfamily, x = 0.55)
plt_left <- plot_grid(
  plt_raw1 + theme(plot.margin = unit(c(0.5,0.5,0,1), 'cm')),
  plt_map1 + theme(plot.margin = unit(c(0,0,0,0), 'cm')),
  ncol = 1, labels = c('A','C'), rel_heights = rh, align = 'v',
  label_x = lx, label_y = ly, label_fontface = label_fontface,
  label_fontfamily = label_fontfamily)
plt_left <- plot_grid(title_left, plt_left,
                      ncol = 1, rel_heights = c(0.05, 0.95))
plt_right <- plot_grid(
  plt_raw2 + theme(plot.margin = unit(c(0.5,0.5,0,1), 'cm')),
  plt_map2 + theme(plot.margin = unit(c(0,0,0,0), 'cm')),
  ncol = 1, labels = c('B','D'), rel_heights = rh, align = 'v',
  label_x = lx, label_y = ly, label_fontface = label_fontface,
  label_fontfamily = label_fontfamily)
plt_right <- plot_grid(title_right, plt_right,
                       ncol = 1, rel_heights = c(0.05, 0.95))
plt_leg <- plot_grid(leg_raw, leg_map,
                     ncol = 1, rel_heights = rh, align = 'v')
plt_biota_combined <- plot_grid(plt_left, plt_right, plt_leg,
                                ncol = 3, rel_widths = c(1,1,0.2)) + 
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))

a4_dim <- c(8.3, 11.7)
pw <- 1.5 * a4_dim[1]
ph <- 1.5 * 0.5 * a4_dim[2]
ggsave('biota_prevalence_ranking.png', plt_biota_combined, 'png', dir_plots,
       width = pw, height = ph, units = 'in')

# SST and pH
plt_raw1 <- plt_sst_raw + theme(legend.position = 'none')
plt_raw2 <- plt_pH_raw + theme(legend.position = 'none')
leg_raw <- ggdraw(get_legend(plt_sst_raw))
plt_map1 <- plt_sst$plot +
  theme(legend.position = 'none') + labs(title = element_blank())
plt_map2 <- plt_pH$plot +
  theme(legend.position = 'none') + labs(title = element_blank())
leg_map <- plt_sst$legend

title_left <- ggdraw() +
  draw_label('SST', fontface = label_fontface,
             fontfamily = label_fontfamily, x = 0.55)
title_right <- ggdraw() +
  draw_label('pH', fontface = label_fontface,
             fontfamily = label_fontfamily, x = 0.55)
plt_left <- plot_grid(
  plt_raw1 + theme(plot.margin = unit(c(0.5,0.5,0,1), 'cm')),
  plt_map1 + theme(plot.margin = unit(c(0,0,0,0), 'cm')),
  ncol = 1, labels = c('A','C'), rel_heights = rh, align = 'v',
  label_x = lx, label_y = ly, label_fontface = label_fontface,
  label_fontfamily = label_fontfamily)
plt_left <- plot_grid(title_left, plt_left,
                      ncol = 1, rel_heights = c(0.05, 0.95))
plt_right <- plot_grid(
  plt_raw2 + theme(plot.margin = unit(c(0.5,0.5,0,1), 'cm')),
  plt_map2 + theme(plot.margin = unit(c(0,0,0,0), 'cm')),
  ncol = 1, labels = c('B','D'), rel_heights = rh, align = 'v',
  label_x = lx, label_y = ly, label_fontface = label_fontface,
  label_fontfamily = label_fontfamily)
plt_right <- plot_grid(title_right, plt_right,
                       ncol = 1, rel_heights = c(0.05, 0.95))
plt_leg <- plot_grid(leg_raw, leg_map,
                     ncol = 1, rel_heights = rh, align = 'v')
plt_stress_combined <- plot_grid(plt_left, plt_right, plt_leg,
                                 ncol = 3, rel_widths = c(1,1,0.2)) + 
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))

ggsave('physical_stress_ranking.png', plt_stress_combined, 'png', dir_plots,
       width = pw, height = ph, units = 'in')

# shipping and facilities
plt_raw1 <- plt_ship_raw + theme(legend.position = 'none')
plt_raw2 <- plt_stations_raw + theme(legend.position = 'none')
leg_raw <- ggdraw(get_legend(plt_ship_raw))
plt_map1 <- plt_ship$plot +
  theme(legend.position = 'none') + labs(title = element_blank())
plt_map2 <- pltStations$plot +
  theme(legend.position = 'none') + labs(title = element_blank())
leg_map <- pltStations$legend

title_left <- ggdraw() +
  draw_label('Ship traffic', fontface = label_fontface,
             fontfamily = label_fontfamily, x = 0.55)
title_right <- ggdraw() +
  draw_label('Facilities', fontface = label_fontface,
             fontfamily = label_fontfamily, x = 0.55)
plt_left <- plot_grid(
  plt_raw1 + theme(plot.margin = unit(c(0.5,0.5,0,1), 'cm')),
  plt_map1 + theme(plot.margin = unit(c(0,0,0,0), 'cm')),
  ncol = 1, labels = c('A','C'), rel_heights = rh, align = 'v',
  label_x = lx, label_y = ly, label_fontface = label_fontface,
  label_fontfamily = label_fontfamily)
plt_left <- plot_grid(title_left, plt_left,
                      ncol = 1, rel_heights = c(0.05, 0.95))
plt_right <- plot_grid(
  plt_raw2 + theme(plot.margin = unit(c(0.5,0.5,0,1), 'cm')),
  plt_map2 + theme(plot.margin = unit(c(0,0,0,0), 'cm')),
  ncol = 1, labels = c('B','D'), rel_heights = rh, align = 'v',
  label_x = lx, label_y = ly, label_fontface = label_fontface,
  label_fontfamily = label_fontfamily)
plt_right <- plot_grid(title_right, plt_right,
                       ncol = 1, rel_heights = c(0.05, 0.95))
plt_leg <- plot_grid(leg_raw, leg_map,
                     ncol = 1, rel_heights = rh, align = 'v')
plt_activity_combined <- plot_grid(plt_left, plt_right, plt_leg,
                                   ncol = 3, rel_widths = c(1,1,0.2)) + 
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))

ggsave('human_disturbance_ranking.png', plt_activity_combined, 'png', dir_plots,
       width = pw, height = ph, units = 'in')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge data sets to create risk rankings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method of combining rankings: mean of score, or sum of score
combineMethods <- c('mean', 'sum')
combineMethod <- combineMethods[1]
combineTitle <- paste0(toupper(substr(combineMethod, 1, 1)),
                       substr(combineMethod, 2, nchar(combineMethod)))

# Biota

# plt_biota

chl$type <- 'chl'
krill$type <- 'krill'
columns <- c('cell_index', 'type', 'score')
# Merge the data frames
dat <- rbind(chl[columns], krill[columns])
all_polygons <- sort(unique(dat$cell_index)) # find all mapped polygons
npolys <- length(all_polygons)
score_poly <- switch(
  combineMethod,
  mean = sapply(all_polygons,
                function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
  sum = pmin(nranks,
             sapply(all_polygons,
                    function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
# Rescale combined ranks to 0-5
score_poly <- nranks * {score_poly - min(score_poly)} / diff(range(score_poly))
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_biota <- dat

plt_chl_krill <- plot_fun(dat_biota, Title = combineTitle,
                          legend_title = 'prevalence')
plt_biota <- ggdraw(arrangeGrob(
  plt_chl$plot + theme(legend.position = 'none'),
  plt_krill$plot + theme(legend.position = 'none'),
  plt_chl_krill$plot + theme(legend.position = 'none'),
  plt_chl$legend,
  nrow = 1, widths = c(0.3, 0.3, 0.3, 0.1), top = 'Biota'))
plt_biota

# Stressors

# plt_stress

sst$type <- 'sst'
pH$type <- 'pH'
# Merge the data frames
dat <- rbind(sst[columns], pH[columns])
all_polygons <- sort(unique(dat$cell_index)) # find all mapped polygons
npolys <- length(all_polygons)
score_poly <- switch(
  combineMethod,
  mean = sapply(all_polygons,
                function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
  sum = pmin(nranks,
             sapply(all_polygons,
                    function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
score_poly <- nranks * {score_poly - min(score_poly)} / diff(range(score_poly))
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_stress <- dat

plt_sst_pH <- plot_fun(dat, Title = combineTitle, legend_title = 'stress')
# plt_sst_pH$plot

plt_stress <- ggdraw(arrangeGrob(
  plt_sst$plot + theme(legend.position = 'none'),
  plt_pH$plot + theme(legend.position = 'none'),
  plt_sst_pH$plot + theme(legend.position = 'none'),
  plt_sst$legend,
  nrow = 1, widths = c(0.3, 0.3, 0.3, 0.1), top = 'Stressors'))
plt_stress

# Activity

# plt_activity

ship$type <- 'ship'
grid_cells_reduced$type <- 'station'

# Merge the data frames
dat <- rbind(ship[columns], grid_cells_reduced[columns])
all_polygons <- sort(unique(dat$cell_index)) # find all mapped polygons
npolys <- length(all_polygons)
score_poly <- switch(
  combineMethod,
  mean = sapply(all_polygons,
                function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
  sum = pmin(nranks,
             sapply(all_polygons,
                    function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
score_poly <- nranks * {score_poly - min(score_poly)} / diff(range(score_poly))
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_activity <- dat

plt_ship_station <- plot_fun(dat, Title = combineTitle, legend_title = 'presence')
# plt_ship_station$plot

plt_activity <- ggdraw(arrangeGrob(
  plt_ship$plot + theme(legend.position = 'none'),
  pltStations$plot + theme(legend.position = 'none'),
  plt_ship_station$plot + theme(legend.position = 'none'),
  pltStations$legend,
  nrow = 1, widths = c(0.3, 0.3, 0.3, 0.1), top = 'Activity'))
plt_activity


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Where does activity overlap with biota & stressors?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Biota
plt_biota_activity <- ggdraw(arrangeGrob(
  plt_chl_krill$plot + theme(legend.position = 'none') + labs(title = 'Biota'),
  plt_chl_krill$legend,
  plt_ship_station$plot + theme(legend.position = 'none') + labs(title = 'Activity'),
  plt_ship_station$legend,
  nrow = 1, widths = c(0.4, 0.1, 0.4, 0.1)))
plt_biota_activity


dat_biota$type <- 'biota'
dat_activity$type <- 'activity'
dat <- rbind(dat_biota, dat_activity)

all_polygons <- sort(unique(dat$cell_index))
joint_polygons <- all_polygons[all_polygons %in% dat_biota$cell_index & 
                                 all_polygons %in% dat_activity$cell_index]
npolys <- length(joint_polygons)
score_poly <- switch(
  combineMethod,
  mean = sapply(joint_polygons,
                function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
  sum = pmin(nranks,
             sapply(all_polygons,
                    function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
score_poly <- switch(
  combineMethod, # rescale to 0-5 if needed
  sum = nranks * {score_poly - min(score_poly)} / diff(range(score_poly)),
  mean = score_poly)
dat <- grid_cells_all[grid_cells_all$cell_index %in% joint_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_interaction <- dat

plt_biota_activity <- plot_fun(dat_interaction, Title = combineTitle, legend_title = 'risk')
# plt_biota_activity$plot

plt_interaction <- ggdraw(arrangeGrob(
  plt_chl_krill$plot + labs(title = 'Biota') + theme(legend.position = 'none'),
  plt_chl_krill$legend,
  plt_ship_station$plot + labs(title = 'Activity') + theme(legend.position = 'none'),
  plt_ship_station$legend,
  plt_biota_activity$plot + labs(title = 'Interaction') + theme(legend.position = 'none'),
  plt_biota_activity$legend,
  nrow = 1, widths = c(0.25, 1/12, 0.25, 1/12, 0.25, 1/12), top = 'Interaction with biota'))
plt_interaction


# Stressors
plt_stessor_activity <- ggdraw(arrangeGrob(
  plt_sst_pH$plot + theme(legend.position = 'none') + labs(title = 'Stressors'),
  plt_sst_pH$legend,
  plt_ship_station$plot + theme(legend.position = 'none') + labs(title = 'Activity'),
  plt_ship_station$legend,
  nrow = 1, widths = c(0.4, 0.1, 0.4, 0.1)))
plt_stessor_activity

dat_stress$type <- 'stress'
dat_activity$type <- 'activity'
dat <- rbind(dat_stress, dat_activity)

all_polygons <- sort(unique(dat$cell_index))
joint_polygons <- all_polygons[all_polygons %in% dat_stress$cell_index & 
                                 all_polygons %in% dat_activity$cell_index]
npolys <- length(joint_polygons)
score_poly <- switch(
  combineMethod,
  mean = sapply(joint_polygons,
                function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
  sum = pmin(nranks,
             sapply(all_polygons,
                    function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
score_poly <- switch(
  combineMethod, # rescale to 0-5 if needed
                     sum = nranks * {score_poly - min(score_poly)} / diff(range(score_poly)),
                     mean = score_poly)
dat <- grid_cells_all[grid_cells_all$cell_index %in% joint_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_multistressor <- dat

plt_stress_activity <- plot_fun(dat_multistressor, Title = combineTitle, legend_title = 'risk')
# plt_stress_activity$plot

plt_multistressor <- ggdraw(arrangeGrob(
  plt_sst_pH$plot + labs(title = 'Stressors') + theme(legend.position = 'none'),
  plt_sst_pH$legend,
  plt_ship_station$plot + labs(title = 'Activity') + theme(legend.position = 'none'),
  plt_ship_station$legend,
  plt_stress_activity$plot + labs(title = 'Interaction') + theme(legend.position = 'none'),
  plt_stress_activity$legend,
  nrow = 1, widths = c(0.25, 1/12, 0.25, 1/12, 0.25, 1/12), top = 'Contribution to stressors'))
plt_multistressor


# Total interaction -- activity with biota & stressors
dat_interaction$type <- 'interaction'
dat_multistressor$type <- 'multistressor'
dat <- rbind(dat_interaction, dat_multistressor)
all_polygons <- sort(unique(dat$cell_index))
npolys <- length(all_polygons)
score_poly <- switch(
  combineMethod,
  mean = sapply(all_polygons,
                function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
  sum = pmin(nranks,
             sapply(all_polygons,
                    function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
score_poly <- nranks * {score_poly - min(score_poly)} / diff(range(score_poly)) # rescale to 0-5
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_hotspot <- dat

plt_hotspot <- plot_fun(dat_hotspot, Title = combineTitle, legend_title = 'risk')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine and save these plots for paper
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We need 6 plots: mean ranks for biota, stressors & activity on left column, and
# interaction ranks for biota-activity & stressors-activity, and overall interaction
# rank on the right column.

lx <- 0.05 # label adjustment factors
ly <- 0.9

plt_map1 <- plt_chl_krill$plot + theme(plot.title = element_blank()) +
  guides(fill_new = guide_legend('biota\nprevalence'))
leg1 <- ggdraw(get_legend(plt_map1))
plt_map1 <- plt_map1 + theme(legend.position = 'none')

plt_map2 <- plt_sst_pH$plot + theme(plot.title = element_blank()) +
  guides(fill_new = guide_legend('abiotic\nstress')) 
leg2 <- ggdraw(get_legend(plt_map2))
plt_map2 <- plt_map2 + theme(legend.position = 'none')

plt_map3 <- plt_ship_station$plot + theme(plot.title = element_blank()) +
  guides(fill_new = guide_legend('human\nactivity'))
leg3 <- ggdraw(get_legend(plt_map3))
plt_map3 <- plt_map3 + theme(legend.position = 'none')

plt_map4 <- plt_biota_activity$plot + theme(plot.title = element_blank()) +
  guides(fill_new = guide_legend('interaction\nwith biota'))
leg4 <- ggdraw(get_legend(plt_map4))
plt_map4 <- plt_map4 + theme(legend.position = 'none')

plt_map5 <- plt_stress_activity$plot + theme(plot.title = element_blank()) +
  guides(fill_new = guide_legend('interaction\nwith stressors'))
leg5 <- ggdraw(get_legend(plt_map5))
plt_map5 <- plt_map5 + theme(legend.position = 'none')

plt_map6 <- plt_hotspot$plot + theme(plot.title = element_blank()) +
  guides(fill_new = guide_legend('interaction\noverall'))
leg6 <- ggdraw(get_legend(plt_map6))
plt_map6 <- plt_map6 + theme(legend.position = 'none')

rw <- c(0.8, 0.2)
plt_left <- plot_grid(
  plt_map1,
  plt_map2,
  plt_map3,
  ncol = 1, labels = c('A','C','E'), label_x = lx, label_y = ly,
  label_fontface = label_fontface, label_fontfamily = label_fontfamily)
leg_left <- plot_grid(
  leg1,
  leg2,
  leg3,
  ncol = 1)
plt_left <- plot_grid(
  plt_left, leg_left,
  ncol = 2, rel_widths = rw)

plt_right <- plot_grid(
  plt_map4,
  plt_map5,
  plt_map6,
  ncol = 1, labels = c('B','D','F'), label_x = lx, label_y = ly,
  label_fontface = label_fontface, label_fontfamily = label_fontfamily)
leg_right <- plot_grid(
  leg4,
  leg5,
  leg6,
  ncol = 1)
plt_right <- plot_grid(
  plt_right, leg_right,
  ncol = 2, rel_widths = rw)


plt_results <- plot_grid(
  plt_left, plt_right,
  ncol = 2) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))


A4_w <- 8.3 # A4 dimensions
A4_h <- 11.7
sc_w <- 1.5 # scaling factors
sc_h <- 1.2
pw <- sc_w * A4_w
ph <- sc_h * A4_h
ggsave('main_results_ranking.png', plt_results, 'png', dir_plots, width = pw,
       height = ph, units = 'in')


