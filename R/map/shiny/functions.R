# Functions used in interactive map shiny app

# Utility function to make map grid polygons
createPolygons <- function(dat){
  # Generate list of polygons describing lon/lat bounding box for each row of input
  # data frame. The input data must contain columns called 'lonmin' ,'lonmax', 'latmin',
  # and 'latmax'that define the area corresponding the measurement variable(s).
  # (I think this function may be too complicated -- perhaps we only need call to
  # 'Polygon' function, and not the 'Polygons' function. But this is not important...)
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

# Load and organise data required for map
get_data <- function(res, baseDirectory, shinyDirectory, sstType = 'trend', pHType = 'trend',
                     shipSummaryDataOrRaw = 'summary', sstTrend_significantOnly = TRUE,
                     pHTrend_significantOnly = TRUE, significanceLevel = 0.05, significanceContours = c(0.05, 0.25, 0.5),
                     SST_overallTrend = TRUE, pH_overallTrend = TRUE,
                     roundShipTime = FALSE, indexGridCells = TRUE){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load data with spatial resolution given by res = '3x1' or '9x3'.
  # Required outputs are assigned to the function calling environment.
  # baseDirectory is the project directory where everything is stored.
  # shinyDirectory is where the Shiny app is stored.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load packages
  library(reshape2)
  library(mapdata)
  library(sp)
  library(sf)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  res_map <- setNames(as.numeric(strsplit(res, 'x')[[1]]), c('lon','lat'))
  assign('res_map', res_map, envir = parent.frame())
  assign('sstType', sstType, envir = parent.frame())
  assign('pHType', pHType, envir = parent.frame())
  assign('shipSummaryDataOrRaw', shipSummaryDataOrRaw, envir = parent.frame())
  assign('significanceContours', significanceContours, envir = parent.frame())
  
  assign('SST_overallTrend', SST_overallTrend, envir = parent.frame())
  assign('pH_overallTrend', pH_overallTrend, envir = parent.frame())
  
  # Coordinate reference system
  crs_world <- 4326 # standard longitude-latitude coordinates for world map
  crs_use <- 3031 # Antarctic projection
  assign('crs_world', crs_world, envir = parent.frame())
  assign('crs_use', crs_use, envir = parent.frame())
  

  # Map grid cells ----------------------------------------------------------
  if(indexGridCells){
    
    lonrange <- c(-180, 180)
    latrange <- c(-90, -45)
    longrid <- seq(lonrange[1], lonrange[2], res_map['lon'])
    latgrid <- seq(latrange[1], latrange[2], res_map['lat'])
    nlon <- length(longrid) - 1
    nlat <- length(latgrid) - 1
    lonmin <- head(longrid, -1)
    lonmax <- tail(longrid, -1)
    latmin <- head(latgrid, -1)
    latmax <- tail(latgrid, -1)
    lonmin <- as.numeric(matrix(lonmin, ncol = 1)[,rep(1,nlat)])
    lonmax <- as.numeric(matrix(lonmax, ncol = 1)[,rep(1,nlat)])
    latmin <- as.numeric(matrix(latmin, nrow = 1)[rep(1,nlon),])
    latmax <- as.numeric(matrix(latmax, nrow = 1)[rep(1,nlon),])
    grid_cells_all <- data.frame(lonmin = lonmin, lonmax = lonmax, latmin = latmin, latmax = latmax)
    # Generate spatial polygons
    grid_cells_all <- createPolygons(grid_cells_all)
    grid_cells_all <- SpatialPolygons(grid_cells_all, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
    grid_cells_all <- st_as_sf(grid_cells_all)
    grid_cells_all$cell_index <- 1:nrow(grid_cells_all)
    grid_cell_areas <- st_area(grid_cells_all) # get grid cell areas
    if(all(attr(grid_cell_areas, 'units')$numerator == 'm')) grid_cell_areas <- grid_cell_areas * 1e-6  # convert m^2 to km^2
    grid_cells_all$area <- as.numeric(grid_cell_areas)
    grid_centroids <- st_centroid(grid_cells_all)
    centroid_coords <- st_coordinates(grid_centroids)
    grid_cells_all$centroid_lon <- centroid_coords[,1]
    grid_cells_all$centroid_lat <- centroid_coords[,2]
    grid_centroids_xy <- st_transform(grid_centroids, crs_use)
    centroid_coords <- st_coordinates(grid_centroids_xy)
    grid_cells_all$centroid_x <- centroid_coords[,1]
    grid_cells_all$centroid_y <- centroid_coords[,2]
    
  }
  
  # Load map shape files -----------------------------------------------------
  verbose <- FALSE
  f1 <- 'data/map/Antarctic coastline polygons/medium res/add_coastline_medium_res_polygon_v7_6.shp'
  f2 <- 'data/map/sub Antarctic coastline polygons/sub_antarctic_coastline_high_res_polygon_v1.0.shp'
  f1 <- paste(baseDirectory, f1, sep = '/')
  f2 <- paste(baseDirectory, f2, sep = '/')
  nc1 <- st_read(f1, quiet = !verbose)
  nc2 <- st_read(f2, quiet = !verbose)
  
  # Merge coastline data sets
  nc1$location = 'Antarctica'
  nc2$source <- NULL
  nc2$surface <- 'land'
  nc1 <- nc1[,c('location','surface','geometry')]
  nc2 <- nc2[,c('location','surface','geometry')]
  nc <- rbind(nc1, nc2)
  
  # attr(crs_world, 'class') <- 'crs'
  
  # Reset coastline bounding box to map limits -- defined by most northern latitude
  latmax <- -50
  mapLimits <- data.frame(Longitude = -180:180, 
                          Latitude = rep(latmax, 361))
  mapLimits <- st_as_sf(mapLimits, coords = 1:2, crs = crs_world)
  mapLimits <- st_transform(mapLimits, crs_use)
  attr(st_geometry(nc), 'bbox') <- st_bbox(mapLimits)
  
  nc$surface = factor(nc$surface, levels = unique(nc$surface))
  
  assign('nc', nc, envir = parent.frame())
  
  # Load ecoregion shape files -----------------------------------------------------
  # verbose = TRUE
  f <- 'data/marine ecoregions/data/commondata/data0/meow_ecos_expl_clipped_expl.shp'
  f <- paste(baseDirectory, f, sep = '/')
  eco <- st_read(f, quiet = !verbose)
  
  # Filter ecoregions by province
  keepProvinces <- c(60, 61, 48, 59, 62) # 60=Scotia Sea, 61=Continental High Antarctic, 48=Magellenic, 59=Subantarctic Islands, 62=Subantarctic New Zealand
  eco <- eco[eco$PROV_CODE %in% keepProvinces,]
  
  # Change coordinate reference system
  eco <- st_transform(eco, crs_use)
  
  # Match bounding box to coastline data
  attr(st_geometry(eco), 'bbox') <- st_bbox(nc)
  
  # Crop ecoregions to map limits
  eco <- st_crop(eco, st_bbox(eco))
  
  assign('eco', eco, envir = parent.frame())
  
  
  
  # Load shipping data  -----------------------------------------------------
  
  switch(shipSummaryDataOrRaw,
         summary = {
           
           f <- 'data/shipping/McCarthy_2022/so_ecoregion_combined.csv' # data table with all ship types
           f <- paste(baseDirectory, f, sep = '/')
           ship <- read.table(f, sep = ',', header = TRUE)
           
           # Set time units to days
           ship$total_time <- ship$total_time / 86400 # seconds to days
           ship$mean_time <- ship$mean_time / 86400
           ship$median_time <- ship$median_time / 86400
           
           ship_classes <- unique(ship$ship_class)
           n_ship_classes <- length(ship_classes)
           
           # Match to spatial data on ecoregions
           ship_poly_ <- eco[eco$ECO_CODE_X %in% ship$ecoregion_number,]
           measure_vars <- c("n_voyages", "n_ships", "total_time", "median_time", "mean_time", "n_time" )
           ship_poly <- do.call('rbind',
                                lapply(1:n_ship_classes, function(z){
                                  sc <- ship_classes[z]
                                  d <- ship[ship$ship_class == sc,]
                                  sp <- ship_poly_
                                  sp$ship_class <- sapply(sp$ECO_CODE_X, function(z) d$ship_class[d$ecoregion_number == z])
                                  sp$n_voyages <- sapply(sp$ECO_CODE_X, function(z) d$n_voyages[d$ecoregion_number == z])
                                  sp$n_ships <- sapply(sp$ECO_CODE_X, function(z) d$n_ships[d$ecoregion_number == z])
                                  sp$total_time <- sapply(sp$ECO_CODE_X, function(z) d$total_time[d$ecoregion_number == z])
                                  sp$median_time <- sapply(sp$ECO_CODE_X, function(z) d$median_time[d$ecoregion_number == z])
                                  sp$mean_time <- sapply(sp$ECO_CODE_X, function(z) d$mean_time[d$ecoregion_number == z])
                                  sp$n_time <- sapply(sp$ECO_CODE_X, function(z) d$n_time[d$ecoregion_number == z])
                                  # handle missing values - regions not visited by certain ship types
                                  if(is.list(sp$ship_class)){
                                    j <- sapply(sp$ship_class, function(y) length(y) == 0)
                                    sp$ship_class[j] <- as.vector(rep(sc, sum(j)), 'list')
                                    sp[j,measure_vars] <- NA
                                    for(i in c('ship_class', measure_vars)) sp[[i]] <- unlist(sp[[i]])
                                  }
                                  sp
                                })
           )
           
           
           # # Match to spatial data on ecoregions
           # ship_poly <- eco[eco$ECO_CODE_X %in% ship$Ecoregion.number,]
           # ship_poly$Number.of.visits <- sapply(ship_poly$ECO_CODE_X, function(z) ship$Number.of.visits[ship$Ecoregion.number == z])
           # ship_poly$Mean.time.in.ecoregion <- sapply(ship_poly$ECO_CODE_X, function(z) ship$Mean.time.in.ecoregion..days.[ship$Ecoregion.number == z])
           # ship_poly$time_total <- sapply(ship_poly$ECO_CODE_X, function(z) ship$time_total[ship$Ecoregion.number == z])
           
           
           assign('ship_poly', ship_poly, envir = parent.frame())
           
         },
         
         raw = {
           
           f <- 'data/shipping/McCarthy_2022/ship_time_res_9x3.csv'
           # f <- 'data/shipping/McCarthy_2022/ship_time_res_3x1.csv'
           f <- paste(baseDirectory, f, sep = '/')
           ship <- read.table(f, sep = ',', header = TRUE)
           
           fm <- 'data/shipping/McCarthy_2022/vessel_metadata.csv'
           fm <- paste(baseDirectory, fm, sep = '/')
           ship_metadata <- read.table(fm , sep = ',', header = TRUE)
           
           # ship_yrs <- sort(unique(ship$year))
           # n_ship_yrs <- length(ship_yrs)
           
           # head(ship)
           # head(ship_metadata)
           
           ship <- merge(ship, ship_metadata, by = 'vessel_id')
           
           # For each year and vessel activity, sum ship time over all vessels
           ship <- aggregate(ship_time ~ lonmin + lonmax + latmin + latmax + year + activity_type, data = ship, FUN = sum)

           # Sum ship time over vessel activity
           ship_ <- aggregate(ship_time ~ lonmin + lonmax + latmin + latmax + year, data = ship, FUN = sum)
           ship_$activity_type <- 'all'
           ship_ <- ship_[names(ship)]
           ship <- rbind(ship, ship_)

           # Average ship_time over years
           ship_ <- aggregate(ship_time ~ lonmin + lonmax + latmin + latmax + activity_type, data = ship, FUN = mean)
           ship_$year <- 'all'
           ship_ <- ship_[names(ship)]
           ship <- rbind(ship, ship_)
           
           # Use the annual average only, discard data pertaining to specific years
           ship <- subset(ship, year == 'all')
           
           if(roundShipTime){
             # Round ship_time to nearest day -- this removes grid cells with negligible
             # ship presence, creating a more useful colour scale with greater focus
             # on areas of high traffic
             # ship$ship_time <- round(ship$ship_time)
             # Hmmm, I only want to crop the small measures, so do this instead of
             # rounding all data.
             lowTraffic <- ship$ship_time < 1
             ship$ship_time[lowTraffic] <- round(ship$ship_time[lowTraffic])
           }

           ship <- subset(ship, ship_time != 0) # remove zero measures
           
           ship$activity_type <- factor(ship$activity_type, 
                                        levels = c('all','fishing','tourism','supply','research','other'),
                                        labels = c('all','fishing','tourism','supply','research','other'))

           # Generate spatial polygons
           ship_poly <- createPolygons(ship)
           ship_poly <- SpatialPolygons(ship_poly, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
           ship_poly <- st_as_sf(ship_poly)
           
           if(indexGridCells){
             # Match each ship grid cell to one listed in grid_cells_all
             ship_centroids <- st_centroid(ship_poly)
             distBetweenCells <- st_distance(grid_centroids, ship_centroids)
             ship_poly$cell_index <- grid_cells_all$cell_index[apply(distBetweenCells, 2, which.min)]
           }
           ship_poly <- st_transform(ship_poly, crs = crs_use) # convert coordinates
           ship_poly$year <- ship$year
           ship_poly$activity <- ship$activity_type
           ship_poly$value <- ship$ship_time
           
           # assign('ship_poly', ship_poly, envir = parent.frame())
           
         }
  )
  
  weblink_shipping <- 'https://doi.org/10.1073/pnas.2110303118'
  assign('weblink_shipping', weblink_shipping, envir = parent.frame())
  

  # Load plastic data -----------------------------------------------------
  
  # filepath = 'data/plastic_quantity'
  # # filename = 'plastic_quantity.csv'
  # filename = 'plastic_quantity_new2.csv'
  # f <- paste(wd_base, filepath, filename, sep = '/')
  # DATA = read.csv(f, stringsAsFactors = TRUE)
  
  filename = 'plastic_quantity_new3.csv'
  f <- paste(shinyDirectory, filename, sep = '/')
  DATA = read.csv(f, stringsAsFactors = TRUE)
  
  loadTooltipFromFile <- TRUE
  saveTooltipData <- TRUE
  # Does data with tooltip variable exist?
  fw <- paste0(paste(sub('.csv', '', f), 'tooltip', 'wide', sep = '_'), '.Rds')
  fl <- paste0(paste(sub('.csv', '', f), 'tooltip', 'long', sep = '_'), '.Rds')
  dat_tooltip_saved <- file.exists(fw) & file.exists(fl)
  if(loadTooltipFromFile){
    if(dat_tooltip_saved){
      DATA_wide <- readRDS(fw)
      DATA_long <- readRDS(fl)
    }else{loadTooltipFromFile <- FALSE}
  }
  
  # Data pre-processing ----
  
  # Include spaces in names of data for better display
  DATA_names_orig <- names(DATA)
  DATA_names_spaces <- DATA_names_orig
  DATA_names_spaces[c(2:9, 11:12)] <- c('Sample type','Sample gear','Litter ID method','Litter category','Litter scale','Plastic form','Plastic size','Sample at station','Site category','Sample ID')
  DATA_name_swap <- as.data.frame(cbind(original = DATA_names_orig, spaces = DATA_names_spaces))
  # names(DATA) <- DATA_name_swap$spaces
  assign('DATA_name_swap', DATA_name_swap, envir = parent.frame())
  
    
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
  DATA$Depth[DATA$Depth == 'na'] <- ''
  # DATA$Depth[DATA$Depth == ''] <- NA
  
  dep <- DATA$Depth
  surf <- suppressWarnings(as.numeric(dep) <= 5)
  surf[is.na(surf)] <- FALSE
  surf <- surf | dep %in% c('5m','<1m', '5 m', '<1 m', 'surface')
  subsurf <- suppressWarnings(as.numeric(dep) > 5)
  subsurf[is.na(subsurf)] <- FALSE
  subsurf <- subsurf | dep %in% c('subsurface')
  surf <- marine & surf
  subsurf <- marine & subsurf
  s[surf] <- 'marine surface'
  s[subsurf] <- 'marine subsurface'
  DATA$SampleType <- factor(s, levels = unique(s))
  
  # Make a dummy grouping variable for SampleType
  allSampleTypes <- unique(DATA$SampleType)
  seawater <- allSampleTypes[grepl('sea', allSampleTypes) | grepl('marine', allSampleTypes)]
  freshwater <- allSampleTypes[grepl('fresh', allSampleTypes) | grepl('waste', allSampleTypes)]
  sediment <- allSampleTypes[grepl('sediment', allSampleTypes)]
  air <- allSampleTypes[grepl('air', allSampleTypes)]
  ice <- allSampleTypes[grepl('ice', allSampleTypes)]
  DATA$SampleType_grouped[DATA$SampleType %in% seawater] <- 'seawater'
  DATA$SampleType_grouped[DATA$SampleType %in% freshwater] <- 'freshwater/wastewater'
  DATA$SampleType_grouped[DATA$SampleType %in% sediment] <- 'sediment'
  DATA$SampleType_grouped[DATA$SampleType %in% air] <- 'air'
  DATA$SampleType_grouped[DATA$SampleType %in% ice] <- 'ice'
  DATA$SampleType_grouped <- factor(DATA$SampleType_grouped, levels = c('seawater','freshwater/wastewater','sediment','ice','air'))
  
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
  DATA$Variable <- factor(DATA$Variable, levels = c('concentration', 'density', 'mass concentration', 'mass density', 'presence/absence'))
  
  
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
  DATA$Statistic <- factor(DATA$Statistic, levels = c('raw','mean','s.d.','min','max','none'))
  
  # Set non-numeric Values to -Inf to distinguish them from NAs -- this for presence/absence data
  DATA$Value[is.na(DATA$Value) & DATA$Observation != ''] <- -Inf
  
  #DATA$Observation <- as.character(DATA$Observation)
  
  # Order the data
  DATA <- DATA[order(DATA$Source, DATA$SampleType, DATA$SampleID, DATA$LitterScale, DATA$PlasticSize, 
                     DATA$PlasticForm, DATA$Variable, DATA$Statistic, DATA$Replicate),]
  # DATA <- DATA[order(DATA$Source, DATA$SampleID, DATA$`Litter scale`, DATA$`Plastic size`, 
  #                    DATA$`Plastic form`, DATA$Variable, DATA$Statistic, DATA$Replicate),]
  
  
  
  # Include a data_id variable unique to each sample
  #head(DATA)
  DATA$order <- 1:nrow(DATA)
  d <- unique(DATA[c('Source', 'SampleType', 'SampleAtStation',  'SampleID')])
  d$data_id <- paste('sample', 1:nrow(d))
  
  DATA <- merge(DATA, d, sort = FALSE)
  DATA <- DATA[order(DATA$order),]
  DATA <- DATA[,names(DATA) != 'order']
  
  # for(i in 1:nrow(DATA_name_swap)){
  #   names(DATA)[names(DATA) == DATA_name_swap$original[i]] = DATA_name_swap$spaces[i]
  # }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Pre-allocate/load tooltip
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Interactive tooltip showing plastic sample info as table is a code bottle-neck
  # To avoid slow map rendering, save data with tooltip variable to file.
  
  if(loadTooltipFromFile){
    if(!all(dim(DATA) == dim(DATA_long) - c(0,1))) loadTooltipFromFile <- FALSE
    if(loadTooltipFromFile){
      dataMatched <- as.vector(DATA == DATA_long[,1:ncol(DATA)])
      dataMatched <- all(dataMatched[!is.na(dataMatched)])
      if(!dataMatched) loadTooltipFromFile <- FALSE
    }
  }
  
  if(!loadTooltipFromFile){
    # d <- DATA
    dn <- names(DATA)
    # Change to wide-form -- choose variables to produce one row per data_id
    dummyGroupingVariables <- c('PlasticForm_grouped')
    yvars <- c('LitterScale', 'PlasticSize', 'PlasticForm', 'Variable', 'Statistic', 'Replicate', 'Unit') # variables to change to wide format
    xvars <- dn[!dn %in% c(yvars, 'Value', dummyGroupingVariables)] # everything else, apart from Value and dummy grouping variables, stays in long format
    # Handle spaces in data column names
    xvars_ <- xvars
    xvars_s <- grepl(' ', xvars)
    xvars_[xvars_s] <- paste0('`', xvars_[xvars_s], '`')
    yvars_ <- yvars
    yvars_s <- grepl(' ', yvars)
    yvars_[yvars_s] <- paste0('`', yvars_[yvars_s], '`')
    DATA_wide <- suppressMessages(
      dcast(DATA, list(xvars_, yvars_), value.var = 'Value'))
    wellOrganised <- nrow(DATA_wide) == length(unique(DATA$data_id))
    if(!wellOrganised) warning('Data are not well organised! Some samples are associated with multiple rows of data: we require a single row per sample in the wide-form data set.')
    # Define interactive tooltip using flextable...
    # The yvars (now stretched into wide form) are displayed (long form) in
    # tooltip. Now also choose other variables to display in a single row.
    keepVars <- c('SampleGear', 'Depth', 'Sample date', 'Coordinates', 'Coordinates (start)', 'Coordinates (end)', 'SampleAtStation', 'Observation') # SampleAtStation and Observation are needed to choose appropriate coordinates display
    omitVars <- xvars[!xvars %in% keepVars]
    temp_dat <- as.data.frame(DATA_wide)[!names(DATA_wide) %in% omitVars]
    temp_dat <- split(temp_dat, seq_len(nrow(temp_dat)))
    sType <- as.character(DATA_wide$SampleType) # units depend on SampleType
    sType[grepl('water', sType) | grepl('marine', sType)] <- 'water'
    DATA_wide$tooltip <- sapply(1:length(temp_dat), function(z) fun_flextable(x = temp_dat[[z]], longVars = yvars, singleRowVars = keepVars, sampleType = sType[z], DATA_name_swap)) # this is time-consuming and could maybe be moved to pre-amble/saved data set...
    xvars <- c(xvars, 'tooltip')
    # Melt data back to original shape
    DATA_long <- melt(DATA_wide, id.vars = xvars, value.name = 'Value')#, na.rm = TRUE)
    l <- strsplit(as.character(DATA_long$variable), '_')
    dl <- setNames(as.data.frame(do.call('rbind', l)), yvars)
    DATA_long <- subset(DATA_long, select = - variable)
    DATA_long <- cbind(DATA_long, dl)
    # Remove NAs
    mv <- unique(DATA$data_id[is.na(DATA$Value)]) # samples with missing values in the data
    j <- (DATA_long$data_id %in% mv) # index melted data with missing samples
    ds <- DATA[DATA$data_id %in% mv, c('data_id', yvars)] # discard rows not appearing in original data
    dms <- DATA_long[j, c('data_id', yvars)]
    k <- sapply(1:nrow(dms), function(z){
      m1 <- as.matrix(ds)
      m2 <- as.matrix(dms[z,])[rep(1, nrow(m1)),]
      t <- m1 == m2
      any(apply(t, 1, all))
    })
    j[which(j)[!k]] <- FALSE
    j <- j | !is.na(DATA_long$Value)
    DATA_long <- DATA_long[j,]
    # Merge dummyGroupingVariables
    DATA_long <- merge(DATA_long, cbind(DATA, order = 1:nrow(DATA)))
    # Reorder
    DATA_long <- DATA_long[order(DATA_long$order),]
    DATA_long <- DATA_long[,names(DATA_long) != 'order']
    DATA_long <- DATA_long[,c(dn, 'tooltip')]
    
    if(saveTooltipData){
      fn <- paste0(paste(sub('.csv', '', f), 'tooltip', 'wide', sep = '_'), '.Rds')
      saveRDS(DATA_wide, fn)
      # fp <- paste(filepath, fn, sep = '/')
      # saveRDS(DATA_wide, fp)
      fn <- paste0(paste(sub('.csv', '', f), 'tooltip', 'long', sep = '_'), '.Rds')
      saveRDS(DATA_long, fn)
      # fp <- paste(filepath, fn, sep = '/')
      # saveRDS(DATA_long, fp)
    }
    
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # # Transform to correct CRS
  # dat_plastic_recast <- st_as_sf(dat_plastic_recast, coords = c('Longitude', 'Latitude'), crs = crs_world)
  # dat_plastic_recast <- st_transform(dat_plastic_recast, crs_use)
  # xvars <- c(xvars, 'geometry')
  
  
  # Mapping coordinates
  # DATA_sf = st_as_sf(DATA, coords = c('Longitude', 'Latitude'), crs = crs_world)
  # DATA_sf = st_transform(DATA_sf, crs_use)
  DATA_sf = st_as_sf(DATA_wide, coords = c('Longitude', 'Latitude'), crs = crs_world)
  DATA_sf = st_transform(DATA_sf, crs_use)
  
  assign('DATA_sf', DATA_sf, envir = parent.frame())
  assign('DATA_long', DATA_long, envir = parent.frame())
  
  # Load research station data ----------------------------------------------
  filename <- 'COMNAP_Antarctic_Facilities_Master.csv'
  filepath <- 'data/research stations/COMNAP'
  f <- paste(baseDirectory, filepath, filename, sep = '/')
  STATIONS <- read.csv(f, stringsAsFactors = TRUE)
  
  weblink_facilities <- 'https://github.com/PolarGeospatialCenter/comnap-antarctic-facilities'
  assign('weblink_facilities', weblink_facilities, envir = parent.frame())

  # Data pre-processing ----
  STATIONS$Type = reorder(STATIONS$Type, rep(1, nrow(STATIONS)), sum, decreasing = TRUE)
  
  # Define tooltip variable for interactive info...
  
  # Attribute values cannot contain single quotes "'", so remove them... I want to
  # find a way to include these...
  STATIONS$Official_Name <- gsub("'", "", STATIONS$Official_Name)
  stationTypes <- levels(STATIONS$Type)
  # STATIONS$Type <- as.character(STATIONS$Type)
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
  
  # STATIONS$tooltip <- paste0(STATIONS$Official_Name, ' ', STATIONS$Type, '\n ',
  #                            'Primary operator: ', STATIONS$Operator_primary, '\n ',
  #                            'Peak population size = ', STATIONS$Peak_Population)
  
  STATIONS$Record_ID <- paste('facility', STATIONS$Record_ID)
  
  # Mapping coordinates
  STATIONS_sf = st_as_sf(STATIONS, coords = c('Longitude_DD', 'Latitude_DD'), crs = crs_world)
  
  
  # Match each station to a grid cell listed in grid_cells_all
  # station_centroids <- st_centroid(STATIONS_sf)
  # distBetweenCells <- st_distance(grid_centroids, STATIONS_sf)
  # STATIONS_sf$cell_index <- grid_cells_all$cell_index[apply(distBetweenCells, 2, which.min)]

  STATIONS_sf = st_transform(STATIONS_sf, crs_use)
  xy = matrix(unlist(STATIONS_sf$geometry), 2, nrow(STATIONS))
  STATIONS$x = xy[1,]
  STATIONS$y = xy[2,]
  
  # Calculate each facility's distance from the coast. The resulting distance 
  # estimates are quite crude (mostly) because the coastline data is not very
  # detailed (despite being called hi-res!)
  coastline <- map('worldHires', 'antarctica', plot = FALSE)
  coastPoints <- data.frame(X = coastline$x, Y = coastline$y)
  coastPoints <- coastPoints[!{is.na(coastPoints$X) | is.na(coastPoints$Y)},]
  
  # Use haversine formula for distance calculation
  hav <- function(theta){
    y <- 0.5 * {1 - cos(theta)}
    y[y > 1] <- 1 # control for floating point error
    y}
  greatCircDist <- function(lon1, lon2, lat1, lat2, r){
    dlon <- lon2 - lon1
    dlat <- lat2 - lat1
    2 * r * asin({hav(dlat) + cos(lat1)*cos(lat2)*hav(dlon)} ^ 0.5)}
  
  # Convert degrees into radians
  coastPoints <- pi / 180 * coastPoints
  station_coords <- STATIONS[c('Longitude_DD','Latitude_DD')]
  station_coords[c('Longitude_rad','Latitude_rad')] <- pi / 180 * station_coords # convert degrees to radians
  
  STATIONS_sf$Distance_From_Coast <- sapply(1:nrow(station_coords), function(z)
    min(greatCircDist(coastPoints[,1], station_coords$Longitude_rad[z], coastPoints[,2], station_coords$Latitude_rad[z], r = 6371))
  )
  
  
  # There is something wrong with this distance calculation... need to figure it out... maybe the polygon points are not all coastline...
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  STATIONS_sf$tooltip <- paste0(STATIONS_sf$Official_Name, ' ', STATIONS_sf$Type, '\n ',
                                'Primary operator: ', STATIONS_sf$Operator_primary, '\n ',
                                'Peak population size = ', STATIONS_sf$Peak_Population, '\n ',
                                'Distance from coast = ', round(STATIONS_sf$Distance_From_Coast, 1), ' km')
  
  assign('STATIONS_sf', STATIONS_sf, envir = parent.frame())
  
  
  # Load krill data ---------------------------------------------------------
  filename = 'krill_data_mapped.csv'
  filepath = 'MatLab/temp'
  f <- paste(baseDirectory, filepath, filename, sep = '/')
  KRILL = read.csv(f, stringsAsFactors = TRUE)
  
  weblink_krill <- 'https://www.doi.org/10.5285/8b00a915-94e3-4a04-a903-dd4956346439'
  assign('weblink_krill', weblink_krill, envir = parent.frame())

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
  
  # # I think this function is too complicated -- maybe only need call to Polygon function, and not the Polygons function...
  # createPolygons <- function(dat){
  #   # Returns list of polygons describing lon/lat bounding box for each row of dat
  #   n = nrow(dat)
  #   lapply(1:n, FUN = function(z){
  #     x = dat[z,]
  #     Polygons(
  #       list(
  #         Polygon(
  #           cbind(
  #             c(x$lonmin, x$lonmax, x$lonmax, x$lonmin, x$lonmin),
  #             c(x$latmin, x$latmin, x$latmax, x$latmax, x$latmin)
  #           ))), paste0('row', z)
  #     )}
  #   )
  # }
  
  krill_poly <- createPolygons(KRILL) # generate polygons
  krill_poly <- SpatialPolygons(krill_poly, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
  krill_poly <- st_as_sf(krill_poly)
  if(indexGridCells){
    # Match polygons to those listed in grid_cells_all
    krill_centroids <- st_centroid(krill_poly)
    distBetweenCells <- st_distance(grid_centroids, krill_centroids)
    krill_poly$cell_index <- grid_cells_all$cell_index[apply(distBetweenCells, 2, which.min)]
  }
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
  
  # assign('krill_poly', krill_poly, envir = parent.frame())
  
  
  # Load chlorophyll data ---------------------------------------------------------
  
  filename <- paste0('chl_data_mapped_res_', res, '.csv')
  filepath <- 'MatLab/temp'
  f <- paste(baseDirectory, filepath, filename, sep = '/')
  
  CHL <- read.csv(f, stringsAsFactors = TRUE)
  
  weblink_chlorophyll <- 'https://catalogue.ceda.ac.uk/uuid/a91ae89c4c6d011f84cdfbc3d41a73c8'
  assign('weblink_chlorophyll', weblink_chlorophyll, envir = parent.frame())
  
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
  if(indexGridCells){
    # Match polygons to those listed in grid_cells_all
    chl_centroids <- st_centroid(chl_poly)
    distBetweenCells <- st_distance(grid_centroids, chl_centroids)
    chl_poly$cell_index <- grid_cells_all$cell_index[apply(distBetweenCells, 2, which.min)]
  }
  chl_poly <- st_transform(chl_poly, crs = crs_use) # convert coordinates
  chl_poly$month <- CHL$month
  chl_poly$value <- CHL$value
  
  # assign('chl_poly', chl_poly, envir = parent.frame())
  
  
  # Load SST trend (or anomaly) data --------------------------------------
  
  # sstType = 'trend' # 'trend' or 'anomaly
  # res <- '3x1' # choose data resolution: either  '3x1' or '9x3'
  
  if(SST_overallTrend){
    filename <- paste0('sst_time_series_trend_res_', res, '.csv') # my month-by-month trend calculation    
  }else{
    filename <- paste0('sst_', sstType, '_res_', res, '.csv') # my month-by-month trend calculation
  }
  filepath <- 'MatLab/temp'
  f <- paste(baseDirectory, filepath, filename, sep = '/')
  
  SST <- read.csv(f)
  
  weblink_sst <- 'https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/download?dataset=ESACCI-GLO-SST-L4-REP-OBS-SST'
  assign('weblink_sst', weblink_sst, envir = parent.frame())
  
  # Data pre-processing -----------------------------------------------------     
  
  if(SST_overallTrend){

    coords <- c('lonmin','lonmax','latmin','latmax')
    SST <- cbind(SST[,coords], metric = 'trend', value = SST$linearTrend, pval = SST$pValue)
    
    if(sstTrend_significantOnly){
      # Set to zero any non-significant linear trends
      insignificant <- SST$pval > significanceLevel
      insignificant[is.na(insignificant)] <- FALSE
      SST$value[insignificant] <- 0
    }
    
    SST <- subset(SST, !is.nan(value))  # omit missing data
    
    # To make it easier to plot according to p-value, or to somehow highlight
    # p-values, define p-level as a category variable
    nsiglevels <- length(significanceContours) + 1
    SST$plevel <- apply(outer(SST$pval, significanceContours, '>'), 1, sum) + 1
    siglabels <- c(paste('p <=', significanceContours), paste('p >', tail(significanceContours, 1)))
    SST$plevel <- factor(SST$plevel, levels = 1:nsiglevels, labels = siglabels)
    
  }else{
    
    if(sstType == 'anomaly'){
      sst_ymax <- max(SST$year)
      SST <- subset(SST, year == sst_ymax, select = -c(year, value))
      names(SST)[ncol(SST)] = 'value'
      SST <- subset(SST, !is.nan(value))  # omit missing data
      mon <- 1:3
      SST <- subset(SST, month %in% mon)
      # Mean SST anomaly across selected (summer) months
      SST_allMonths <- aggregate(value ~ lonmin + lonmax + latmin + latmax, SST, 'mean')
      SST_allMonths$month <- 'all'
      SST_allMonths <- SST_allMonths[,c(1:4,6,5)]
      SST <- rbind(SST, SST_allMonths)
    }
    
    if(sstType == 'trend'){
      coords <- c('lonmin','lonmax','latmin','latmax')
      pvals <- SST[,c(names(SST) %in% coords | grepl('pValue', names(SST)))]
      names(pvals) <- gsub('pValue_', '', names(pvals))
      SST <- SST[,!grepl('pValue', names(SST))]
      SST <- melt(SST, id.vars = coords, variable.name = 'month', value.name = 'value')
      pvals <- melt(pvals, id.vars = coords, variable.name = 'month', value.name = 'value')
      
      SST$metric <- sstType
      SST$pval <- pvals$value
      
      if(sstTrend_significantOnly){
        # Set to zero any non-significant linear trends
        insignificant <- SST$pval > significanceLevel
        insignificant[is.na(insignificant)] <- FALSE
        SST$value[insignificant] <- 0
      }
      
      SST <- subset(SST, !is.nan(value))  # omit missing data
      
      # To make it easier to plot according to p-value, or to somehow highlight
      # p-values, define p-level as a category variable
      nsiglevels <- length(significanceContours) + 1
      SST$plevel <- apply(outer(SST$pval, significanceContours, '>'), 1, sum) + 1
      siglabels <- c(paste('p <=', significanceContours), paste('p >', tail(significanceContours, 1)))
      SST$plevel <- factor(SST$plevel, levels = 1:nsiglevels, labels = siglabels)
      
      # Include summer months only
      mon <- month.abb[1:3]
      SST <- subset(SST, month %in% mon | month == 'Summer')
      SST$month <- as.character(SST$month)
      SST$month[SST$month == 'Summer'] <- 'all'
      for(i in 1:length(mon)) SST$month[SST$month == mon[i]] <- i
      SST$month <- as.factor(SST$month)
    }
    
  }

  # Generate spatial polygons
  sst_poly <- createPolygons(SST)
  sst_poly <- SpatialPolygons(sst_poly, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
  sst_poly <- st_as_sf(sst_poly)
  
  if(indexGridCells){
    # Match polygons to those listed in grid_cells_all
    sst_centroids <- st_centroid(sst_poly)
    distBetweenCells <- st_distance(grid_centroids, sst_centroids)
    sst_poly$cell_index <- grid_cells_all$cell_index[apply(distBetweenCells, 2, which.min)]
  }
  sst_poly <- st_transform(sst_poly, crs = crs_use) # convert coordinates
  sst_poly$month <- SST$month
  sst_poly$metric <- SST$metric
  sst_poly$value <- SST$value
  
  # Get p-values for linear trends if returning all trends regardless of their
  # statistical significance
  if(sstType != 'anomaly'){
    sst_poly$pval <- SST$pval
    sst_poly$plevel <- SST$plevel
  }
  
  # assign('sst_poly', sst_poly, envir = parent.frame())
  
  # Load pH trend (or anomaly) data -----------------------------------------
  
  # pHType = 'trend' # 'trend' or 'anomaly
  # res <- '3x1' # choose data resolution: either  '3x1' or '9x3'

  if(pH_overallTrend){
    filename <- paste0('pH_time_series_trend_res_', res, '.csv') # my month-by-month trend calculation    
  }else{
    filename <- paste0('pH_', pHType, '_res_', res, '.csv') # my month-by-month trend calculation
  }
  filepath <- 'MatLab/temp'
  f <- paste(baseDirectory, filepath, filename, sep = '/')
  
  pH <- read.csv(f)
  
  weblink_pH <- 'https://data.marine.copernicus.eu/product/MULTIOBS_GLO_BIO_CARBON_SURFACE_REP_015_008/description'
  assign('weblink_pH', weblink_pH, envir = parent.frame())
  
  
  # Data pre-processing -----------------------------------------------------     
  
  if(pH_overallTrend){
    
    coords <- c('lonmin','lonmax','latmin','latmax')
    pH <- cbind(pH[,coords], metric = 'trend', value = pH$linearTrend, pval = pH$pValue)
    
    if(pHTrend_significantOnly){
      # Set to zero any non-significant linear trends
      insignificant <- pH$pval > significanceLevel
      insignificant[is.na(insignificant)] <- FALSE
      pH$value[insignificant] <- 0
    }
    
    pH <- subset(pH, !is.nan(value))  # omit missing data
    
    # To make it easier to plot according to p-value, or to somehow highlight
    # p-values, define p-level as a category variable
    nsiglevels <- length(significanceContours) + 1
    pH$plevel <- apply(outer(pH$pval, significanceContours, '>'), 1, sum) + 1
    siglabels <- c(paste('p <=', significanceContours), paste('p >', tail(significanceContours, 1)))
    pH$plevel <- factor(pH$plevel, levels = 1:nsiglevels, labels = siglabels)
    
  }else{
    
    if(pHType == 'anomaly'){
      pH_ymax <- max(pH$year)
      pH <- subset(pH, year == pH_ymax, select = -c(year, pH, pH_sd))
      names(pH)[ncol(pH)] = 'value'
      pH <- subset(pH, !is.nan(value))  # omit missing data
      mon <- 1:3
      pH <- subset(pH, month %in% mon)
      # Mean pH anomaly across selected (summer) months
      pH_allMonths <- aggregate(value ~ lonmin + lonmax + latmin + latmax, pH, 'mean')
      pH_allMonths$month <- 'all'
      pH_allMonths <- pH_allMonths[,c(1:4,6,5)]
      pH <- rbind(pH, pH_allMonths)
      pH$month <- as.factor(pH$month)
    }
    
    if(pHType == 'trend'){
      coords <- c('lonmin','lonmax','latmin','latmax')
      pvals <- pH[,c(names(pH) %in% coords | grepl('pValue', names(pH)))]
      names(pvals) <- gsub('pValue_', '', names(pvals))
      pH <- pH[,!grepl('pValue', names(pH))]
      pH <- melt(pH, id.vars = coords, variable.name = 'month', value.name = 'value')
      pvals <- melt(pvals, id.vars = coords, variable.name = 'month', value.name = 'value')
      
      pH$metric <- pHType
      pH$pval <- pvals$value
      
      # A small number of grid cells have p-value = NA due to few (only 2?) measurements -- set these p-values to 1
      pH$pval[!is.na(pH$value) & is.na(pH$pval)] <- 1
      
      if(pHTrend_significantOnly){
        # Set to zero any non-significant linear trends
        insignificant <- pH$pval > significanceLevel
        insignificant[is.na(insignificant)] <- FALSE
        pH$value[insignificant] <- 0
      }
      
      pH <- subset(pH, !is.nan(value))  # omit missing data
      
      # To make it easier to plot according to p-value, or to somehow highlight
      # p-values, define p-level as a category variable
      nsiglevels <- length(significanceContours) + 1
      pH$plevel <- apply(outer(pH$pval, significanceContours, '>'), 1, sum) + 1
      siglabels <- c(paste('p <=', significanceContours), paste('p >', tail(significanceContours, 1)))
      pH$plevel <- factor(pH$plevel, levels = 1:nsiglevels, labels = siglabels)
      
      # Include summer months only
      mon <- month.abb[1:3]
      pH <- subset(pH, month %in% mon | month == 'Summer')
      pH$month <- as.character(pH$month)
      pH$month[pH$month == 'Summer'] <- 'all'
      for(i in 1:length(mon)) pH$month[pH$month == mon[i]] <- i
      pH$month <- as.factor(pH$month)
    }
    
  }
  
  # Generate spatial polygons
  pH_poly <- createPolygons(pH)
  pH_poly <- SpatialPolygons(pH_poly, proj4string = CRS(paste0('+init=epsg:', as.character(crs_world)))) # create spatial object
  pH_poly <- st_as_sf(pH_poly)
  
  if(indexGridCells){
    # Match polygons to those listed in grid_cells_all
    pH_centroids <- st_centroid(pH_poly)
    distBetweenCells <- st_distance(grid_centroids, pH_centroids)
    pH_poly$cell_index <- grid_cells_all$cell_index[apply(distBetweenCells, 2, which.min)]
  }
  pH_poly <- st_transform(pH_poly, crs = crs_use) # convert coordinates
  pH_poly$month <- pH$month
  pH_poly$metric <- pH$metric
  pH_poly$value <- pH$value
  
  # Get p-values for linear trends if returning all trends regardless of their
  # statistical significance
  if(pHType != 'anomaly'){
    pH_poly$pval <- pH$pval
    pH_poly$plevel <- pH$plevel
  }
  
  # assign('pH_poly', pH_poly, envir = parent.frame())
  
  # Tidy up indexing --------------------------------------------------------

  if(indexGridCells){
    
    gridCellsWithData <- unique(c(
      ship_poly$cell_index,
      # STATIONS_sf$cell_index,
      chl_poly$cell_index,
      krill_poly$cell_index,
      sst_poly$cell_index,
      pH_poly$cell_index
    ))
    
    i <- grid_cells_all$cell_index %in% gridCellsWithData
    grid_cells_all <- grid_cells_all[i,]
    relabelCells <- data.frame(old = grid_cells_all$cell_index)
    relabelCells$new <- 1:nrow(relabelCells)
    
    relabelFun <- function(dat) sapply(1:nrow(dat), function(z) relabelCells$new[relabelCells$old == dat$cell_index[z]])
    ship_poly$cell_index <- relabelFun(ship_poly)
    # STATIONS_sf$cell_index <- relabelFun(STATIONS_sf)
    chl_poly$cell_index <- relabelFun(chl_poly)
    krill_poly$cell_index <- relabelFun(krill_poly)
    sst_poly$cell_index <- relabelFun(sst_poly)
    pH_poly$cell_index <- relabelFun(pH_poly)
    grid_cells_all$cell_index <- relabelFun(grid_cells_all)
    
    grid_cells_all <- st_transform(grid_cells_all, crs = crs_use) # convert coordinates
    
    assign('grid_cells_all', grid_cells_all, envir = parent.frame())
    
    plt_map_cells <-
      ggplot(grid_cells_all) +
      geom_sf(data = nc,
              aes(fill = surface),
              show.legend = FALSE) +
      scale_fill_manual(values = c('grey','skyblue','skyblue','grey')) +
      geom_sf(alpha = 0) +
      geom_sf_text(aes(label = cell_index), size = 3) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.background = element_rect(fill = 'white', colour = 'white')
      )
    
    assign('plt_map_cells', plt_map_cells, envir = parent.frame())
    
  }
  
  assign('ship_poly', ship_poly, envir = parent.frame())
  assign('krill_poly', krill_poly, envir = parent.frame())
  assign('chl_poly', chl_poly, envir = parent.frame())
  assign('sst_poly', sst_poly, envir = parent.frame())
  assign('pH_poly', pH_poly, envir = parent.frame())

}

# Define interactive tooltip for plastic data
fun_flextable <- function(x, longVars, singleRowVars, sampleType, DATA_name_swap){
  # Function defining tables to interactively display when mouse hovers over mapped (plastic) data points.
  # Inputs: x = data, pu = plastic_units, singleRowVars = variables only taking one row of table
  # (Some nice examples, worth a look, presented here https://davidgohel.github.io/ggiraph/articles/offcran/shiny.html)
  anyMeasures <- any(!is.na(x[,!names(x) %in% singleRowVars]))
  if(anyMeasures){
    # Arrange data frame
    x <-suppressMessages(melt(x, id.vars = singleRowVars)) # change to long-form
    x <- x[!is.na(x$value),]
    # Disentangle variables...
    l <- strsplit(as.character(x$variable), '_')
    # y <- as.data.frame(do.call('rbind', l))
    # names(y) <- longVars[1:ncol(y)]
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
      x$Value[x$Statistic == 's.d.'] <- v
      x$Statistic[x$Statistic == 's.d.'] <- 'mean +- s.d.'
      x <- x[x$Statistic != 'mean',]
    }
    if(any(x$Statistic == 'none')){
      longVars <- longVars[longVars != 'Statistic']}
    # Sample replicates
    if(length(unique(x$Replicate)) == 1){
      x <- x[,-which(names(x) == 'Replicate')]
      longVars <- longVars[longVars != 'Replicate']}
    # Depth
    if(x$Depth[1] == ''){
      singleRowVars <- singleRowVars[!singleRowVars == 'Depth']
    }else{
      # Include unit if Depth is numeric
      depth <- suppressWarnings(as.numeric(x$Depth))
      if(!any(is.na(depth)))
        x$Depth <- paste(depth, 'm')
    }
    x$Measure <- x$Variable # MAYBE RENAME THE Variable COLUMN AS Measure IN THE MAIN DATA SET...
    n <- nrow(x)
    if(x$SampleAtStation[1]){
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates (start)', 'Coordinates (end)', 'SampleAtStation')]
    }else{
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates', 'SampleAtStation')]}
    if(any(x$Unit == 'none')){
      longVars <- longVars[longVars != 'Unit']}
    if(any(x$Value == -Inf)){
      x$Value <- x$Observation
    }else{
      x$Value <- paste(x$Value, x$Unit)}
    singleRowVars <- singleRowVars[singleRowVars != 'Observation']
    
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
    # Create flextable
    ft <- flextable(x, col_keys = c(singleRowVars, longVars, 'Value'))
    ft <- bold(ft, part = 'header', bold = TRUE)
    ft <- set_table_properties(ft, layout = 'autofit')
    as.character(htmltools_value(ft, ft.shadow = FALSE))
    
  }else{
    
    if(x$SampleAtStation[1]){
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates (start)', 'Coordinates (end)', 'SampleAtStation')]
    }else{
      singleRowVars <- singleRowVars[!singleRowVars %in% c('Coordinates', 'SampleAtStation')]}
    if(!any(x$Observation != '')) singleRowVars <- singleRowVars[singleRowVars != 'Observation']
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

# Set the plotting parameters
set_plot_params <- function(base_map, plastics, stations){
  
  # Plot symbols - store in separate data frame
  # Research stations use (mostly) unfilled symbols
  pltShapes <- c(1, 0, 5, 2, 6, 19)
  stationTypes <- levels(stations$Type)
  nstationTypes <- length(stationTypes)
  pltSymbols <- data.frame(Class = rep('ResearchStation', nstationTypes), Type = stationTypes, symbol = pltShapes[1:nstationTypes])
  # Plastic samples use filled plot symbols (21:25) that differ according to sample type.
  pltShapes <- 21:25
  sampleTypes <- levels(plastics$SampleType_grouped)
  nsampleTypes <- length(sampleTypes)
  pltSymbols <- rbind(pltSymbols,
                      data.frame(Class = rep('PlasticSample', nsampleTypes), Type = sampleTypes, symbol = pltShapes[1:nsampleTypes]))
  pltSymbols$Type <- factor(pltSymbols$Type, levels = pltSymbols$Type)
  
  # Plot colours
  # Use a qualitative palette for plastic sources -- the default is limited to 12 colours
  ncol <- length(unique(plastics$Source))
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
  
  sources <- levels(plastics$Source)
  nsources <- length(sources)
  pltColours <- data.frame(Source = sources, colour = pltColours[1:nsources])
  
  du <- unique(data.frame(Source = plastics$Source, URL = plastics$URL))
  pltColours <- merge(pltColours, du, by = 'Source')
  
  # Plot bounding box
  bbox_map <- st_bbox(base_map)
  bbox_dat <- st_bbox(plastics)
  # bbox_krill <- st_bbox(krill_poly)
  # bbox_chl <- st_bbox(chl_poly)
  BBox <- bbox_map
  aspectRatio = unname(diff(BBox[c(1,3)]) / diff(BBox[c(2,4)]))
  linebreaks <- function(n){HTML(strrep(br(), n))} # convenience function
  
  assign('pltSymbols', pltSymbols, envir = parent.frame())
  assign('pltColours', pltColours, envir = parent.frame())
  assign('BBox', BBox, envir = parent.frame())
  assign('aspectRatio', aspectRatio, envir = parent.frame())
  assign('linebreaks', linebreaks, envir = parent.frame())
  
  
}


# The main plotting function
make_plot <- function(dat, background = 'none', displayEcoregions = FALSE, backgroundOnly = FALSE, plasticOnly = FALSE, ptSize = 6, legPtSize = 4, alpha = 0.6, polyLineWidth = 0.75, contourLineWidth = 0.75, legWidth = 0.3, mapAspectRatio = aspectRatio, plotSignificanceContours = FALSE, significanceContours = NULL, sst_sigContours = significanceContours, pH_sigContours = significanceContours){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to generate map plot, called from inside the server function after
  # data have been filtered by user selected inputs.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Inputs:
  # dat               = list of all required (reactively filtered) data sets
  # background        = character specifying which background layer to plot
  # displayEcoregions = TRUE/FALSE: plot distinct ecoregions -- required for shipping data from McCarthy
  # alpha             = point transparency
  # Output:
  # A ggplot object containing interactive geoms compatible with girafe
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nc <- dat$nc
  dat_background <- dat$background
  dat_plastic <- dat$plastic
  dat_stations <- dat$stations
  symbols <- dat$symbols
  # pltColours <- dat$pltColours
  
  anyPlastic <- nrow(dat_plastic) > 0
  anyStations <- nrow(dat_stations) > 0
  anyBackground <- background != 'none'
  
  #~~~~~~~~~~~~~~~~~
  # Background layer
  #~~~~~~~~~~~~~~~~~
  switch(background,
         
         none = {
           plt_background <- ggplot()
         },
         
         krill = {
           legTitle <- 'Krill\n(number/m\u00B2)'
           if(exists('weblink_krill')){
             plt_background <- 
               ggplot() +
               geom_sf(data = dat_background, aes(fill = colourgroup)) +
               scale_fill_viridis_d_interactive(option = 'plasma',
                                                name = label_interactive(legTitle,
                                                                         data_id = 'legend.title',
                                                                         hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                         onclick = paste0("window.open(`", as.character(weblink_krill), "`);"),
                                                                         tooltip = weblink_krill))
           }else{
             plt_background <- 
               ggplot() +
               geom_sf(data = dat_background, aes(fill = colourgroup)) +
               scale_fill_viridis_d(option = 'plasma', name = legTitle)
           }
           
         },
         
         chl = {
           legTitle <- 'Chlorophyll\n(mg/m\u00B3)'
           # legTitle <- 'Chlorophyll\u1D44E\n(mg/m\u00B3)'
           # legTitle <- 'Chlorophyll 𝑎\n(mg/m\u00B3)'

           if(exists('weblink_chlorophyll')){
             plt_background <-
               ggplot() +
               geom_sf(data = dat_background, aes(fill = value)) +
               scale_fill_viridis_c_interactive(option = 'viridis', trans = 'log10',
                                                name = label_interactive(legTitle,
                                                                         data_id = 'legend.title',
                                                                         hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                         onclick = paste0("window.open(`", as.character(weblink_chlorophyll), "`);"),
                                                                         tooltip = weblink_chlorophyll))
           }else{
             plt_background <-
               ggplot() +
               geom_sf(data = dat_background, aes(fill = value)) +
               scale_fill_viridis_c(option = 'viridis', trans = 'log10',
                                    name = legTitle)
           }
           
         },
         
         sst = {
           if(exists('sstType')){
             legLabel <- switch(sstType,
                                trend = 'SST trend\n(\u00B0C/year)',# 'SST trend\n(\u2103/year)',
                                anomaly = 'SST anomaly\n(\u00B0C)')# 'SST anomaly\n(\u2103)')
           }else{
             legLabel <- 'SST'
           }
           weblink_exists <- exists('weblink_sst')
           nv <- dat_background$value[dat_background$value < 0]
           pv <- dat_background$value[dat_background$value > 0]
           any_positive <- length(pv) > 0
           any_negative <- length(nv) > 0
           ncol_sst <- 11
           ncol_sst_ <- 0.5 * {ncol_sst - 1}
           colfunc <- colorRampPalette(c('dodgerblue4', 'white', 'red3'))
           Cols <- colfunc(ncol_sst)
           cval <- c(quantile(nv, seq(0, 1, length = ncol_sst_ + 1))[1:ncol_sst_], 0, quantile(pv, seq(0, 1, length = ncol_sst_ + 1))[2:{ncol_sst_ + 1}])
           if(any_positive & any_negative){
             cval <- {cval - min(cval)} / diff(range(cval))
             if(weblink_exists){
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn_interactive(colours = Cols, values = cval,
                                                  name = label_interactive(legLabel,
                                                                           data_id = 'legend.title',
                                                                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                           onclick = paste0("window.open(`", as.character(weblink_sst), "`);"),
                                                                           tooltip = weblink_sst))
             }else{
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn(colours = Cols, values = cval,
                                      name = legLabel)
             }
           }else{if(any_negative){
             cval <- {cval[1:{ncol_sst_+1}] - min(cval[1:{ncol_sst_+1}])} / diff(range(cval[1:{ncol_sst_+1}]))
             cval <- {cval - min(cval)} / diff(range(cval))
             if(weblink_exists){
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn_interactive(colours = Cols[1:{ncol_sst_+1}], values = cval,
                                                  name = label_interactive(legLabel,
                                                                           data_id = 'legend.title',
                                                                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                           onclick = paste0("window.open(`", as.character(weblink_sst), "`);"),
                                                                           tooltip = weblink_sst))
             }else{
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn(colours = Cols[1:{ncol_sst_+1}], values = cval,
                                      name = legLabel)
             }
           }else{if(any_positive){
             cval <- {cval[{ncol_sst_+1}:ncol_sst] - min(cval[{ncol_sst_+1}:ncol_sst])} / diff(range(cval[{ncol_sst_+1}:ncol_sst]))
             cval <- {cval - min(cval)} / diff(range(cval))
             if(weblink_exists){
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn_interactive(colours = Cols[{ncol_sst_+1}:ncol_sst], values = cval,
                                                  name = label_interactive(legLabel,
                                                                           data_id = 'legend.title',
                                                                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                           onclick = paste0("window.open(`", as.character(weblink_sst), "`);"),
                                                                           tooltip = weblink_sst))
             }else{
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn(colours = Cols[{ncol_sst_+1}:ncol_sst], values = cval,
                                      name = legLabel)
             }
           }}}
         },
         
         pH = {
           
           if(exists('pHtype')){
             legLabel <- switch(pHType,
                                trend = 'pH trend\n(year\u207B\u00B9)',# bquote(atop(pH ~ trend, (1/year))),
                                anomaly = 'pH anomaly\n(dimensionless)')# bquote(pH ~ anomaly))
           }else{
             legLabel <- 'pH'# bquote(pH)
           }
           weblink_exists <- exists('weblink_pH')
           nv <- dat_background$value[dat_background$value < 0]
           pv <- dat_background$value[dat_background$value > 0]
           any_positive <- length(pv) > 0
           any_negative <- length(nv) > 0
           ncol_pH <- 11
           ncol_pH_ <- 0.5 * {ncol_pH - 1}
           colfunc <- colorRampPalette(c('darkorange2', 'white', 'darkorchid4'))
           Cols <- colfunc(ncol_pH)
           cval <- c(quantile(nv, seq(0, 1, length = ncol_pH_ + 1))[1:ncol_pH_], 0, quantile(pv, seq(0, 1, length = ncol_pH_ + 1))[2:{ncol_pH_ + 1}])
           # Colour scale selection depends on range of values
           if(any_positive & any_negative){
             cval <- {cval - min(cval)} / diff(range(cval))
             if(weblink_exists){
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn_interactive(colours = Cols, values = cval,
                                                  name = label_interactive(legLabel,
                                                                           data_id = 'legend.title',
                                                                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                           onclick = paste0("window.open(`", as.character(weblink_pH), "`);"),
                                                                           tooltip = weblink_pH))
             }else{
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn(colours = Cols, values = cval,
                                      name = legLabel)
             }
           }else{if(any_negative){
             cval <- {cval[1:{ncol_pH_+1}] - min(cval[1:{ncol_pH_+1}])} / diff(range(cval[1:{ncol_pH_+1}]))
             cval <- {cval - min(cval)} / diff(range(cval))
             if(weblink_exists){
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn_interactive(colours = Cols[1:{ncol_pH_+1}], values = cval,
                                                  name = label_interactive(legLabel,
                                                                           data_id = 'legend.title',
                                                                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                           onclick = paste0("window.open(`", as.character(weblink_pH), "`);"),
                                                                           tooltip = weblink_pH))
             }else{
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn(colours = Cols[1:{ncol_pH_+1}], values = cval,
                                      name = legLabel)
             }
           }else{if(any_positive){
             cval <- {cval[{ncol_pH_+1}:ncol_pH] - min(cval[{ncol_pH_+1}:ncol_pH])} / diff(range(cval[{ncol_pH_+1}:ncol_pH]))
             cval <- {cval - min(cval)} / diff(range(cval))
             
             if(weblink_exists){
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn_interactive(colours = Cols[{ncol_pH_+1}:ncol_pH], values = cval,
                                                  name = label_interactive(legLabel,
                                                                           data_id = 'legend.title',
                                                                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                           onclick = paste0("window.open(`", as.character(weblink_pH), "`);"),
                                                                           tooltip = weblink_pH))
             }else{
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background, aes(fill = value)) +
                 scale_fill_gradientn(colours = Cols[{ncol_pH_+1}:ncol_pH], values = cval,
                                      name = legLabel)  
             }
           }}}
         },
         
         ship = {
           
           if(exists(shipSummaryDataOrRaw)){
             
             weblink_exists <- exists('weblink_shipping')
             
             plt_background <-switch(shipSummaryDataOrRaw,
                                     
                                     summary = {
                                       shipType <- unique(dat_background$ship_class)
                                       ShipType <- paste0(toupper(substr(shipType,1,1)), substr(shipType,2,nchar(shipType)))
                                       leg_lab <- bquote(.(paste0('Ship time:', '\n', shipType, ' vessels', '\n', '(days)')))
                                       vesselsPresent <- !is.na(dat_background$total_time)
                                       if(weblink_exists){
                                         plt <- 
                                           ggplot() + 
                                           geom_sf(data = dat_background[vesselsPresent,], aes(fill = total_time), colour = 'black', linewidth = polyLineWidth) +
                                           scale_fill_viridis_c_interactive(option = 'mako', trans = 'log10', direction = -1,
                                                                            name = label_interactive(leg_lab,
                                                                                                     data_id = 'legend.title',
                                                                                                     hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                                                     onclick = paste0("window.open(`", weblink_shipping, "`);"),
                                                                                                     tooltip = weblink_shipping)) +
                                           geom_sf(data = dat_background[!vesselsPresent,], fill = 'white', colour = 'black', linewidth = polyLineWidth)
                                         plt
                                       }else{
                                         plt <- 
                                           ggplot() + 
                                           geom_sf(data = dat_background[vesselsPresent,], aes(fill = total_time), colour = 'black', linewidth = polyLineWidth) +
                                           scale_fill_viridis_c(option = 'mako', trans = 'log10', direction = -1,
                                                                name = leg_lab) +
                                           geom_sf(data = dat_background[!vesselsPresent,], fill = 'white', colour = 'black', linewidth = polyLineWidth)
                                         plt
                                       }
                                     },
                                     
                                     raw = {
                                       shipActivity <- unique(dat_background$activity)
                                       leg_lab <- paste0('Ship time\n', shipActivity, ' vessels\n', '(days/year)')
                                       # leg_lab <- bquote(.(paste0('Ship time:', '\n', shipActivity, ' vessels', '\n', '(days/yr)')))
                                       if(weblink_exists){
                                         plt <- 
                                           ggplot() + 
                                           geom_sf(data = dat_background, aes(fill = value), colour = 'black', linewidth = polyLineWidth) +
                                           scale_fill_viridis_c_interactive(option = 'mako', trans = 'log10', direction = -1,
                                                                            name = label_interactive(leg_lab,
                                                                                                     data_id = 'legend.title',
                                                                                                     hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                                                     onclick = paste0("window.open(`", weblink_shipping, "`);"),
                                                                                                     tooltip = weblink_shipping))
                                         plt
                                       }else{
                                         plt <- 
                                           ggplot() + 
                                           geom_sf(data = dat_background, aes(fill = value), colour = 'black', linewidth = polyLineWidth) +
                                           scale_fill_viridis_c(option = 'mako', trans = 'log10', direction = -1,
                                                                name = leg_lab) # +
                                         # geom_sf(data = dat_background[!vesselsPresent,], fill = 'white', colour = 'black', linewidth = polyLineWidth)
                                         plt
                                       }
                                     }
             )
             
           }else{
             
             plt <- ggplot()
             plt
             
           }
         }
  )
  
  if(plotSignificanceContours){
    if(background %in% c('sst', 'pH')){
      # Create contours in standard lat-lon coordinates
      dat_ <- st_transform(dat_background, crs_world)
      centroids <- as.data.frame(st_coordinates(st_centroid(dat_$geometry))) # get coordinates of polygon centroids
      # Account for rounding error after coordinate transform
      centroids <- round(centroids, 4)
      centroids$pval <- dat_$pval
      centroids$plevel <- dat_$plevel
      xgrid <- sort(unique(centroids$X))
      ygrid <- sort(unique(centroids$Y))
      nx <- length(xgrid)
      ny <- length(ygrid)
      xygrid <- expand.grid(X = xgrid, Y = ygrid)
      centroids <- merge(xygrid, centroids, all.x = TRUE)
      centroids <- centroids[order(centroids$Y, centroids$X),]
      zgrid <- matrix(centroids$pval, nx, ny)
      
      if(background == 'pH'){
        if(is.null(pH_sigContours)){
          log_pval <- log10(dat_$pval)
          seq_pval <- seq(min(log_pval), max(log_pval), length = 5)
          seq_pval <- pretty(seq_pval)
          significanceContours <- 10 ^ seq_pval[2:4]
        }
      }
      if(background == 'sst'){
        if(is.null(sst_sigContours)){
          significanceContours <- c(0.05, 0.25, 0.5)
        }
      }
      cL <- contourLines(xgrid, ygrid, zgrid, levels = significanceContours)      
      # Transform contour lines to spatial object
      # library(maptools)
      crs_world_ <- paste0('+init=epsg:', crs_world)
      cLdf <- ContourLines2SLDF(cL, proj4string = CRS(crs_world_))
      # Transform to sf
      cLdf_sf <- st_as_sf(cLdf)
      # Tranform to stereogrpahic coordinates
      cLdf_sf <- st_transform(cLdf_sf, crs_use)
      names(cLdf_sf)[names(cLdf_sf) == 'level'] <- 'p-value'
      plt_background <- plt_background +
        geom_sf(data = cLdf_sf, aes(linetype = `p-value`), linewidth = contourLineWidth)#, inherit.aes = FALSE)# +
    }
  }
  
  plt_background <- plt_background +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.background = element_rect(fill = 'white', colour = 'white')
    )
  
  leg_background <- get_legend(plt_background)
  
  #~~~~~~~~~~~~~~~~~~~~~~
  # Main map -- no legend
  #~~~~~~~~~~~~~~~~~~~~~~
  
  plt_map <- plt_background + guides(fill = 'none') + guides(linetype = 'none')# + theme(legend.position = 'right')
  
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
  
  plt_map <- plt_map + guides(fill_new = 'none')
  
  
  if(backgroundOnly){
    return(
      list(plot = ggdraw(plot_grid(plt_map, leg_background, ncol = 2, rel_widths = c(0.8, 0.2))))
    )
  }
  
  if(anyStations | anyPlastic){
    plt_map <- plt_map +
      scale_shape_manual(values = setNames(symbols$symbol, symbols$Type))
  }
  
  # Research stations
  if(anyStations){
    plt_map <- plt_map +
      new_scale_colour() +
      geom_sf_interactive(data = dat_stations,
                          aes(shape = Type, colour = Seasonality, data_id = Record_ID, tooltip = tooltip),
                          alpha = alpha, size = ptSize, stroke = 2, show.legend = FALSE) +
      scale_colour_manual(values = c('forestgreen','purple'))
    # scale_colour_manual(values = c('forestgreen','firebrick'))
  }
  
  # Plastic samples
  if(anyPlastic){
    plt_map <- plt_map +
      new_scale_fill() +
      geom_sf_interactive(data = dat_plastic,
                          aes(fill = Source, shape = SampleType_grouped, data_id = data_id, tooltip = tooltip),
                          alpha = alpha, size = ptSize, show.legend = FALSE) +
      scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
  }
  
  plt_map <- plt_map + guides(fill_new_new = 'none')
  
  # plt_map <- plt_map +
  #   coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax']))# +
  plt_map <- plt_map +
    coord_sf(xlim = c(BBox['xmin'], BBox['xmax']), ylim = c(BBox['ymin'], BBox['ymax']))# +
  
  #~~~~~~~
  # Legend
  #~~~~~~~
  
  # Research stations
  if(anyStations){
    weblink_exists <- exists('weblink_facilities')
    symbols_ <- symbols[symbols$Class == 'ResearchStation',]
    plt_stations <-
      ggplot() +
      geom_sf(data = dat_stations,
              aes(shape = Type, colour = Seasonality),
              alpha = 1, size = ptSize, stroke = 2) +
      scale_colour_manual(values = c('forestgreen','purple'))
    # scale_colour_manual(values = c('forestgreen','firebrick'))
    if(weblink_exists){
      plt_stations <- plt_stations +
        scale_shape_manual_interactive(values = setNames(symbols_$symbol, symbols_$Type),
                                       name = label_interactive(label = 'Facility',
                                                                data_id = 'legend.title',
                                                                hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                                                                onclick = paste0("window.open(`", weblink_facilities, "`);"),
                                                                tooltip = weblink_facilities)) +
        guides(shape = guide_legend_interactive(override.aes = list(size = legPtSize, stroke = 1)),
               colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
    }else{
      plt_stations <- plt_stations +
        scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
        guides(shape = guide_legend(title = 'Facility', override.aes = list(size = legPtSize)),
               colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
    }
    plt_stations <- plt_stations +
      theme(legend.key = element_blank())# +
      # guides(shape = guide_legend(override.aes = list(size = legPtSize)),
      #        colour = guide_legend(override.aes = list(size = legPtSize, shape = 1)))
    leg_stations <- get_legend(plt_stations)
    # get legend dimensions
    lw = (leg_stations$widths)
    lh = (leg_stations$heights)
    lw = lw[grepl('cm', lw)]
    lh = lh[grepl('cm', lh)]
    leg_stations_width = sum(as.numeric(gsub('cm', '', lw)))
    leg_stations_height = sum(as.numeric(gsub('cm', '', lh)))
  }else{
    leg_stations <- NULL
    leg_stations_width = 0.5
    leg_stations_height = 0.5
  }
  
  if(anyPlastic){
    symbols_ = symbols[symbols$Class == 'PlasticSample',]
    plt_plastic_samples <-
      ggplot() +
      geom_sf(data = dat_plastic,
              aes(fill = Source, shape = SampleType_grouped),
              alpha = 1, size = ptSize) +
      guides(
        shape = guide_legend(
          title = 'Sample type',
          override.aes = list(size = legPtSize)),
        fill = guide_legend_interactive(
          override.aes = list(
            shape = 21,
            size = legPtSize))
      ) +
      scale_fill_manual_interactive(
        values = setNames(pltColours$colour, pltColours$Source),
        data_id = setNames(pltColours$Source, pltColours$Source),
        tooltip = setNames(pltColours$URL, pltColours$Source),
        labels = function(Source) {
          lapply(Source, function(s) {
            u <- unique(dat_plastic$URL[dat_plastic$Source == s])
            label_interactive(
              as.character(s),
              data_id = as.character(s),
              hover_css = 'fill:blue;font-size:13px;font-weight:bold',
              onclick = paste0("window.open(`", as.character(u), "`);"),
              tooltip = as.character(u)
            )
          })
        }
      ) +
      scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
      theme(legend.key = element_blank())
    leg_plastic <- get_legend(plt_plastic_samples)
    # get legend dimensions
    lw = (leg_plastic$widths)
    lh = (leg_plastic$heights)
    lw = lw[grepl('cm', lw)]
    lh = lh[grepl('cm', lh)]
    leg_plastic_width = sum(as.numeric(gsub('cm', '', lw)))
    leg_plastic_height = sum(as.numeric(gsub('cm', '', lh)))
  }else{
    leg_plastic <- NULL
    leg_plastic_width = 0.5
    leg_plastic_height = 0.5
  }
  
  # Background
  if(anyBackground){
    # leg_background <- get_legend(plt_background)
    # get legend dimensions
    lw = (leg_background$widths)
    lh = (leg_background$heights)
    lw = lw[grepl('cm', lw)]
    lh = lh[grepl('cm', lh)]
    leg_background_width = sum(as.numeric(gsub('cm', '', lw)))
    leg_background_height = sum(as.numeric(gsub('cm', '', lh)))
  }else{
    leg_background <- NULL
    leg_background_width = 0.5
    leg_background_height = 0.5
  }
  
  if(plasticOnly){
    return(
      list(plot = ggdraw(plot_grid(plt_map, leg_plastic, leg_stations, ncol = 3, rel_widths = c(0.7, 0.15, 0.15))))
    )
  }
  
  # Find size (cm) of combined legend -- include spacings
  leg_width <- {2 * 1} + 3 * max(c(leg_plastic_width, leg_stations_width, leg_background_width))
  leg_height <- {3 * 2} + max(c(leg_plastic_height, leg_stations_height, leg_background_height))
  # Size of map & complete plot
  tot_width = leg_width / legWidth
  map_width = tot_width - leg_width
  map_height <- map_width / mapAspectRatio
  tot_height <- max(map_height, leg_height)
  
  # Combine all map components
  leg_complete <- plot_grid(leg_background, leg_plastic, leg_stations,
                            ncol = 3)
  
  plt_complete <- ggdraw(
    plot_grid(
      plt_map,
      leg_complete,
      ncol = 2, rel_widths = c(1-legWidth, legWidth), align = 'h', axis = 't'))
  
  return(
    list(
      plot = plt_complete, width = tot_width, height = tot_height
    )
  )
}

