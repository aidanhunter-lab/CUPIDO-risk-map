#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Custom functions for CUPIDO-risk-map project
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data organisation functions ---------------------------------------------

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


# Get southern hemisphere map data sourced from Natural Earth
getMapData <- function(
    dataDirectory, mapSource = 'Natural Earth',
    fileNames = list(coastline = 'ne_10m_coastline.shp',
                     land = 'ne_10m_land.shp',
                     ocean = 'ne_10m_ocean.shp',
                     iceshelf = 'ne_10m_antarctic_ice_shelves_polys.shp'),
    lon_lim = c(-180, 180), lat_lim = c(-90, 90),
    crs_base = 4326, crs_use = NULL,
    createGrid = FALSE, res = NULL, removeFrac = 0.05,
    verbose = FALSE, loadFromFile = TRUE, autoSave = TRUE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load and return map data as a spatial data frame.
  # Default arguments return entire southern hemisphere.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Load packages & set directories -----------------------------------------
  # Required packages
  library(sf)
  library(sp)
  
  # Specify directories of stored shape files
  if(!is.null(mapSource)) dataDirectory <- paste(dataDirectory, mapSource, sep = '/')
  if(!dir.exists(dataDirectory)) stop(paste0('Specified dataDirectory (', dataDirectory, ') does not exist!'))
  subDirectories <- basename(list.dirs(dataDirectory))
  subDirectories <- subDirectories[subDirectories != basename(dataDirectory)]
  anySubDirectories <- length(subDirectories) > 0
  if(!anySubDirectories){
    coastlineFile <- paste(dataDirectory, fileNames$coastline, sep = '/')
    landFile      <- paste(dataDirectory, fileNames$land,      sep = '/')
    oceanFile     <- paste(dataDirectory, fileNames$ocean,     sep = '/')
    iceshelfFile  <- paste(dataDirectory, fileNames$iceshelf,  sep = '/')
  }else{
    coastDir <- subDirectories[grep('coast', subDirectories)]
    landDir  <- subDirectories[grep('land', subDirectories)]
    oceanDir <- subDirectories[grep('ocean', subDirectories)]
    iceDir   <- subDirectories[grep('ice', subDirectories)]
    coastlineFile <- paste(dataDirectory, coastDir, fileNames$coastline, sep = '/')
    landFile      <- paste(dataDirectory, landDir,  fileNames$land,      sep = '/')
    oceanFile     <- paste(dataDirectory, oceanDir, fileNames$ocean,     sep = '/')
    iceshelfFile  <- paste(dataDirectory, iceDir,   fileNames$iceshelf,  sep = '/')
  }
  
  # Processed map data is saved into '/data/map/temp' directory
  filePath <- paste(dirname(dataDirectory), 'temp', sep = '/')
  if(!filePath %in% list.dirs(dirname(dataDirectory))) dir.create(filePath)
  nMaps <- 0.25 * length(lon_lim) * length(unlist(lat_lim))
  if(is.numeric(lat_lim) & length(lat_lim) == 2) lat_lim <- list(lat_lim)
  fileName_map <- sapply(1:nMaps, function(z){
    paste('map data',
          paste('lon', paste(lon_lim, collapse = ' ')),
          paste('lat', paste(lat_lim[[z]], collapse = ' ')), sep = '_')})
  
  fullFile_map <- paste0(paste(filePath, fileName_map, sep = '/'), '.rds')

  # Load compiled map from file ---------------------------------------------
  # Load from file to save running this script
  mapIsLoaded <- FALSE
  if(loadFromFile){
    if(all(file.exists(fullFile_map))){
      dat <- lapply(fullFile_map, readRDS)
      # dat <- readRDS(fullFile_map)
      # crs <- st_crs(dat)
      # crs_base <- st_crs(st_read(coastlineFile, quiet = !verbose))
      mapIsLoaded <- TRUE}}
  
  # Create map --------------------------------------------------------------
  if(!mapIsLoaded){
    # Load the data
    coastline <- st_read(coastlineFile, quiet = !verbose)
    land      <- st_read(landFile,      quiet = !verbose)
    iceshelf  <- st_read(iceshelfFile,  quiet = !verbose)
    
    coastline <- coastline[!is.na(coastline$featurecla),]
    land <- land[!is.na(land$featurecla),]
    iceshelf <- iceshelf[!is.na(iceshelf$featurecla),]
    
    # Merge data into single data frame
    land      <- land[c('featurecla')]
    iceshelf  <- iceshelf[c('featurecla')]
    coastline <- coastline[c('featurecla')]
    land$featurecla      <- 'land'
    iceshelf$featurecla  <- 'ice'
    coastline$featurecla <- 'coastline'
    dat <- rbind(land, iceshelf, coastline)
    names(dat)[1] <- 'feature'
    
    dat <- st_make_valid(dat)
    
    # Set the coordinate reference system (CRS)
    if(is.numeric(crs_base)) crs_base <- st_crs(crs_base)
    if(is.numeric(crs_use)) crs_use <- st_crs(crs_use)

    # Create map bounding line
    mapLimit_line <- lapply(1:nMaps, function(z){
      x <- matrix(c(Longitude = seq(lon_lim[1], lon_lim[2], length.out = 10 * diff(lon_lim)),
                    Latitude = rep(lat_lim[[z]][2], 10 * diff(lon_lim))), ncol = 2)
      x <- st_linestring(x, dim = 'XY')
      x <- st_sf(st_as_sfc(list(x), crs = crs_base))
      st_geometry(x) <- 'geometry'
      st_transform(x, crs = crs_use)})
    mapLimit_poly <- lapply(1:nMaps, function(z){
      st_cast(mapLimit_line[[z]], 'POLYGON')}) 
    
    # Crop map data to within boundary
    dat <- st_transform(dat, crs = crs_use)
    if(!verbose){
      dat <- lapply(1:nMaps, function(z)
        suppressWarnings(st_intersection(mapLimit_poly[[z]], dat)))}else{
          dat <- lapply(1:nMaps, function(z)
            st_intersection(mapLimit_poly[[z]], dat))}
    #   z=1
    # ggplot() + geom_sf(data = subset(dat[[z]], feature == 'land'), fill = 'darkgreen') + geom_sf(data = subset(dat[[z]], feature == 'ice'), fill = 'lightblue') + geom_sf(data = subset(dat[[z]], feature == 'coastline'), colour = 'black')
    
    # Include map boundary line in dat
    mapLimit_line <- lapply(1:nMaps, function(z){
      x <- mapLimit_line[[z]]
      x$feature <- 'boundary'
      x[c('feature', 'geometry')]})
    dat <- lapply(1:nMaps, function(z) rbind(dat[[z]], mapLimit_line[[z]]))
    
    # Create an ocean polygon -- it's tricky to eliminate narrow edges between adjacent polygons...
    # The land and ice polygons are not perfectly aligned so ocean can creep into
    # the gaps. Solve this by expanding the coastline points into tiny polygons then
    # merging all polygons to mask all dry regions.
    coastline <- lapply(1:nMaps, function(z)
      st_cast(subset(dat[[z]], feature == 'coastline'), 'MULTILINESTRING'))
    xy <- lapply(1:nMaps, function(z) st_coordinates(coastline[[z]]))
    if(!verbose){
      coastline <- lapply(1:nMaps, function(z)
        as.data.frame(suppressWarnings(st_cast(coastline[[z]], 'POINT'))))}else{
          coastline <- lapply(1:nMaps, function(z)
            as.data.frame(st_cast(coastline, 'POINT')))}
    coastline <- lapply(1:nMaps, function(z)
      coastline[[z]][,names(coastline[[z]]) != 'geometry', drop = FALSE])
    bb <- lapply(1:nMaps, function(z) st_bbox(dat[[z]]))
    x <- lapply(1:nMaps, function(z) outer(xy[[z]][,'X'], c(-1,1) * diff(bb[[z]][c('xmin','xmax')]) * 1e-3, '+'))
    y <- lapply(1:nMaps, function(z) outer(xy[[z]][,'Y'], c(-1,1) * diff(bb[[z]][c('ymin','ymax')]) * 1e-3, '+'))
    nxy <- lapply(1:nMaps, function(z) nrow(xy[[z]]))
    ids <- lapply(1:nMaps, function(z) rownames(coastline[[z]]))
    xy <- lapply(1:nMaps, function(w){
      xw <- x[[w]]
      yw <- y[[w]]
      lapply(1:nxy[[w]], function(z){
        Polygons(list(
          Polygon(matrix(c(xw[z,1], yw[z,1],
                           xw[z,2], yw[z,1],
                           xw[z,2], yw[z,2],
                           xw[z,1], yw[z,2],
                           xw[z,1], yw[z,1]), ncol = 2, byrow = TRUE))), ids[[w]][z])
      })
    })
    xy <- lapply(1:nMaps, function(z){
      w <- SpatialPolygons(xy[[z]])
      w <- SpatialPolygonsDataFrame(w, data = coastline[[z]])
      w <- st_as_sf(w)
      st_crs(w) <- crs_use
      return(w)})
    ocean <- lapply(1:nMaps, function(z){
      w <- st_difference(mapLimit_poly[[z]], st_make_valid(st_union(rbind(dat[[z]], xy[[z]]))))
      w$feature <- 'ocean'
      w[c('feature','geometry')]})
    dat <- lapply(1:nMaps, function(z) rbind(dat[[z]], ocean[[z]]))

    # z = 2
    # ggplot() + geom_sf(data = subset(dat[[z]], feature == 'ocean'), fill = 'blue') + geom_sf(data = subset(dat[[z]], feature == 'land'), fill = 'darkgreen') + geom_sf(data = subset(dat[[z]], feature == 'ice'), fill = 'lightblue') + geom_sf(data = subset(dat[[z]], feature == 'coastline'), colour = 'black') + geom_sf(data = subset(dat[[z]], feature == 'boundary'), colour = 'black', linewidth = 1)
    
  }
  
  # Create grid -------------------------------------------------------------
  # Define map grid using resolution 'res'
  made_grid <- FALSE
  # if(!gridIsLoaded){
  #   if(createGrid){
  #     if(is.null(res)){
  #       warning("Must specify 'res' to create grid.")}else{
  #         made_grid <- TRUE
  #         bb <- matrix(c(lon_lim[1], lat_lim[1],
  #                        lon_lim[2], lat_lim[1],
  #                        lon_lim[2], lat_lim[2],
  #                        lon_lim[1], lat_lim[2],
  #                        lon_lim[1], lat_lim[1]
  #         ), ncol = 2, byrow = TRUE)
  #         bb <- st_linestring(bb, dim = 'XY')
  #         bb <- st_sf(st_as_sfc(list(bb), crs = crs_base))
  #         st_geometry(bb) <- 'geometry'
  #         bb <- st_cast(bb, 'POLYGON')
  #         map_grid <- st_make_grid(bb, cellsize = res)
  #         map_grid <- st_sf(map_grid)
  #         st_geometry(map_grid) <- 'geometry'
  #         # Get grid cell mid points (not centroids, but midway between lat-lon bounds)
  #         cell_centroids <- st_centroid(map_grid)
  #         centroid_coords <- st_coordinates(cell_centroids)
  #         lon_mids <- seq(lon_lim[1] + 0.5 * res['lon'], lon_lim[2] - 0.5 * res['lon'], res['lon'])
  #         map_grid$mid_lon_box <- lon_mids[apply(outer(centroid_coords[,'X'], lon_mids, FUN = function(x,y) abs(x-y)), 1, which.min)]
  #         lat_mids <- seq(lat_lim[1] + 0.5 * res['lat'], lat_lim[2] - 0.5 * res['lat'], res['lat'])
  #         map_grid$mid_lat_box <- lat_mids[apply(outer(centroid_coords[,'Y'], lat_mids, FUN = function(x,y) abs(x-y)), 1, which.min)]
  #         # map_grid$centroid_lon_box <- centroid_coords[,1]
  #         # map_grid$centroid_lat_box <- centroid_coords[,2]
  #         cell_centroids <- st_transform(cell_centroids, crs)
  #         centroid_coords <- st_coordinates(cell_centroids)
  #         map_grid$centroid_x_box <- centroid_coords[,1]
  #         map_grid$centroid_y_box <- centroid_coords[,2]
  #         map_grid <- st_transform(map_grid, crs = crs) # grid covers entire map domain
  #         
  #         # Adjust grid to remove (portions of) cells masked by land/ice
  #         if(!verbose){
  #           mask <- suppressWarnings(st_difference(mapLimit_poly, ocean))
  #           map_grid <- suppressWarnings(st_difference(map_grid, mask))}else{
  #             mask <- st_difference(mapLimit_poly, ocean)
  #             map_grid <- st_difference(map_grid, mask)}
  #         
  #         # Get true grid cell centroids (post cropping) and grid cell areas
  #         if(!verbose){
  #           cell_centroids <- suppressWarnings(st_centroid(map_grid))}else{
  #             cell_centroids <- st_centroid(map_grid)}
  #         centroid_coords <- st_coordinates(cell_centroids)
  #         map_grid$centroid_x_true <- centroid_coords[,1]
  #         map_grid$centroid_y_true <- centroid_coords[,2]
  #         map_grid <- st_transform(map_grid, crs = crs_base)
  #         if(!verbose){
  #           cell_centroids <- suppressWarnings(st_centroid(map_grid))}else{
  #             cell_centroids <- st_centroid(map_grid)}
  #         centroid_coords <- st_coordinates(cell_centroids)
  #         map_grid$centroid_lon_true <- centroid_coords[,1]
  #         map_grid$centroid_lat_true <- centroid_coords[,2]
  #         cell_areas <- st_area(map_grid)
  #         if(all(attr(cell_areas, 'units')$numerator == 'm')) cell_areas <- cell_areas * 1e-6  # convert m^2 to km^2
  #         map_grid$area_km2 <- as.numeric(cell_areas)
  #         map_grid <- st_transform(map_grid, crs = crs)
  #         
  #         # Omit tiny cells that may appear within intricate coastline. Remove cells that
  #         # are smaller than 0<removeFrac<1 of their expected (uncropped) area.
  #         lat_include <- lat_mids %in% map_grid$mid_lat_box
  #         j <- sum(lat_include) # number of expected areas to calculate
  #         y_ <- lat_mids[lat_include]
  #         x_ <- rep(0, length(y))
  #         x <- outer(x_, 0.5 * c(-1,1) * res['lon'], '+')
  #         y <- outer(y_, 0.5 * c(-1,1) * res['lat'], '+')
  #         pol <- lapply(1:j, function(z){
  #           Polygons(list(Polygon(matrix(c(x[z,1],y[z,1],
  #                                          x[z,2],y[z,1],
  #                                          x[z,2],y[z,2],
  #                                          x[z,1],y[z,2],
  #                                          x[z,1],y[z,1]), ncol = 2, byrow = TRUE))), as.character(z))
  #         })
  #         pol <- SpatialPolygons(pol)
  #         pol <- SpatialPolygonsDataFrame(pol, data = data.frame(pol = seq_along(pol)))
  #         pol <- st_as_sf(pol)
  #         st_crs(pol) <- crs_base
  #         a <- st_area(pol)
  #         if(all(attr(a, 'units')$numerator == 'm')) a <- a * 1e-6  # convert m^2 to km^2
  #         a <- as.numeric(a) # expected areas
  #         j <- apply(outer(map_grid$mid_lat_box, lat_mids[lat_include], '=='), 1, which)
  #         remove <- map_grid$area_km2 / a[j] < removeFrac
  #         map_grid <- map_grid[!remove,]
  #         
  #         # Put grid cells in order
  #         rowOrder <- order(map_grid$mid_lat_box, map_grid$mid_lon_box, decreasing = FALSE)
  #         colOrder <- c('feature','mid_lon_box','mid_lat_box','centroid_x_box','centroid_y_box', 'centroid_lon_true','centroid_lat_true','centroid_x_true','centroid_y_true','area_km2','geometry')
  #         map_grid <- map_grid[rowOrder,colOrder]
  #         
  #         # ggplot() + geom_sf(data = subset(dat, feature == 'ocean'), fill = 'skyblue') + geom_sf(data = subset(dat, feature == 'land'), fill = 'darkgreen') + geom_sf(data = subset(dat, feature == 'ice'), fill = 'white') + geom_sf(data = subset(dat, feature == 'coastline'), colour = 'black') + geom_sf(data = subset(dat, feature == 'boundary'), colour = 'black', linewidth = 1) + geom_sf(data = st_cast(map_grid, 'LINESTRING'), linewidth = 0.1, colour = 'red')
  #       }
  #   }
  # }
  
  # Store output ------------------------------------------------------------
  output <- list()
  output$lon_lim <- lon_lim
  output$lat_lim <- lat_lim
  output$map_data <- setNames(dat, paste('lat', sapply(1:nMaps, function(z) paste(lat_lim[[z]], collapse = '_')), sep = '_'))
  if(createGrid){
    output$map_grid <- map_grid
    output$grid_res <- res}
  output$crs_use <- crs_use
  output$crs_base <- crs_base

  
  # Save --------------------------------------------------------------------
  if(autoSave){
    # Save map data into 'temp' directory within current directory
    if(!mapIsLoaded){
      for(i in 1:nMaps) saveRDS(dat[[i]], fullFile_map[i])}
    # Save the grid
    if(made_grid){
      if(!gridIsLoaded) saveRDS(map_grid, fullFile_grid)}
  }
  
  # End ---------------------------------------------------------------------
  return(output)
}


# Load and organise data required for maps
get_data <- function(
    baseDirectory, dataDirectory, mapDirectory, res,
    # range_plastic_or_cells = 'cells',
    allLitterTypes = FALSE,
    sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
    sstTrend_significantOnly = TRUE, pHTrend_significantOnly = TRUE,
    significanceLevel = 0.05, significanceContours = c(0.05, 0.25, 0.5),
    SST_overallTrend = TRUE, pH_overallTrend = TRUE, sigDigits = 2, roundShipTime = FALSE,
    shipOrPersonTime = NULL, indexGridCells = TRUE, loadTooltipFromFile = TRUE,
    verbose = FALSE, theseData = NULL){
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
  
  crs_world <- st_crs(crs_world)
  crs_use <- st_crs(crs_use)
  
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
    grid_cells_all <- SpatialPolygons(grid_cells_all, proj4string = CRS(st_as_text(crs_world))) # create spatial object
    grid_cells_all <- st_as_sf(grid_cells_all)
    grid_cells_all$cell_index <- 1:nrow(grid_cells_all)
    grid_cell_areas <- st_area(grid_cells_all) # get grid cell areas
    if(all(attr(grid_cell_areas, 'units')$numerator == 'm')) grid_cell_areas <- grid_cell_areas * 1e-6  # convert m^2 to km^2
    grid_cells_all$area <- as.numeric(grid_cell_areas)
    grid_centroids <- suppressWarnings(st_centroid(grid_cells_all))
    centroid_coords <- st_coordinates(grid_centroids)
    grid_cells_all$centroid_lon <- centroid_coords[,1]
    grid_cells_all$centroid_lat <- centroid_coords[,2]
    grid_centroids_xy <- st_transform(grid_centroids, crs_use)
    centroid_coords <- st_coordinates(grid_centroids_xy)
    grid_cells_all$centroid_x <- centroid_coords[,1]
    grid_cells_all$centroid_y <- centroid_coords[,2]
    
  }
  
  # Load ecoregion shape files -----------------------------------------------------
  
  if(is.null(theseData) || 'ecoregions' %in% theseData){
  
  # verbose = TRUE
  f <- 'data/commondata/data0/meow_ecos_expl_clipped_expl.shp'
  f <- paste(dataDirectory, 'marine ecoregions', f, sep = '/')
  eco <- st_read(f, quiet = !verbose)
  
  # Filter ecoregions by province
  keepProvinces <- c(60, 61, 48, 59, 62) # 60=Scotia Sea, 61=Continental High Antarctic, 48=Magellenic, 59=Subantarctic Islands, 62=Subantarctic New Zealand
  eco <- eco[eco$PROV_CODE %in% keepProvinces,]
  
  # Change coordinate reference system
  eco <- st_transform(eco, crs_use)
  
# #  adjust from here down
#   
#   # Match bounding box to coastline data
#   attr(st_geometry(eco), 'bbox') <- st_bbox(nc)
#   
#   # Crop ecoregions to map limits
#   eco <- st_crop(eco, st_bbox(eco))
#   
#   assign('eco', eco, envir = parent.frame())
  }  
  
  # Load shipping data  -----------------------------------------------------
  
  if(is.null(theseData) || 'ship' %in% theseData){
    
  switch(shipSummaryDataOrRaw,
         summary = {
           # Summary section is redundant
           
           f <- 'McCarthy_2022/so_ecoregion_combined.csv' # data table with all ship types
           f <- paste(dataDirectory, 'shipping', f, sep = '/')
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
           
           assign('ship_poly', ship_poly, envir = parent.frame())
           
         },
         
         raw = {
           # It doesn't matter what value is set for 'res', just use 9x3 data
           f <- 'McCarthy_2022/compiled/ship_time_res_9x3.csv'
           f <- paste(dataDirectory, 'shipping', f, sep = '/')
           ship <- read.table(f, sep = ',', header = TRUE)
           
           fm <- 'McCarthy_2022/compiled/vessel_metadata.csv'
           fm <- paste(dataDirectory, 'shipping', fm, sep = '/')
           ship_metadata <- read.table(fm , sep = ',', header = TRUE)
           
           ship <- merge(ship, ship_metadata, by = 'vessel_id')
           
           ship$person_capacity[is.na(ship$person_capacity)] <- NA
           ship$person_time[is.na(ship$person_capacity)] <- 0

           # For each year and vessel activity, sum ship time over all vessels
           ship <- aggregate(cbind(ship_time, person_time) ~ lonmin + lonmax + latmin + latmax + year + activity_type, data = ship, FUN = sum)
           
           # Sum ship time over vessel activity
           ship_ <- aggregate(cbind(ship_time, person_time) ~ lonmin + lonmax + latmin + latmax + year, data = ship, FUN = sum)
           ship_$activity_type <- 'all'
           ship_ <- ship_[names(ship)]
           ship <- rbind(ship, ship_)
           
           # Average ship_time over years
           ship_ <- aggregate(cbind(ship_time, person_time) ~ lonmin + lonmax + latmin + latmax + activity_type, data = ship, FUN = mean)
           ship_$year <- 'all'
           ship_ <- ship_[names(ship)]
           ship <- rbind(ship, ship_)
           
           # Use the annual average only, discard data pertaining to specific years
           ship <- subset(ship, year == 'all')
           
           if(roundShipTime){
             # Round ship_time to nearest day -- this removes grid cells with negligible
             # ship presence, creating a more useful colour scale with greater focus
             # on areas of high traffic.
             # ship$ship_time <- round(ship$ship_time)
             # Hmmm, I only want to crop the small measures, so only round numbers
             # less than one.
             lowTraffic <- ship$ship_time < 1
             ship$ship_time[lowTraffic] <- round(ship$ship_time[lowTraffic])
             lowTraffic <- ship$person_time < 1
             ship$person_time[lowTraffic] <- round(ship$person_time[lowTraffic])
           }
           
           if(is.null(shipOrPersonTime)){
             ship <- melt(ship, measure.vars = c('ship_time', 'person_time'))
           }else{
             if(!shipOrPersonTime %in% c('ship', 'person')) error("'shipOrPersonTime' must be either 'ship' or 'person'")
             ship <- switch(shipOrPersonTime,
                            ship = ship[,names(ship) != 'person_time'],
                            person = ship[,names(ship) != 'ship_time'])
             ship <- melt(ship, measure.vars = paste0(shipOrPersonTime, '_time'))
           }
           ship$variable <- gsub('_', ' ', ship$variable)
           
           # Zero measurements are grid cells where no AIS pings were recorded -- 
           # these could be omitted, but it's probably best to leave them in.
           
           # ship <- subset(ship, value != 0) # remove zero measures
           
           ship$activity_type <- factor(ship$activity_type, 
                                        levels = c('all','fishing','tourism','supply','research','other'),
                                        labels = c('all','fishing','tourism','supply','research','other'))
           
          lat_max_ship <- max(ship$latmax)
           
           # Generate spatial polygons
           ship_poly <- createPolygons(ship)
           ship_poly <- SpatialPolygons(ship_poly, proj4string = CRS(st_as_text(crs_world))) # create spatial object
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
           ship_poly$variable <- ship$variable
           ship_poly$value <- ship$value
           
         }
  )
  
  weblink_shipping <- 'https://doi.org/10.1073/pnas.2110303118'
  assign('weblink_shipping', weblink_shipping, envir = parent.frame())
  
  ship_poly_local <- ship_poly
  
  }
  
  # Load plastic data -----------------------------------------------------
  
  if(is.null(theseData) || 'plastic' %in% theseData){
  
  # filepath = 'data/plastic_quantity'
  # # filename = 'plastic_quantity.csv'
  # filename = 'plastic_quantity_new2.csv'
  # f <- paste(wd_base, filepath, filename, sep = '/')
  # DATA = read.csv(f, stringsAsFactors = TRUE)
  
  filename <- 'plastic_quantity.csv'
  f <- paste(dataDirectory, 'plastic samples', 'collated', filename, sep = '/')
  # f <- paste(shinyDirectory, filename, sep = '/')
  DATA <- read.csv(f, stringsAsFactors = TRUE)
  
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
  
  }
  
  # Data pre-processing ----
  
  if(is.null(theseData) || 'plastic' %in% theseData){
    
  # Include spaces in names of data for better display
  DATA_names_orig <- names(DATA)
  DATA_names_spaces <- DATA_names_orig
  DATA_names_spaces[c(2:9, 11:12)] <- c('Sample type','Sample gear','Litter ID method','Litter category','Litter scale','Plastic form','Plastic size','Sample at station','Site category','Sample ID')
  DATA_name_swap <- as.data.frame(cbind(original = DATA_names_orig, spaces = DATA_names_spaces))
  assign('DATA_name_swap', DATA_name_swap, envir = parent.frame())
  
  # Replace any (MatLab style) NaN values with (R style) NA values
  NaN2NA <- function(z, numericOnly = TRUE){
    if(numericOnly){
      if(is.numeric(z)){
        z[is.nan(z)] <- NA
        as.numeric(z)}else{
          z}}else{
            z[is.nan(z)] <- NA
            z}}
  for(j in 1:ncol(DATA)) DATA[,j] <- NaN2NA(DATA[,j])
  
  # Convert logical variables saved as zeros and one into logical class
  DATA$SampleAtStation <- as.logical(DATA$SampleAtStation)
  
  # Date format
  DATA$Date_1 <- as.character(as.Date(DATA$Date_1, '%d-%b-%Y'))
  DATA$Date_2 <- as.character(as.Date(DATA$Date_2, '%d-%b-%Y'))
  DATA$Date_1[is.na(DATA$Date_1)] <- ''
  DATA$Date_2[is.na(DATA$Date_2)] <- ''
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
  DATA$Year <- as.integer(substr(DATA$Date_1, 1, 4))
  
  # Sample type -- sub/surface water, sediment, etc
  s <- as.character(DATA$SampleType)
  marine <- s == 'seawater'
  DATA$Depth <- as.character(DATA$Depth)
  DATA$Depth[DATA$Depth == 'na'] <- ''
  
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
  DATA$Coordinates[as.logical(DATA$SampleAtStation)] <- coord_lab[as.logical(DATA$SampleAtStation)]
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
  DATA$Variable <- factor(DATA$Variable, levels = c('concentration', 'density', 'mass concentration', 'mass density', 'flux', 'presence/absence'))
  
  
  plastic_units <- unique(DATA[,c('SampleType', 'Variable','Unit')])
  plastic_units$SampleType <- as.character(plastic_units$SampleType)
  plastic_units$SampleType[grepl('water', plastic_units$SampleType) | 
                             grepl('marine', plastic_units$SampleType)] <- 'water'
  plastic_units <- unique(plastic_units)
  plastic_units$SampleType <- factor(plastic_units$SampleType, levels = unique(plastic_units$SampleType))
  
  
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
  
  
  # Plastic scale -- SCALES NEED SPECIFIED BETTER IN DATA SET SO THE BELOW WILL NEED ADJUSTED...
  DATA$LitterScale <- as.character(DATA$LitterScale)
  DATA$LitterScale[grepl('micro and meso', DATA$LitterScale)] <- 'micro and meso'
  DATA$LitterScale[DATA$LitterScale %in% c('see size data','')] <- 'unspecified'
  DATA$LitterScale <- factor(DATA$LitterScale, levels = c('micro', 'meso', 'micro and meso', 'macro', 'unspecified'))
  
  
  # Plastic size
  DATA$PlasticSize <- as.character(DATA$PlasticSize)
  DATA$PlasticSize[DATA$PlasticSize == ''] <- 'unspecified'
  DATA$PlasticSize <- factor(DATA$PlasticSize, levels = unique(DATA$PlasticSize))
  
  
  # Litter category -- only include plastic litter
  DATA$LitterCategory <- as.character(DATA$LitterCategory)
  if(!allLitterTypes){
    DATA <- subset(DATA, LitterCategory == 'plastic')
  }
  
  
  # Measurement statistic
  DATA$Statistic <- as.character(DATA$Statistic)
  DATA$Statistic[DATA$Statistic == 'stddev'] <- 's.d.'
  DATA$Statistic <- factor(DATA$Statistic, levels = c('raw','mean','s.d.','min','max','none'))
  
  # Set non-numeric Values to -Inf to distinguish them from NAs -- this for presence/absence data
  DATA$Value[is.na(DATA$Value) & DATA$Variable == 'presence/absence'] <- -Inf
  
  
  # Order the data
  DATA <- DATA[order(DATA$Source, DATA$SampleType, DATA$SampleID, DATA$LitterScale, DATA$PlasticSize, 
                     DATA$PlasticForm, DATA$Variable, DATA$Statistic, DATA$Replicate),]
  
  
  # Include a data_id variable unique to each sample
  DATA$order <- 1:nrow(DATA)
  d <- unique(DATA[c('Source', 'SampleType', 'SampleAtStation',  'SampleID')])
  d$data_id <- paste('sample', 1:nrow(d))
  
  DATA <- merge(DATA, d, sort = FALSE)
  DATA <- DATA[order(DATA$order),]
  DATA <- DATA[,names(DATA) != 'order']
  
  
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
    dn <- names(DATA)
    # Change to wide-form -- choose variables to produce one row per data_id
    dummyGroupingVariables <- c('PlasticForm_grouped')
    yvars <- c('LitterCategory', 'LitterScale', 'PlasticSize', 'PlasticForm', 'Variable', 'Statistic', 'Replicate', 'Unit') # variables to change to wide format -- displayed over multiple rows
    xvars <- dn[!dn %in% c(yvars, 'Value', dummyGroupingVariables)] # everything else, apart from Value and dummy grouping variables, stays in long format
    # Choose other variables to display in a single row.
    keepVars <- c('Source', 'SampleGear', 'Depth', 'Sample date', 'Coordinates', 'Coordinates (start)', 'Coordinates (end)', 'SampleAtStation', 'Observation')
    # Some of these variables may require multiple rows depending on the data source/sample
    expandRows <- sapply(1:length(keepVars), function(u){
      any(sapply(1:length(sources), function(v){
        d <- DATA[DATA$Source == sources[v],]
        atStation <- unique(d$SampleAtStation)
        any(sapply(1:length(atStation), function(w){
          d2 <- d[d$SampleAtStation == atStation[w],]
          sampleTypes <- unique(d2$SampleType)
          any(sapply(1:length(sampleTypes), function(x){
            allSamples <- unique(d2$SampleID[d2$SampleType == sampleTypes[x]])
            any(sapply(1:length(allSamples), function(y){
              length(unique(d2[[keepVars[u]]][d2$SampleID == allSamples[y] & d2$SampleType == sampleTypes[x]])) > 1
            }))
          }))
        }))
      }))
    })
    xvars <- xvars[!xvars %in% keepVars[expandRows]]
    yvars <- c(keepVars[expandRows], yvars)
    keepVars <- keepVars[!expandRows]
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
    # tooltip.
    omitVars <- xvars[!xvars %in% keepVars]
    temp_dat <- as.data.frame(DATA_wide)[!names(DATA_wide) %in% omitVars]
    temp_dat <- split(temp_dat, seq_len(nrow(temp_dat)))
    sType <- as.character(DATA_wide$SampleType) # units depend on SampleType
    sType[grepl('water', sType) | grepl('marine', sType)] <- 'water'
    DATA_wide$tooltip <- sapply(1:length(temp_dat), function(z) fun_flextable(x = temp_dat[[z]], longVars = yvars, singleRowVars = keepVars, sampleType = sType[z], DATA_name_swap)) # this is time-consuming and could maybe be moved to pre-amble/saved data set...
    xvars <- c(xvars, 'tooltip')
    # Melt data back to original shape
    DATA_long <- melt(DATA_wide, id.vars = xvars, value.name = 'Value')
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
      fn <- paste0(paste(sub('.csv', '', f), 'tooltip', 'long', sep = '_'), '.Rds')
      saveRDS(DATA_long, fn)
    }
    
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Mapping coordinates
  DATA_sf <- st_as_sf(DATA_wide, coords = c('Longitude', 'Latitude'), crs = crs_world)
  DATA_sf <- st_transform(DATA_sf, crs_use)
  
  assign('DATA_sf', DATA_sf, envir = parent.frame())
  assign('DATA_long', DATA_long, envir = parent.frame())
  
  }
  
  # Load research station data ----------------------------------------------
  
  if(is.null(theseData) || 'facility' %in% theseData){
    
  
  filename <- 'COMNAP_Antarctic_Facilities_Master.csv'
  f <- paste(dataDirectory, 'research stations', 'COMNAP', filename, sep = '/')
  STATIONS <- read.csv(f, stringsAsFactors = TRUE)
  
  weblink_facilities <- 'https://github.com/PolarGeospatialCenter/comnap-antarctic-facilities'
  assign('weblink_facilities', weblink_facilities, envir = parent.frame())
  
  }
  
  # Data pre-processing ----

  if(is.null(theseData) || 'facility' %in% theseData){
    
  STATIONS$Type <- reorder(STATIONS$Type, rep(1, nrow(STATIONS)), sum, decreasing = TRUE)
  
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
  
  STATIONS$Record_ID <- paste('facility', STATIONS$Record_ID)
  
  # Mapping coordinates
  STATIONS_sf <- st_as_sf(STATIONS, coords = c('Longitude_DD', 'Latitude_DD'), crs = crs_world)
  
  STATIONS_sf <- st_transform(STATIONS_sf, crs_use)
  xy <- matrix(unlist(STATIONS_sf$geometry), 2, nrow(STATIONS))
  STATIONS$x <- xy[1,]
  STATIONS$y <- xy[2,]
  
  # Calculate each facility's distance from the coast. The resulting distance 
  # estimates are quite crude (mostly) because the coastline data is not very
  # detailed (despite being called hi-res!)
  coastline <- map('worldHires', 'antarctica', plot = FALSE)
  coastPoints <- data.frame(X = coastline$x, Y = coastline$y)
  coastPoints <- coastPoints[!{is.na(coastPoints$X) | is.na(coastPoints$Y)},]
  
  # Use haversine formula for distance calculation.
  # Should use the great circle distance built into a loaded package instead of this...
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
  
  
  STATIONS_sf$tooltip <- paste0(STATIONS_sf$Official_Name, ' ', STATIONS_sf$Type, '\n ',
                                'Primary operator: ', STATIONS_sf$Operator_primary, '\n ',
                                'Peak population size = ', STATIONS_sf$Peak_Population, '\n ',
                                'Distance from coast = ', round(STATIONS_sf$Distance_From_Coast, 1), ' km')
  
  assign('STATIONS_sf', STATIONS_sf, envir = parent.frame())
  
  }
  
  # Load krill data ---------------------------------------------------------
  
  if(is.null(theseData) || 'krill' %in% theseData){
    
  filename <- 'krill_data_mapped.csv'
  f <- paste(dataDirectory, 'KRILLBASE', 'compiled', filename, sep = '/')
  KRILL <- read.csv(f, stringsAsFactors = TRUE)
  
  weblink_krill <- 'https://www.doi.org/10.5285/8b00a915-94e3-4a04-a903-dd4956346439'
  assign('weblink_krill', weblink_krill, envir = parent.frame())
  
  }
  
  # Data pre-processing -----------------------------------------------------
  
  if(is.null(theseData) || 'krill' %in% theseData){

  KRILL <- KRILL[!is.nan(KRILL$value),] # omit missing data
  # Mean abundance across months for each grid cell
  KRILL_allMonths <- aggregate(value ~ lonmin + lonmax + latmin + latmax, KRILL, 'mean')
  KRILL_allMonths$month <- 'all'
  KRILL_allMonths <- KRILL_allMonths[,c(1:4,6,5)]
  KRILL <- rbind(KRILL, KRILL_allMonths)
  # Round to sigDigits significant figures
  KRILL$value <- signif(KRILL$value, sigDigits)
  
  lat_max_krill <- max(KRILL$latmax)
  
  krill_poly <- createPolygons(KRILL) # generate polygons
  krill_poly <- SpatialPolygons(krill_poly, proj4string = CRS(st_as_text(crs_world))) # create spatial object
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
  clims[ncols+1] <- Inf
  for(i in 1:ncols){
    j <- clims[i] <=  krill_poly$value & krill_poly$value < clims[i+1]
    krill_poly$colourgroup[j] = colgroup[i]
  }
  krill_poly$colourgroup <- factor(krill_poly$colourgroup, levels = colgroup)
  
  krill_poly_local <- krill_poly
  
  }
  
  # Load chlorophyll data ---------------------------------------------------------
  
  if(is.null(theseData) || 'chl' %in% theseData){
    
  filename <- paste0('chl_data_mapped_res_', res, '.csv')
  f <- paste(dataDirectory, 'chlorophyll', 'SeaWiFS', filename, sep = '/')
  
  CHL <- read.csv(f, stringsAsFactors = TRUE)
  
  weblink_chlorophyll <- 'https://catalogue.ceda.ac.uk/uuid/a91ae89c4c6d011f84cdfbc3d41a73c8'
  assign('weblink_chlorophyll', weblink_chlorophyll, envir = parent.frame())
  
  }
  
  # Data pre-processing -----------------------------------------------------
  
  if(is.null(theseData) || 'chl' %in% theseData){

  CHL <- CHL[!is.nan(CHL$value),] # omit missing data
  # Mean abundance across months for each grid cell
  CHL_allMonths <- aggregate(value ~ lonmin + lonmax + latmin + latmax, CHL, 'mean')
  CHL_allMonths$month <- 'all'
  CHL_allMonths <- CHL_allMonths[,c(1:4,6,5)]
  CHL <- rbind(CHL, CHL_allMonths)
  
  # Round to 2 significant figures
  CHL$value <- signif(CHL$value, sigDigits)
  
  lat_max_chl <- max(CHL$latmax)
  
  chl_poly <- createPolygons(CHL) # generate polygons
  chl_poly <- SpatialPolygons(chl_poly, proj4string = CRS(st_as_text(crs_world))) # create spatial object
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
  
  chl_poly_local <- chl_poly
  
  }  
  
  # Load SST trend (or anomaly) data --------------------------------------
  
  if(is.null(theseData) || 'sst' %in% theseData){
    
  if(SST_overallTrend){
    filename <- paste0('sst_time_series_trend_res_', res, '.csv') # my month-by-month trend calculation    
  }else{
    filename <- paste0('sst_', sstType, '_res_', res, '.csv') # my month-by-month trend calculation
  }
  f <- paste(dataDirectory, 'sst', 'ESA', 'compiled', filename, sep = '/')
  
  SST <- read.csv(f)
  
  weblink_sst <- 'https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/download?dataset=ESACCI-GLO-SST-L4-REP-OBS-SST'
  assign('weblink_sst', weblink_sst, envir = parent.frame())
  
  }
  
  # Data pre-processing -----------------------------------------------------     
  
  if(is.null(theseData) || 'sst' %in% theseData){
    
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
  
  lat_max_sst <- max(SST$latmax)
  
  # Generate spatial polygons
  sst_poly <- createPolygons(SST)
  sst_poly <- SpatialPolygons(sst_poly, proj4string = CRS(st_as_text(crs_world))) # create spatial object
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

  sst_poly_local <- sst_poly
  
  }
  
  # Load pH trend (or anomaly) data -----------------------------------------
  
  if(is.null(theseData) || 'pH' %in% theseData){
    
  if(pH_overallTrend){
    filename <- paste0('pH_time_series_trend_res_', res, '.csv') # my month-by-month trend calculation    
  }else{
    filename <- paste0('pH_', pHType, '_res_', res, '.csv') # my month-by-month trend calculation
  }
  f <- paste(dataDirectory, 'pH', 'SOCAT', 'compiled', filename, sep = '/')
  
  pH <- read.csv(f)
  
  weblink_pH <- 'https://data.marine.copernicus.eu/product/MULTIOBS_GLO_BIO_CARBON_SURFACE_REP_015_008/description'
  assign('weblink_pH', weblink_pH, envir = parent.frame())
  
  }
  
  # Data pre-processing -----------------------------------------------------     
  
  if(is.null(theseData) || 'pH' %in% theseData){
    
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
  
  lat_max_pH <- max(pH$latmax)
  
  # Generate spatial polygons
  pH_poly <- createPolygons(pH)
  pH_poly <- SpatialPolygons(pH_poly, proj4string = CRS(st_as_text(crs_world))) # create spatial object
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
  
  pH_poly_local <- pH_poly
  
  }

  # Load map shape files -----------------------------------------------------
  
  if(is.null(theseData) || 'map' %in% theseData){
    
  lat_lim_plastic <- ceiling(max(DATA_long$Latitude) / 5) * 5
  lat_lim_plastic <- c(-90, lat_lim_plastic)
  lat_lim_cells <- c(-90, max(c(#lat_max_chl,
                                #lat_max_krill,
                                lat_max_sst#,
                                #lat_max_pH,
                                #lat_max_ship
                                )))
  
  assign('lat_lim_plastic', lat_lim_plastic, envir = parent.frame())
  assign('lat_lim_cells', lat_lim_cells, envir = parent.frame())
  
  nc <- getMapData(
    dataDirectory = mapDirectory,
    lat_lim = list(lat_lim_plastic, lat_lim_cells),
    crs_base = crs_world, crs_use = crs_use)
  
  nc <- nc$map_data
  
  # coerce new data to form of previous data -- adjust later if wanted...
  nc <- lapply(1:length(nc), function(z){
    x <- nc[[z]]
    x$surface <- x$feature
    x <- x[x$surface %in% c('land', 'ice'),]
    x$surface[x$surface == 'ice'] <- 'ice shelf'
    x$location <- 'Antarctica'
    x <- x[,c('location', 'surface', 'geometry')]
    x$surface <- factor(x$surface, levels = unique(x$surface))
    x
  })
  
  assign('nc_plastic', nc[[1]], envir = parent.frame())
  assign('nc_cells', nc[[2]], envir = parent.frame())
  
  # Trim ecoregions to coastline boundaries
  eco <- lapply(1:length(nc), function(z){
    x <- eco
    attr(st_geometry(x), 'bbox') <- st_bbox(nc[[z]])
    suppressWarnings(st_crop(x, st_bbox(x)))})
  
  assign('eco_plastic', eco[[1]], envir = parent.frame())
  assign('eco_cells', eco[[2]], envir = parent.frame())
  
  # Below is code for loading old data

  # verbose <- FALSE
  # f1 <- 'Antarctic coastline polygons/medium res/add_coastline_medium_res_polygon_v7_6.shp'
  # f2 <- 'sub Antarctic coastline polygons/sub_antarctic_coastline_high_res_polygon_v1.0.shp'
  # f1 <- paste(dataDirectory, 'map', f1, sep = '/')
  # f2 <- paste(dataDirectory, 'map', f2, sep = '/')
  # nc1 <- st_read(f1, quiet = !verbose)
  # nc2 <- st_read(f2, quiet = !verbose)
  # 
  # # Merge coastline data sets
  # nc1$location = 'Antarctica'
  # nc2$source <- NULL
  # nc2$surface <- 'land'
  # nc1 <- nc1[,c('location','surface','geometry')]
  # nc2 <- nc2[,c('location','surface','geometry')]
  # nc <- rbind(nc1, nc2)
  # 
  # # Reset coastline bounding box to map limits -- defined by most northern latitude
  # latmax <- -50
  # mapLimits <- data.frame(Longitude = -180:180, 
  #                         Latitude = rep(latmax, 361))
  # mapLimits <- st_as_sf(mapLimits, coords = 1:2, crs = crs_world)
  # crs_world_full <- st_crs(mapLimits)
  # mapLimits <- st_transform(mapLimits, crs_use)
  # crs_use_full <- st_crs(mapLimits)
  # attr(st_geometry(nc), 'bbox') <- st_bbox(mapLimits)
  # 
  # nc$surface = factor(nc$surface, levels = unique(nc$surface))
  # 
  # assign('nc', nc, envir = parent.frame())
  
  }
  
  # Tidy up indexing --------------------------------------------------------
  
  x <- paste(c('ship','chl','krill','sst','pH'),'poly', sep = '_')
  xl <- paste(x, 'local', sep = '_')
  xi <- length(x)
  e <- sapply(xl, function(z) exists(z))
  
  if(!indexGridCells){
    if(is.null(theseData)){
      assign('ship_poly', ship_poly, envir = parent.frame())
      assign('krill_poly', krill_poly, envir = parent.frame())
      assign('chl_poly', chl_poly, envir = parent.frame())
      assign('sst_poly', sst_poly, envir = parent.frame())
      assign('pH_poly', pH_poly, envir = parent.frame())
    }else{
      for(i in 1:xi){
        if(e[i]) assign(x[i], get(x[i]), envir = parent.frame())
      }
    }
  }else{
    if(is.null(theseData)){
      gridCellsWithData <- unique(c(
        ship_poly$cell_index,
        chl_poly$cell_index,
        krill_poly$cell_index,
        sst_poly$cell_index,
        pH_poly$cell_index
      ))
    }else{
      gridCellsWithData <- c()
      for(i in 1:xi){
        if(e[i]){
          y <- get(x[i])
          gridCellsWithData <- c(gridCellsWithData, y$cell_index)
        }
      }
      gridCellsWithData <- unique(gridCellsWithData)
    }
    i <- grid_cells_all$cell_index %in% gridCellsWithData
    grid_cells_all <- grid_cells_all[i,]
    relabelCells <- data.frame(old = grid_cells_all$cell_index)
    relabelCells$new <- 1:nrow(relabelCells)
    relabelFun <- function(dat) sapply(1:nrow(dat), function(z) relabelCells$new[relabelCells$old == dat$cell_index[z]])
    
    if(is.null(theseData)){
      ship_poly$cell_index <- relabelFun(ship_poly)
      chl_poly$cell_index <- relabelFun(chl_poly)
      krill_poly$cell_index <- relabelFun(krill_poly)
      sst_poly$cell_index <- relabelFun(sst_poly)
      pH_poly$cell_index <- relabelFun(pH_poly)
      grid_cells_all$cell_index <- relabelFun(grid_cells_all)
    }else{
      for(i in 1:xi){
        if(e[i]) assign(x[i], relabelFun(get(x[i])))
      }
    }
    
    grid_cells_all <- st_transform(grid_cells_all, crs = crs_use) # convert coordinates
    
    assign('grid_cells_all', grid_cells_all, envir = parent.frame())
    
    plt_map_cells <-
      ggplot(grid_cells_all) +
      geom_sf(data = nc_cells,
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
    
    if(is.null(theseData)){
      assign('ship_poly', ship_poly, envir = parent.frame())
      assign('krill_poly', krill_poly, envir = parent.frame())
      assign('chl_poly', chl_poly, envir = parent.frame())
      assign('sst_poly', sst_poly, envir = parent.frame())
      assign('pH_poly', pH_poly, envir = parent.frame())
    }else{
      for(i in 1:xi){
        if(e[i]) assign(x[i], get(x[i]), envir = parent.frame())
      }
    }
  }
  
  
  
  
  
  # if(indexGridCells){
  #   
  #   if(is.null(theseData)){
  #     gridCellsWithData <- unique(c(
  #       ship_poly$cell_index,
  #       chl_poly$cell_index,
  #       krill_poly$cell_index,
  #       sst_poly$cell_index,
  #       pH_poly$cell_index
  #     ))
  #   }else{
  #     x <- paste(c('ship','chl','krill','sst','pH'),'poly','local',sep='_')
  #     xi <- length(x)
  #     e <- rep(FALSE, xi)
  #     gridCellsWithData <- c()
  #     for(i in 1:xi){
  #       if(exists(x[i])){
  #         e[i] <- TRUE
  #         y <- get(x[i])
  #         gridCellsWithData <- c(gridCellsWithData, y$cell_index) 
  #       }
  #     }
  #     gridCellsWithData <- unique(gridCellsWithData)
  #   }
  #   
  #   i <- grid_cells_all$cell_index %in% gridCellsWithData
  #   grid_cells_all <- grid_cells_all[i,]
  #   relabelCells <- data.frame(old = grid_cells_all$cell_index)
  #   relabelCells$new <- 1:nrow(relabelCells)
  #   
  #   relabelFun <- function(dat) sapply(1:nrow(dat), function(z) relabelCells$new[relabelCells$old == dat$cell_index[z]])
  #   
  #   if(is.null(theseData)){
  #     ship_poly$cell_index <- relabelFun(ship_poly)
  #     chl_poly$cell_index <- relabelFun(chl_poly)
  #     krill_poly$cell_index <- relabelFun(krill_poly)
  #     sst_poly$cell_index <- relabelFun(sst_poly)
  #     pH_poly$cell_index <- relabelFun(pH_poly)
  #     grid_cells_all$cell_index <- relabelFun(grid_cells_all)
  #   }else{
  #     for(i in 1:xi){
  #       if(e[i]){
  #         assign(gsub('_local','',x[i]), relabelFun(get(x[i])))
  #       }
  #     }
  #   }
  #   
  #   grid_cells_all <- st_transform(grid_cells_all, crs = crs_use) # convert coordinates
  #   
  #   assign('grid_cells_all', grid_cells_all, envir = parent.frame())
  #   
  #   plt_map_cells <-
  #     ggplot(grid_cells_all) +
  #     geom_sf(data = nc_cells,
  #             aes(fill = surface),
  #             show.legend = FALSE) +
  #     scale_fill_manual(values = c('grey','skyblue','skyblue','grey')) +
  #     geom_sf(alpha = 0) +
  #     geom_sf_text(aes(label = cell_index), size = 3) +
  #     theme(
  #       axis.title = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank(),
  #       panel.background = element_rect(fill = 'white', colour = 'white'),
  #       plot.background = element_rect(fill = 'white', colour = 'white')
  #     )
  #   
  #   assign('plt_map_cells', plt_map_cells, envir = parent.frame())
  #   
  # }
  # 
  # assign('ship_poly', ship_poly, envir = parent.frame())
  # assign('krill_poly', krill_poly, envir = parent.frame())
  # assign('chl_poly', chl_poly, envir = parent.frame())
  # assign('sst_poly', sst_poly, envir = parent.frame())
  # assign('pH_poly', pH_poly, envir = parent.frame())
  
}



# Interactive (shiny app) map functions -----------------------------------

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
    longVars <- longVars[longVars != 'Unit']
    
    for(j in 1:length(singleRowVars)) x[[singleRowVars[j]]] <- as.character(x[[singleRowVars[j]]])
    
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
  
  # Set 6 colour scales: blue, green, red, orange, purple and brown
  
  # blues <- c('skyblue', '#00009c') # choose starting colours then adjust to get very light and very dark shades
  # blues_ <- rgb2hsv(col2rgb(blues))
  # blues_[2,1] <- 0.35
  # blues_[3,1] <- 1
  # blues_[1,2] <- 230/360
  # blues_[3,2] <- 1
  # hsv(blues_[1,1],blues_[2,1],blues_[3,1])
  # hsv(blues_[1,2],blues_[2,2],blues_[3,2])
  blues <- c('#A6E5FF', '#002BFF') # colours chosen by altering (commented) hsv values above
  
  # greens <- c('palegreen', 'darkgreen')
  # greens_ <- rgb2hsv(col2rgb(greens))
  # greens_[1,1] <- 150/360
  # greens_[2,1] <- 0.3
  # greens_[1,2] <- 150/360
  # greens_[3,2] <- 0.3
  # hsv(greens_[1,1],greens_[2,1],greens_[3,1])
  # hsv(greens_[1,2],greens_[2,2],greens_[3,2])
  greens <- c('#B0FBD5', '#004D26') # '#004020'
  
  # reds <- c('pink', 'red3')
  # reds_ <- rgb2hsv(col2rgb(reds))
  # reds_[2,1] <- 0.075
  # reds_[3,2] <- 0.6
  # hsv(reds_[1,1],reds_[2,1],reds_[3,1])
  # hsv(reds_[1,2],reds_[2,2],reds_[3,2])
  reds <- c('#FFECEF', '#990000')# '#800000' '#990000')
  
  # oranges <- c('orange', '#fe5a1d')
  # oranges_ <- rgb2hsv(col2rgb(oranges))
  # oranges_[2,1] <- 0.4
  # oranges_[2,2] <- 1
  # oranges_[3,2] <- 0.9
  # hsv(oranges_[1,1],oranges_[2,1],oranges_[3,1])
  # hsv(oranges_[1,2],oranges_[2,2],oranges_[3,2])
  oranges <- c('#FFDB99', '#E63E00')
  
  # purples <- c('#bf94e4', '#4b0082')
  # purples_ <- rgb2hsv(col2rgb(purples))
  # purples_[2,1] <- 0.15
  # purples_[3,1] <- 1
  # purples_[3,2] <- 0.6
  # hsv(purples_[1,1],purples_[2,1],purples_[3,1])
  # hsv(purples_[1,2],purples_[2,2],purples_[3,2])
  purples <- c('#EDD9FF', '#580099')
  
  # browns <- c('#efcc00', '#7b3f00')
  # browns_ <- rgb2hsv(col2rgb(browns))
  # browns_[2,1] <- 0.6
  # browns_[3,1] <- 1
  # browns_[3,2] <- 0.4
  # hsv(browns_[1,1],browns_[2,1],browns_[3,1])
  # hsv(browns_[1,2],browns_[2,2],browns_[3,2])
  browns <- c('#FFE966', '#663400')
  
  # # create a custom colour palette based on the 'Pairs' palette...
  # blues <- c('skyblue', '#00009c')
  # greens <- c('palegreen', 'darkgreen')
  # reds <- c('#ffc1cc', '#c80815')
  # oranges <- c('#f8b878', '#fe5a1d')
  # purples <- c('#dcd0ff', '#663399')
  # browns <- c('lightgoldenrod1', '#7b3f00')
  
  pltColours <- c(blues, reds, purples, oranges, greens, browns)
  # pltColours <- c(blues, greens, reds, oranges, purples, browns)
  
  # pltColours <- brewer.pal(12, 'Paired') # the Set3 and Paired palettes has a maximum of 12 colours
  # "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
  
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
      as.vector(t(mc[seq(ne, 1, -1),])),
      pltColours[seq(1, 11, 2)]
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

# Coordinate axes for map
overlay_coordinate_grid <- function(
    plt, latlim = NULL, lat_increments = NULL, lat_step = 10, lat_step_minor = 0.5 * lat_step,
    lon_vals = NULL, singleLatAxis = FALSE, textSize = 4, minorTicks = TRUE){
  # Round latlim to nearest lat_step
  if(!is.null(latlim)) latlim[2] <- ceiling(latlim[2] / lat_step_minor) * lat_step_minor
  # Set default lat/lon values for axes
  if(is.null(latlim)) latlim <- c(-90, -50)
  if(is.null(lon_vals)) lon_vals <- c(0, 90, 180, 270)
  if(is.null(lat_increments)) lat_increments <- seq(latlim[1], latlim[2], lat_step)
  if(minorTicks) lat_increments_minor <- seq(latlim[1] + lat_step_minor, latlim[2], lat_step)
  #~~~~~~~~~
  # Latitude
  #~~~~~~~~~
  # axis lines
  if(singleLatAxis) lons <- 0 else lons <- lon_vals
  nlon <- length(lons)
  lat_lines <- lapply(1:nlon, function(z) data.frame(lon = rep(lons[z], 2), lat = latlim))
  lat_lines <- lapply(1:nlon, function(z){
    lat_line <- lat_lines[[z]]
    lat_line <- st_linestring(x = as.matrix(lat_line), dim = 'XY')
    lat_line <- st_sfc(lat_line, crs = crs_world)
    lat_line <- st_transform(lat_line, crs_use)
  })
  for(j in 1:nlon) plt <- plt + geom_sf(data = lat_lines[[j]]) # axis line
  # axis ticks
  nlat <- length(lat_increments)
  lr <- 1 / {{lat_increments - min(lat_increments)} / diff(range(lat_increments))}
  lr[1] <- 0
  lon_width0 <- 1.4
  lon_width <- lr * lon_width0
  gridPos <- lapply(1:nlon, function(w){
    lon <- lons[w]
    lapply(1:nlat,
           function(z){
             y <- data.frame(lon = lon + c(-1,1) * lon_width[z],
                             lat = rep(lat_increments[z], 2))
             y <- st_linestring(x = as.matrix(y), dim = 'XY')
             y <- st_sfc(y, crs = crs_world)
             y <- st_transform(y, crs = crs_use)
             y})})
  for(i in 1:nlon)
    for(j in 1:nlat)
      plt <- plt + geom_sf(data = gridPos[[i]][[j]]) # axis ticks
  # Repeat for minor ticks
  if(minorTicks){
    nlat_m <- length(lat_increments_minor)
    lr <- 1 / {{lat_increments_minor - min(lat_increments)} / diff(range(lat_increments))}
    # lr[1] <- 0
    lon_width <- 0.4 * lon_width0
    lon_width <- lr * lon_width
    gridPos_m <- lapply(1:nlon, function(w){
      lon <- lons[w]
      lapply(1:nlat_m,
             function(z){
               y <- data.frame(lon = lon + c(-1,1) * lon_width[z],
                               lat = rep(lat_increments_minor[z], 2))
               y <- st_linestring(x = as.matrix(y), dim = 'XY')
               y <- st_sfc(y, crs = crs_world)
               y <- st_transform(y, crs = crs_use)
               y})})
    for(i in 1:nlon)
      for(j in 1:nlat_m)
        plt <- plt + geom_sf(data = gridPos_m[[i]][[j]]) # axis ticks
  }
  
  
  tick_labs <- lapply(1:nlat, function(z) paste(abs(lat_increments[z]), '* degree ~ S'))
  xnudge <- 0.025 * diff(BBox[c('xmin','xmax')])
  # if(singleLatAxis) i <- 1:nlat else i <- 2:nlat
  for(j in 2:nlat){
    plt <- plt +
      geom_sf_text(
        data = gridPos[[which(lons == 0)]][[j]], label = tick_labs[[j]],
        parse = TRUE, hjust = 0, nudge_x = xnudge, size = textSize)
  }
  
  #~~~~~~~~~~
  # Longitude
  #~~~~~~~~~~
  nlon <- length(lon_vals)
  lon_points <- data.frame(lon = lon_vals, lat = rep(max(latlim), nlon))
  gridPos <- lapply(1:nlon,
                    function(z){
                      y <- lon_points[z,]
                      y <- st_point(as.matrix(y), dim = 'XY')
                      y <- st_sfc(y, crs = crs_world)
                      y <- st_transform(y, crs = crs_use)
                      y})
  lon_labs <- lapply(1:nlon, function(z) paste(lon_vals[z] %% 360, '* degree ~ E')) 
  ynudge <- 0.0175 * diff(BBox[c('ymin','ymax')])
  for(j in 1:nlon){
    vj <- 0
    lv <- lon_vals[j] %% 360
    if({0 <= lv & lv <= 90} | 270 <= lv) inv <- FALSE else inv <- TRUE
    r <- {360 - lv} %% 360 # rotation angle for longitude labels
    xn <- xnudge * sin(lv * pi / 180)
    yn <- ynudge * cos(lv * pi / 180)
    if(inv){
      r <- {r - 180} %% 360
      vj <- 1
    }
    plt <- plt + 
      geom_sf_text(
        data = gridPos[[j]], label = lon_labs[[j]], parse = TRUE, vjust = vj,
        nudge_x = xn, nudge_y = yn, size = textSize, angle = r)
  }
  
  plt <- plt + theme(axis.title = element_blank())
  return(plt)
}

# The main plotting function
make_plot <- function(
    dat, background = 'none', displayEcoregions = FALSE, backgroundOnly = FALSE,
    stationsOnly = FALSE, plasticOnly = FALSE, ptSize = 6, legPtSize = 4,
    legPtSize_fill = 6, ptStroke = 2, alpha = 0.85, polyLineWidth = 0.75,
    contourLineWidth = 0.75, legWidth = 0.25, mapAspectRatio = aspectRatio,
    plotSignificanceContours = FALSE, significanceContours = NULL,
    sst_sigContours = significanceContours,
    pH_sigContours = significanceContours, stationPopSize = FALSE,
    expandBorder = FALSE, showCoordGrid = TRUE, singleLatAxis = FALSE,
    latlim = NULL, axisTextSize = 4, discreteColourScheme = TRUE,
    nColours = NULL, manual_legend_breaks = FALSE, cellBorderCol = 'grey',
    cellBorderWidth = 0.1, overlay_labels = NULL, na.colour = 'grey40',
    background.na.remove = FALSE, na.in.legend = TRUE, ice.colour = 'skyblue',
    land.colour = 'grey85'){
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
  
  library(viridisLite)
  
  nc <- dat$nc
  dat_background <- dat$background
  dat_plastic <- dat$plastic
  dat_stations <- dat$stations
  symbols <- dat$symbols
  # pltColours <- dat$pltColours
  
  if(background.na.remove) dat_background <- dat_background[!is.na(dat_background$value),]
  
  if(is.null(dat_background)) dat_background <- data.frame()
  if(is.null(dat_plastic)) dat_plastic <- data.frame()
  if(is.null(dat_stations)) dat_stations <- data.frame()
  
  anyPlastic <- nrow(dat_plastic) > 0
  anyStations <- nrow(dat_stations) > 0
  anyBackground <- background != 'none'
  
  if(discreteColourScheme & is.null(nColours))
    nColours <- 7 # default number of discrete colours
  
  isNA <- is.na(dat_background$value)
  anyNA <- any(isNA)
  
  e2sci <- function(x, sig.dig = 2, lim = 4){
    # Convert scientific notation using e to standard x10 format.
    # This only seems to work for small numbers -- don't know why...
    oscipen <- options()$scipen
    options(scipen = 1e6)
    isNeg <- x < 0
    x <- abs(x)
    large <- x >= 10 ^ lim
    small <- x <= 10 ^ -{lim-sig.dig+1}
    y <- as.character(x)
    uniSuperDigit <- c('\u2070','\u00B9','\u00B2','\u00B3','\u2074','\u2075',
                       '\u2076','\u2077','\u2078','\u2079') # superscripts 0-9
    uniSuperNeg <- '\u207B'
    if(any(large)){
      y_ <- substr(y[large], 1, 1)
      if(sig.dig > 1) y_ <- paste(y_, substr(y[large], 2, sig.dig), sep = '.')
      y_ <- paste0(y_, 'x10', uniSuperDigit[nchar(y[large])])
      y[large] <- y_
    }
    if(any(small)){
      y_ <- sapply(strsplit(y[small], '\\.'), function(z) z[2])
      countZeros <- gregexpr('0', y_)
      nZeros <- sapply(countZeros, function(z){
        w <- rbind(1:length(z), z)
        sum(w[1,] == w[2,]) + 1})
      y_ <- substr(y_, nZeros, nchar(y_))
      if(sig.dig == 1) y_ <- substr(y_, 1, 1)
      if(sig.dig > 1) y_ <- paste(substr(y_, 1, 1), substr(y_, 2, sig.dig), sep = '.')
      y_ <- paste0(y_, 'x10', uniSuperNeg, uniSuperDigit[nZeros+1])
      y[small] <- y_
    }
    y[isNeg] <- paste0('-', y[isNeg])
    options(scipen = oscipen)
    return(y)
  }
  #~~~~~~~~~~~~~~~~~
  # Background layer
  #~~~~~~~~~~~~~~~~~
  background.legend.id <- 'background.legend.title'
  switch(background,
         
         none = {
           plt_background <- ggplot()
         },
         
         krill = {
           legTitle <- 'Krill\n(number m\u207B\u00B2)'
           
           if(!discreteColourScheme){
             # Continuous colour scheme
             
             # To use log-scale and distinguish zeros from NAs, replace zeros
             # with half the minimum positive value
             dat_background$value[dat_background$value == 0] <- 0.5 * min(dat_background$value[dat_background$value > 0], na.rm = TRUE)
             
             # # handle zeros in log-scale using NA functionality
             # dat_background$value[dat_background$value == 0] <- NA
             # zero_col <- viridis_pal(begin = 0, end = 0, option = 'plasma')(1)
             
             plt_background <- 
               ggplot() +
               geom_sf(data = dat_background,
                       aes(fill = value), colour = cellBorderCol,
                       linewidth = cellBorderWidth)
             
             if(!na.in.legend & anyNA){
               plt_background_leg <- ggplot() +
                 geom_sf(data = subset(dat_background, !isNA),
                         aes(fill = value), colour = cellBorderCol,
                         linewidth = cellBorderWidth)}
             

             if(exists('weblink_krill')){
               plt_background <- plt_background +
                 scale_fill_viridis_c_interactive(
                   option = 'plasma', trans = 'log10',
                   labels = comma, name = label_interactive(
                     legTitle, data_id = background.legend.id,
                     hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                     onclick = paste0("window.open(`",
                                      as.character(weblink_krill), "`);"),
                     tooltip = weblink_krill),
                   na.value = na.colour)
               if(!na.in.legend & anyNA){
                 plt_background_leg <- plt_background_leg +
                   scale_fill_viridis_c_interactive(
                     option = 'plasma', trans = 'log10',
                     labels = comma, name = label_interactive(
                       legTitle, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_krill), "`);"),
                       tooltip = weblink_krill),
                     na.value = na.colour)
               }
               }else{
               plt_background <- plt_background +
                 scale_fill_viridis_c(
                   option = 'plasma', name = legTitle, trans = 'log10',
                   na.value = na.colour, labels = comma)
               if(!na.in.legend & anyNA){
                 plt_background_leg <- plt_background_leg +
                   scale_fill_viridis_c(
                     option = 'plasma', name = legTitle, trans = 'log10',
                     na.value = na.colour, labels = comma)
               }
             }
             
             
           }else{
             # Discrete colour scheme
             # hard-code a colour scale for krill -- set this to TRUE to create
             # colour scale automatically
             autoScaleColours <- FALSE
             if(autoScaleColours){
               # handle zeros in log-scale by setting to minimum positive value
               dat_background$value[dat_background$value == 0] <- min(
                 dat_background$value[dat_background$value > 0])
               dat_background$value_log10 <- log10(dat_background$value)
               
               # Create break points
               dat_background$breaks <- cut(dat_background$value_log10, nColours)
               # Adjust breaks for pretty legend labels
               brks <- gsub(']', '', gsub('\\(', '',
                                          levels(dat_background$breaks)))
               brks <- unique(as.numeric(unlist(strsplit(brks, ','))))
               mag <- 10 ^ floor(brks)
               mag_0.5 <- 0.5 * mag
               brks <- 10 ^ brks
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_0.5[i]) * mag_0.5[i]
               brks[1] <- floor(brks[1] / mag_0.5[1]) * mag_0.5[1]
               brks[nbrks] <- ceiling(brks[nbrks] / mag_0.5[nbrks]) * mag_0.5[nbrks]
               brks_ <- e2sci(brks)
               legLabs <- sapply(
                 1:nColours, function(z) paste0('[', brks_[z], ', ', brks_[z+1], ')'))
               legLabs[nbrks-1] <- gsub('\\)', ']', legLabs[nbrks-1])
               brks <- log10(brks)
               dat_background$breaks <- cut(dat_background$value_log10, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               
               levelsInData <- levels(dat_background$breaks) %in% dat_background$breaks
               
               dat_background$breaks <- factor(dat_background$breaks,
                                               levels = c(rev(levels(dat_background$breaks))))
               
               plt_background <- 
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth)
               
               if(!na.in.legend & anyNA){
                 plt_background_leg <- ggplot() +
                   geom_sf(data = subset(dat_background, !is.na(value)),
                           aes(fill = breaks), colour = cellBorderCol,
                           linewidth = cellBorderWidth)}
               
               if(exists('weblink_krill')){
                 plt_background <- plt_background + 
                   scale_fill_viridis_d_interactive(
                     option = 'plasma',
                     direction = -1,
                     name = label_interactive(
                       legTitle, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_krill), "`);"),
                       tooltip = weblink_krill),
                     na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
                 if(!na.in.legend & anyNA){
                   plt_background_leg <- plt_background_leg +
                     scale_fill_viridis_d_interactive(
                       option = 'plasma',
                       direction = -1,
                       name = label_interactive(
                         legTitle, data_id = background.legend.id,
                         hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                         onclick = paste0("window.open(`",
                                          as.character(weblink_krill), "`);"),
                         tooltip = weblink_krill),
                       na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
                 }
                 
               }else{
                 plt_background <- plt_background + 
                   scale_fill_viridis_d(
                     option = 'plasma', direction = -1, name = legTitle,
                     na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
                 if(!na.in.legend & anyNA){
                   plt_background_leg <- plt_background_leg +
                     scale_fill_viridis_d(
                       option = 'plasma', direction = -1, name = legTitle,
                       na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
                 }
               }
               
             }else{

               # Manually specify the colour scale break points for krill so that
               # we avoid too many fractional values
               brks <- c(0,0.1,1,5,10,25,50,200)
               nb <- length(brks)
               legLabs <- sapply(
                 1:{nb-1}, function(z) paste0('[', brks[z], ', ', brks[z+1], ')'))
               legLabs[nb-1] <- gsub('\\)', ']', legLabs[nb-1])
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)

               levelsInData <- levels(dat_background$breaks) %in% dat_background$breaks

               dat_background$breaks <- factor(dat_background$breaks,
                                               levels = rev(levels(dat_background$breaks)))

               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth)

               if(!na.in.legend & anyNA){
                 plt_background_leg <- ggplot() +
                   geom_sf(data = subset(dat_background, !is.na(value)),
                           aes(fill = breaks), colour = cellBorderCol,
                           linewidth = cellBorderWidth)}

               if(exists('weblink_krill')){
                 plt_background <- plt_background +
                   scale_fill_viridis_d_interactive(
                     option = 'plasma',
                     direction = -1,
                     name = label_interactive(
                       legTitle, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_krill), "`);"),
                       tooltip = weblink_krill),
                     na.value = na.colour,
                     labels = c(levels(dat_background$breaks)[levelsInData], 'no data')
                     # drop = FALSE
                     )
                 if(!na.in.legend & anyNA){
                   plt_background_leg <- plt_background_leg +
                     scale_fill_viridis_d_interactive(
                       option = 'plasma',
                       direction = -1,
                       name = label_interactive(
                         legTitle, data_id = background.legend.id,
                         hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                         onclick = paste0("window.open(`",
                                          as.character(weblink_krill), "`);"),
                         tooltip = weblink_krill),
                       na.value = na.colour,
                       labels = c(levels(dat_background$breaks)[levelsInData], 'no data')
                       # drop = FALSE
                       )
                 }
               }else{
                 plt_background <- plt_background +
                   scale_fill_viridis_d(
                     option = 'plasma', direction = -1, name = legTitle,
                     na.value = na.colour,
                     labels = c(levels(dat_background$breaks)[levelsInData], 'no data')
                     # drop = FALSE
                     )
                 if(!na.in.legend & anyNA){
                   plt_background_leg <- plt_background_leg +
                     scale_fill_viridis_d(
                       option = 'plasma', direction = -1, name = legTitle,
                       na.value = na.colour,
                       labels = c(levels(dat_background$breaks)[levelsInData], 'no data')
                       # drop = FALSE
                       )
                 }
               }
               
             }
           }
         },
         
         chl = {
           legTitle <- 'Chlorophyll \uD835\uDC4E\n(mg m\u207B\u00B3)'
           
           if(!discreteColourScheme){
             # Continuous colour scheme
             plt_background <-
               ggplot() +
               geom_sf(data = dat_background,
                       aes(fill = value), colour = cellBorderCol,
                       linewidth = cellBorderWidth)
             if(exists('weblink_chlorophyll')){
               plt_background <- plt_background + 
                 scale_fill_viridis_c_interactive(
                   option = 'viridis', trans = 'log10', name = label_interactive(
                     legTitle, data_id = background.legend.id,
                     hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                     onclick = paste0("window.open(`", as.character(
                       weblink_chlorophyll), "`);"),
                     tooltip = weblink_chlorophyll))
             }else{
               plt_background <- plt_background +
                 scale_fill_viridis_c(option = 'viridis', trans = 'log10',
                                      name = legTitle)
             }
           }else{
             # Discrete colour scheme
             dat_background$value_log10 <- log10(dat_background$value)
             # Create break points
             dat_background$breaks <- cut(dat_background$value_log10, nColours)
             # Adjust breaks for pretty legend labels
             brks <- gsub(']', '', gsub('\\(', '',
                                        levels(dat_background$breaks)))
             brks <- unique(as.numeric(unlist(strsplit(brks, ','))))
             mag <- 10 ^ floor(brks)
             mag_0.5 <- 0.5 * mag
             brks <- 10 ^ brks
             nbrks <- length(brks)
             i <- 2:{nbrks-1}
             brks[i] <- round(brks[i] / mag_0.5[i]) * mag_0.5[i]
             brks[1] <- floor(brks[1] / mag_0.5[1]) * mag_0.5[1]
             brks[nbrks] <- ceiling(brks[nbrks] / mag_0.5[nbrks]) * mag_0.5[nbrks]
             legLabs <- sapply(
               1:nColours, function(z) paste0('[', brks[z], ', ', brks[z+1], ')'))
             legLabs[nbrks-1] <- gsub('\\)', ']', legLabs[nbrks-1])
             brks <- log10(brks)
             dat_background$breaks <- cut(dat_background$value_log10, brks,
                                          labels = legLabs,
                                          include.lowest = TRUE)
             
             dat_background$breaks <- factor(dat_background$breaks,
                                             levels = rev(levels(dat_background$breaks)))

             plt_background <-
               ggplot() +
               geom_sf(data = dat_background,
                       aes(fill = breaks), colour = cellBorderCol,
                       linewidth = cellBorderWidth) +
               guides(fill = guide_legend(reverse = FALSE))
             if(exists('weblink_chlorophyll')){
               plt_background <- plt_background + 
                 scale_fill_viridis_d_interactive(
                   option = 'viridis', direction = -1, name = label_interactive(
                     legTitle, data_id = background.legend.id,
                     hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                     onclick = paste0("window.open(`", as.character(
                       weblink_chlorophyll), "`);"),
                     tooltip = weblink_chlorophyll),
                   na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
             }else{
               plt_background <- plt_background +
                 scale_fill_viridis_d(
                   option = 'viridis', direction = -1, name = legTitle,
                   na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
             }
           }
         },
         
         sst = {
           if(exists('sstType')){
             legLabel <- switch(sstType,
                                trend = 'SST trend\n(\u00B0C year\u207B\u00B9)',
                                anomaly = 'SST anomaly\n(\u00B0C)')
           }else{
             legLabel <- 'SST'
           }
           weblink_exists <- exists('weblink_sst')
           nv <- dat_background$value[dat_background$value < 0]
           pv <- dat_background$value[dat_background$value > 0]
           any_positive <- length(pv) > 0
           any_negative <- length(nv) > 0
           ncol_sst <- 11
           if(!is.null(nColours)) ncol_sst <- nColours
           if(ncol_sst %% 2 != 1) ncol_sst + 1 # ensure odd number of colours
           ncol_sst_ <- 0.5 * {ncol_sst - 1}
           colfunc <- colorRampPalette(c('dodgerblue4', 'white', 'red3'))
           colfunc_neg <- colorRampPalette(c('white', 'dodgerblue4'))
           colfunc_pos <- colorRampPalette(c('red3', 'white'))
           Cols <- colfunc(ncol_sst)
           Cols_neg <- colfunc_neg(ncol_sst)
           Cols_pos <- colfunc_pos(ncol_sst)
           cval <- c(
             quantile(nv, seq(0, 1, length = ncol_sst_ + 1))[1:ncol_sst_],
             0,
             quantile(pv, seq(0, 1, length = ncol_sst_ + 1))[2:{ncol_sst_ + 1}])
           if(!discreteColourScheme){
             # Continuous colour scale
             if(any_positive & any_negative){
               cval <- {cval - min(cval)} / diff(range(cval))
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = value), colour = cellBorderCol,
                         linewidth = cellBorderWidth)
               if(weblink_exists){
                 plt_background <- plt_background +
                   scale_fill_gradientn_interactive(
                     colours = Cols, values = cval, name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_sst), "`);"),
                       tooltip = weblink_sst), 
                     na.value = na.colour)
               }else{
                 plt_background <- plt_background +
                   scale_fill_gradientn(colours = Cols, values = cval,
                                        name = legLabel, na.value = na.colour)
               }
             }else{if(any_negative){
               cval <- {cval[1:{ncol_sst_+1}] - min(cval[1:{ncol_sst_+1}])} / 
                 diff(range(cval[1:{ncol_sst_+1}]))
               cval <- {cval - min(cval)} / diff(range(cval))
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = value, colour = value)) +
                 scale_colour_gradientn(colours = Cols[1:{ncol_sst_+1}],
                                        values = cval) +
                 guides(colour = 'none')
               if(weblink_exists){
                 plt_background <- plt_background +
                   scale_fill_gradientn_interactive(
                     colours = Cols[1:{ncol_sst_+1}], values = cval,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_sst), "`);"),
                       tooltip = weblink_sst))
               }else{
                 plt_background <- plt_background +
                   scale_fill_gradientn(colours = Cols[1:{ncol_sst_+1}],
                                      values = cval, name = legLabel)
               }
             }else{if(any_positive){
               cval <- {cval[{ncol_sst_+1}:ncol_sst] - 
                   min(cval[{ncol_sst_+1}:ncol_sst])} / 
                 diff(range(cval[{ncol_sst_+1}:ncol_sst]))
               cval <- {cval - min(cval)} / diff(range(cval))
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = value), colour = cellBorderCol)
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_gradientn_interactive(
                     colours = Cols[{ncol_sst_+1}:ncol_sst], values = cval,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_sst), "`);"),
                       tooltip = weblink_sst),
                     na.value = na.colour)
               }else{
                 plt_background <- plt_background + 
                   scale_fill_gradientn(colours = Cols[{ncol_sst_+1}:ncol_sst],
                                        values = cval, name = legLabel,
                                        na.value = na.colour)
               }
             }}}
           }else{
             # Discrete colour scheme
             if(any_positive & any_negative){
               # Create break points (all non-significant trends have value 0)
               nq <- quantile(nv, seq(0, 1, length = ncol_sst_ + 1))
               pq <- quantile(pv, seq(0, 1, length = ncol_sst_ + 1))
               brks <- unname(c(nq, pq))
               # Adjust breaks for pretty legend labels
               mag <- 10 ^ floor(log10(abs(brks)))
               mag_0.5 <- 0.5 * mag
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_0.5[i]) * mag_0.5[i]
               brks[1] <- floor(brks[1] / mag_0.5[1]) * mag_0.5[1]
               brks[nbrks] <- ceiling(brks[nbrks] / mag_0.5[nbrks]) * mag_0.5[nbrks]
               brks_ <- e2sci(brks)
               legLabs1 <- sapply(1:ncol_sst_, function(z)
                 paste0('(', brks_[z+1], ', ', brks_[z], ']'))
               legLabs2 <- sapply({ncol_sst_ + 2}:ncol_sst, function(z)
                 paste0('(', brks_[z], ', ', brks_[z+1], ']'))
               legLabs1[ncol_sst_] <- gsub('\\(', '[', legLabs1[ncol_sst_])
               legLabs2[1] <- gsub('\\(', '[', legLabs2[1])
               legLabs <- c(legLabs1, '0', legLabs2)
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               
               dat_background$breaks <- factor(dat_background$breaks,
                                               levels = rev(levels(dat_background$breaks)))
               
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth) +
                 guides(fill = guide_legend(reverse = FALSE))
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_manual_interactive(
                     values = rev(Cols),
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`", as.character(
                         weblink_sst), "`);"),
                       tooltip = weblink_sst),
                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }else{
                 plt_background <- plt_background +
                   scale_fill_manual(values = rev(Cols), name = legLabel,
                                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }
             }else{if(any_negative){
               # Create break points (all non-significant trends have value 0)
               nq <- quantile(nv, seq(0, 1, length = ncol_sst))
               brks <- unname(c(nq, 0))
               # Adjust breaks for pretty legend labels
               mag <- 10 ^ floor(log10(abs(brks)))
               mag_0.5 <- 0.5 * mag
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_0.5[i]) * mag_0.5[i]
               brks[1] <- floor(brks[1] / mag_0.5[1]) * mag_0.5[1]
               # brks[nbrks] <- ceiling(brks[nbrks] / mag_0.5[nbrks]) * mag_0.5[nbrks]
               brks_ <- e2sci(brks)
               legLabs <- sapply(1:{ncol_sst - 1}, function(z)
                 paste0('(', brks_[z+1], ', ', brks_[z], ']'))
               legLabs[nbrks - 2] <- gsub('\\(', '[', legLabs[nbrks - 2])
               legLabs[nbrks - 1] <- '0'
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               
               dat_background$breaks <- factor(dat_background$breaks,
                                               levels = rev(levels(dat_background$breaks)))
               
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth) +
                 guides(fill = guide_legend(reverse = TRUE))
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_manual_interactive(
                     values = Cols_neg,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`", as.character(
                         weblink_sst), "`);"),
                       tooltip = weblink_sst),
                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }else{
                 plt_background <- plt_background +
                   scale_fill_manual(values = Cols_neg, name = legLabel,
                                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }
             }else{if(any_positive){
               # Create break points (all non-significant trends have value 0)
               pq <- quantile(pv, seq(0, 1, length = ncol_sst))
               brks <- unname(c(0, pq))
               # Adjust breaks for pretty legend labels
               mag <- 10 ^ floor(log10(abs(brks)))
               mag_0.5 <- 0.5 * mag
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_0.5[i]) * mag_0.5[i]
               brks[nbrks] <- ceiling(brks[nbrks] / mag_0.5[nbrks]) * mag_0.5[nbrks]
               brks_ <- e2sci(brks)
               legLabs <- sapply(2:ncol_sst, function(z)
                 paste0('(', brks_[z], ', ', brks_[z+1], ']'))
               legLabs[1] <- gsub('\\(', '[', legLabs[1])
               legLabs <- c('0', legLabs)
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               
               dat_background$breaks <- factor(dat_background$breaks,
                                               levels = rev(levels(dat_background$breaks)))
               
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth) +
                 guides(fill = guide_legend(reverse = FALSE))
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_manual_interactive(
                     values = Cols_pos,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`", as.character(
                         weblink_sst), "`);"),
                       tooltip = weblink_sst),
                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }else{
                 plt_background <- plt_background +
                   scale_fill_manual(values = Cols_pos, name = legLabel,
                                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }
             }}}
           }
         },
         
         pH = {
           if(exists('pHType')){
             legLabel <- switch(pHType,
                                trend = 'pH trend\n(year\u207B\u00B9)',
                                anomaly = 'pH anomaly\n(dimensionless)')
           }else{
             legLabel <- 'pH'
           }
           weblink_exists <- exists('weblink_pH')
           notNA <- !is.na(dat_background$value)
           nv <- dat_background$value[notNA & dat_background$value < 0]
           pv <- dat_background$value[notNA & dat_background$value > 0]
           any_positive <- length(pv) > 0
           any_negative <- length(nv) > 0
           any_zero <- any(notNA & dat_background$value == 0)
           ncol_pH <- 11
           if(!is.null(nColours)) ncol_pH <- nColours
           if(ncol_pH %% 2 != 1) ncol_pH + 1 # ensure odd number of colours
           ncol_pH_ <- 0.5 * {ncol_pH - 1}
           colfunc <- colorRampPalette(c('darkorange2', 'white', 'darkorchid4'))
           colfunc_neg <- colorRampPalette(c('white', 'darkorange2'))
           colfunc_pos <- colorRampPalette(c('darkorchid4', 'white'))
           Cols <- colfunc(ncol_pH)
           Cols_neg <- colfunc_neg(ncol_pH)
           Cols_pos <- colfunc_pos(ncol_pH)
           cval <- c(
             quantile(nv, seq(0, 1, length = ncol_pH_ + 1))[1:ncol_pH_],
             0,
             quantile(pv, seq(0, 1, length = ncol_pH_ + 1))[2:{ncol_pH_ + 1}])
           if(!discreteColourScheme){
             # Continuous colour scale
             if(any_positive & any_negative){
               cval <- {cval - min(cval)} / diff(range(cval))
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = value), colour = cellBorderCol)               
               if(weblink_exists){
                 plt_background <- plt_background +
                   scale_fill_gradientn_interactive(
                     colours = Cols, values = cval, name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_pH), "`);"),
                       tooltip = weblink_pH),
                     na.value = na.colour)
               }else{
                 plt_background <- plt_background +
                   scale_fill_gradientn(colours = Cols, values = cval,
                                        name = legLabel)
               }
             }else{if(any_negative){
               cval <- {cval[1:{ncol_pH_+1}] - min(cval[1:{ncol_pH_+1}])} / 
                 diff(range(cval[1:{ncol_pH_+1}]))
               cval <- {cval - min(cval)} / diff(range(cval))
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = value), colour = cellBorderCol)
               if(weblink_exists){
                 plt_background <- plt_background +
                   scale_fill_gradientn_interactive(
                     colours = rev(Cols[1:{ncol_pH_+1}]), values = cval,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_pH), "`);"),
                       tooltip = weblink_pH),
                     na.value = na.colour) +
                   guides(fill = guide_colorbar(reverse = TRUE))
               }else{
                 plt_background <- plt_background +
                   scale_fill_gradientn(colours = rev(Cols[1:{ncol_pH_+1}]),
                                        values = cval, name = legLabel) +
                   guides(fill = guide_colorbar(reverse = TRUE))
                 
               }
             }else{if(any_positive){
               cval <- {cval[{ncol_pH_+1}:ncol_pH] - 
                   min(cval[{ncol_pH_+1}:ncol_pH])} / 
                 diff(range(cval[{ncol_pH_+1}:ncol_pH]))
               cval <- {cval - min(cval)} / diff(range(cval))
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = value), colour = cellBorderCol)
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_gradientn_interactive(
                     colours = rev(Cols[{ncol_pH_+1}:ncol_pH]), values = cval,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        as.character(weblink_pH), "`);"),
                       tooltip = weblink_pH),
                     na.value = na.colour)
               }else{
                 plt_background <- plt_background + 
                   scale_fill_gradientn(colours = rev(Cols[{ncol_pH_+1}:ncol_pH]),
                                        values = cval, name = legLabel,
                                        na.value = na.colour)
               }
             }}}
           }else{
             # Discrete colour scheme
             if(any_positive & any_negative){
               # Create break points (all non-significant trends have value 0)
               nq <- quantile(nv, seq(0, 1, length = ncol_pH_ + 1))
               pq <- quantile(pv, seq(0, 1, length = ncol_pH_ + 1))
               brks <- unname(c(nq, pq))
               # Adjust breaks for pretty legend labels
               mag <- 10 ^ floor(log10(abs(brks)))
               mag_ <- 0.1 * mag
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_[i]) * mag_[i]
               brks[1] <- floor(brks[1] / mag_[1]) * mag_[1]
               brks[nbrks] <- ceiling(brks[nbrks] / mag_[nbrks]) * mag_[nbrks]
               brks_ <- e2sci(brks)
               legLabs1 <- sapply(1:ncol_pH_, function(z)
                 paste0('(', brks_[z+1], ', ', brks_[z], ']'))
               legLabs2 <- sapply({ncol_pH_ + 2}:ncol_pH, function(z)
                 paste0('(', brks_[z], ', ', brks_[z+1], ']'))
               legLabs1[ncol_pH_] <- gsub('\\(', '[', legLabs1[ncol_pH_])
               legLabs2[1] <- gsub('\\(', '[', legLabs2[1])
               legLabs <- c(legLabs1, '0', legLabs2)
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               
               dat_background$breaks <- factor(dat_background$breaks, 
                                               levels = levels(dat_background$breaks))
               # dat_background$breaks <- factor(dat_background$breaks, 
               #                                 levels = rev(levels(dat_background$breaks)))
               
               levelsInData <- levels(dat_background$breaks) %in% dat_background$breaks
               
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth) +
                 guides(fill = guide_legend(reverse = FALSE))
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_manual_interactive(
                     values = Cols[levelsInData],
                     # values = rev(Cols[levelsInData]),
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`", as.character(
                         weblink_pH), "`);"),
                       tooltip = weblink_pH),
                     na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data')#,
#                     drop = FALSE
                     )
               }else{
                 plt_background <- plt_background +
                   scale_fill_manual(values = Cols[levelsInData], name = legLabel,
                                     na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
               }
             }else{if(any_negative){
               # Create break points (all non-significant trends have value 0)
               brks <- unname(quantile(nv, seq(0, 1, length = ncol_pH)))
               if(any_zero) brks <- c(brks, 0)
               # Adjust breaks for pretty legend labels
               mag <- 10 ^ floor(log10(abs(brks)))
               mag_ <- 0.1 * mag
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_[i]) * mag_[i]
               brks[1] <- floor(brks[1] / mag_[1]) * mag_[1]
               brks[nbrks] <- ceiling(brks[nbrks] / mag_[nbrks]) * mag_[nbrks]
               brks <- unique(brks)
               nbrks <- length(brks)
               brks_ <- e2sci(brks)
               legLabs <- sapply(1:{nbrks - 1}, function(z)
                 paste0('(', brks_[z+1], ', ', brks_[z], ']'))
               j <- ifelse(any_zero, nbrks - 2, nbrks - 1)
               legLabs[j] <- gsub('\\(', '[', legLabs[j])
               if(any_zero) legLabs[nbrks - 1] <- '0'
               ncol <- length(legLabs)
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               if(any_zero){
                 Cols_neg <- colfunc_neg(ncol)}else{
                   Cols_neg <- colfunc_neg(ncol+1)[2:{ncol+1}]}
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth) +
                 guides(fill = guide_legend(reverse = FALSE))
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_manual_interactive(
                     values = Cols_neg,
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`", as.character(
                         weblink_pH), "`);"),
                       tooltip = weblink_pH),
                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }else{
                 plt_background <- plt_background +
                   scale_fill_manual(values = Cols_neg, name = legLabel,
                                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }
             }else{if(any_positive){
               # Create break points (all non-significant trends have value 0)
               brks <- unname(quantile(pv, seq(0, 1, length = ncol_pH)))
               if(any_zero) brks <- c(0, brks)
               # Adjust breaks for pretty legend labels
               mag <- 10 ^ floor(log10(abs(brks)))
               mag_ <- 0.1 * mag
               nbrks <- length(brks)
               i <- 2:{nbrks-1}
               brks[i] <- round(brks[i] / mag_[i]) * mag_[i]
               brks[1] <- floor(brks[1] / mag_[1]) * mag_[1]
               brks[nbrks] <- ceiling(brks[nbrks] / mag_[nbrks]) * mag_[nbrks]
               brks <- unique(brks)
               nbrks <- length(brks)
               brks_ <- e2sci(brks)
               legLabs <- sapply(1:{nbrks-1}, function(z)
                 paste0('(', brks_[z], ', ', brks_[z+1], ']'))
               legLabs[1] <- gsub('\\(', '[', legLabs[1])
               if(any_zero) legLabs <- c('0', legLabs)
               ncol <- length(legLabs)
               dat_background$breaks <- cut(dat_background$value, brks,
                                            labels = legLabs,
                                            include.lowest = TRUE)
               dat_background$breaks <- factor(dat_background$breaks, 
                                               levels = rev(levels(dat_background$breaks)))
               
               
               if(any_zero){
                 Cols_pos <- colfunc_pos(ncol)}else{
                   Cols_pos <- colfunc_pos(ncol+1)[1:ncol]}
               plt_background <-
                 ggplot() +
                 geom_sf(data = dat_background,
                         aes(fill = breaks), colour = cellBorderCol,
                         linewidth = cellBorderWidth) +
                 guides(fill = guide_legend(reverse = FALSE))
               if(weblink_exists){
                 plt_background <- plt_background + 
                   scale_fill_manual_interactive(
                     values = rev(Cols_pos),
                     name = label_interactive(
                       legLabel, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`", as.character(
                         weblink_pH), "`);"),
                       tooltip = weblink_pH),
                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }else{
                 plt_background <- plt_background +
                   scale_fill_manual(values = rev(Cols_pos), name = legLabel,
                                     na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
               }
             }}}
           }
         },
         
         ship = {
           
           if(exists(shipSummaryDataOrRaw)){
             
             weblink_exists <- exists('weblink_shipping')
             
             if(shipSummaryDataOrRaw == 'summary'){
               # Summary data were not used. See option for raw data below...
               shipType <- unique(dat_background$ship_class)
               ShipType <- paste0(toupper(substr(shipType,1,1)),
                                  substr(shipType,2,nchar(shipType)))
               leg_lab <- paste0('Ship time:', '\n', shipType,
                                 ' vessels', '\n', '(days year\u207B\u00B9)')
               vesselsPresent <- !is.na(dat_background$total_time)
               if(weblink_exists){
                 plt_background <- 
                   ggplot() + 
                   geom_sf(data = dat_background[vesselsPresent,],
                           aes(fill = total_time),
                           colour = cellBorderCol,
                           linewidth = polyLineWidth) +
                   scale_fill_viridis_c_interactive(
                     option = 'mako', trans = 'log10', direction = -1,
                     name = label_interactive(
                       leg_lab, data_id = background.legend.id,
                       hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                       onclick = paste0("window.open(`",
                                        weblink_shipping, "`);"),
                       tooltip = weblink_shipping)) +
                   geom_sf(data = dat_background[!vesselsPresent,],
                           fill = 'white', colour = cellBorderCol,
                           linewidth = polyLineWidth)
               }else{
                 plt_background <- 
                   ggplot() + 
                   geom_sf(data = dat_background[vesselsPresent,],
                           aes(fill = total_time),
                           colour = cellBorderCol,
                           linewidth = polyLineWidth) +
                   scale_fill_viridis_c(option = 'mako', trans = 'log10',
                                        direction = -1, name = leg_lab) +
                   geom_sf(data = dat_background[!vesselsPresent,],
                           fill = 'white', colour = cellBorderCol,
                           linewidth = polyLineWidth)
               }
             }
             
             if(shipSummaryDataOrRaw == 'raw'){
               
               shipActivity <- unique(dat_background$activity)
               
               v <- unique(dat_background$variable)
               if(length(v) != 1) error("Ship traffic metric must be specified")
               leg_lab <- switch(v,
                                 `ship time` = paste0('Ship traffic:', '\n', shipActivity,
                                                      ' vessels', '\n', '(ship days year\u207B\u00B9)'),
                                 `person time` = paste0('Ship traffic:', '\n', shipActivity,
                                                        ' vessels', '\n', '(person days year\u207B\u00B9)'))
               
               if(!discreteColourScheme){
                 # Continuous colour scheme
                 plt_background <- 
                   ggplot() + 
                   geom_sf(data = dat_background,
                           aes(fill = value), colour = cellBorderCol)
                 
                 if(!na.in.legend & anyNA){
                   plt_background_leg <- ggplot() +
                     geom_sf(data = subset(dat_background, !is.na(value)),
                             aes(fill = value), colour = cellBorderCol)}

                 if(weblink_exists){
                   plt_background <- plt_background + 
                     scale_fill_viridis_c_interactive(
                       option = 'mako', trans = 'log10', direction = -1,
                       name = label_interactive(
                         leg_lab, data_id = background.legend.id,
                         hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                         onclick = paste0("window.open(`",
                                          weblink_shipping, "`);"),
                         tooltip = weblink_shipping),
                       na.value = na.colour)
                   if(!na.in.legend & anyNA){
                     plt_background_leg <- plt_background_leg +
                       scale_fill_viridis_c_interactive(
                         option = 'mako', trans = 'log10', direction = -1,
                         name = label_interactive(
                           leg_lab, data_id = background.legend.id,
                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                           onclick = paste0("window.open(`",
                                            weblink_shipping, "`);"),
                           tooltip = weblink_shipping),
                         na.value = na.colour)}
                 }else{
                   plt_background <- plt_background +
                     scale_fill_viridis_c(option = 'mako', trans = 'log10',
                                          direction = -1, name = leg_lab,
                                          na.value = na.colour)
                   if(!na.in.legend & anyNA){
                     plt_background_leg <- plt_background_leg +
                       scale_fill_viridis_c(option = 'mako', trans = 'log10',
                                            direction = -1, name = leg_lab,
                                            na.value = na.colour)}
                 }
               }else{
                 # Discrete colour scheme
                 isNA <- is.na(dat_background$value)
                 isZero <- !isNA & dat_background$value == 0
                 anyZeros <- any(isZero)
                 dat_background$value_log10 <- log10(dat_background$value)
                 mv <- round(dat_background$value[!isNA & !isZero])
                 mv <- min(mv[mv > 0]) # minimum positive integer
                 
                 # Create break points
                 if(all(manual_legend_breaks == FALSE)){
                   if(!anyZeros){
                     i <- !isNA & dat_background$value > mv
                     brks <- seq(min(dat_background$value_log10[i]), max(dat_background$value_log10[i]), length.out = nColours)
                     mag <- 10 ^ floor(brks)
                     mag_ <- 0.5 * mag
                     brks <- 10 ^ brks
                     i <- 1:{nColours-1}
                     brks[i] <- round(brks[i] / mag[i]) * mag[i]
                     i <- nColours
                     brks[i] <- ceiling(brks[i] / mag_[i]) * mag_[i]
                     brks <- c(0, brks)
                     # Adjust breaks for pretty legend labels
                     brks_ <- brks
                     brks_[brks > 0] <- e2sci(brks[brks > 0])
                     legLabs <- sapply(
                       1:nColours, function(z) paste0('[', brks_[z], ', ', brks_[z+1], ')'))
                     legLabs[nColours-1] <- gsub('\\)', ']', legLabs[nColours-1])
                     dat_background$breaks <- cut(dat_background$value, brks, 
                                                  labels = legLabs, include.lowest = TRUE)
                   }else{
                     i <- !isNA & dat_background$value > mv
                     brks <- seq(min(dat_background$value_log10[i]), max(dat_background$value_log10[i]), length.out = nColours-1)
                     mag <- 10 ^ floor(brks)
                     mag_ <- 0.5 * mag
                     brks <- 10 ^ brks
                     i <- 1:{nColours-2}
                     brks[i] <- round(brks[i] / mag[i]) * mag[i]
                     i <- nColours-1
                     brks[i] <- ceiling(brks[i] / mag_[i]) * mag_[i]
                     brks <- c(0, brks)
                     # Adjust breaks for pretty legend labels
                     brks_ <- brks
                     brks_[brks > 0] <- e2sci(brks[brks > 0])
                     legLabs <- sapply(
                       1:{nColours-1}, function(z) paste0('[', brks_[z], ', ', brks_[z+1], ')'))
                     legLabs[nColours-1] <- gsub('\\)', ']', legLabs[nColours-1])
                     legLabs <- c('0', legLabs)
                     
                     i <- !isZero & !isNA
                     v <- as.character(cut(dat_background$value[i], brks, 
                                           labels = legLabs[-1], include.lowest = TRUE))
                     dat_background$breaks[i] <- v
                     dat_background$breaks[isZero] <- legLabs[1]
                     dat_background$breaks <- factor(dat_background$breaks, levels = legLabs)
                   }
                 }else{
                   # break points are input as function argument
                   brks <- manual_legend_breaks
                   brks_ <- brks
                   brks_[brks > 0] <- e2sci(brks[brks > 0])
                   if(!anyZeros){
                     legLabs <- sapply(
                       1:nColours, function(z) paste0('[', brks_[z], ', ', brks_[z+1], ')'))
                     legLabs[nColours-1] <- gsub('\\)', ']', legLabs[nColours-1])
                     dat_background$breaks <- cut(dat_background$value, brks, 
                                                  labels = legLabs, include.lowest = TRUE)                   
                   }else{
                     legLabs <- sapply(
                       1:{nColours-1}, function(z) paste0('[', brks_[z], ', ', brks_[z+1], ')'))
                     legLabs[nColours-1] <- gsub('\\)', ']', legLabs[nColours-1])
                     legLabs <- c('0', legLabs)
                     
                     i <- !isZero & !isNA
                     v <- as.character(cut(dat_background$value[i], brks, 
                                           labels = legLabs[-1], include.lowest = TRUE))
                     dat_background$breaks[i] <- v
                     dat_background$breaks[isZero] <- legLabs[1]
                     dat_background$breaks <- factor(dat_background$breaks, levels = legLabs)
                   }
                 }
                 
                 dat_background$breaks <- factor(dat_background$breaks,
                                                 levels = rev(levels(dat_background$breaks)))
                 
                 levelsInData <- levels(dat_background$breaks) %in% dat_background$breaks
                 
                 plt_background <-
                   ggplot() +
                   geom_sf(data = dat_background,
                           aes(fill = breaks), colour = cellBorderCol,
                           linewidth = cellBorderWidth) +
                   guides(fill = guide_legend(reverse = FALSE))
                 
                 if(!na.in.legend & anyNA){
                   plt_background_leg <- ggplot() +
                     geom_sf(data = subset(dat_background, !is.na(value)),
                             aes(fill = breaks), colour = cellBorderCol,
                             linewidth = cellBorderWidth) +
                     guides(fill = guide_legend(reverse = FALSE))}
                 
                 if(weblink_exists){
                   plt_background <- plt_background + 
                     scale_fill_viridis_d_interactive(
                       option = 'mako', direction = -1,
                       # drop = FALSE,
                       name = label_interactive(
                         leg_lab, data_id = background.legend.id,
                         hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                         onclick = paste0("window.open(`", as.character(
                           weblink_shipping), "`);"),
                         tooltip = weblink_shipping),
                       na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
                   if(!na.in.legend & anyNA){
                     plt_background_leg <- plt_background_leg +
                       scale_fill_viridis_d_interactive(
                         option = 'mako', direction = -1,
                         # drop = FALSE,
                         name = label_interactive(
                           leg_lab, data_id = background.legend.id,
                           hover_css = 'fill:blue;font-size:13px;font-weight:bold',
                           onclick = paste0("window.open(`", as.character(
                             weblink_shipping), "`);"),
                           tooltip = weblink_shipping),
                         na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData]))}
                   
                 }else{
                   plt_background <- plt_background +
                     scale_fill_viridis_d(option = 'mako', direction = -1,
                                          # drop = FALSE,
                                          name = leg_lab,
                                          na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))
                   if(!na.in.legend & anyNA){
                     plt_background_leg <- plt_background_leg +
                       scale_fill_viridis_d(option = 'mako', direction = -1,
                                            # drop = FALSE,
                                            name = leg_lab,
                                            na.value = na.colour, labels = c(levels(dat_background$breaks)[levelsInData], 'no data'))}
                 }
               }
             }
             
             
             # plt_background <-switch(shipSummaryDataOrRaw,
             #   summary = {},
             # 
             #   raw = {
             #     shipActivity <- unique(dat_background$activity)
             #     
             #     v <- unique(dat_background$variable)
             #     if(length(v) != 1) error("Ship traffic metric must be specified")
             #     leg_lab <- switch(v,
             #                       `ship time` = paste0('Ship traffic:', '\n', shipActivity,
             #                                            ' vessels', '\n', '(ship days year\u207B\u00B9)'),
             #                       `person time` = paste0('Ship traffic:', '\n', shipActivity,
             #                                              ' vessels', '\n', '(person days year\u207B\u00B9)'))
             #     
             #     if(!discreteColourScheme){
             #       # Continuous colour scheme
             #       plt <- 
             #         ggplot() + 
             #         geom_sf(data = dat_background,
             #                 aes(fill = value), colour = cellBorderCol)
             #       if(weblink_exists){
             #         plt <- plt + 
             #           scale_fill_viridis_c_interactive(
             #             option = 'mako', trans = 'log10', direction = -1,
             #             name = label_interactive(
             #               leg_lab, data_id = background.legend.id,
             #               hover_css = 'fill:blue;font-size:13px;font-weight:bold',
             #               onclick = paste0("window.open(`",
             #                                weblink_shipping, "`);"),
             #               tooltip = weblink_shipping),
             #             na.value = na.colour)
             #       }else{
             #         plt <- plt +
             #           scale_fill_viridis_c(option = 'mako', trans = 'log10',
             #                                direction = -1, name = leg_lab,
             #                                na.value = na.colour)
             #       }
             #     }else{
             #       # Discrete colour scheme
             #       dat_background$value_log10 <- log10(dat_background$value)
             #       # Create break points
             #       dat_background$breaks <- cut(dat_background$value_log10, nColours)
             #       # Adjust breaks for pretty legend labels
             #       brks <- gsub(']', '', gsub('\\(', '',
             #                                  levels(dat_background$breaks)))
             #       brks <- unique(as.numeric(unlist(strsplit(brks, ','))))
             #       mag <- 10 ^ floor(brks)
             #       mag_ <- 0.5 * mag
             #       # mag_ <- mag
             #       brks <- 10 ^ brks
             #       nbrks <- length(brks)
             #       i <- 2:{nbrks-1}
             #       brks[i] <- round(brks[i] / mag[i]) * mag[i]
             #       brks[nbrks] <- ceiling(brks[nbrks] / mag_[nbrks]) * mag_[nbrks]
             #       brks[1] <- round(brks[1])
             #       brks_ <- e2sci(brks)
             #       legLabs <- sapply(
             #         1:nColours, function(z) paste0('[', brks_[z], ', ', brks_[z+1], ')'))
             #       legLabs[nbrks-1] <- gsub('\\)', ']', legLabs[nbrks-1])
             #       brks <- log10(brks)
             #       dat_background$breaks <- cut(dat_background$value_log10, brks,
             #                                    labels = legLabs,
             #                                    include.lowest = TRUE)
             #       
             #       dat_background$breaks <- factor(dat_background$breaks,
             #                                       levels = rev(levels(dat_background$breaks)))
             #       
             #       plt <-
             #         ggplot() +
             #         geom_sf(data = dat_background,
             #                 aes(fill = breaks), colour = cellBorderCol,
             #                 linewidth = cellBorderWidth) +
             #         guides(fill = guide_legend(reverse = FALSE))
             #       if(weblink_exists){
             #         plt <- plt + 
             #           scale_fill_viridis_d_interactive(
             #             option = 'mako', direction = -1, drop = FALSE,
             #             name = label_interactive(
             #               leg_lab, data_id = background.legend.id,
             #               hover_css = 'fill:blue;font-size:13px;font-weight:bold',
             #               onclick = paste0("window.open(`", as.character(
             #                 weblink_shipping), "`);"),
             #               tooltip = weblink_shipping),
             #             na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
             #       }else{
             #         plt <- plt +
             #           scale_fill_viridis_d(option = 'mako', direction = -1, drop = FALSE,
             #                                name = leg_lab,
             #                                na.value = na.colour, labels = c(levels(dat_background$breaks), 'no data'))
             #       }
             #     }
             #   }
             # )
             
             
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
        geom_sf(data = cLdf_sf, aes(linetype = `p-value`),
                linewidth = contourLineWidth)#, inherit.aes = FALSE)# +
    }
  }
  
  plt_background <- plt_background +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.background = element_rect(fill = 'white', colour = 'white')
    )
  
  
  if(na.in.legend | !anyNA){
    leg_background <- get_plot_component(
      plt_background + 
        theme(legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0)),
      'guide-box', return_all = TRUE)[[1]] # this is equivalent to commented line below, but doesn't throw a warning    
    # leg_background <- get_legend(plt_background)
  }else{
    leg_background <- get_plot_component(
      plt_background_leg + 
        theme(legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0)),
      'guide-box', return_all = TRUE)[[1]]
    # leg_background <- get_legend(plt_background_leg)
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~
  # Main map -- no legend
  #~~~~~~~~~~~~~~~~~~~~~~
  
  plt_map <- plt_background + guides(fill = 'none') + guides(linetype = 'none')# + theme(legend.position = 'right')
  
  # Ecoregions
  if(displayEcoregions){
    plt_map <- plt_map +
      geom_sf(data = eco, 
              linewidth = polyLineWidth, colour = 'black',
              fill = alpha('white', 0))
  }
  
  # Coastline
  plt_map <- plt_map +
    new_scale_fill() +
    geom_sf(data = nc,
            aes(fill = surface),
            show.legend = FALSE) +
    scale_fill_manual(values = c(land.colour, ice.colour))
  # scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
  
  plt_map <- plt_map + guides(fill_new = 'none')
  
  plt_map <- plt_map + coord_sf(expand = expandBorder)
  
  include_overlay_labels <- !is.null(overlay_labels) & 
    all(class(overlay_labels) == c('sf','data.frame'))
  
  if(include_overlay_labels){
    overlay_labels_plot <- function(m, x){
      n <- nrow(x)
      r <- x$lon %% 360
      inv <- rep(FALSE, n)
      inv[90 < r & r < 270] <- TRUE
      r <- {360 - r} %% 360
      r[inv] <- {r[inv] - 180} %% 360
      Size <- x$size
      Colour <- x$colour
      hj <- x$hjust
      vj <- x$vjust
      if(is.null(x$family)) x$family <- 'serif'
      if(is.null(x$fontface)) x$family <- 1
      # fam <- ifelse(!is.null(x$family), x$family, 'serif')
      # ff <- ifelse(!is.null(x$fontface), x$fontface, 1)
      plt <- m + 
        geom_sf_text(
          data = x,
          aes(label = Name, family = family, fontface = fontface),
          size = Size, colour = Colour, angle = r, hjust = hj, vjust = vj,
          lineheight = 0.9#, 
          # family = fam, fontface = ff
        )
      plt}
  }
  
  if(backgroundOnly){
    if(showCoordGrid) plt_map <- overlay_coordinate_grid(
      plt_map, singleLatAxis = singleLatAxis, latlim = latlim, textSize = axisTextSize)
    if(include_overlay_labels) plt_map <- overlay_labels_plot(plt_map, overlay_labels)
    return(
      list(plot = ggdraw(plt_map),
           legend = ggdraw(leg_background),
           plot_complete = ggdraw(plot_grid(
             plt_map, leg_background, ncol = 2, rel_widths = c(0.8, 0.2))))
    )
  }
  
  if(anyStations | anyPlastic){
    plt_map <- plt_map +
      suppressWarnings(
      scale_shape_manual(values = setNames(symbols$symbol, symbols$Type))
    )
  }
  
  # Research stations
  if(anyStations){
    if(!stationsOnly) showStationLegend <- FALSE else showStationLegend <- TRUE
    if(!stationPopSize){
      plt_map <- plt_map +
        new_scale_colour() +
        geom_sf_interactive(data = dat_stations,
                            aes(shape = Type, colour = Seasonality,
                                data_id = Record_ID, tooltip = tooltip),
                            alpha = alpha, size = ptSize, stroke = ptStroke,
                            show.legend = showStationLegend) +
        scale_colour_manual(values = c('forestgreen','purple'),
                            guide = guide_legend(order = 2))
    }else{
      # Infill missing values using means
      i <- is.na(dat_stations$Peak_Population)
      typesMissing <- unique(dat_stations$Type[i])
      for(j in 1:length(typesMissing)){
        m <- dat_stations$Type == typesMissing[j]
        p <- dat_stations$Peak_Population[m & !i]
        infill <- round(prod(p) ^ {1/length(p)})
        dat_stations$Peak_Population[m & i] <- infill
      }
      popRanks <- c(0, 10, 25, 50, 100, Inf)
      nranks <- length(popRanks) - 1
      tt <- outer(dat_stations$Peak_Population, popRanks[1:nranks], '>=') &
        outer(dat_stations$Peak_Population, popRanks[2:{nranks+1}], '<')
      dat_stations$Population <- apply(tt, 1 , which)
      dat_stations$Population <- factor(
        dat_stations$Population, levels = 1:nranks,
        labels = c(paste('\u2265', popRanks[1:nranks])))
      plt_map <- plt_map +
        new_scale_colour() +
        geom_sf_interactive(
          data = dat_stations,
          aes(shape = Type, colour = Seasonality, size = Population,
              data_id = Record_ID, tooltip = tooltip),
          alpha = alpha, stroke = ptStroke, show.legend = showStationLegend) +
        scale_colour_manual(values = c('forestgreen','purple'))
    }
  }
  
  if(stationsOnly){
    plt_map <- plt_map +
      guides(
        shape = guide_legend(override.aes = list(
          size = legPtSize, stroke = 1),
          order = 1),
        colour = guide_legend(override.aes = list(
          size = legPtSize,
          shape = symbols$symbol[symbols$Type == 'Station']),
          order = 2))
    if(stationPopSize){
      plt_map <- plt_map + 
        guides(size = guide_legend(override.aes = list(
          shape = symbols$symbol[symbols$Type == 'Station'],
          stroke = 1),
          order = 3))
    }
    plt_map <- plt_map + theme(legend.key = element_blank())
    leg_stations <- get_legend(plt_map)
    plt_map <- plt_map + theme(legend.position = 'none')
    if(showCoordGrid) plt_map <- overlay_coordinate_grid(
      plt_map, singleLatAxis = singleLatAxis, latlim = latlim, textSize = axisTextSize)
    if(include_overlay_labels) plt_map <- overlay_labels_plot(plt_map, overlay_labels)
    return(
      list(plot = ggdraw(plt_map),
           legend = ggdraw(leg_stations),
           plot_complete = ggdraw(plot_grid(
             plt_map, leg_stations, ncol = 2, rel_widths = c(0.8, 0.2))))
    )
  }
  
  # Plastic samples
  if(anyPlastic){
    
    # To maximise visibility of samples, order the plot layers so that sources with
    # the most samples in distinct locations are plotted first and those with fewer
    # samples are layered on top. Try to also account for multiple samples collected
    # very close together.
    # In fact, it's easier and better to simply use the spread of points.
    sources <- levels(dat_plastic$Source)
    nsources <- length(sources)
    # ndistinct <- setNames(sapply(
    #   1:nsources,
    #   function(z) nrow(unique(st_coordinates(
    #     dat_plastic[dat_plastic$Source == sources[z],])))),
    #   sources)
    spread <- setNames(sapply(1:nsources, function(z){
      x <- unique(st_coordinates(
        dat_plastic[dat_plastic$Source == sources[z],]))
      prod(apply(x, 2, function(a) diff(range(a))))
    }), sources)
    # spread <- {spread - min(spread)} / diff(range(spread))
    # ord <- ndistinct * spread
    
    dat_plastic$spread <- unname(spread[as.character(dat_plastic$Source)])
    dat_plastic_ <- dat_plastic[order(dat_plastic$spread, decreasing = TRUE),]
    dat_plastic_$Source <- factor(dat_plastic_$Source,
                                  levels = unique(dat_plastic_$Source))
    
     plt_map <- plt_map +
      new_scale_fill() +
      geom_sf_interactive(data = dat_plastic_,
                          aes(fill = Source, shape = SampleType_grouped,
                              data_id = data_id, tooltip = tooltip),
                          alpha = alpha, size = ptSize, show.legend = FALSE) +
      scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
  }
  
  plt_map <- plt_map + guides(fill_new_new = 'none', colour_new = 'none')
  
  plt_map <- plt_map +
    coord_sf(xlim = c(BBox['xmin'], BBox['xmax']),
             ylim = c(BBox['ymin'], BBox['ymax']))
  
  if(showCoordGrid) plt_map <- overlay_coordinate_grid(
    plt_map, singleLatAxis = singleLatAxis, latlim = latlim, textSize = axisTextSize)
  if(include_overlay_labels) plt_map <- overlay_labels_plot(plt_map, overlay_labels)
  
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
              alpha = 1, size = ptSize, stroke = ptStroke) +
      scale_colour_manual(values = c('forestgreen','purple'))
    if(weblink_exists){
      plt_stations <- plt_stations +
        scale_shape_manual_interactive(
          values = setNames(symbols_$symbol, symbols_$Type),
          name = label_interactive(
            label = 'Facility',
            data_id = 'stations.legend.title',
            hover_css = 'fill:blue;font-size:13px;font-weight:bold',
            onclick = paste0("window.open(`", weblink_facilities, "`);"),
            tooltip = weblink_facilities)) +
        guides(
          shape = guide_legend_interactive(
            override.aes = list(size = legPtSize, stroke = 1)),
          colour = guide_legend(
            override.aes = list(size = legPtSize, shape = 1)))
    }else{
      plt_stations <- plt_stations +
        scale_shape_manual(values = setNames(symbols_$symbol, symbols_$Type)) +
        guides(
          shape = guide_legend(
            title = 'Facility', override.aes = list(size = legPtSize)),
          colour = guide_legend(
            override.aes = list(size = legPtSize, shape = 1)))
    }
    plt_stations <- plt_stations +
      theme(legend.key = element_blank())
    leg_stations <- get_plot_component(
      plt_stations + 
        theme(legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0)),
      'guide-box', return_all = TRUE)[[1]]
    # leg_stations <- get_legend(plt_stations)
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
          order = 1,
          override.aes = list(size = legPtSize)),
        fill = guide_legend_interactive(
          order = 2,
          override.aes = list(
            shape = 21,
            size = legPtSize_fill))
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
    leg_plastic <- get_plot_component(plt_plastic_samples,# + theme(legend.margin = margin(0,0,0,0)),
                                      'guide-box', return_all = TRUE)[[1]]
    # leg_plastic <- get_legend(plt_plastic_samples)
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
      list(plot = plt_map,
           legend = ggdraw(leg_plastic),
           plot_complete = ggdraw(plot_grid(
             plt_map, leg_plastic, ncol = 2, rel_widths = c(0.8, 0.2)))
      )
    )
  }
  
  # Find size (cm) of combined legend -- include spacings
  leg_width <- {1 * 1} + 2 * max(c(leg_plastic_width, leg_stations_width,
                                   leg_background_width))
  leg_height <- {3 * 2} + max(c(leg_plastic_height, leg_stations_height + leg_background_height))
  # Size of map & complete plot
  tot_width = leg_width / legWidth
  map_width = tot_width - leg_width
  map_height <- map_width / mapAspectRatio
  tot_height <- max(map_height, leg_height)
  
  leg_blank <- patchwork::plot_spacer() + theme_void()# + theme(plot.margin = unit(c(0,0,0,0),'cm'))
  
  l1 <- plot_grid(
    leg_blank, leg_plastic, leg_blank,
    nrow = 3, rel_heights = c(
    0.5*{tot_height-leg_plastic_height}, leg_plastic_height, 0.5*{tot_height-leg_plastic_height}))
  l2 <- plot_grid(
    leg_blank, leg_stations, leg_blank, leg_background, leg_blank,
    nrow = 5, rel_heights = c(
      0.45*{tot_height - {leg_stations_height + leg_background_height}},
      leg_stations_height,
      0.1*{tot_height - {leg_stations_height + leg_background_height}},
      leg_background_height,
      0.45*{tot_height - {leg_stations_height + leg_background_height}}
      ), align = 'v', axis = 'l')
  leg_complete <- plot_grid(l1, l2, ncol=2)
  
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


# Risk map functions ------------------------------------------------------

# geometric mean of vector x
geomean <- function(x, na.rm = FALSE){
  if(na.rm) x <- x[!is.na(x)]
  return(exp(sum(log(x)) / length(x)))
}


# function to plot risk levels
plot_fun <- function(data, n = nranks, v_option = 'viridis', v_direction = -1,
                     legend_title = 'risk', Title = NULL, show_legend = TRUE,
                     showAxis = TRUE, latlim = NULL, axisTextSize = 4,
                     cellBorderColour = 'grey', cellBorderWidth = 0.1,
                     na.colour = 'grey40', background.na.remove = FALSE,
                     ice.colour = 'skyblue', land.colour = 'grey85'){
  data$rank <- factor(data$rank, levels = n:1, labels = n:1)
  # get the legend
  plt <-
    ggplot() +
    geom_sf(data = data, aes(fill = rank), colour = cellBorderColour,
            linewidth = cellBorderWidth, show.legend = TRUE) +
    scale_fill_viridis_d(option = v_option, direction = v_direction,
                         drop = FALSE, na.translate = FALSE, name = legend_title)
  leg <- ggdraw(get_legend(plt))
  # make full plot
  plt <-
    ggplot() +
    geom_sf(data = data, aes(fill = rank), colour = cellBorderColour,
            linewidth = cellBorderWidth, show.legend = show_legend) +
    scale_fill_viridis_d(option = v_option, direction = v_direction,
                         drop = FALSE, na.translate = TRUE, name = legend_title,
                         na.value = na.colour, labels = c(levels(data$rank), 'no data')) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.background = element_rect(fill = 'white', colour = 'white')
    ) +
    new_scale('fill') +
    geom_sf(data = nc,
            aes(fill = surface)) +
    scale_fill_manual(values = c(land.colour,ice.colour,ice.colour,land.colour), guide = 'none')
  
  if(showAxis){
    assign('BBox', st_bbox(nc), envir = .GlobalEnv)
    plt <- overlay_coordinate_grid(plt, textSize = axisTextSize, latlim = latlim)
  }
  
  if(!is.null(Title)) plt <- plt + labs(title = Title)
  return(list(plot = plt, legend = leg))
}

# function to plot stations
plot_fun_stations <- function(data_stations, data_grid_cells, rank_type = 'additive', n = nranks,
                              v_option = 'viridis', v_direction = -1, alpha = 0.75,
                              stroke = 1, rank_limits_PopSize = NULL, Title = NULL, show_legend = TRUE,
                              showAxis = TRUE, latlim = NULL, axisTextSize = 4){
  
  data_stations$rank <- factor(data_stations$rank, levels = n:1, labels = n:1)
  
  if(is.null(rank_limits_PopSize)){
    data_stations$rank_pop_labels <- data_stations$rank_pop
  }else{
    rl <- paste('\u2265', rank_limits_PopSize[1:n])
    data_stations$rank_pop_labels <- factor(data_stations$rank_pop, levels = 1:n, labels = rl)
  }
  
  seasons <- tolower(levels(data_stations$Seasonality))
  data_stations$Seasonality <- factor(data_stations$Seasonality, labels = seasons)
  
  data_grid_cells$rank <- factor(data_grid_cells$rank, levels = n:1, labels = n:1)
  
  plt <-
    ggplot(data = data_stations) + 
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.background = element_rect(fill = 'white', colour = 'white'),
      legend.key = element_blank(),
      legend.background = element_blank()
    ) +
    geom_sf(data = data_grid_cells,
            aes(fill = rank)) +
    scale_fill_viridis_d(option = v_option, direction = v_direction, name = 'activity',
                         drop = FALSE, na.translate = FALSE, guide = guide_legend(order = 3)) +
    new_scale('fill') +
    geom_sf(data = nc,
            aes(fill = surface)) +
    scale_fill_manual(values = c('grey','skyblue','skyblue','grey'), guide = 'none') +
    new_scale('fill') +
    geom_sf(aes(shape = Seasonality, size = rank_pop_labels, fill = rank),
            alpha = alpha, stroke = stroke) +
    scale_fill_viridis_d(option = v_option, direction = v_direction, guide = 'none') +
    scale_fill_viridis_d(aesthetics = 'fill_new_new', option = v_option, direction = v_direction, name = 'activity', drop = FALSE, na.translate = FALSE, guide = guide_legend(override.aes = list(size = NA, shape = NA))) +
    scale_size_manual(values = 1:n, guide = guide_legend(title = 'population', order = 1, override.aes = list(fill = NA, linetype = 0))) +
    scale_shape_manual(values = c(21, 22), guide = guide_legend(title = 'seasonality', order = 2, override.aes = list(size = n, fill = NA, linetype = 0)))
  
  if(showAxis){
    assign('BBox', st_bbox(nc), envir = .GlobalEnv)
    plt <- overlay_coordinate_grid(plt, latlim = latlim, textSize = axisTextSize)
  }
  
  if(!is.null(Title)) plt <- plt + labs(title = Title)
  leg <- suppressWarnings(ggdraw(get_legend(plt)))
  if(!show_legend) plt <- plt + theme(legend.position = 'none')
  output <- list(plot = plt, legend = leg)
  return(output)
}


# # The 'subset' function doesn't work with sf objects. Use this simple function to subset sf data
# subset_sf <- function(dat, variable, value){
#   x <- dat[dat[[variable]] %in% value,]
#   if(is.factor(x[[variable]])) x[[variable]] <- as.factor(as.character(x[[variable]]))
#   return(x)}

# Use a continuous 'score' to assess risk. This behaves better than integer ranks
# under averaging and may be transformed into an integer 'rank' by taking the ceiling
# (or floor).
score2rank <- function(x, base = 0){
  y <- ceiling(x)
  if(is.numeric(base)) y[y == base] <- base + 1
  y
}

# functions for rounding -- used in automated process for selecting ranking limit values
zerosAfterDecimal <- function(x){
  if(x == 0) out <- 0 else{
    options(scipen = 999)
    y <- as.character(abs(x))
    y <- strsplit(y, '\\.')
    if(as.numeric(y[[1]][1]) > 0) out <- 0 else{
      y <- y[[1]][2]
      y <- as.numeric(strsplit(y, '')[[1]]) == 0
      out <- which(!y)[1] - 1
    }
    options(scipen = 0)
  }
  return(out)
}
orderOfMagnitude <- Vectorize(function(x){
  if(x == 0) out <- 1 else{
    z <- zerosAfterDecimal(x)
    if(z != 0) out <- 10 ^ -{z+1} else{
      y <- strsplit(as.character(abs(x)), '\\.')[[1]][1]
      if(y == '0') out <- 0.1 else{
        out <- 10 ^ {nchar(y) - 1}
      }}}
  return(out)
})
# orderOfMagnitude <- Vectorize(orderOfMagnitude)

