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

# getwd()
setwd('~/Documents/Git Repos/CUPIDO-risk-map')

# Load map shape file -----------------------------------------------------
verbose <- FALSE
nc <- st_read('data/map/add_coastline_medium_res_polygon_v7.3.shp', quiet = !verbose)
nc$surface = factor(nc$surface, levels = unique(nc$surface))
# Coordinate reference system
crs_world = 4326
crs_use = 3031



# Load plastic data -----------------------------------------------------
filename = 'plastic_quantity.csv'
filepath = 'data/plastic_quantity/'
DATA = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing ----
# Date format
DATA$Date_1 = as.Date(DATA$Date_1, '%d-%b-%Y')
DATA$Date_2 = as.Date(DATA$Date_2, '%d-%b-%Y')
# Include Year column
DATA$Year = as.integer(substr(DATA$Date_1, 1, 4))

# Create sample type data column
DATA$SampleType[DATA$Depth %in% c('5m','<1m','surface')] = 'surface'
DATA$SampleType[DATA$Depth %in% c('subsurface')] = 'subsurface'
DATA$SampleType = factor(DATA$SampleType, levels = c('surface', 'subsurface'))

# Reorder data sources by publication date
sources <- levels(DATA$Source)
nsources <- length(sources)
ncs <- nchar(sources)
sourceYear <- suppressWarnings(
  as.integer(sapply(1:nsources, FUN = function(z) substr(sources[z], ncs[z]-4, ncs[z]-1)))
)
DATA$Source <- factor(DATA$Source, levels = sources[order(sourceYear)])


# Mapping coordinates
DATA_sf = st_as_sf(DATA, coords = c("Longitude", "Latitude"), crs = crs_world)
DATA_sf = st_transform(DATA_sf, crs_use)
xy = matrix(unlist(DATA_sf$geometry), 2, nrow(DATA))
DATA$x = xy[1,]
DATA$y = xy[2,]

# # Set (minimum) size of plot bounding box
# bbox_map = st_bbox(nc)
# bbox_dat = st_bbox(DATA_sf)
# bbox = setNames(c(min(bbox_map$xmin, bbox_dat$xmin),
#                   min(bbox_map$ymin, bbox_dat$ymin),
#                   max(bbox_map$xmax, bbox_dat$xmax),
#                   max(bbox_map$ymax, bbox_dat$ymax)), names(bbox_map))
# aspectRatio = diff(bbox[c(1,3)]) / diff(bbox[c(2,4)])


# Load research station data ----------------------------------------------
filename = 'COMNAP_Antarctic_Facilities_Master.csv'
filepath = 'data/research stations/COMNAP/'
STATIONS = read.csv(paste0(filepath, filename), stringsAsFactors = TRUE)

# Data pre-processing ----
STATIONS$Type = reorder(STATIONS$Type, rep(1, nrow(STATIONS)), sum, decreasing = TRUE)

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
pltShapes <- c(3,4,8,7,9,10,12,13,14,11)
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
bbox <- setNames(c(min(bbox_map$xmin, bbox_dat$xmin, bbox_krill$xmin),
                   min(bbox_map$ymin, bbox_dat$ymin, bbox_krill$ymin),
                   max(bbox_map$xmax, bbox_dat$xmax, bbox_krill$xmax),
                   max(bbox_map$ymax, bbox_dat$ymax, bbox_krill$ymax)), names(bbox_map))
aspectRatio = diff(bbox[c(1,3)]) / diff(bbox[c(2,4)])




# Define UI for plastic data map app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Mapping Southern Ocean Plastic Data"),
  
  # fluidRow(
  #   column(12,
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: year range -------------------------------------------------------
      sliderInput("YearRange", "Years:",
                  min = min(DATA$Year, na.rm = TRUE) , max = max(DATA$Year, na.rm = TRUE),
                  value = range(DATA$Year, na.rm = TRUE), step = 1, sep = ''),
      
      # Input: measurement variable ----
      selectInput("Variable", "Measurement:",
                  c("Concentration (pieces/m3)" = "concentration",
                    "Density (pieces/km2)" = "density",
                    "Mass density (g/km2)" = "massDensity"),
                  multiple = TRUE, selected = c("concentration", "density", "massDensity")),
      
      # Input: plastic type ----
      # I NEED TO THINK ABOUT HOW TO HANDLE THE TOTAL COLUMN... OMIT IT FOR NOW
      # BUT MAYBE INCLUDE IT LATER WHEN GRAPHING DATA...
      selectInput("Type", "Plastic type:",
                  c("Fragment" = "fragment",
                    "Fibre" = "fibre",
                    "Film" = "film"),
                  multiple = TRUE, selected = c("fragment", "fibre", "film")),
      
      # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment) ----
      selectInput("SampleType", "Sample type:",
                  c("Near-surface" = "surface",
                    "Subsurface" = "subsurface"),
                  multiple = TRUE, selected = c("surface")),
      
      # Input: sample type (depth, for now, later extended to sub/surface-beach-sediment) ----
      selectInput("StationType", "Facility:",
                  c("Station" = "Station",
                    "Camp" = "Camp",
                    "Refuge" = "Refuge",
                    "Airfield" = "Airfield Camp",
                    "Laboratory" = "Laboratory",
                    "Depot" = "Depot"),
                  multiple = TRUE, selected = c("Station")),

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
                   selected = 'none'),

      
      # # Input: plot value transformation ----
      # radioButtons("tran", "Scale:",
      #              c("Natural" = "norm",
      #                "Log" = "log")),
      
      # # Input: Checkbox for whether outliers should be included ----
      # checkboxInput("outliers", "Show outliers", TRUE),
      
      width = 3
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Plot of the requested variable against mpg ----
      # h3("Data plot"),
      
      plotOutput('blank', width = '100%', height = '1px'),
      plotOutput('plt'), #, inline = T)
      
      width = 9
      
      # br(),
      # 
      # # Output: Verbatim text for data summary ----
      # h3("Data summary"),
      # verbatimTextOutput("summary")
      # 
      # br(),
      # 
      # # Output: HTML table with requested number of observations ----
      # 
      #   h3("Raw data"),
      #   tableOutput("view")
        

    )
    
  )
  #   )
  # )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output, session) {
  
  # Function returning filtered data specified by Shiny inputs
  filter_plastic_data = reactive({
    return(
      subset(DATA,
             input$YearRange[1] <= Year & Year <= input$YearRange[2] &
               DATA$Variable %in% input$Variable &
               DATA$Type %in% input$Type &
               DATA$SampleType %in% input$SampleType
      )
    )
  })
  
  # Subset data required for map overlay and transform to correct CRS
  filtered_plastic_data = reactive({
    # Filter
    d = unique(subset(filter_plastic_data(), select = c(Source, Station, SampleType, Longitude, Latitude)))
    # Transform
    d = st_as_sf(d, coords = c("Longitude", "Latitude"), crs = crs_world)
    d = st_transform(d, crs_use)
    return(d)
  })
  
  # Filter the research station data
  filtered_station_data = reactive({
    return(
      subset(STATIONS_sf,
               STATIONS_sf$Type %in% input$StationType
      )
    )
  })
  
  # Filter background data
  filtered_background_dat <- reactive({
    backgroundOptions <- c('none', 'krill', 'chl')
    t <- sapply(backgroundOptions, FUN = function(z) grepl(z, input$background))
    background <- names(which(t)) # get the chosen background
    background_dat <- switch(background,
                             krill = krill_poly,
                             chl = chl_poly
    ) # get background data
    if(background != 'none'){
      m <- strsplit(input$background, '_')[[1]][2] # get the chosen month
      background_dat <- subset(background_dat, month == m) # filter by month
    }
    return(list(background, background_dat))
  })
  
  # Filter plot symbols based on sample & station types
  Symbols <- reactive({
    return(
      subset(pltSymbols, Type %in% c(input$SampleType, input$StationType))
    )
  })
  
  # Create blank plot to specify sizing/aspect ratio
  output$blank = renderPlot({
   ggplot(data.frame(x = 1), aes(NA, NA)) + geom_blank() + theme_void()
  })
  
  blankwidth = reactive({
    # this is the magic that makes it work
    bw = session$clientData$output_blank_width
    return(0.85 * bw)
#    if (bw > 400) 400 else bw
  })
  
  blankheight <- reactive({
    blankwidth() / aspectRatio
  })
  
  
  # Build up the plot layer by layer, extracting each legend separately, then
  # piece it together.
  
  output$plt <- renderPlot({
    
    dat_plastic <- filtered_plastic_data()
    dat_stations <- filtered_station_data()
    background_info <- filtered_background_dat()
    background <- background_info[[1]]
    dat_background <- background_info[[2]]
    symbols <- Symbols()
    
    # dat_plastic <- DATA_sf
    # dat_stations <- STATIONS_sf
    # dat_krill <- krill_poly
    # dat_chl <- chl_poly
    # background <- 'chl'
    # symbols <- pltSymbols
    
    anyPlastic <- nrow(dat_plastic) > 0
    anyStations <- nrow(dat_stations) > 0
    
    # Background layer
    anyBackground <- input$background != 'none'
    switch(background,
           krill = {
             plt_background <- 
               ggplot() +
               geom_sf(data = dat_background, aes(fill = colourgroup)) +
               scale_fill_viridis_d(option = 'plasma',
                                    name = bquote(atop(Krill, individuals / m^2)))
               # guides(fill = guide_legend(title = bquote(atop(Krill, individuals / m^2))))
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
                alpha = 1, size = 4, stroke = 1) +
        scale_colour_manual(values = c('forestgreen','firebrick')) +
        scale_shape_manual(values = symbols$symbol[symbols$Class == 'ResearchStation']) +
        # scale_shape_manual(values = Symbols()$symbol[Symbols()$Class == 'ResearchStation']) +
        # scale_shape_manual(values = c(3,4,8,7,9,10)) +
        # guides(shape = guide_legend(title = 'Facility', override.aes = list(shape = unique(dat_stations$symbol))),
        #        colour = guide_legend(override.aes = list(colour = c('forestgreen','firebrick'))))
        guides(shape = guide_legend(title = 'Facility'), 
               colour = guide_legend(override.aes = list(colour = c('forestgreen','firebrick'), shape = 3)))
    }
    
    # Plastic samples
    if(anyPlastic){
      plt_plastic_samples <- 
        ggplot() +
        geom_sf(data = dat_plastic,
                aes(fill = Source, shape = SampleType),
                alpha = 1, size = 4) +
        scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source)) +
        #      scale_colour_brewer(type = 'qual', palette = 'Set1', aesthetics = c('fill')) +
        scale_shape_manual(values = symbols$symbol[symbols$Class == 'PlasticSample']) +
        guides(shape = guide_legend(title = 'Sample type'),
               fill = guide_legend(override.aes = list(shape = c(21))))
      # guides(fill = guide_legend(override.aes = list(shape = c(21))))
    }
    
    # Extract legends
    if(anyStations){
      leg_stations <- gtable_filter(ggplot_gtable(ggplot_build(plt_stations)), "guide-box") 
      leg_stations_grob <- grobTree(leg_stations)
    }
    
    if(anyPlastic){
      leg_plastic_samples <- gtable_filter(ggplot_gtable(ggplot_build(plt_plastic_samples)), "guide-box") 
      leg_plastic_samples_grob <- grobTree(leg_plastic_samples)
    }
    
    # Combine plot layers
    plt_no_legend <- 
      
      # Coastline & (optional) background
      plt_map
      
      # ggplot() + 
      
      # # Coastline
      # geom_sf(data = nc,
      #         aes(fill = surface),
      #         show.legend = FALSE) +
      # scale_fill_manual(values = c('grey','skyblue','skyblue','grey')) +
      
    if(anyStations | anyPlastic){
      plt_no_legend <- plt_no_legend +
        scale_shape_manual(values = setNames(symbols$symbol, symbols$Type))
    }
      
    # Research stations
    if(anyStations){
      plt_no_legend <- plt_no_legend +
      geom_sf(data = dat_stations,
              aes(shape = Type, colour = Seasonality),
              alpha = 1, size = 4, stroke = 1, show.legend = FALSE) +
        scale_colour_manual(values = c('forestgreen','firebrick'))
        # scale_shape_manual(values = symbols$symbol) +
        # scale_shape_manual(values = symbols$symbol[symbols$Class == 'ResearchStation']) +
        # scale_shape_manual(values = Symbols()$symbol) + # IT'S CRUCIAL TO ENSURE THAT THE FACTOR ORDER IS CONSISTENT BETWEEN DATA SETS AND THE SYMBOLS DATA SET
        #      scale_shape_manual(values = c(3,4,8,7,9,10,21:25)) + # 1st 6 values for stations, remaining values for plastic samples
    }
      
    # Plastic samples
    if(anyPlastic){
      plt_no_legend <- plt_no_legend +
      new_scale('fill') +
        geom_sf(data = dat_plastic,
                aes(fill = Source, shape = SampleType),
                alpha = 1, size = 4, show.legend = FALSE) +
        scale_fill_manual(values = setNames(pltColours$colour, pltColours$Source))
      # scale_colour_brewer(type = 'qual', palette = 'Set1', aesthetics = c('fill')) +
      # scale_shape_manual(values = symbols$symbol[symbols$Class == 'PlasticSample']) +
    }
      
    plt_no_legend <- plt_no_legend +
      coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), ylim = c(bbox['ymin'], bbox['ymax'])) +
      theme_void()
    
    # Combine legends
    legs <- 
      ggplot(data = data.frame(x = 0:1, y = 0:1)) +
      geom_blank(aes(x = x, y = y)) + 
      theme_minimal() + 
      ylab(NULL) + xlab(NULL) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() 
      )
    if(anyPlastic){
      legs <- legs +
        annotation_custom(leg_plastic_samples_grob, xmin = 0.5, xmax = 0.5, ymin = 0.5, ymax = 1)
    }
    if(anyStations){
      legs <- legs +
        annotation_custom(leg_stations_grob, xmin = 0.35, xmax = 0.35, ymin = 0, ymax = 0.5) 
    }

    # Use cowplot package to arrange plot
    # plt <- ggdraw(plt_no_legend, clip = 'on')# + 
    plt <- ggdraw() + 
      draw_plot(plt_no_legend, 0, 0, 0.7, 1)
    if(anyPlastic){
      plt <- plt +
        draw_plot(leg_plastic_samples_grob, 0.7, 0, 0.15, 1)
    }
    if(anyStations){
      plt <- plt +
        draw_plot(leg_stations_grob, 0.85, 0, 0.15, 1)
    }
   plt
    
    
    # # Use gridExtra package to arrange plot 
    # # Output plot
    # plt <- arrangeGrob(plt_no_legend, legs, ncol=2, widths=c(0.8, 0.2))
    # 
    # grid.newpage()
    # grid.draw(plt)    
    
  },
height = blankheight, width = blankwidth)
  


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
