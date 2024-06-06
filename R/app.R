#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Shiny app for interactive map of Southern Ocean displaying plastic & environmental data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This interactive map, version 1.0, was created using R version 4.4

# Packages ----------------------------------------------------------------
# Required packages
pkg.list <- c('remotes','shiny','sp','sf','mapdata',
              # 'maptools', # deprecated package
              'ggplot2','gtable','grid',
              'gridExtra','cowplot','RColorBrewer','scales','ggiraph','reshape2',
              'flextable','DT','patchwork') #,'devtools')

# It seems that the code for generating multi-layered plots with multiple colour/fill
# scales is not robust to updated package versions, so using specific package is
# required. This is an issue to address in future...

# Package versions used for interactive map version 1.0:
pkg.version <- c(remotes = '2.5.0', shiny = '1.7.4', sp = '2.1.4', sf = '1.0.16', mapdata = '2.3.1',
                 # maptools = '1.1.6',
                 ggplot2 = '3.4.4', gtable = '0.3.3', grid = '4.4.0',
                 gridExtra = '2.3', cowplot = '1.1.1', RColorBrewer = '1.1.3', scales = '1.2.1',
                 ggiraph = '0.8.7', reshape2 = '1.4.4', flextable = '0.9.1',
                 DT = '0.27', patchwork = '1.2.0') #, devtools = '2.4.5')

# Notation for archived versions on CRAN varies slightly -- using either a period or dash for minor versions
pkg.version.archive <- c(remotes = '2.5.0', shiny = '1.7.4', sp = '2.1-4', sf = '1.0-16', mapdata = '2.3-1',
                         # maptools = '1.1.6',
                         ggplot2 = '3.4.4', gtable = '0.3.3', grid = '4.4-0',
                         gridExtra = '2.3', cowplot = '1.1.1', RColorBrewer = '1.1-3', scales = '1.2.1',
                         ggiraph = '0.8.7', reshape2 = '1.4.4', flextable = '0.9.1',
                         DT = '0.27', patchwork = '1.2.0') #, devtools = '2.4.5')

# Any problems arising after a fresh install may be due incompatibilities between
# newer package versions, so if there are problems then try running with the package
# versions listed here.
install.listed.pkg.version <- FALSE # if any package is missing, install listed version (TRUE) or latest version (FALSE)?

# Load packages, installing if necessary.
for(i in 1:length(pkg.list)){
  pkg <- pkg.list[i]
  j <- library(pkg, character.only = TRUE, logical.return = TRUE)
  if(!j){
    if(install.listed.pkg.version & pkg != 'remotes'){
      # Problems may emerge with this method if package dependencies are not
      # installed. There may be an automatic way to install dependencies while
      # repos = NULL, otherwise dependencies may need installed individually...
      # method 1
      # pkg_ <- paste0('https://cran.r-project.org/package=', pkg, '&version=', pkg.version[pkg])
      # install.packages(pkgs = pkg_, repos = NULL)
      # method 2
      # pkg_ <- paste0('http://cran.r-project.org/src/contrib/Archive/', pkg, '/', pkg, '_', pkg.version.archive[pkg], '.tar.gz')
      # install.packages(pkg_, repos = NULL, type = 'source')
      # method 3
      tryCatch(remotes::install_version(
        pkg, version = pkg.version[i], build = TRUE, upgrade = FALSE),
        error = function(e) NA, warning = function(w) NA)
    }else{
      install.packages(pkg, repos = 'https://cloud.r-project.org')
    }
    library(pkg, character.only = TRUE)
  }else{
    if(install.listed.pkg.version){
      desired.version <- pkg.version[i] == packageVersion(pkg)
      if(!desired.version){
        # method 1
        # pkg_ <- paste0('https://cran.r-project.org/package=', pkg, '&version=', pkg.version[pkg])
        # install.packages(pkgs = pkg_, repos = NULL)
        # method 2
        # pkg_ <- paste0('http://cran.r-project.org/src/contrib/Archive/', pkg, '/', pkg, '_', pkg.version.archive[pkg], '.tar.gz')
        # install.packages(pkg_, repos = NULL, type = 'source')
        # method 3
        tryCatch(remotes::install_version(
          pkg, version = pkg.version[i], build = TRUE, upgrade = FALSE),
                 error = function(e) NA, warning = function(w) NA)
      }
    }
  }
}

# We need another package, ggnewscale, that, at the time of writing, needs to be
# installed from GitHub. It also needs to be an old version as there appeared to
# be a bug in the latest version.
# The package developer proposed a solution to the bug here: https://github.com/eliocamp/ggnewscale/issues/45
# but I'm not convinced this helped! The issue seemed to solved by including '_new'
# as a suffix in arguments to the 'guides' function.

# ggnewscale_version <- '0.4.3'
ggnewscale_version <- '0.4.10'
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

dir_base <- dirname(getwd())
dir_data <- paste(dir_base, 'data', sep = '/')
dir_map <- paste(dir_data, 'map', sep = '/')


# Load data ---------------------------------------------------------------

source('functions.R', local = TRUE)

significantTrendsOnly <- TRUE
loadTooltipFromFile <- TRUE # set to FALSE when new plastic data is included or any modifications are made to the data
displayAllLitterTypes <- TRUE # If FALSE then only plastics are displayed. Set to TRUE to show samples of non-plastics -- this is important when trawls bring up items that appear to be plastic (e.g., Kuklinski data)
displayBackgroundNAs <- TRUE

roundShipTime <- FALSE

if(!displayBackgroundNAs){
  get_data(
    baseDirectory = dir_base, dataDirectory = dir_data, mapDirectory = dir_map,
    res = '3x1', 
    allLitterTypes = displayAllLitterTypes,
    sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
    sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly,
    roundShipTime = roundShipTime, shipOrPersonTime = NULL, indexGridCells = FALSE,
    loadTooltipFromFile = loadTooltipFromFile, theseData = NULL)
}else{
  # Include NA values for empty grid cells in all background data sets. Use the
  # sea surface temperature data to define the grid cells because this data set
  # has values for every cell.
  get_data(
    baseDirectory = dir_base, dataDirectory = dir_data, mapDirectory = dir_map,
    res = '9x3', 
    allLitterTypes = displayAllLitterTypes,
    sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
    sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly,
    roundShipTime = roundShipTime, shipOrPersonTime = NULL, indexGridCells = FALSE,
    loadTooltipFromFile = loadTooltipFromFile,
    theseData = c('map', 'plastic', 'facility', 'ecoregions', 'krill','ship','sst'))
  
  sst_poly_9x3 <- sst_poly
  
  get_data(
    baseDirectory = dir_base, dataDirectory = dir_data, mapDirectory = dir_map,
    res = '3x1', 
    allLitterTypes = displayAllLitterTypes,
    sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
    sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly,
    roundShipTime = roundShipTime, shipOrPersonTime = NULL, indexGridCells = FALSE,
    loadTooltipFromFile = loadTooltipFromFile,
    theseData = c('chl','sst','pH'))
  
  all_cells_9x3 <- unique(sst_poly_9x3[,'geometry'])
  all_cells_3x1 <- unique(sst_poly[,'geometry'])
  rm(sst_poly_9x3)
  # Chlorophyll
  for(l in unique(chl_poly$month)){
    d <- chl_poly[chl_poly$month == l,]
    i <- !all_cells_3x1$geometry %in% d$geometry
    if(sum(i) == 0) next
    j <- all_cells_3x1[i,]
    j$month <- l
    j$value <- NA
    chl_poly <- rbind(chl_poly, j)
    chl_poly <- chl_poly[order(chl_poly$month),]
  }
  # Krill
  for(l in unique(krill_poly$month)){
    d <- krill_poly[krill_poly$month == l,]
    i <- !all_cells_9x3$geometry %in% d$geometry
    if(sum(i) == 0) next
    j <- all_cells_9x3[i,]
    j$month <- l
    j$value <- NA
    j$colourgroup <- NA
    krill_poly <- rbind(krill_poly, j)
    krill_poly <- krill_poly[order(krill_poly$month),]
  }
  # Sea surface temperature
  for(l in unique(sst_poly$metric)){
    d <- sst_poly[sst_poly$metric == l,]
    i <- !all_cells_3x1$geometry %in% d$geometry
    if(sum(i) == 0) next
    j <- all_cells_3x1[i,]
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
    i <- !all_cells_3x1$geometry %in% d$geometry
    if(sum(i) == 0) next
    j <- all_cells_3x1[i,]
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
      i <- !all_cells_9x3$geometry %in% d$geometry
      if(sum(i) == 0) next
      j <- all_cells_9x3[i,]
      j$year <- 'all'
      j$activity <- l
      j$variable <- m
      j$value <- NA
      ship_poly <- rbind(ship_poly, j)
      ship_poly <- ship_poly[order(ship_poly$activity, ship_poly$variable),]
    }
  }
}

nc <- nc_plastic


# It will be nice to able to select between trend data and anomalies from within
# the Shiny app. This will involve returning both data sets from the get_data
# function, then including a selection option in the Shiny dashboard that chooses
# between them. This is not critical to include, but would be a nice extension
# to what I've already produced

# Flextable options ------------------------------------------------------

set_flextable_defaults(
  font.size = 10,
  font.family = 'Roboto',
  theme_fun = 'theme_vanilla'
)

# Plotting parameters -----------------------------------------------------

set_plot_params(nc, DATA_sf, STATIONS_sf)

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
  'Chlorophyll (Mar)'
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
  'chl_3'
)

if(SST_overallTrend){
  backgroundData_choiceNames <- c(backgroundData_choiceNames, 'SST trend')
  backgroundData_choiceValues <- c(backgroundData_choiceValues, 'sst')
}else{
  backgroundData_choiceNames <- c(backgroundData_choiceNames, c(paste('SST', sstType, '(Jan--Mar)', sep = ' '),
                                                               paste('SST', sstType, '(Jan)', sep = ' '),
                                                               paste('SST', sstType, '(Feb)', sep = ' '),
                                                               paste('SST', sstType, '(Mar)', sep = ' ')))
  backgroundData_choiceValues <- c(backgroundData_choiceValues, c('sst_all',
                                                                  'sst_1',
                                                                  'sst_2',
                                                                  'sst_3'))
}

if(pH_overallTrend){
  backgroundData_choiceNames <- c(backgroundData_choiceNames, 'pH trend')
  backgroundData_choiceValues <- c(backgroundData_choiceValues, 'pH')
}else{
  backgroundData_choiceNames <- c(backgroundData_choiceNames, c(paste('pH', pHType, '(Jan--Mar)', sep = ' '),
                                                                paste('pH', pHType, '(Jan)', sep = ' '),
                                                                paste('pH', pHType, '(Feb)', sep = ' '),
                                                                paste('pH', pHType, '(Mar)', sep = ' ')))
  backgroundData_choiceValues <- c(backgroundData_choiceValues, c('pH_all',
                                                                  'pH_1',
                                                                  'pH_2',
                                                                  'pH_3'))
}

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
                           'Downward flux (pieces/m2/day)' = 'flux',
                           'Presence/absence' = 'presence/absence'),
                         multiple = TRUE, selected = c('concentration', 'mass concentration', 'density', 'mass density', 'flux', 'presence/absence')),
             
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
             
             # Input: shipping metric (ship or person time)
             radioButtons('ShipMetric', 'Ship traffic metric:',
                          choiceNames = c('ship time', 'person time'),
                          choiceValues = c('ship time', 'person time')),

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
           plotOutput('blank', width = paste0(9/12*100,'%'), height = '10px'),
           # plotOutput('blank', width = '100%', height = '10px'),
           girafeOutput('plt')
    ),
    
  ),
  
  fluidRow(
    column(width = 12, #offset = 3,
           
           # wellPanel(
           h4('Selected plastic data'),
           DT::DTOutput('datatab_plastic')
           # shiny::dataTableOutput('datatab_plastic')
           # )
           
    )
  ),
  
  fluidRow(
    column(width = 12,
           
           # wellPanel(
           h4('Selected facilities'),
           DT::DTOutput('datatab_stations')
           # shiny::dataTableOutput('datatab_stations')
           # )
           
    )
  )
)

# Define server logic ----
server <- function(input, output, session){
  
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
        StationType = input$StationType,
        ShipMetric = input$ShipMetric
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
  filtered_plastic_data <- reactive({
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
        }# else warning('Background data cannot be filtered by month: check that R Shiny ui input options match the data.')
      }else{
        if(background == 'ship'){
          y <- listInputs()
          background_dat <- switch(shipSummaryDataOrRaw,
                                   summary = {
                                     x <- strsplit(input$background, '_')[[1]]
                                     if(x[2] %in% background_dat$ship_class){
                                       background_dat <- subset(background_dat, ship_class == x[2])
                                     } else{warning('Background data cannot be filtered by ship type: check that R Shiny input options match the data.')}
                                   },
                                   raw = {
                                     z <- subset(background_dat, variable == y$ShipMetric)
                                     x <- strsplit(input$background, '_')[[1]]
                                     if(x[2] %in% background_dat$activity){
                                       background_dat <- subset(z, activity == x[2])
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
    ggplot(data.frame()) + geom_point() + xlim(0,1) + ylim(0,1) + theme_void()
    # ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point() + theme_void()
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
                     plotSignificanceContours = display_trend_pvals(), latlim = lat_lim_plastic)
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
                        # opts_tooltip(css = 'background:white;'),
                        opts_tooltip(css = 'background:white;', offy = 30, opacity = 0.9, placement = 'auto'),
                        opts_zoom(min = 0.5, max = 10),
                        opts_hover(css = 'opacity:1.0;stroke-width:4;cursor:pointer;', reactive = TRUE),
                        opts_hover_inv(css = 'opacity:0.2;cursor:pointer;'),
                        opts_selection(type = 'multiple', css = 'opacity:1.0;stroke-width:4;'),
                        opts_selection_key(css = girafe_css("stroke:none;",
                                                            text = "stroke:none;fill:black;"))#,
                        # opts_selection_key(css = girafe_css("stroke:red; stroke-width:2px",
                        #                                     text = "stroke:none;fill:red;font-size:12px"))#,
                        # opts_hover(css = 'opacity:1.0;stroke-width:4;cursor:pointer;', reactive = TRUE)
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

  displayed_data_table_plastic <- reactive({
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
    # d$URL <- paste0("<a href='", d$URL, "'>", d$URL,"</a>")
    d$URL <- paste0("<a href='", d$URL, "' target=_blank", ">", d$URL,"</a>")
    row.names(d) <- NULL
    d <- datatable(d,
                   extensions = 'Buttons',
                   options = list(
                     paging = FALSE, searching = TRUE, autoWidth = TRUE, scrollX = TRUE,
                     # fillContainer = TRUE,
                     # bPaginate = FALSE,
                     # autoWidth = TRUE,
                     dom = 'Bftp',
                     buttons = c('copy', 'csv')
                   ), escape = FALSE)
    d
  })
  
  output$datatab_plastic <- DT::renderDT(displayed_data_table_plastic())
  # output$datatab_plastic <- shiny::renderDataTable(displayed_data_table_plastic())#,
                                            # extensions = 'Buttons',
                                            # options = list(dom = 'Bfrtip',
                                            #                buttons = c('copy', 'csv')))

  display_data_table_stations <- reactive({
    d <- filtered_station_data()
    d <- as.data.frame(d)
    d <- d[d$Record_ID %in% selected_data(),]
    if(nrow(d) < 1) return(NULL)
    d <- subset(d, select = -c(Record_ID, Elevation_Datum, tooltip, geometry, order, Distance_From_Coast))#, Webcam_URL))
    d$Photo_URL <- paste0("<a href='", d$Photo_URL, "' target=_blank", ">", d$Photo_URL,"</a>")
    d$Webcam_URL <- paste0("<a href='", d$Webcam_URL, "' target=_blank", ">", d$Webcam_URL,"</a>")
    row.names(d) <- NULL
    d <- datatable(d,
                   extensions = 'Buttons',
                   options = list(
                     paging = TRUE, searching = TRUE, autoWidth = TRUE, scrollX = TRUE,
                     dom = 'Bftp',
                     buttons = c('copy', 'csv')
                   ), escape = FALSE)
    d
  })
  
  output$datatab_stations <- DT::renderDT(display_data_table_stations())
  # output$datatab_stations <- shiny::renderDataTable(display_data_table_stations())

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
shinyApp(ui, server, options = list(launch.browser = TRUE))
