# individual map plots for paper

# 1. background data

make_background_map <- function(background, group, eco = FALSE, trendType = 'unspecified'){
  d <- get(paste(background, 'poly', sep = '_'))
  d <- switch(background,
              chl = subset(d, month == group),
              krill = subset(d, month == group),
              sst = subset(d, month == group),
              pH = subset(d, month == group),
              ship = subset(d, ship_class == group))
  dat <- list(background = d, plastic = DATA_sf[0,],  stations = STATIONS_sf[0,], pltSymbols = pltSymbols)
  mp <- make_plot(dat = dat, background = background, backgroundOnly = TRUE, 
                  displayEcoregions = eco, sstType = trendType, pHType = trendType)
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
ship_class <- 'all'
p_ship <- make_background_map(background, ship_class, eco = FALSE)# + theme_nothing()
p_ship

# Sea surface temperature
background <- 'sst'
month <- 'all'
type = 'trend'
p_sst <- make_background_map(background, month, trendType = type)# + theme_nothing()
p_sst

# pH
background <- 'pH'
month <- 'all'
type = 'trend'
p_pH <- make_background_map(background, month, trendType = type)# + theme_nothing()
p_pH

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










