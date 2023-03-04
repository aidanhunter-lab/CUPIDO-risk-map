# Create map of risk

# Preamble: packages/functions --------------------------------------------

# Set directories
wd_orig <- getwd()
wd_base <- '~/Documents/Git Repos/CUPIDO-risk-map'

# Load packages
library(MASS) # for Box-Cox transform
library(car) # for Yeo-Johnson
# library(VGAM)
library(ggplot2)
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


# Load/create functions
source('functions.R') # function (get_data) to load/organise data

significantTrendsOnly <- TRUE

get_data(res = '9x3', baseDirectory = wd_base, shinyDirectory = wd_orig, 
         sstType = 'trend', pHType = 'trend',
         sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly,
         shipSummaryDataOrRaw = 'raw')



sst_poly <- sst_poly[sst_poly$metric != 'p-value',]
pH_poly <- pH_poly[pH_poly$metric != 'p-value',]

plot_fun <- function(data, n = nranks, clrs = Cols){
  # function to plot risk levels
  data$rank <- factor(data$rank, levels = n:1, labels = n:1)
  plt <-
    ggplot() +
    geom_sf(data = data, aes(fill = rank)) +
    scale_fill_manual(values = clrs, breaks = n:1, name = 'risk') +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.background = element_rect(fill = 'white', colour = 'white')
    ) +
    new_scale('fill') +
    geom_sf(data = nc,
            aes(fill = surface),
            show.legend = FALSE) +
    scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
  plt
}



# plot_fun <- function(data, n = nranks){
#   # function to plot risk levels
#   data$rank <- factor(data$rank, levels = n:1, labels = n:1)
#   plt <-
#     ggplot() +
#     geom_sf(data = data, aes(fill = rank)) +
#     scale_fill_brewer(type = 'seq',  palette = 'RdYlGn', direction = 1,
#                       name = 'risk') +
#     theme(
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       panel.background = element_rect(fill = 'white', colour = 'white'),
#       plot.background = element_rect(fill = 'white', colour = 'white')
#     ) +
#     new_scale('fill') +
#     geom_sf(data = nc,
#             aes(fill = surface),
#             show.legend = FALSE) +
#     scale_fill_manual(values = c('grey','skyblue','skyblue','grey'))
#   plt
# }

# The 'subset' function doesn't work with sf objects. Use this simple function to subset sf data
subset_sf <- function(dat, variable, value){
  x <- dat[dat[[variable]] %in% value,]
  if(is.factor(x[[variable]])) x[[variable]] <- as.factor(as.character(x[[variable]]))
  return(x)}

# Use a continuous 'score' to assess risk. This behaves better than integer ranks
# under averaging and may be transformed into an integer 'rank' by taking the ceiling
# (or floor).
score2rank <- function(x, base = 0){
  y <- ceiling(x)
  y[y == base] <- base + 1
  y
}


# Filter by month(s) -- summer months only.
month <- 'all' # select month: 1, 2, 3, and/or 'all'

chl <- subset_sf(chl_poly, 'month', month)
krill <- subset_sf(krill_poly, 'month', month)
sst <- subset_sf(sst_poly, 'month', month)
pH <- subset_sf(pH_poly, 'month', month)

nranks <- 5 # number of risk categories

colfunc <- colorRampPalette(c('firebrick3', 'yellow', 'forestgreen'))
Cols <- setNames(colfunc(nranks), nranks:1)


# # Take a simple approach of using statistical quantiles to rank each map grid
# # cell for each data type. This is not appropriate for all the data types...
# # In fact, it may not be appropriate at all!
# qvals <- seq(0, 1, length = nranks + 1) # evenly distributed quantiles
# qgroup_data <- function(x, q) quantile(x, q)


# Chlorophyll -------------------------------------------------------------

# High productivity regions are associated with greater potential for environmental
# risk, so risk increases with chlorophyll concentration.

# Manually specify values for data rankings by observation of data distribution.

x <- {1:nrow(chl)} / nrow(chl)
y <- sort(chl$value)
plot(x, y)

r <- rep(NA, nranks + 1)
r[nranks:{nranks+1}] <- c(1, max(y)) # values greater than 1 have maximum rank
y_ <- y[y < r[nranks]]
x_ <- {1:length(y_)} / length(y_)
plot(x_, y_)
r[1:2] <- c(min(y), 0.25) # values less than 0.25 have lowest rank

y_ <- y[r[2] < y & y < r[nranks]]
x_ <- {1:length(y_)} / length(y_)
plot(x_, y_)
# Set the remaining ranks evenly across the intermediate values
rna <- is.na(r)
nrna <- sum(rna)
r[rna] <- seq(r[2], r[nranks], length = nrna + 2)[2:{1+nrna}]
plot(x, y)
abline(h = r)

rank_limits_chl <- r

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

plt_chl <- plot_fun(chl)
plt_chl


# Krill -------------------------------------------------------------------

# Krill is a Southern Ocean key-stone species. Regions with high krill density 
# have greater potential for environmental risk.
# The quantile approach will need to be modified for krill because there are lots
# of grid cells with measurement zeros. Where we have zeros the risk must be assigned
# the lowest rank.

# Manually specify values for data rankings by observation of data distribution.

x <- {1:nrow(krill)} / nrow(krill)
y <- sort(krill$value)
par(mfrow = c(1,1))
plot(x, y)
plot(x, log10(y))

zero_krill <- krill$value == 0
observed_krill <- !zero_krill

krill_p <- krill[observed_krill,] # positive values

x_p <- {1:nrow(krill_p)} / nrow(krill_p)
y_p <- sort(krill_p$value)
plot(x_p, y_p)
plot(x_p, log10(y_p))

r <- rep(NA, nranks + 1)
r[nranks:{nranks+1}] <- c(50, max(y_p)) # values greater than 1 have maximum rank
y_ <- y_p[y_p < r[nranks]]
x_ <- {1:length(y_)} / length(y_)
plot(x_, y_)
r[nranks - 1] <- 10
y_ <- y_[y_ < min(r, na.rm = TRUE)]
x_ <- {1:length(y_)} / length(y_)
plot(x_, y_)
r[1:3] <- c(0,1,5)

plot(x, log10(y))
abline(h = log10(r[-1]))

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

plt_krill <- plot_fun(krill)
plt_krill


# Temperature -------------------------------------------------------------

# Potential environmental risk due to temperature change is typically associated
# rising temperature, but the SST data record regions of both rising and falling
# temperatures. I think it must be assumed that a zero trend be associated with
# zero risk, and that regions of rising OR falling temperatures be associated with
# greater risk. A problem with this approach is distinguishing between rising and
# falling temperatures on the resulting risk map. I'm not sure how to address this.

# As above, specify risk rankings by inspection of the data distributions.
ind_sst_zero <- sst$value == 0
ind_sst_neg <- sst$value < 0
ind_sst_pos <- sst$value > 0
any_neg <- any(ind_sst_neg)
any_pos <- any(ind_sst_pos)
sst_neg <- sst$value[ind_sst_neg]
sst_pos <- sst$value[ind_sst_pos]

# Negative values
if(any_neg){
  hist(sst_neg)
  y <- sort(sst_neg)
  ny <- length(y)
  x <- {1:ny} / ny
  plot(x,y)
  # In this case I think that simply splitting the data into evenly spaced groups
  # will produce an appropriate risk distribution.
  if(significantTrendsOnly){
    outlier <- rep(FALSE, ny)
    outlier[c(1:3, ny)] <- TRUE
    ry <- range(y[!outlier])
    rneg <- seq(ry[1], ry[2], length = nranks + 1)
    rneg[c(1, nranks+1)] <- range(y) # adjust to account outlying values
    abline(h = rneg)
  }else{
    outlier <- rep(FALSE, ny)
    outlier[c(1:4, ny)] <- TRUE
    ry <- range(y[!outlier])
    rneg <- seq(ry[1], ry[2], length = nranks + 1)
    rneg[c(1, nranks+1)] <- range(y) # adjust to account outlying values
    abline(h = rneg)
  }
}

# Positive values
if(any_pos){
  # For regions with positive sst trend the same approach will be fine provided the
  # single outlying point is handled separately
  hist(sst_pos)
  y <- sort(sst_pos)
  ny <- length(y)
  x <- {1:ny} / ny
  plot(x,y)
  if(significantTrendsOnly){
    outlier <- rep(FALSE, ny)
    outlier[c(1, ny)] <- TRUE
    ry <- range(y[!outlier])
    rpos <- seq(ry[1], ry[2], length = nranks + 1)
    rpos[c(1, nranks+1)] <- range(y) # adjust to account outlying values
    abline(h = rpos)
  }else{
    outlier <- rep(FALSE, ny)
    outlier[c({ny-2}:ny)] <- TRUE
    ry <- range(y[!outlier])
    rpos <- seq(ry[1], ry[2], length = nranks + 1)
    rpos[c(1, nranks+1)] <- range(y) # adjust to account outlying values
    abline(h = rpos)
  }
}

rank_limits_sst <- rep(0, 2*nranks)
if(any_neg) rank_limits_sst[1:nranks] <- rneg[1:nranks]
if(any_pos) rank_limits_sst[{nranks+1}:{2*nranks}] <- rpos[2:{nranks+1}]

sst$score[ind_sst_zero] <- 0
for(i in 1:nranks){
  # negative values
  rneg <- rank_limits_sst[i:{i+1}]
  if(i > 2){
    j <- rneg[1] < sst_neg & sst_neg <= rneg[2]} else{
      j <- rneg[1] <= sst_neg & sst_neg <= rneg[2]}
  v <- sst_neg[j]
  sst$score[ind_sst_neg][j] <- {nranks - i} + {rneg[2] - v} / diff(rneg)
  # positive values
  rpos <- rank_limits_sst[{nranks+i-1}:{nranks+i}]
  if(i < nranks){
    j <- rpos[1] <= sst_pos & sst_pos < rpos[2]} else{
      j <- rpos[1] <= sst_pos & sst_pos <= rpos[2]}
  v <- sst_pos[j]
  sst$score[ind_sst_pos][j] <- {i - 1} + {v - rpos[1]} / diff(rpos)
}

sst$rank <- score2rank(sst$score)

plt_sst <- plot_fun(sst)
plt_sst


# pH ----------------------------------------------------------------------

# Handle pH similarly to temperature. Bear in mind that this assumption makes less
# sense because pH trends are less evenly distributed around zero. This will
# need changed later...
# Ask Clara about the nature of risk regarding pH changes... Are increases in pH,
# in any way, risky? Probably not for calcifying organisms, but what about more
# generally? Is it reasonable to assign regions exhibiting an increase in pH the
# minimum risk rating? This would make things easier, but I just don;t know if it makes ecological sense...

ind_pH_zero <- pH$value == 0
ind_pH_neg <- pH$value < 0
ind_pH_pos <- pH$value > 0
any_neg <- any(ind_pH_neg)
any_pos <- any(ind_pH_pos)
pH_neg <- pH$value[ind_pH_neg]
pH_pos <- pH$value[ind_pH_pos]

# Negative values
if(any_neg){
  hist(pH_neg)
  y <- sort(pH_neg)
  ny <- length(y)
  x <- {1:ny} / ny
  plot(x,y)
  # In this case I think that simply splitting the data into evenly spaced groups
  # will produce an appropriate risk distribution but, as above, handle outlying
  # points separately.
  if(significantTrendsOnly){
    outlier <- rep(FALSE, ny)
    outlier[c(1:4, ny)] <- TRUE
    ry <- range(y[!outlier])
    rneg <- seq(ry[1], ry[2], length = nranks + 1)
    rneg[c(1, nranks+1)] <- range(y) # adjust to account for outlying values
  }else{
    outlier <- rep(FALSE, ny)
    outlier[c(1:4, {ny-1}:ny)] <- TRUE
    ry <- range(y[!outlier])
    rneg <- seq(ry[1], ry[2], length = nranks + 1)
    rneg[c(1, nranks+1)] <- range(y) # adjust to account for outlying value
  }
  abline(h = rneg)
}

# Positive values
if(any_pos){
  hist(pH_pos)
  y <- sort(pH_pos)
  ny <- length(y)
  x <- {1:ny} / ny
  plot(x,y)
  if(significantTrendsOnly){
    outlier <- rep(FALSE, ny)
    ry <- range(y[!outlier])
    rpos <- seq(ry[1], ry[2], length = nranks + 1)
    rpos[c(1, nranks+1)] <- range(y) # adjust to account for outlying values
  }else{
    outlier <- rep(FALSE, ny)
    ry <- range(y[!outlier])
    rpos <- seq(ry[1], ry[2], length = nranks + 1)
    rpos[c(1, nranks+1)] <- range(y) # adjust to account for outlying values
  }
  abline(h = rpos)
}

rank_limits_pH <- rep(0, 2*nranks)
if(any_neg) rank_limits_pH[1:nranks] <- rneg[1:nranks]
if(any_pos) rank_limits_pH[{nranks+1}:{2*nranks}] <- rpos[2:{nranks+1}]

pH$score[ind_pH_zero] <- 0
for(i in 1:nranks){
  # negative values
  rneg <- rank_limits_pH[i:{i+1}]
  if(i > 2){
    j <- rneg[1] < pH_neg & pH_neg <= rneg[2]} else{
      j <- rneg[1] <= pH_neg & pH_neg <= rneg[2]}
  v <- pH_neg[j]
  pH$score[ind_pH_neg][j] <- {nranks - i} + {rneg[2] - v} / diff(rneg)
  # positive values
  rpos <- rank_limits_pH[{nranks+i-1}:{nranks+i}]
  if(i < nranks){
    j <- rpos[1] <= pH_pos & pH_pos < rpos[2]} else{
      j <- rpos[1] <= pH_pos & pH_pos <= rpos[2]}
  v <- pH_pos[j]
  pH$score[ind_pH_pos][j] <- {i - 1} + {v - rpos[1]} / diff(rpos)
}

pH$rank <- score2rank(pH$score)

plt_pH <- plot_fun(pH)
plt_pH



# Combine data ------------------------------------------------------------

# Combine data sets to generate single risk metric

chl$type <- 'chl'
krill$type <- 'krill'
sst$type <- 'sst'
pH$type <- 'pH'

cols <- c('type','score')

# Merging & manipulating data frames of class 'sf' seems a bit troublesome...
# Maybe there's a slick way of doing this, but for now I've just cobbled something together!
dat <- rbind(chl[cols], krill[cols], sst[cols], pH[cols]) # all data
dat_ <- as.data.frame(dat$geometry)
dat_ <- do.call('rbind', lapply(1:nrow(dat_), function(z) unlist(dat_$geometry[z])))

all_polygons <- unique(dat$geometry)
npolys <- length(all_polygons)

ind_poly <- lapply(1:npolys, function(z){
  i <- dat_ == matrix(unlist(all_polygons[[z]]), nrow = 1)[rep(1, nrow(dat_)),]
  which(apply(i, 1, all))
})
score_poly <- sapply(1:npolys, function(z){
  i <- ind_poly[[z]]
  mean(dat$score[i], na.rm = TRUE)
})

dat <- unique(dat['geometry'])
dat$score <- score_poly

dat$rank <- score2rank(dat$score)

# dat$rank <- factor(dat$rank, levels = nranks:1, labels = nranks:1)


plt_all <- plot_fun(dat)
plt_all

plt_chl <- plt_chl + labs(title = 'Chlorophyll')
plt_krill <- plt_krill + labs(title = 'Krill')
plt_sst <- plt_sst + labs(title = 'Sea surface temperature')
plt_pH <- plt_pH + labs(title = 'pH')
plt_all <- plt_all + labs(title = 'Combined (mean) risk')

plt_combined <- gridExtra::grid.arrange(plt_chl, plt_krill, plt_sst, plt_pH, plt_all,
                                        nrow = 3)

plot(plt_combined)


# Make some summary plots....

# Chlorophyll
x <- {1:nrow(chl)} / nrow(chl)
y <- sort(chl$value)
plot(x, y)
abline(h = rank_limits_chl)

chl$rank <- factor(chl$rank, levels = nranks:1, labels = nranks:1)

plt_dist_chl <- ggplot(chl) + 
  geom_histogram(aes(x = value, fill = rank), binwidth = 0.01) +
  xlab(expression(Chlorophyll~a~(mg~m^{-3}))) +
  labs(title = 'Data distribution displaying cut-off values used for risk rankings') +
  scale_fill_manual(values = Cols, name = 'risk') +
  theme_bw()
plt_dist_chl

# Krill
x <- {1:nrow(krill)} / nrow(krill)
y <- sort(krill$value)
plot(x, y)
abline(h = rank_limits_krill)
plot(x, y, log = 'y')
abline(h = log10(rank_limits_krill[-1]))

krill$rank <- factor(krill$rank, levels = nranks:1, labels = nranks:1)

plt_dist_krill <- ggplot(krill) + 
  geom_histogram(aes(x = value, fill = rank), binwidth = 0.05) +
  scale_x_log10() + 
  xlab(expression(Krill~(number~m^{-2}))) +
  labs(title = 'Data distribution displaying cut-off values used for risk rankings', subtitle = 'Note the log scale omits measurement zeros') +
  scale_fill_manual(values = Cols, name = 'risk') +
  theme_bw()
plt_dist_krill


# SST
x <- {1:nrow(sst)} / nrow(sst)
y <- sort(sst$value)
plot(x, y)
abline(h = rank_limits_sst)

sst$rank <- factor(sst$rank, levels = nranks:1, labels = nranks:1)

plt_dist_sst <- ggplot(subset(sst, value != 0)) + 
  geom_histogram(aes(x = value, fill = rank), binwidth = 0.0005) +
  xlab(expression(SST~linear~trend~(degree*C~year^{-1}))) +
  labs(title = 'Data distribution displaying cut-off values used for risk rankings',
       subtitle = 'Note that both increases & decreases in SST are deemed risky') +
  scale_fill_manual(values = Cols, name = 'risk') +
  theme_bw()
plt_dist_sst


# pH
x <- {1:nrow(pH)} / nrow(pH)
y <- sort(pH$value)
plot(x, y)
abline(h = rank_limits_pH)

pH$rank <- factor(pH$rank, levels = nranks:1, labels = nranks:1)

plt_dist_pH <- ggplot(subset(pH, value != 0)) + 
  geom_histogram(aes(x = value, fill = rank), binwidth = 0.000025) +
  xlab(expression(pH~linear~trend~(year^{-1}))) +
  labs(title = 'Data distribution displaying cut-off values used for risk rankings',
       subtitle = 'Note that both increases & decreases in pH are deemed risky...\nshould positive trends be ascribed the lowest risk?') +
  scale_fill_manual(values = Cols, name = 'risk') +
  theme_bw()
plt_dist_pH


chl_combined <- gridExtra::grid.arrange(plt_dist_chl, plt_chl, nrow = 2)
krill_combined <- gridExtra::grid.arrange(plt_dist_krill, plt_krill, nrow = 2)
sst_combined <- gridExtra::grid.arrange(plt_dist_sst, plt_sst, nrow = 2)
pH_combined <- gridExtra::grid.arrange(plt_dist_pH, plt_pH, nrow = 2)

# save plots

ggsave('chl_risk.png', plot = chl_combined, device = 'png', width = 6, height = 12, units = 'in')
ggsave('krill_risk.png', plot = krill_combined, device = 'png', width = 6, height = 12, units = 'in')
ggsave('sst_risk.png', plot = sst_combined, device = 'png', width = 6, height = 12, units = 'in')
ggsave('pH_risk.png', plot = pH_combined, device = 'png', width = 6, height = 12, units = 'in')

ggsave('combined_risk_map.png', plot = plt_combined, device = 'png', width = 8, height = 12, units = 'in')



# 2nd attempt -- after meeting with Clara & Sally -------------------------

# The novelty of the work is assessing the potential of plastic to contribute to
# environmental stress. Thus, at each stage of risk mapping I need to account for
# impact of plastic.
# Consider two kinds of 'risk':
# (1) Risk of interaction between plastic and biota.
# (2) Plastic contribution to multi-stressor risk.
# For (1) risk of interaction we use data on spatial distributions of primary productivity
# (chlorophyll) and of a keystone zooplankton species (krill) to assess where the
# distributions of these biota are likely to overlap with areas of relatively high
# plastic emissions.
# For (2) multi-stressor risk we spatial distribution data on long-term trends of
# sea surface temperature and pH to investigate where existing areas of concern 
# overlap with areas of relatively high plastics emissions. These variables were
# selected as many regard their long term variability to pose a risk to the existing
# ecosystem structure and function.

# Note that although temperature and pH may decrease or increase over time and that
# there are regions where each of these variables display increases/decreases, we
# concentrate only on the risk of increasing temperature and decreasing pH because
# these trends are currently the focus of much research effort, whereas the impacts
# of decreasing temperature and rising pH in Antarctic waters are lesser known due
# to much less research effort/funding. This is important to note with respect to
# temperature because there is a greater sea surface area of temperature decline
# than there is of temperature rise, whereas almost all Antarctic waters exhibit
# a decline in pH.

# Take roughly the same approach as my, above, 1st attempt...

# Preamble: packages/functions --------------------------------------------

# Set directories
wd_orig <- getwd()
wd_base <- '~/Documents/Git Repos/CUPIDO-risk-map'

# Load packages
library(MASS) # for Box-Cox transform
library(car) # for Yeo-Johnson
# library(VGAM)
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

# Load/create functions
source('functions.R') # function (get_data) to load/organise data

significantTrendsOnly <- TRUE

get_data(res = '9x3', baseDirectory = wd_base, shinyDirectory = wd_orig, 
         sstType = 'trend', pHType = 'trend', shipSummaryDataOrRaw = 'raw',
         sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly)

plot_fun <- function(data, n = nranks, v_option = 'viridis', v_direction = -1,
                     legend_title = 'risk', Title = NULL, show_legend = TRUE){
  # function to plot risk levels
  data$rank <- factor(data$rank, levels = n:1, labels = n:1)
  # get the legend
  plt <-
    ggplot() +
    geom_sf(data = data, aes(fill = rank), show.legend = TRUE) +
    scale_fill_viridis_d(option = v_option, direction = v_direction,
                         drop = FALSE, na.translate = FALSE, name = legend_title)
  leg <- ggdraw(get_legend(plt))
  # make full plot
  plt <-
    ggplot() +
    geom_sf(data = data, aes(fill = rank), show.legend = show_legend) +
    scale_fill_viridis_d(option = v_option, direction = v_direction,
                         drop = FALSE, na.translate = FALSE, name = legend_title) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      plot.background = element_rect(fill = 'white', colour = 'white')
    ) +
    new_scale('fill') +
    geom_sf(data = nc,
            aes(fill = surface)) +
    scale_fill_manual(values = c('grey','skyblue','skyblue','grey'), guide = 'none')
  if(!is.null(Title)) plt <- plt + labs(title = Title)
  return(list(plot = plt, legend = leg))
}



# We need a separate function to plot stations
plot_fun_stations <- function(data_stations, data_grid_cells, rank_type = 'additive', n = nranks,
                              v_option = 'viridis', v_direction = -1, alpha = 0.75,
                              stroke = 1, rank_limits_PopSize = NULL, Title = NULL, show_legend = TRUE){

  # if(!{'rank' %in% names(data_stations)}){
  #   dn <- names(data_stations)
  #   dr <- dn[grepl(rank_type, dn) & grepl('rank', dn)]
  #   data_stations$rank <- data_stations[[dr]]
  # }
  data_stations$rank <- factor(data_stations$rank, levels = n:1, labels = n:1)

  if(is.null(rank_limits_PopSize)){
    data_stations$rank_pop_labels <- data_stations$rank_pop
  }else{
    # rl <- paste(rank_limits_PopSize[1:n], '\u2264')
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
    scale_fill_viridis_d(option = v_option, direction = v_direction, name = 'disturbance',
                         drop = FALSE, na.translate = FALSE, guide = guide_legend(order = 3)) +
    new_scale('fill') +
    geom_sf(data = nc,
            aes(fill = surface)) +
    scale_fill_manual(values = c('grey','skyblue','skyblue','grey'), guide = 'none') +
    new_scale('fill') +
    geom_sf(aes(shape = Seasonality, size = rank_pop_labels, fill = rank),
            alpha = alpha, stroke = stroke) +
    scale_fill_viridis_d(option = v_option, direction = v_direction, guide = 'none') +
    scale_fill_viridis_d(aesthetics = 'fill_new_new', option = v_option, direction = v_direction, name = 'disturbance', drop = FALSE, na.translate = FALSE, guide = guide_legend(override.aes = list(size = NA, shape = NA))) +
    scale_size_manual(values = 1:n, guide = guide_legend(title = 'population', order = 1, override.aes = list(fill = NA, linetype = 0))) +
    scale_shape_manual(values = c(21, 22), guide = guide_legend(title = 'seasonality', order = 2, override.aes = list(size = n, fill = NA, linetype = 0)))
  
  if(!is.null(Title)) plt <- plt + labs(title = Title)
  leg <- suppressWarnings(ggdraw(get_legend(plt)))
  if(!show_legend) plt <- plt + theme(legend.position = 'none')
  output <- list(plot = plt, legend = leg)
  return(output)
}


# The 'subset' function doesn't work with sf objects. Use this simple function to subset sf data
subset_sf <- function(dat, variable, value){
  x <- dat[dat[[variable]] %in% value,]
  if(is.factor(x[[variable]])) x[[variable]] <- as.factor(as.character(x[[variable]]))
  return(x)}

# Use a continuous 'score' to assess risk. This behaves better than integer ranks
# under averaging and may be transformed into an integer 'rank' by taking the ceiling
# (or floor).
score2rank <- function(x, base = 0){
  y <- ceiling(x)
  if(is.numeric(base)) y[y == base] <- base + 1
  y
}

# Useful functions for rounding -- used in automated process for selecting ranking limit values
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
orderOfMagnitude <- function(x){
  if(x == 0) out <- 1 else{
    z <- zerosAfterDecimal(x)
    if(z != 0) out <- 10 ^ -{z+1} else{
      y <- strsplit(as.character(abs(x)), '\\.')[[1]][1]
      if(y == '0') out <- 0.1 else{
        out <- 10 ^ {nchar(y) - 1}
      }}}
  return(out)
}
orderOfMagnitude <- Vectorize(orderOfMagnitude)


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
# ship <- subset_sf(ship_poly, 'activity', ship_type)


nranks <- 5 # number of risk categories

# Chlorophyll -------------------------------------------------------------

# High productivity regions are associated with greater potential for environmental
# risk, so risk increases with chlorophyll concentration.

# Specify values for data rankings by observation of data distribution.

x <- {1:nrow(chl)} / nrow(chl)
y <- sort(chl$value)
plot(x, y)
plot(x, y, log = 'y')

# r <- rep(NA, nranks + 1)
# r[nranks:{nranks+1}] <- c(1, max(y)) # values greater than 1 have maximum rank
# y_ <- y[y < r[nranks]]
# x_ <- {1:length(y_)} / length(y_)
# plot(x_, y_)
# r[1:2] <- c(min(y), 0.25) # values less than 0.25 have lowest rank
# 
# y_ <- y[r[2] < y & y < r[nranks]]
# x_ <- {1:length(y_)} / length(y_)
# plot(x_, y_)
# # Set the remaining ranks evenly across the intermediate values
# rna <- is.na(r)
# nrna <- sum(rna)
# r[rna] <- seq(r[2], r[nranks], length = nrna + 2)[2:{1+nrna}]
# plot(x, y)
# abline(h = r)
# 
# rank_limits_chl <- r


yl <- log10(y)
r <- 10 ^ seq(min(yl), max(yl), length = nranks + 1) # even spread on log scale
r <- round(r * 10) / 10 # round to nearest 0.1
# r[2:5] <- c(0.2, 0.4, 0.6, 1) # manually adjust for nice rounded values
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

options(scipen = 999)
plt_chl_raw <- 
  ggplot(chl) + 
  geom_point(aes(x = index, y = value, fill = rank), shape = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis_c(guide = 'legend') +
  guides(fill = guide_legend(title = 'prevalence', reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0) +
  xlab(expression(measurement ~ index)) +
  ylab(expression(chlorophyll ~ abundance ~ (mg / m^3))) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())
options(scipen = 0)

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

r <- rep(NA, nranks + 1)
r[nranks:{nranks+1}] <- c(50, max(y_p)) # values greater than 1 have maximum rank
y_ <- y_p[y_p < r[nranks]]
x_ <- {1:length(y_)} / length(y_)
plot(x_, y_)
r[nranks - 1] <- 10
y_ <- y_[y_ < min(r, na.rm = TRUE)]
x_ <- {1:length(y_)} / length(y_)
plot(x_, y_)
r[1:3] <- c(0,1,5)

plot(x, y, log = 'y')
abline(h = r[-1])

rank_limits_krill <- r

# The above, manually specified, ranks are nice and sensible, but try specifying
# rankings as an even spread on log scale after removing outliers.
plot(x_p, y_p, log = 'y')
y_r <- y_p[-c(1,length(y_p))]
x_r <- 1:length(y_r)
plot(x_r, y_r, log = 'y')
yl <- log10(y_r)

r <- 10 ^ seq(min(yl), max(yl), length = nranks + 1)
r[c(1, nranks+1)] <- c(0, max(y)) # adjust for outliers and measurement zeros
r <- round(r / orderOfMagnitude(r)) * orderOfMagnitude(r) # round to nice values
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

options(scipen = 999)
plt_krill_raw <- 
  ggplot(krill) + 
  geom_point(aes(x = index, y = value, fill = rank), shape = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis_c(guide = 'legend') +
  guides(fill = guide_legend(title = 'prevalence', reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0) +
  xlab(expression(measurement ~ index)) +
  ylab(expression(krill ~ abundance ~ (individuals / m^2))) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())
options(scipen = 0)


# Temperature -------------------------------------------------------------

# Potential environmental risk due to temperature change is typically associated
# rising temperature, but the SST data record regions of both rising and falling
# temperatures. I think it must be assumed that a zero trend be associated with
# zero risk, and that regions of rising OR falling temperatures be associated with
# greater risk. A problem with this approach is distinguishing between rising and
# falling temperatures on the resulting risk map. I'm not sure how to address this.

# Clara and Sally disputed my above statement about falling temperature being equally
# associated with risk of disturbance, saying that I should map the risk of increasing
# temperature only. I disagree and believe it more honest and useful to map regions 
# of both increasing and decreasing temperatures, but let's try it their way...

sst_pos <- sst$value > 0 # index grid cells with trend of increasing temperature (all p-values)
sst_neg <- !sst_pos # cells with negative values are ignored(!) and assigned either minimum risk rank or omitted from map

# # index cells based on significance of trend -- ignoring negative trends :(
# sig_0.05 <- sst_pos & sst$pval <= 0.05 # all significant positive trends
# sig_0.1 <- sst_pos &  0.05 < sst$pval & sst$pval <= 0.1 # insignificant positive trends
# sig_0.2 <- sst_pos &  0.1 < sst$pval & sst$pval <= 0.2 # even more insignificant positive trends
# # Redefine the categorical variable, plevel
# sst$plevel <- NA
# sst$plevel[sig_0.05] <- 0.05
# sst$plevel[sig_0.1] <- 0.1
# sst$plevel[sig_0.2] <- 0.2
# sst$plevel <- factor(sst$plevel, levels = c(0.05, 0.1, 0.2), labels = c(paste('p', '\u2264', 0.05),
#                                                           paste(0.05, '\u2264', 'p', '\u2264', 0.1),
#                                                           paste(0.1, '\u2264', 'p', '\u2264', 0.2)))

sig_0.05 <- sst_pos & sst$pval <= 0.05 # all significant positive trends
sig_0.25 <- sst_pos &  0.05 < sst$pval & sst$pval <= 0.25 # insignificant positive trends
sig_0.5 <- sst_pos &  0.25 < sst$pval & sst$pval <= 0.5 # insignificant positive trends
sig_1 <- sst_pos &  0.5 < sst$pval & sst$pval <= 1 # even more insignificant positive trends
# Redefine the categorical variable, plevel
sst$plevel <- NA
sst$plevel[sig_0.05] <- 0.05
sst$plevel[sig_0.25] <- 0.25
sst$plevel[sig_0.5] <- 0.5
sst$plevel[sig_1] <- 1
sst$plevel <- factor(sst$plevel, levels = c(0.05, 0.25, 0.5, 1), labels = c(paste('p', '\u2264', 0.05),
                                                                        paste(0.05, '\u2264', 'p', '\u2264', 0.25),
                                                                        paste(0.25, '\u2264', 'p', '\u2264', 0.5),
                                                                        paste(0.5, '\u2264', 'p', '\u2264', 1)))

sst_selected_cells <- !is.na(sst$plevel)
sst_selected <- sst[sst_selected_cells,]
o <- order(sst_selected$value)

y <- sst_selected$value[o]
pv <- sst_selected$plevel[o]
x <- 1:length(y)
plot(x, y, col = pv)
legend('topleft', legend = levels(pv), pch = 1, col = 1:length(levels(pv)))

# trim outliers then group
# ol <- c(1:2, length(y))
# ol <- c(1:4, length(y))
ol <- c({length(y)-3}:length(y))
y_r <- y[-ol]
pv_r <- pv[-ol]
x_r <- 1:length(y_r)
# plot(x_r, y_r, col = pv_r)
# legend('topleft', legend = levels(pv), pch = 1, col = 1:length(levels(pv_r)))

yrange <- range(y_r)
# rlo <- round(yrange[1] / orderOfMagnitude(yrange[1])) * orderOfMagnitude(yrange[1]) # round to nearest 0.001
rlo <- floor(yrange[1] / orderOfMagnitude(yrange[1])) * orderOfMagnitude(yrange[1]) # round to nearest 0.001
# rhi <- ceiling(yrange[2] / {orderOfMagnitude(yrange[2]) / 10 * 5}) * orderOfMagnitude(yrange[2]) / 10 * 5 # round up to nearest 0.005
rhi <- ceiling(yrange[2] / orderOfMagnitude(yrange[2]) * 10) * orderOfMagnitude(yrange[2]) / 10 # round up
# r <- c(0, seq(rlo, rhi, length = nranks))
r <- c(seq(rlo, rhi, length = nranks + 1))
r <- round(r / orderOfMagnitude(r) * 1) * orderOfMagnitude(r) / 1
r[c(1, nranks + 1)] <- range(y)
abline(h = r[2:nranks], lty = 3)

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
# sst$score[!sst_selected_cells] <- NA
sst <- sst[sst_selected_cells,]

sst$rank <- score2rank(sst$score)

plt_sst <- plot_fun(sst, Title = 'SST rising trends', legend_title = 'stress')
plt_sst$plot

# Nice plot of raw data
sst$order <- 1:nrow(sst)
sst <- sst[order(sst$value),]
sst$index <- 1:nrow(sst)
sst <- sst[order(sst$order),]

options(scipen = 999)
plt_sst_raw <- 
  ggplot(sst) + 
  geom_point(aes(x = index, y = value, fill = rank), shape = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis_c(guide = 'legend') +
  guides(fill = guide_legend(title = 'stress', reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0) +
  xlab(expression(measurement ~ index)) +
  ylab(expression(SST ~ trend ~ (degree * C / year))) +
#  scale_y_log10() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())
options(scipen = 0)



# pH -------------------------------------------------------------

pH_neg <- pH$value < 0 # index grid cells with trend of decreasing pH (all p-values)
pH_pos <- !pH_neg # cells with negative values are ignored(!) and assigned either minimum risk rank or omitted from map

pH_selected_cells <- pH_neg
pH_selected <- pH[pH_selected_cells,]

y <- sort(pH_selected$value)
x <- 1:length(y)
plot(x, y)

# trim outliers then group
# ol <- c(1:2, {length(y)-1}:length(y))
ol <- c(1, length(y))
y_r <- y[-ol]
x_r <- 1:length(y_r)
# plot(x_r, y_r)

yrange <- range(y_r)
rlo <- floor(yrange[1] / orderOfMagnitude(yrange[1]) * 10) * orderOfMagnitude(yrange[1]) / 10 # round to nearest 0.0001
rhi <- ceiling(yrange[2] / orderOfMagnitude(yrange[2]) *10) * orderOfMagnitude(yrange[2]) / 10 # round up to nearest 0.0001
r <- seq(rlo, rhi, length = nranks+1)
r <- round(r / orderOfMagnitude(r) * 10) * orderOfMagnitude(r) / 10
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
# pH$score[!pH_selected_cells] <- NA
pH <- pH[pH_selected_cells,]
pH$rank <- score2rank(pH$score)

plt_pH <- plot_fun(pH, Title = 'pH declining trends', legend_title = 'stress')
plt_pH$plot

# Nice plot of raw data
pH$order <- 1:nrow(pH)
pH <- pH[order(pH$value),]
pH$index <- 1:nrow(pH)
pH <- pH[order(pH$order),]

options(scipen = 999)
plt_pH_raw <- 
  ggplot(pH) + 
  geom_point(aes(x = index, y = value, fill = rank), shape = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis_c(guide = 'legend') +
  guides(fill = guide_legend(title = 'stress', reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0) +
  xlab(expression(measurement ~ index)) +
  ylab(expression(pH ~ trend ~ (1 / year))) +
  # scale_y_log10() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())
options(scipen = 0)


# Shipping ----------------------------------------------------------------

x <- {1:nrow(ship)} / nrow(ship)
y <- sort(ship$value)
par(mfrow = c(1,1))
plot(x, y)
plot(x, y, log = 'y')

# For now, just manually choose the cut-offs for risk ranks...
r <- c(min(y),1,10,30,100,max(y))
abline(h = r)
r_orig <- r

rank_limits_ship <- r

# Again, try taking a more automated approach.
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
r <- round(r / orderOfMagnitude(r)) * orderOfMagnitude(r) # round to nice values
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

plt_ship <- plot_fun(ship, Title = 'Ship traffic', legend_title = 'disturbance')
plt_ship$plot

# Nice plot of raw data
ship$order <- 1:nrow(ship)
ship <- ship[order(ship$value),]
ship$index <- 1:nrow(ship)
ship <- ship[order(ship$order),]

options(scipen = 999)
plt_ship_raw <- 
  ggplot(ship) + 
  geom_point(aes(x = index, y = value, fill = rank), shape = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis_c(guide = 'legend') +
  guides(fill = guide_legend(title = 'prevalence', reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0) +
  xlab(expression(measurement ~ index)) +
  ylab(expression(ship ~ traffic ~ (ship ~ days / year))) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())
options(scipen = 0)



# Facilities --------------------------------------------------------------

# Each facility is assigned a rank determined by its population size, then mapped
# grid cells are assigned a disturbance rank determined by their population density
# (people/km^2) calculated as a function of distance to facilities.

stations <- STATIONS_sf

stations$Seasonality <- factor(stations$Seasonality, levels = c('Year-Round','Seasonal'), labels = c('year-round', 'seasonal'))

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

# For now, just manually choose the cut-offs for risk ranks...
r <- c(min(y),10,25,50,100,max(y))
abline(h = r)

# Account for seasonality by reducing population size for seasonal facilities.
# Leave population size of year-round facilities unchanged.
stations$Peak_Population_adjusted <- stations$Peak_Population * 0.5 * {{stations$Seasonality == 'year-round'} + 1}
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

# # Make data frame of mapped grid cells
# grid_cells <- subset(chl, select = geometry) # get mapped grid cells
# # Get grid cell centroids
# grid_cells <- st_transform(grid_cells, crs_world)
# grid_centroids <- st_centroid(grid_cells)
# grid_centroids_ <- as.data.frame(st_coordinates(grid_centroids))
# grid_cells$centroid_x <- grid_centroids_$X
# grid_cells$centroid_y <- grid_centroids_$Y

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
# stations$grid_cell[ind] <- nearest_cell[ind]


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
if(all(attr(grid_cell_area, 'units')$numerator == 'm')) grid_cell_area <- grid_cell_area * 1e-6 # convert m^2 to km^2
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
# distanceLimit <- NULL
distanceLimit <- 1000 # set a 1000 km distance limit, beyond which facility populations do not have impact

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
  # a <- grid_cells$area[cell]
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

rank_limits_stations[2:{nranks}] <- round(rank_limits_stations[2:{nranks}] / orderOfMagnitude(rank_limits_stations[2:{nranks}])) * orderOfMagnitude(rank_limits_stations[2:{nranks}])

rank_limits_stations[nranks+1] <- max(y) # account for outliers
# plot(x_, y_, log = 'y')
abline(h = rank_limits_stations)

# Manually specify ranking limits
if(is.null(distanceLimit))
  rank_limits_stations <- c(0.001, 0.01, 0.1, 1, 10, 200)
# abline(h = rank_limits_stations)

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

# plt_ship <- plot_fun(ship, legend_title = 'disturbance', show_legend = FALSE)
# plt_ship$plot <- plt_ship$plot + labs(title = 'Ship traffic')

# pltStations <- plot_fun_stations(data_stations = stations, data_grid_cells = grid_cells_reduced, rank_limits_PopSize = rank_limits_stationPopSize, Title = 'Facilities', show_legend = FALSE)

grid.arrange(plt_ship$plot + theme(legend.position = 'none'), pltStations$plot + theme(legend.position = 'none'), pltStations$legend,
             ncol = 3, widths = c(0.425,0.425,0.15))


# Nice plot of raw data
grid_cells_reduced$order <- 1:nrow(grid_cells_reduced)
grid_cells_reduced <- grid_cells_reduced[order(grid_cells_reduced$popDensity),]
grid_cells_reduced$index <- 1:nrow(grid_cells_reduced)
grid_cells_reduced <- grid_cells_reduced[order(grid_cells_reduced$order),]

options(scipen = 999)
plt_stations_raw <- 
  ggplot(grid_cells_reduced) + 
  geom_point(aes(x = index, y = popDensity, fill = rank), shape = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis_c(guide = 'legend') +
  guides(fill = guide_legend(title = 'disturbance', reverse = TRUE)) +
  geom_hline(yintercept = r[2:nranks], linetype = 3) +
  annotate('text', x = 1, y = r[2:nranks], label = as.character(r[2:nranks]), vjust = 0) +
  xlab(expression(measurement ~ index)) +
  # ylab(expression(population ~ density ~ 'in' ~ grid ~ cells ~ (people / 10000*km^2))) +
  ylab(expression(pop. ~ density ~ (people / 10000 ~ km^2))) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank())
options(scipen = 0)





# Combine -----------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~
# Combine all plots
#~~~~~~~~~~~~~~~~~~

plt_chl_combined <- ggdraw(arrangeGrob(plt_chl_raw,
                                       plt_chl$plot + labs(title = element_blank()),
                                       nrow = 2, heights = c(0.35, 0.65), top = plt_chl$plot$labels$title))
plt_krill_combined <- ggdraw(arrangeGrob(plt_krill_raw,
                                         plt_krill$plot + labs(title = element_blank()),
                                         nrow = 2, heights = c(0.35, 0.65), top = plt_krill$plot$labels$title))
plt_sst_combined <- ggdraw(arrangeGrob(plt_sst_raw,
                                       plt_sst$plot + labs(title = element_blank()),
                                       nrow = 2, heights = c(0.35, 0.65), top = plt_sst$plot$labels$title))
plt_pH_combined <- ggdraw(arrangeGrob(plt_pH_raw,
                                      plt_pH$plot + labs(title = element_blank()),
                                      nrow = 2, heights = c(0.35, 0.65), top = plt_pH$plot$labels$title))
plt_ship_combined <- ggdraw(arrangeGrob(plt_ship_raw,
                                        plt_ship$plot + labs(title = element_blank()),
                                        nrow = 2, heights = c(0.35, 0.65), top = plt_ship$plot$labels$title))
plt_stations_combined <- ggdraw(arrangeGrob(plt_stations_raw,
                                            pltStations$plot + labs(title = element_blank()),
                                            nrow = 2, heights = c(0.35, 0.65), top = pltStations$plot$labels$title))

plt_biota <- ggdraw(arrangeGrob(plt_chl$plot + theme(legend.position = 'none'),
                                plt_krill$plot + theme(legend.position = 'none'),
                                plt_chl$legend,
                                nrow = 1, widths = c(0.425, 0.425, 0.15), top = 'Biota'))
plt_stress <- ggdraw(arrangeGrob(plt_sst$plot + theme(legend.position = 'none'),
                                 plt_pH$plot + theme(legend.position = 'none'),
                                 plt_pH$legend,
                                 nrow = 1, widths = c(0.425, 0.425, 0.15), top = 'Stressors'))
plt_activity <- ggdraw(arrangeGrob(plt_ship$plot + theme(legend.position = 'none'),
                                   pltStations$plot + theme(legend.position = 'none'),
                                   pltStations$legend,
                                   nrow = 1, widths = c(0.425, 0.425, 0.15), top = 'Activity'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge data sets to create risk rankings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Method of combining rankings: mean of score, or sum of score
combineMethods <- c('mean', 'sum')
combineMethod <- combineMethods[1]
combineTitle <- paste0(toupper(substr(combineMethod, 1, 1)), substr(combineMethod, 2, nchar(combineMethod)))

# Biota

plt_biota

chl$type <- 'chl'
krill$type <- 'krill'
columns <- c('cell_index', 'type', 'score')
# Merge the data frames
dat <- rbind(chl[columns], krill[columns])
all_polygons <- sort(unique(dat$cell_index)) # find all mapped polygons
npolys <- length(all_polygons)
score_poly <- switch(combineMethod,
                     mean = sapply(all_polygons, function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
                     sum = pmin(nranks, sapply(all_polygons, function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
                     ) # calculate score for each polygon
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_biota <- dat

plt_chl_krill <- plot_fun(dat_biota, Title = combineTitle, legend_title = 'prevalence')
plt_chl_krill$plot

plt_biota <- ggdraw(arrangeGrob(plt_chl$plot + theme(legend.position = 'none'),
                                plt_krill$plot + theme(legend.position = 'none'),
                                plt_chl_krill$plot + theme(legend.position = 'none'),
                                plt_chl$legend,
                                nrow = 1, widths = c(0.3, 0.3, 0.3, 0.1), top = 'Biota'))
plt_biota

# Stressors

plt_stress

sst$type <- 'sst'
pH$type <- 'pH'
# Merge the data frames
dat <- rbind(sst[columns], pH[columns])
all_polygons <- sort(unique(dat$cell_index)) # find all mapped polygons
npolys <- length(all_polygons)
score_poly <- switch(combineMethod,
                     mean = sapply(all_polygons, function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
                     sum = pmin(nranks, sapply(all_polygons, function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_stress <- dat

plt_sst_pH <- plot_fun(dat, Title = combineTitle, legend_title = 'stress')
plt_sst_pH$plot

plt_stress <- ggdraw(arrangeGrob(plt_sst$plot + theme(legend.position = 'none'),
                                 plt_pH$plot + theme(legend.position = 'none'),
                                 plt_sst_pH$plot + theme(legend.position = 'none'),
                                 plt_sst$legend,
                                 nrow = 1, widths = c(0.3, 0.3, 0.3, 0.1), top = 'Stressors'))
plt_stress

# Activity

plt_activity

ship$type <- 'ship'
grid_cells_reduced$type <- 'station'

# Merge the data frames
dat <- rbind(ship[columns], grid_cells_reduced[columns])
all_polygons <- sort(unique(dat$cell_index)) # find all mapped polygons
npolys <- length(all_polygons)
score_poly <- switch(combineMethod,
                     mean = sapply(all_polygons, function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
                     sum = pmin(nranks, sapply(all_polygons, function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon
dat <- grid_cells_all[grid_cells_all$cell_index %in% all_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_activity <- dat

plt_ship_station <- plot_fun(dat, Title = combineTitle, legend_title = 'presence')
plt_ship_station$plot

plt_activity <- ggdraw(arrangeGrob(plt_ship$plot + theme(legend.position = 'none'),
                                   pltStations$plot + theme(legend.position = 'none'),
                                   plt_ship_station$plot + theme(legend.position = 'none'),
                                   pltStations$legend,
                                   nrow = 1, widths = c(0.3, 0.3, 0.3, 0.1), top = 'Activity'))
plt_activity


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Where does activity overlap with biota & stressors?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Biota
plt_biota_activity <- ggdraw(arrangeGrob(plt_chl_krill$plot + theme(legend.position = 'none') + labs(title = 'Biota'),
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
score_poly <- switch(combineMethod,
                     mean = sapply(joint_polygons, function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
                     sum = pmin(nranks, sapply(all_polygons, function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon

dat <- grid_cells_all[grid_cells_all$cell_index %in% joint_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_interaction <- dat

plt_biota_activity <- plot_fun(dat_interaction, Title = combineTitle, legend_title = 'risk')
plt_biota_activity$plot

plt_interaction <- ggdraw(arrangeGrob(plt_chl_krill$plot + labs(title = 'Biota') + theme(legend.position = 'none'),
                                      plt_chl_krill$legend,
                                      plt_ship_station$plot + labs(title = 'Activity') + theme(legend.position = 'none'),
                                      plt_ship_station$legend,
                                      plt_biota_activity$plot + labs(title = 'Interaction') + theme(legend.position = 'none'),
                                      plt_biota_activity$legend,
                                      nrow = 1, widths = c(0.25, 1/12, 0.25, 1/12, 0.25, 1/12), top = 'Interaction with biota'))
plt_interaction


# Stressors
plt_stessor_activity <- ggdraw(arrangeGrob(plt_sst_pH$plot + theme(legend.position = 'none') + labs(title = 'Stressors'),
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
score_poly <- switch(combineMethod,
                     mean = sapply(joint_polygons, function(z) mean(dat$score[dat$cell_index == z], na.rm = TRUE)),
                     sum = pmin(nranks, sapply(all_polygons, function(z) sum(dat$score[dat$cell_index == z], na.rm = TRUE)))
) # calculate score for each polygon

dat <- grid_cells_all[grid_cells_all$cell_index %in% joint_polygons,] # create new data frame
dat$score <- score_poly
dat$rank <- score2rank(dat$score)
dat_multistressor <- dat

plt_stress_activity <- plot_fun(dat_multistressor, Title = combineTitle, legend_title = 'risk')
plt_stress_activity$plot

plt_multistressor <- ggdraw(arrangeGrob(plt_sst_pH$plot + labs(title = 'Stressors') + theme(legend.position = 'none'),
                                      plt_sst_pH$legend,
                                      plt_ship_station$plot + labs(title = 'Activity') + theme(legend.position = 'none'),
                                      plt_ship_station$legend,
                                      plt_stress_activity$plot + labs(title = 'Interaction') + theme(legend.position = 'none'),
                                      plt_stress_activity$legend,
                                      nrow = 1, widths = c(0.25, 1/12, 0.25, 1/12, 0.25, 1/12), top = 'Contribution to stressors'))
plt_multistressor





