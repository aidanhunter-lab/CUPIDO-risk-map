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
library(ggnewscale)

# Load/create functions
source('import_data.R') # function (get_data) to load/organise data

significantTrendsOnly <- FALSE

get_data(res = '9x3', baseDirectory = wd_base, shinyDirectory = wd_orig, 
         sstType = 'trend', pHType = 'trend',
         sstTrend_significantOnly = significantTrendsOnly, pHTrend_significantOnly = significantTrendsOnly)

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

