suppressMessages({
  library("ggplot2")
  library("dplyr")
  library("plyr")
  library("ggpubr")
  library("viridis") # viridis colour scale
  library("reshape2")

  ## Load ismev library for marginal GPD fits:
  library("ismev")
  library("evd")

  # Removed these dependencies as they were causing problems. They are just
  # needed to make the plots pretty (e.g., showing degrees on the axis).
  #library("ggspatial")
  #library("ggOceanMapsData")
  #library("ggOceanMaps")

  library("scales") # plotting against time
})

# Base map for the Red Sea
#dt <- data.frame(lon = c(37.1, 37.1,  43,  43),
#                 lat = c(15.5, 20.4, 20.4, 15.5))
#RedSea <- basemap(data = dt)
RedSea <- ggplot() # no base map

# ---- Helper functions ----

## Function to transform to Laplace margins:
laplaceMargins <- function(univariateData, thrQ = 0.95){

  size <- length(univariateData)
  empiricalMargins <- rank(univariateData)/(size+1)

  empiricalMargins[univariateData>quantile(univariateData, thrQ)] <- NA

  invisible(capture.output(fit <- gpd.fit(univariateData, threshold=quantile(univariateData,thrQ))$mle))
  scale <- fit[1]
  shape <- fit[2]

  empiricalMargins[is.na(empiricalMargins)] <- 1 - pgpd(univariateData[univariateData>quantile(univariateData, thrQ)], quantile(univariateData, thrQ), scale, shape, lower.tail=FALSE)*(1-thrQ)

  laplaceMargins <- rep(NA, length(empiricalMargins))
  laplaceMargins[empiricalMargins <= 0.5] <- log(2*empiricalMargins[empiricalMargins<=0.5])
  laplaceMargins[empiricalMargins > 0.5]  <- -log(2*(1-empiricalMargins[empiricalMargins>0.5]))

  return(laplaceMargins)
}


## Saving the location matrix, S, and choosing and saving the conditioning site, s0
choose_s0 <- function(S, path) {

  # Save the locations as a separate matrix:
  save(file = paste0(path, "/S.rda"), S)

  # Save the conditioning site, s0, here chosen to lie roughly in the centre of
  # the spatial domain:
  s0 <- apply(S, 2, mean)
  s0 <- matrix(s0, nrow = 1)
  # s0 as defined above does not necessarily correspond to an exact location in
  # the data set. Here, we find the location that is as close as possible to the
  # average lon/lat point.
  H <- apply(S, 1, function(s) s - s0)
  H_norm <- apply(H, 2, function(h) sqrt(sum(h^2)))
  s0_idx <- which.min(H_norm)
  s0 <- S[s0_idx, ]
  s0 <- matrix(s0, nrow = 1)
  save(file = paste0(path, "/s0.rda"), s0)

  ## Move s0 to the global environment for access in later functions
  s0_idx <<- s0_idx

  return(s0)
}

## Finding fields where the temperature at s0 is "extreme"
standard_Laplace_quantile <- function(p) if (p < 0.5) log(2 * p) else -log(2 - 2 * p)
find_extreme_fields <- function(data_L, s0_idx, path) {

  # First, define the threshold, u, for an observation to be considered extreme.
  # Here, we set u to be a quantile of the standard Laplace distribution:
  u <- standard_Laplace_quantile(0.95)
  save(file = paste0(path, "/u.rda"), u)

  # Find the indices of the extreme observations:
  s0_data <- data_L[, s0_idx]
  extreme_idx <- which(s0_data > u)

  df <- data.frame(time = time[extreme_idx],
                   Z0   = s0_data[extreme_idx],
                   year = year[extreme_idx]) %>%
    mutate(decade =  round_any(year, 10, f = floor))

  # Subset the fields where Z(s0) is extreme:
  extreme_data_L <- data_L[extreme_idx, ]
  extreme_data_L <- t(extreme_data_L) # transpose so that it's in the correct format
  m_e <- ncol(extreme_data_L)
  cat("Number of spatial fields available with Z(s0) > u:", m_e, "\n")
  save(file = paste0(path, "/extreme_data_subset_LaplaceScale.rda"), extreme_data_L)
  save(file = paste0(path, "/m_e.rda"), m_e)

  # Plot fields in the year 2015 to show temporal dependence
  time_ext <- time[extreme_idx]
  idx <- as.numeric(format(time_ext, "%Y")) == 2015
  dfZ <- as.data.frame(extreme_data_L[, idx])
  names(dfZ) <- time_ext[idx]
  dfZ <- reshape2::melt(dfZ, measure.vars = names(dfZ))
  dfZ <- cbind(dfZ, S)
  dfZ <- dplyr::rename(dfZ, time = variable, Z = value)

  # Keep only the first 24 fields so that we can display in a neat 4 x 6 grid.
  # NB: Could also do 5 x 5, but 4 x 6 saves a bit of space.
  dfZ <- dfZ %>% filter(time %in% levels(time)[1:24])

  # Also transform the longitude and latitude back to the original scale for
  # consistency with the other plots.
  dfZ <- dfZ %>% mutate(lon = lon / 1.04, lat = lat / 1.11)

  g <- RedSea
  if (data_type == "regular") {
    g <- g +
      geom_tile(data = dfZ,
                aes(x = lon, y = lat, fill = Z^(1/3)),
                width = 0.19, height = 0.19) +
      scale_fill_viridis_c(option = "magma", name = "")#expression(Z^{1/3}))
  } else {
    g <- ggplot(dfZ, aes(x = lon, y = lat, colour = Z^(1/3))) +
      geom_point() +
      scale_colour_viridis_c(option = "magma", name = "")#expression(Z^{1/3}))
  }


  g <- g +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.background = element_blank())

  g <- g + facet_wrap(~time, nrow = 6)

  ggsave(g, file = "temporal_spread_2015_fields_long.pdf", device = "pdf",
         width = 8.2, height = 10, path = img_path)

  g <- g + facet_wrap(~time, nrow = 4)

  ggsave(g, file = "temporal_spread_2015_fields_wide.pdf", device = "pdf",
         width = 8.2, height = 5.5, path = img_path)

  # Save the years for use as blocks:
  blocks <- year[extreme_idx] - min(year[extreme_idx]) + 1
  B <- length(unique(blocks))
  blocks <- apply(sapply(1:B, function(i) {
    b <- unique(blocks)[i]
    i * (blocks == b)
  }), 1, sum)
  save(file = paste0(path, "/blocks.rda"), blocks)

  # Define and save the regions used for threshold exceedances:
  H <- apply(S, 1, function(s) s - s0)
  H_norm <- apply(H, 2, function(h) sqrt(sum(h^2)))
  num_regions <- 17
  region <- cut(H_norm, num_regions)
  save(file = paste0(path, "/region.rda"), region)
  region_id <- factor(region, labels = 1:num_regions)
  region_id <- as.numeric(region_id)
  save(file = paste0(path, "/region_id.rda"), region_id)

  return(extreme_data_L)
}

plot_data <- function(
  extreme_data_L, S, img_path, geom = c("point", "tile"),
  n_plots = 5, type = c("random", "most_extreme", "least_extreme")) {

  n <- nrow(extreme_data_L)
  n_fields <- ncol(extreme_data_L)

  if (type == "random") {
    plot_idx <- sample(1:n_fields, n_plots)
  } else if (type == "most_extreme") {
    plot_idx <- (n_fields - n_plots + 1):n_fields
  } else if (type == "least_extreme") {
    plot_idx <- 1:n_plots
  } else {
    stop("argument type not recognised")
  }

  Z <- c(extreme_data_L[, plot_idx])
  group <- rep(1:n_plots, each = n)
  df_L  <- data.frame(
    Z = Z, group = group,
    lon = S[, "lon"], lat = S[, "lat"]
  )

  # Don't want to use a shared colour scale, so we can't use facet_wrap()...
  # Instead, we'll create a list of individual plots and combine them.
  plotlist <- lapply(unique(group), function(g) {

    if (geom == "point") {
      gg <- RedSea +
        geom_point(data = filter(df_L, group == g),
                   aes(x = lon / 1.04, y = lat / 1.11, colour = Z)) +
        scale_colour_viridis_c(option = "magma")
    } else if (geom == "tile") {
      gg <- RedSea +
        geom_tile(data = filter(df_L, group == g),
                  aes(x = lon / 1.04, y = lat / 1.11, fill = Z)) +
        scale_fill_viridis_c(option = "magma")
    }
    gg + labs(x = "Longitude", y = "Latitude", colour = "") + theme_bw()
  })

  figure <- ggarrange(plotlist = plotlist, nrow = 1, legend = "top")

  ggsave(figure,
         file = paste0("extreme_fields_", type,".pdf"), device = "pdf",
         width = 15, height = 4, path = img_path)
}


# ---- Load the full Red Sea data set ----

# See this link for a description of this data set and the objects it contains:
# https://hpc.niasra.uow.edu.au/ckan/dataset/red_sea_temperature
load("data/RedSea/redseatemperature.rdata")

## Extract observations from the Summer (July, August, September):
summer_months <- month %in% c(7,8,9)
data  <- data[summer_months, ]
month <- month[summer_months]
year  <- year[summer_months]
time  <- time[summer_months]

## Extract observations from the southern part of interest:
idx  <- loc[, "lat"] < 20.2 & loc[, "lat"] > 15.75
data <- data[, idx]
loc  <- loc[idx, ]

## Transform locations:
loc[, "lon"] <- loc[, "lon"] * 1.04
loc[, "lat"] <- loc[, "lat"] * 1.11

## Plot the full data set
figure <- ggplot() +
  geom_point(aes(x = loc[, "lon"], y = loc[, "lat"]), size = 0.3) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + coord_fixed()
img_path <- "img/RedSea"
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)
ggsave(figure,
       file = "full_data.pdf", device = "pdf",
       width = 6, height = 4, path = img_path)

# Rename the full data so that we can access it in each of the following sections:
full_data <- data
full_loc  <- loc


# ---- Irregular subset of data ----

data_type <- "irregular"
path <- paste0("data/RedSea/", data_type)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

img_path <- paste0("img/RedSea/", data_type)
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)

## Subset the data: Randomly select 678 locations
set.seed(1)
idx <- sample(1:ncol(full_data), 678)
data <- full_data[, idx]
loc  <- full_loc[idx, ]
cat("Number of observations in the irregular Red Sea data set:", length(idx))

## Transform to Laplace margins:
suppressWarnings(data_L <- apply(data, 2, laplaceMargins))

## Save new datasets:
save(data, loc, file = paste0(path, "/data_subset_OriginalScale.RData"))
save(data_L, loc, file = paste0(path, "/data_subset_LaplaceScale.RData"))

## Choose the conditioning site
S <- loc
s0 <- choose_s0(S, path)

# Save the distance matrix
D <- fields::rdist(S)
save(file = paste0(path, "/D.rda"), D)

## Find the extreme data
extreme_data_L <- find_extreme_fields(data_L, s0_idx, path)

## Plot some of the fields
set.seed(1)
plot_data(extreme_data_L, S, img_path = img_path, type = "random", geom = "point")
plot_data(extreme_data_L, S, img_path = img_path, type = "most_extreme", geom = "point")
plot_data(extreme_data_L, S, img_path = img_path, type = "least_extreme", geom = "point")


# ---- Regular subset of data ----

data_type <- "regular"
path <- paste0("./data/RedSea/", data_type)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

img_path <- paste0("img/RedSea/", data_type)
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)

## Subset the spatial locations to every third longitude and latitude value:
## NB: We start the sequence from 3 because there are a couple of missing
## observations in the spatial domain, and starting at 3 ensures that the
## subsetted grid is fully observed. We get n = 678, which is the same value
## as Simpson et al., so they probably subsetted in this manner too.
lon_values <- sort(unique(full_loc[, "lon"]))
lon_values <- lon_values[seq(3, length(lon_values), by = 3)]
lat_values <- sort(unique(full_loc[, "lat"]))
lat_values <- lat_values[seq(3, length(lat_values), by = 3)]
idx  <- which(full_loc[, "lon"] %in% lon_values & full_loc[, "lat"] %in% lat_values)
data <- full_data[, idx]
loc  <- full_loc[idx, ]
cat("Number of locations in each field in the (regular) Red Sea data set:", length(idx))
# Sanity check: ggplot() + geom_point(aes(x = loc[, "lon"], y = loc[, "lat"]))

## Transform the data to Laplace margins:
suppressWarnings(data_L <- apply(data, 2, laplaceMargins))

## Save new data sets:
save(data, loc, file = paste0(path, "/data_subset_OriginalScale.RData"))
save(data_L, loc, file = paste0(path, "/data_subset_LaplaceScale.RData"))

## Choose the conditioning site
S <- loc
s0 <- choose_s0(S, path)

# Save the distance matrix
D <- fields::rdist(S)
save(file = paste0(path, "/D.rda"), D)

## Find the extreme data
extreme_data_L <- find_extreme_fields(data_L, s0_idx, path)

## Plot some of the fields
set.seed(1)
plot_data(extreme_data_L, S, img_path = img_path, type = "random", geom = "tile")
plot_data(extreme_data_L, S, img_path = img_path, type = "most_extreme", geom = "tile")
plot_data(extreme_data_L, S, img_path = img_path, type = "least_extreme", geom = "tile")

## Padding the spatial domain to make a full rectangular grid
lon_values <- sort(unique(S[, "lon"]))
lat_values <- sort(unique(S[, "lat"]))
## Full grid of longitude-latitude pairs:
## Since Julia is column major, we need latitude to run faster than longitude.
full_grid <- expand.grid(lat = sort(lat_values, decreasing = TRUE), lon = lon_values)
# Sanity check: ggplot(full_grid) + geom_point(aes(x = lon, y = lat))
# Save the width and height of the full grid
height <- length(lat_values)
width <- length(lon_values)
save(file = paste0(path, "/height.rda"), height)
save(file = paste0(path, "/width.rda"), width)

# Linear index of the full grid, respecting the column major ordering of Julia:
# See https://docs.julialang.org/en/v1/manual/arrays/
# for each row, l, in loc, find which row of full_grid is equal to l. Then,
# element i of data_idx, that is, data_idx[i], gives the linear index wrt to full_grid
# for the spatial location, loc[i, ].
# Logical: which row in the matrix M equals the vector v?
whichRowMatch <- function(M, v) which(apply(M, 1, function(row) all(row == v)))
data_idx <- apply(loc[, c("lat", "lon")], 1, function(v) whichRowMatch(full_grid, v))
pad_idx  <- (1:nrow(full_grid))[-data_idx]
save(file = paste0(path, "/data_idx.rda"), data_idx)
save(file = paste0(path, "/pad_idx.rda"), pad_idx)

## Plot the data region and the padded region:
df1 <- full_grid[data_idx, ] %>% mutate(region = "data")
df2 <- full_grid[pad_idx, ]  %>% mutate(region = "padded")
df <- rbind(df1, df2)
figure <- RedSea +
  geom_point(data = df,
             aes(x = lon / 1.04, y = lat / 1.11, colour = region)) +
  labs(x = "Longitude", y = "Latitude", colour = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
ggsave(figure,
       file = "data_and_padding.pdf", device = "pdf",
       width = 6, height = 4, path = img_path)
