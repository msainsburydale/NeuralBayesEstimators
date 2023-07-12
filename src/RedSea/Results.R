library("optparse")

option_list = list(
  make_option(c("-d", "--data_type"), type="character", default=NULL,
              help="regular or irregular data?", metavar="character"),
  make_option(c("-n", "--arch"), type="character", default=NULL,
              help="Fully-connected ('DNN') or Convolutional Neural Network ('CNN').", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)
data_type  <- opt$data_type
arch       <- opt$arch

suppressMessages({
library("ggplot2")
library("dplyr")
library("tibble")
library("viridis")  # viridis colour scale
library("reshape2") # melt()
library("ggpubr")
library("gridExtra")
library("xtable")
library("stringr") # str_interp for string interpolation
# library("ggOceanMaps")
options(dplyr.summarise.inform = FALSE) # Suppress summarise info
})
source("src/PlottingFunctions.R")

# If we're using a fully-connected arch, the data could be regular or irregular,
# and the folder names reflect this.
if(arch == "DNN") arch = paste0(arch, data_type)

intermediates_path  <- paste0("intermediates/RedSea/", arch, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")
img_path            <- paste0("img/RedSea/", arch)
data_path           <- paste0("data/RedSea/", data_type, "/")
results_path        <- paste0("results/RedSea/", arch, "/")
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
dir.create(img_path,     showWarnings = FALSE, recursive = TRUE)

# NB: Order of this vector gives the order of the facets in the MAD plot
param_labels <- c(
  "κ"  = expression(kappa),
  "λ"  = expression(lambda),
  "β"    = expression(beta),
  "ρ"  = expression(rho),
  "ν"  = expression(nu),
  "μ"    = expression(mu),
  "τ"  = expression(tau),
  "δ1" = expression(delta[1])
)

# NB: This is a list; definition above is a vector. I need a list for one plot below.
param_labels_list <- list(
  "κ"  = expression(kappa),
  "λ"  = expression(lambda),
  "β"    = expression(beta),
  "ρ"  = expression(rho),
  "ν"  = expression(nu),
  "μ"    = expression(mu),
  "τ"  = expression(tau),
  "δ1" = expression(delta[1])
)

estimator_colours <- c("ND_m150"  = "orange")
estimator_labels <- c("ND_m150" =  bquote(hat(.(boldtheta))[DS]("·")))


# Base map for the Red Sea
# dt <- data.frame(lon = c(36.9, 36.9,  43.1,  43.1),
#                  lat = c(15.5, 20.4, 20.4, 15.5))
# RedSea <- basemap(data = dt)
RedSea <- ggplot() # no base map

alpha <- 0.05 # significance level for confidence intervals

renamedelta <- function(df) {
  colnames(df) <- gsub("δ.", "δ1", colnames(df))
  df
}


# ---- Estimates and bootstrap confidence intervals ----

caption <- "Parameter estimates and confidence intervals for the Red Sea data set."

theta_hat <- estimates_path %>% paste0("real_data_estimates.csv") %>% read.csv %>% renamedelta
rownames(theta_hat) <- "Estimate"

theta_tilde <- estimates_path %>% paste0("bootstrap_samples_nonparametric.csv") %>% read.csv %>% renamedelta
CI <- apply(theta_tilde, 2, function(x) quantile(x, c(alpha/2, 1 - alpha/2)))
write.csv(CI, paste0(results_path, "confidence_interval_nonparametric.csv"))
CI_time <- estimates_path %>% paste0("bootstrap_time_nonparametric.csv") %>% read.csv(header = FALSE)
theta_hat %>%
  rbind(CI) %>%
  xtable(type = "latex",
         caption = caption) %>%
  print(file = paste0(results_path, "confidence_interval_nonparametric.tex"))

# ---- Variable sample size ----

# Load in estimates + true parameters
df     <- estimates_path %>% paste0("merged_test.csv") %>% read.csv
likelihood_path <- paste0(estimates_path, "merged_likelihood_test.csv")
if (file.exists(likelihood_path)) df <- rbind(df, read.csv(likelihood_path))
df$parameter <- gsub("δ₁", "δ1", df$parameter)

x <- lapply(c("MAE", "RMSE", "MSE", "MAD", "zeroone"), function(loss) {

  gg <- df %>%
    diagnosticplot(param_labels, loss = eval(parse(text = loss))) +
    theme(
      strip.text = element_text(size = 17),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 15)
    )

  gg$facet$params$nrow <- 2

  # Change elements of the plots to make it easier to read
  gg <- gg +
    theme(strip.text = element_text(size = 17),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 17),
          legend.key.height = unit(1, "cm"))

  ggsave(
    gg, width = 12, height = 6, device = "pdf", path = img_path,
    file = paste0(loss, "_vs_m.pdf"),
  )
})


# ---- Joint distribution of the estimators ----

# Load in estimates + true parameters
df     <- estimates_path %>% paste0("merged_scenarios.csv") %>% read.csv
likelihood_path <- paste0(estimates_path, "merged_likelihood_scenarios.csv")
if (file.exists(likelihood_path)) df <- rbind(df, read.csv(likelihood_path))
df$parameter <- gsub("δ₁", "δ1", df$parameter)

# Load model realisations
fields <- paste0(estimates_path, "fields_scenarios.csv") %>% read.csv
fields$Z <- fields$Z^3 # FIXME better to do this in Julia

# Sanity check:
if (!all(names(param_labels) %in% unique(df$parameter)) ) {
  stop("Not all parameters have been given labels")
}

n_params <- length(param_labels)
fields   <- fields %>% filter(replicate <= choose(n_params, 2))
num_scenarios <- length(unique(df$k))
all_m         <- unique(df$m)
plotlist <- lapply(all_m, function(j) {
  plotlist <- lapply(1:num_scenarios, function(i) {

    scatterplots <- scatterplotlong(
      df %>% filter(k == i, m == j),
      param_labels
    )

    # Modify legend and remove the axis labels
    scatterplots <- lapply(scatterplots, function(gg) {
      gg +
        guides(colour = guide_legend(override.aes = list(size = 4))) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "left",
          legend.text.align = 0,
          legend.text = element_text(size = 15),
          legend.key.height = unit(1, "cm")
        )
    })

    # Extract legend so that it can be placed in the final plot
    scatterplot_legend.grob <<- get_legend(scatterplots[[1]])
    scatterplots <- lapply(scatterplots, function(gg) gg + theme(legend.position = "none"))

    param_names <- names(param_labels)
    layout      <- matrix(rep(NA, n_params^2), nrow = n_params)

    # Lower-diagonal part of the plot: Parameter estimates
    layout[lower.tri(layout)] <- 1:choose(n_params, 2)
    plotlist <- scatterplots


    # Diagonal part of the plot (parameter names)
    diag(layout) <- (choose(n_params, 2) + 1):(choose(n_params, 2) + n_params)
    diag_plots <- lapply(param_labels_list, function(p) {
      ggplot() +
        annotate("text", x = 0, y = 0, label = p, size = 8) +
        theme_void()
    })


    plotlist <- c(plotlist, diag_plots)

    # Upper-diagonal part of the plot: Field realisations
    layout[upper.tri(layout)] <- (choose(n_params, 2) + n_params + 1):n_params^2
    fields <- fields %>% filter(scenario == i)
    # scalelims <- range(fields$Z)
    scalelims <- quantile(fields$Z, c(0.01, 0.99))
    fields <- fields %>% mutate(Z = pmin(pmax(Z, scalelims[1]), scalelims[2]))

    fieldplots <- lapply(unique(fields$replicate), function(r) {

      field <- fields %>% filter(replicate == r)

      gg <- field_plot(field) %>% change_legend_limits(scalelims)

      return(gg)
    })

    # Extract the legend (which is common to all ).
    fieldplots_legend.grob <<- get_legend(fieldplots[[1]])
    fieldplots <- lapply(fieldplots, function(gg) gg + theme(legend.position = "none"))

    plotlist  <- c(plotlist, fieldplots)

    # Add a column/row to the layout that will store the scatterplot legend
    legend_part <- rep(NA, n_params)
    legend_part[floor(n_params / 2)] <- n_params^2 + 1
    if (n_params %% 2 == 0) legend_part[floor(n_params / 2) + 1] <- n_params^2 + 1
    layout <- cbind(legend_part, layout)
    plotlist <- c(plotlist, list(scatterplot_legend.grob))

    # Add a column/row to the layout that will store the data legend
    legend_part <- rep(NA, n_params)
    legend_part[floor(n_params / 2)] <- n_params^2 + 2
    if (n_params %% 2 == 0) legend_part[floor(n_params / 2) + 1] <- n_params^2 + 2
    layout <- cbind(layout, legend_part)
    plotlist <- c(plotlist, list(fieldplots_legend.grob))

    suppressWarnings(
      figure <- grid.arrange(grobs = plotlist, layout_matrix = layout, widths = c(0.75, rep(1, n_params), 0.75))
    )

    ggsave(figure,
           file = str_interp("Scatterplot_m${j}_scenario${i}.pdf"), device = "pdf",
           width = 14, height = 12, path = img_path)
  })
})


# ---- Observed fields vs. Simulations from the fitted model ----

load(paste0(data_path, "S.rda"))
n_plots <- 4

# Simulations from the fitted model:
fields   <- paste0(estimates_path, "fitted_model_simulations.csv") %>% read.csv
fields   <- rename(fields, lon = x, lat = y)
fields   <- fields %>% filter(replicate <= n_plots)
fields$scenario <- NULL
fields$Z <- fields$Z^3 # FIXME better to do this in Julia

# Observed fields:
load(paste0(data_path, "extreme_data_subset_LaplaceScale.rda"))
n <- nrow(extreme_data_L)
n_fields <- ncol(extreme_data_L)
plot_idx <- sample(1:n_fields, n_plots)
Z     <- c(extreme_data_L[, plot_idx])
replicate <- rep(1:n_plots, each = n)
df_L  <- data.frame(
  Z = Z, replicate = replicate,
  lon = S[, 1], lat = S[, 2]
)

# Combine into one data frame for plotting:
fields$type <- "simulated"
df_L$type   <- "observed"
df <- rbind(fields, df_L)

scalelims <- quantile(df$Z, c(0.01, 0.99))
if (data_type == "regular") {
  figure <- RedSea +
    geom_tile(
      data = df,
      aes(x = lon / 1.04, y = lat / 1.11, fill = pmin(pmax(Z, scalelims[1]), scalelims[2])),
      width = 0.163, height = 0.163
      ) +
    scale_fill_viridis_c(option = "magma")
} else {
  figure <- RedSea +
    geom_point(
      data = df,
      aes(x = lon / 1.04, y = lat / 1.11, colour = pmin(pmax(Z, scalelims[1]), scalelims[2]))
      ) +
    scale_colour_viridis_c(option = "magma")
}



figure <- figure +
  facet_grid(type ~ replicate) +
  labs(x = "Longitude", y = "Latitude", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text = element_text(size = 17),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    )

ggsave(figure,
       file = "observed_fields_and_fitted_model_simulations.pdf", device = "pdf",
       width = 15, height = 6, path = img_path)




# ---- Threshold exceedances: Plotting the regions ----

# Here we plot the empirical vs fitted number of exceedances in certain
# sub-regions at different distances from the conditioning site.

# This diagnostic plot considers the spatial distance between each observation
# location and s0; hence, we need the matrix of spatial locations, S, and the
# location of the conditioning site, s0.
load(paste0(data_path, "S.rda"))
load(paste0(data_path, "s0.rda"))
s0_df <- data.frame(lon = s0[1, 1], lat = s0[1, 2])

# Load the regions and find which index corresponds to s0:
load(paste0(data_path, "region.rda"))
num_regions <- length(levels(region))
cat("Red Sea threshold exceedances - Number of points within each region:", table(region), "\n")

# Plot the regions:
region_df <- as.data.frame(cbind(S, region = region))
region_df <- mutate(region_df, colour = region %% 3) # alternating colour achieved via mod

plot_df <- region_df %>%
  select(lon, lat, colour) %>%
  rbind(cbind(s0_df, colour = "whatever"))

if (data_type == "regular") {
  p1 <- RedSea +
    geom_tile(data = plot_df,
              aes(x = lon / 1.04, y = lat / 1.11, fill = as.factor(colour)),
              width = 0.16, height = 0.16, colour = "white") +
    scale_fill_manual(values = c("azure2", "gray", "black", "red"))
} else {
  p1 <- RedSea +
    geom_point(data = plot_df,
               aes(x = lon / 1.04, y = lat / 1.11, colour = as.factor(colour)),
               size = 3) +
    scale_colour_manual(values = c("azure2", "gray", "black", "red"))
}

p1 <- p1 +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


# ---- Threshold exceedances ----

empirical_exceedances_bootstrap <- paste0("intermediates/RedSea/", data_type, "/empirical_exceedances_bootstrap.csv") %>% read.csv
model_exceedances_bootstrap <- paste0(estimates_path, "threshold_bootstrap_samples_nonparametric.csv") %>% read.csv

# Compute summary statistics (mean, quantile) of threshold exceedance proportions
summarise_prop <- function(df) {
  df %>%
    group_by(region) %>%
    summarise(
      lower_bound = quantile(prop_exceedances, alpha/2),
      upper_bound = quantile(prop_exceedances, 1 - alpha/2),
      prop_exceedances = mean(prop_exceedances)
    )
}
empirical_exceedances <- empirical_exceedances_bootstrap %>% summarise_prop
model_exceedances <- model_exceedances_bootstrap %>% summarise_prop

# Plot the model-based and empirical exceedances:
empirical_exceedances$estimator <- "Empirical"
model_exceedances$estimator     <- "Model"
df <- rbind(empirical_exceedances, model_exceedances)
df$region_id <- factor(df$region, labels = 1:num_regions)

jitter_width <- 0.3
p2 <- ggplot(df,
             aes(x = region_id, y = prop_exceedances, colour = estimator,
                 ymin = lower_bound, ymax = upper_bound)) +
  # geom_errorbar(size = 0.5, width = 0.25) +
  # geom_point(size = 2) +
  geom_linerange(alpha = 0.8, size = 0.5, position = position_dodge(width = jitter_width)) +
  geom_point(alpha = 0.8, size = 1, position = position_dodge(width = jitter_width)) +
  labs(y = expression(paste("Proportion of threshold exceedances given ", Z[0] > u)),
       x = "Region", colour = "") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.text = element_text(size = 15))


figure <- ggarrange(p1, p2)

ggsave(figure,
       file = paste0("threshold_exceedances.pdf"), device = "pdf",
       width = 12.5, height = 4.5, path = img_path)
