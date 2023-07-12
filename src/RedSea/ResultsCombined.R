data_type <- "regular"
data_path <- paste0("data/RedSea/", data_type, "/")
img_path  <- paste0("img/RedSea/", data_type, "/")
results_path  <- paste0("results/RedSea/", data_type, "/")
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

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

# Base map for the Red Sea
# dt <- data.frame(lon = c(37.1, 37.1,  43,  43),
#                  lat = c(15.5, 20.4, 20.4, 15.5))
# RedSea <- basemap(data = dt)
RedSea <- ggplot() # no base map

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

# Unfortunately, bold() doesn't work for greek letters.. leave for now.
threshold_labels <- c(
  "Empirical"  = "Empirical",
  "CNN" = expression(paste("Model: ", widehat(bold(theta))["CNN"])),
  "DNN" = expression(paste("Model: ", widehat(bold(theta))["DNN"])),
  "ML" = expression(paste("Model: ", widehat(bold(theta))["ML"]))
)
estimator_labels <- c(
  "CNN" = expression(widehat(bold(theta))["CNN"]),
  "DNN" = expression(widehat(bold(theta))["DNN"])
)

alpha <- 0.1 # significance level for confidence intervals

renamedelta <- function(df) {
  colnames(df) <- gsub("δ.", "δ1", colnames(df))
  df
}


# ---- Estimates and bootstrap confidence intervals ----

# Combined table of estimates and confidence intervals from DNN, CNN, and eventually ML

caption <- "Parameter estimates and confidence intervals for the Red Sea data set."

theta_hat_DNN <- read.csv("intermediates/RedSea/DNNregular/Estimates/real_data_estimates.csv") %>% renamedelta %>% round(2)
theta_hat_CNN <- read.csv("intermediates/RedSea/CNN/Estimates/real_data_estimates.csv") %>% renamedelta %>% round(2)
rownames(theta_hat_CNN) <- rownames(theta_hat_DNN) <- "Estimate"

theta_tilde_DNN <- read.csv("intermediates/RedSea/DNNregular/Estimates/bootstrap_samples_nonparametric.csv")  %>% renamedelta
theta_tilde_CNN <- read.csv("intermediates/RedSea/CNN/Estimates/bootstrap_samples_nonparametric.csv")  %>% renamedelta

CI_DNN <- apply(theta_tilde_DNN, 2, function(x) quantile(x, c(alpha/2, 1 - alpha/2))) %>% round(2)
CI_CNN <- apply(theta_tilde_CNN, 2, function(x) quantile(x, c(alpha/2, 1 - alpha/2))) %>% round(2)

DNN <- rbind(theta_hat_DNN, CI_DNN)
DNN <- apply(DNN, 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
CNN <- rbind(theta_hat_CNN, CI_CNN)
CNN <- apply(CNN, 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))

combined <- rbind(CNN, DNN)

# Replace unicode with tex version for direct input to tex file (only need to add \ in the tex file)
# colnames(combined) <- colnames(combined) %>%
#   str_replace("κ", "$kappa$") %>%
#   str_replace("λ", "$lambda$") %>%
#   str_replace( "β", "$beta$") %>%
#   str_replace( "ρ", "$rho$") %>%
#   str_replace("ν", "$nu$") %>%
#   str_replace("μ", "$mu$") %>%
#   str_replace("τ", "$tau$") %>%
#   str_replace("δ1", "$delta_1$")

write.csv(combined, paste0(results_path, "estimates_and_confidence_intervals_wide.csv"))
xtable(combined, type = "latex",
       caption = caption) %>%
  print(file = paste0(results_path, "estimates_and_confidence_intervals_wide.tex"))

# The wide format ended up being way too wide for to fit in the manuscript..
# instead, we'll have to have estimates and confidence intervals on separate rows

write.csv(combined, paste0(results_path, "estimates_and_confidence_intervals_long.csv"))
xtable(t(combined), type = "latex",
       caption = caption) %>%
  print(file = paste0(results_path, "estimates_and_confidence_intervals_long.tex"))

# ---- Variable sample size ----

df <- lapply(c("CNN", "DNNregular", "DNNirregular"), function(x) {
  str_interp("intermediates/RedSea/${x}/Estimates/merged_test.csv") %>%
    read.csv  %>%
    filter(estimator == "ND_m150") %>%
    mutate(estimator = x)
})
df <- bind_rows(df)
df$parameter <- gsub("δ₁", "δ1", df$parameter)


estimator_labels <- c(
  "CNN" = expression(widehat(bold(theta))["CNN"]),
  "DNNregular" = expression(widehat(bold(theta))["DNN: regular"]),
  "DNNirregular" = expression(widehat(bold(theta))["DNN: irregular"])
)

estimator_colours <- c(
  "CNN" = "orange",
  "DNNregular" = "#440154FF",
  "DNNirregular" = "#21908CFF"
)

x <- lapply(c("MAE", "RMSE", "MSE", "MAD", "zeroone"), function(loss) {

  # Omit very small sample sizes
  df <- filter(df, m >= 10)

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
    gg, width = 12, height = 6, device = "pdf", path = "img/RedSea/",
    file = paste0(loss, "_vs_m.pdf"),
  )
})



# ---- Variable sample size (regular data only) ----

df <- lapply(c("CNN", "DNNregular"), function(x) {
  str_interp("intermediates/RedSea/${x}/Estimates/merged_test.csv") %>%
    read.csv  %>%
    filter(estimator == "ND_m150") %>%
    mutate(estimator = x)
})
df <- bind_rows(df)
df$parameter <- gsub("δ₁", "δ1", df$parameter)


estimator_labels <- c(
  "CNN" = expression(widehat(bold(theta))["CNN"]),
  "DNNregular" = expression(widehat(bold(theta))["DNN: regular"])
)

estimator_colours <- c(
  "CNN" = "orange",
  "DNNregular" = "#440154FF"
)

x <- lapply(c("MAE", "RMSE", "MSE", "MAD", "zeroone"), function(loss) {

  # Omit very small sample sizes
  df <- filter(df, m >= 10)

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
    gg, width = 12, height = 6, device = "pdf", path = "img/RedSea/",
    file = paste0(loss, "_vs_m_regularOnly.pdf"),
  )
})
