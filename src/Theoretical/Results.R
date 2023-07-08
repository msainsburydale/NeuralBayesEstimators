library("optparse")
option_list <- list(
  make_option("--model", type="character", default=NULL,
              help="A relative path to the folder of the assumed model.", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
model      <- parse_args(opt_parser)$model

# ---- Prep ----

source("src/PlottingFunctions.R")

intermediates_path  <- paste0("intermediates/Theoretical/", model, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")
img_path            <- paste0("img/Theoretical/", model)
dir.create(img_path, recursive = TRUE, showWarnings = FALSE)

param_labels <- c("θ"  = expression(theta))
estimator_subset <- c("BayesEstimator", "q = 5_vanilla", "q = 5", "q = 200", "q = 200_expert")
estimator_labels["ML"] <- "ML estimator"


estimator_order <- estimator_subset  


diagnostic_plots <- function(xi_hat, xi, min_m = 1, 
                             estimator_subset = unique(xi_hat$estimator), 
                             width = 6.5, height = 4) {
  
  xi_hat <- xi_hat %>% filter(m >= min_m, estimator %in% estimator_subset)
  
  ## MAD
  ggsave(
    MAD_vs_m(xi, xi_hat, param_labels) %>% rm_facet_labels %>% add_breaks(min_m), 
    width = width, height = height, device = "pdf", path = img_path,
    file = paste0("MAD_vs_m_grouped_by_estimator_", min_m, ".pdf")
  )
  
  ## RMSE
  RMSE <- function(x, y) sqrt(mean((x - y)^2))
  gg <- MAD_vs_m(xi, xi_hat, param_labels, loss = RMSE) + labs(x = expression(m), y = "RMSE") 
  gg <- gg %>% rm_facet_labels %>% add_breaks(min_m)
  ggsave(
    gg, width = width, height = height, device = "pdf", path = img_path,
    file = paste0("RMSE_vs_m_grouped_by_estimator_", min_m, ".pdf")
  )
  
  ## rΩ = MAE (since we used MAE loss during training)
  MAE <- function(x, y) mean(abs(x - y))
  gg <- MAD_vs_m(xi, xi_hat, param_labels, loss = MAE) + labs(x = expression(m), y = expression(r[Omega](hat(theta))))
  gg <- gg %>% rm_facet_labels %>% add_breaks(min_m)
  ggsave(
    gg, width = width, height = height, device = "pdf", path = img_path, 
    file = paste0("MAE_vs_m_grouped_by_estimator_", min_m, ".pdf")
  )
  
  return(gg)
}

# ---- loss vs. m, grouped by estimator ----

# Load in estimates + true parameters
xi_hat <- estimates_path %>% paste0("estimates_test.csv")  %>% read.csv 
xi     <- estimates_path %>% paste0("parameters_test.csv") %>% read.csv 


rm_facet_labels <- function(gg) {
  gg + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
}

add_breaks <- function(gg, brks) {
  gg$scales$scales[[1]]$breaks <- c(brks, gg$scales$scales[[1]]$breaks)
  gg
}

MAE_1  <- diagnostic_plots(xi_hat, xi, 1, estimator_subset)

# ---- Density plot for small and large sample case  ----

densityplot <- function(xi_hat, truth) {
  ggplot(xi_hat) +
    geom_line(aes(x = θ, group = estimator, colour = estimator), stat = "density") +
    scale_colour_estimator(xi_hat) +
    geom_vline(aes(xintercept = truth), colour = "gray", linetype = "dashed") +
    theme_bw() +
    labs(x = expression(theta), colour = "") +
    theme(
      legend.text.align = 0
    )
}


# Load in estimates + true parameters
xi_hat <- estimates_path %>% paste0("estimates_scenarios.csv")  %>% read.csv 
xi     <- estimates_path %>% paste0("parameters_scenarios.csv") %>% read.csv 
truth  <- xi$θ[1]
xi_hat$estimator <- sub("ND", "q = ", xi_hat$estimator)

# all_m <- c(5, max(xi_hat$m))
all_m <- c(10, 120)
xi_hat <- xi_hat %>% filter(m %in% all_m, estimator %in% estimator_subset)

# Rather than facetting, create a list of plots so that we can easily combine 
# with the MAE plot. 
densityplots <- lapply(all_m, function(i) {
  
  xi_hat <- xi_hat %>% filter(m == i)
  
  densityplot(xi_hat, truth)
})

# ---- Combine the MAE and density plots ----

increase_font <- function(gg) {
  gg + theme(
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 13)
    )
}

plotlist <- c(list(MAE_1), densityplots)
plotlist <- lapply(plotlist, increase_font)

# The first element of plotlist provides the legend
plotlist[[1]] <- plotlist[[1]] + 
  theme(
    legend.position = "top",
    legend.key.size = unit(1.5, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(margin = margin(r = 1.2, unit = 'cm'))
    )

figure <- ggarrange(plotlist = plotlist, nrow = 1, common.legend = TRUE)

ggsave(
  figure,
  file = str_interp("loss_and_density.pdf"), 
  width = 12, height = 4, path = img_path, device = "pdf"
)

