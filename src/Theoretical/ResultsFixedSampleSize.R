library("optparse")
option_list <- list(
  make_option("--model", type="character", default=NULL,
              help="A relative path to the folder of the assumed model.", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
model      <- parse_args(opt_parser)$model

source("src/PlottingFunctions.R")

# ---- Prep ----

intermediates_path  <- paste0("intermediates/Theoretical/", model, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")
img_path            <- paste0("img/Theoretical/", model)
dir.create(img_path, recursive = TRUE, showWarnings = FALSE)

param_labels <- c("θ"  = expression(theta))

estimator_labels <- c(
  "BayesEstimator" = "Bayes estimator",
  "OAATEstimator" = "One-at-a-time estimator",
  "NN_m10"  = "Neural Bayes estimator",
  "ML" = "ML estimator",
  "MAPEstimator" = "MAP estimator"
)

estimator_colours <- c(
  "BayesEstimator" = "#21908CFF",
  "OAATEstimator" = "orange",
  "NN_m10" = "red",
  "ML"    = "#440154FF",
  "MAPEstimator"    = "pink"
)

estimator_linetypes <- c(
  "BayesEstimator" = "solid",
  "OAATEstimator" = "solid",
  "NN_m10" = "solid",
  "ML" = "solid",
  "MAPEstimator" = "solid"
)


# ---- Density plot  ----

# Load in estimates + true parameters
xi_hat <- estimates_path %>% paste0("estimates_scenarios.csv")  %>% read.csv
xi     <- estimates_path %>% paste0("parameters_scenarios.csv") %>% read.csv
truth  <- xi$θ[1]
theta <- truth
a <- 4
b <- 1
m <- 10
atilde <- a + m
c  <- 2^(-1/atilde)
d <- b/c
lowerbound <- d
f <- function(x) (d < x & x < theta * 2^(1/atilde)) * (c*m/theta) * (c*x/theta)^(m-1)

x_lower <- seq(0.6, d, length.out = 100)
x_bulk <- seq(d+0.0001, atilde/m, length.out = 1000)
x_upper <- seq(theta/c+0.0001, 1.45, length.out = 100)
analytic_lower <- data.frame(x = x_lower, f = f(x_lower))
analytic_bulk  <- data.frame(x = x_bulk, f = f(x_bulk))
analytic_upper <- data.frame(x = x_upper, f = f(x_upper))

densityplot <- function(xi_hat, estimator_subset) {

  xi_hat <- xi_hat %>% filter(m == 10, estimator %in% estimator_subset)

  figure <- ggplot(xi_hat) +
    geom_line(aes(x = θ, colour = estimator, linetype = estimator),
              stat = "density", adjust = 1.5) +
    # geom_line(data = analytic_lower, aes(x = x, y = f), colour = "darkgray") +
    # geom_line(data = analytic_bulk, aes(x = x, y = f), colour = "darkgray") +
    # geom_line(data = analytic_upper, aes(x = x, y = f), colour = "darkgray") +
    # scale_colour_estimator(xi_hat) +
    scale_estimator(xi_hat, scale = "colour", values = estimator_colours) +
    scale_estimator(xi_hat, scale = "linetype", values = estimator_linetypes) +
    geom_vline(aes(xintercept = truth), colour = "gray", linetype = "dashed") +
    scale_x_continuous(
      limits = c(0.6, 1.45),
      breaks = c(0.6, 1.0, truth, 1.4),
      labels = c(0.6, 1.0, expression(theta), 1.4)
    ) +
    theme_bw() +
    labs(x = expression(hat(theta)), y = "density") +
    theme(
      legend.text.align = 0,
      legend.title = element_blank(),
      panel.grid = element_blank()
    )

  # if Zi ~ Unif(0, theta), then Y = max(Zi) has a known distribution; see https://math.stackexchange.com/a/3287934
  # MLE_density <- function(y, theta, m) ifelse(y <= theta, m * y^(m-1) / theta^m, 0)
  # figure <- figure +
  #   stat_function(
  #     fun = MLE_density, args = c(truth, 10), n = 1000,
  #     colour = estimator_colours["ML"],
  #     linetype = "dotted", alpha = 0.5
  #     )

  return(figure)
}

estimator_subset <- c("BayesEstimator", "NN_m10", "OAATEstimator", "ML", "MAPEstimator")
estimator_order  <- estimator_subset

ggsave(
  densityplot(xi_hat, estimator_subset),
  file = "density_all.pdf",
  width = 7, height = 3.5, path = img_path, device = "pdf"
)

estimator_subset <- c("BayesEstimator", "NN_m10", "OAATEstimator", "ML")
estimator_order  <- estimator_subset

ggsave(
  densityplot(xi_hat, estimator_subset),
  file = "density_noMAP.pdf",
  width = 7, height = 3.5, path = img_path, device = "pdf"
)


# ---- Checking with histogram ----

# estimator_subset <- c("BayesEstimator", "NN_m10")
# estimator_order  <- estimator_subset
#
# densityplot(xi_hat, estimator_subset)
#
# ggplot(xi_hat %>% filter(estimator %in% estimator_subset)) +
#   geom_histogram(aes(x = θ, colour = estimator), bins = 100, alpha=0.2, position="identity") +
#   scale_estimator(xi_hat, scale = "colour", values = estimator_colours) +
#   theme_bw()
