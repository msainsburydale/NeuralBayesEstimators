library("optparse")
option_list <- list(
  make_option("--model", type="character", default=NULL,
              help="A relative path to the folder of the assumed model.", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
model      <- parse_args(opt_parser)$model

img_path  <- paste0("img/Theoretical/", model)
dir.create(img_path, recursive = TRUE, showWarnings = FALSE)

source("src/PlottingFunctions.R")

# ---- Estimator labelling and colours ----

param_labels <- c("θ"= expression(theta))

estimator_subset <- c(
  "BayesEstimator",
  "NN_m5",
  "NN_m150",
  "NN_m1to150",
  "NN_m1to150expert"
)

estimator_labels <- c(
  "BayesEstimator" = "Bayes estimator",
  "NN_m5" = expression(hat(theta)("·"~";"~gamma[5])),
  "NN_m150" = expression(hat(theta)("·"~";"~gamma[150])),
  "NN_m1to150" = expression(hat(theta)("·"~";"~gamma[1:150])),
  "NN_m1to150expert" = expression(hat(theta)[1:150]("·") * ", " * S("·") == "|Z|")
)

estimator_colours <- c(
  "BayesEstimator" = "#21908CFF",
  "NN_m5"    = "red",
  "NN_m150"    = "orange",
  "NN_m1to150"    = "#440154FF",
  "NN_m1to150expert"    = "green"
)

estimator_order <- estimator_subset

# ---- Prep ----

all_m <- c("5", "150", "1to150")

intermediates_path  <- paste0("intermediates/Theoretical/", model, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")

# Load in estimates + true parameters
xi_hat <- estimates_path %>% paste0("estimates_test.csv")  %>% read.csv
xi     <- estimates_path %>% paste0("parameters_test.csv") %>% read.csv

# NN estimator titles
# NN_estimators <- paste0("NN_", all_m)

# Separate NN information into different columns. This is done by the underscore,
# so first merge _expert into _m.
# xi_hat$estimator <- sub("m1to150_expert", "m1to150expert", xi_hat$estimator)
# suppressWarnings(
#   xi_hat <- xi_hat %>%
#     separate(estimator, into = c('estimator', 'agg', 'q', 'trainingm'), sep = "_", fill = "right")
# )

# Some more wrangling:
# xi_hat$estimator <- paste0(xi_hat$estimator, xi_hat$trainingm)
# xi_hat$trainingm <- NULL
# xi_hat$estimator <- sub("NN", "", xi_hat$estimator)
# xi_hat$estimator <- sub("NA", "", xi_hat$estimator)
xi_hat <- rename(xi_hat, "estimate" = "θ")
xi <- rename(xi, "truth" = "θ")
df <- cbind(xi_hat, xi)

# ---- average risk plot ----

df <- df %>% filter(estimator %in% estimator_subset)

# Compute average risk, r_Ω, for each combination of interest. Here,
# the average risk is the mean absolute error.
MAE <- function(x, y) mean(abs(x - y))
df <- df %>%
  mutate(residual = estimate - truth) %>%
  group_by(estimator, m) %>%
  summarize(risk = MAE(estimate, truth))

breaks <- unique(df$m)

figure <- ggplot(data = df,
                 aes(x = m, y = risk, colour = estimator, group = estimator)) +
  geom_point() +
  geom_line(alpha = 0.75) +
  labs(colour = "", x = expression(m), y = expression(r[Omega](hat(theta)("· ; ·")))) +
  scale_x_continuous(breaks = breaks) +
  scale_colour_estimator(df) +
  theme_bw() +
  theme(
    legend.text=element_text(size = 14),
    legend.text.align = 0,
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Zoom in on the larger sample sizes.
xmin=110; xmax = 150
ymin=0.15; ymax = 0.25

window <- figure +
  theme(
    # axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(position = "right") +
  coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))


figure <- figure +
  geom_rect(aes(xmin=xmin, xmax=xmax + 2, ymin=ymin, ymax=ymax),
            size = 0.3, colour = "grey50", fill = "transparent")

# Extract the legend and convert it to a ggplot object so that it can be added
# to figure in a custom position
legend_plot <- figure %>% get_legend %>% as_ggplot

# Add some padding around window
window <- window + theme(plot.margin = unit(c(10, 10, 20, 10), "points"))

g <- ggarrange(
  figure,
  ggarrange(legend_plot, window, ncol = 1, legend = "none"),
  legend = "none"
)

ggsave(
  g,
  file = str_interp("risk_variableSetSize_meanPoolOnly.pdf"),
  width = 8, height = 4, path = img_path, device = "pdf"
)
