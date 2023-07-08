source("src/PlottingFunctions.R")
library("forcats") # fct_rev

model <- "GaussianProcess/nuFixed"
intermediates_path  <- paste0("intermediates/", model, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")
img_path            <- paste0("img/", model)

# NB: Order of this list gives the order of the facets in the risk plots
param_labels <- c(
  "σ"  = expression(sigma[epsilon]),
  "ρ"  = expression(rho)
)

## Subset of estimators
estimators <- c("N_m1",  "ND_piecewise", "ML")


# ---- Variable sample size ----

# Load in estimates + true parameters
df     <- estimates_path %>% paste0("merged_test.csv") %>% read.csv
likelihood_path <- paste0(estimates_path, "merged_likelihood_test.csv")
if (file.exists(likelihood_path)) df <- rbind(df, read.csv(likelihood_path))

all_loss <- c("MAE", "RMSE", "MSE", "MAD", "zeroone")

risk_plots <- lapply(all_loss, function(loss) {

  gg <- df %>%
    filter(estimator %in% estimators) %>%
    diagnosticplot(param_labels, loss = eval(parse(text = loss))) +
    theme(
      strip.text = element_text(size = 17),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 15)
    )

  ggsave(
    gg, width = 8, height = 3.5, device = "pdf", path = img_path,
    file = paste0(loss, "_vs_m.pdf"),
  )

  gg
})

names(risk_plots) <- all_loss
risk_plot <- risk_plots[["MSE"]]

# ---- Joint distribution ----

# Load in estimates + true parameters
df     <- estimates_path %>% paste0("merged_scenarios.csv") %>% read.csv
likelihood_path <- paste0(estimates_path, "merged_likelihood_scenarios.csv")
if (file.exists(likelihood_path)) df <- rbind(df, read.csv(likelihood_path))

all_m  <- unique(df$m)

x <- lapply(all_m, function(i) {

  # filter the estimates to only a subset of m and estimators
  df <- df %>% filter(m == i, estimator %in% estimators)

  figure <- scatterplotlong(df, param_labels) + facet_wrap(fct_rev(as.factor(k)) ~ ., scales = "free")

  ggsave(figure,
         file = str_interp("Scatterplot_m${i}.pdf"), device = "pdf",
         width = 6.7, height = 5, path = img_path)
})


# ---- Single joint distribution ----

estimates_plot <- df %>%
  filter(m == 150, k == 1, estimator %in% estimators) %>%
  scatterplotlong(param_labels)


# ---- Combined risk and scatter plot ----

# To align the top border of the panels, we use a blank title and change its size appropriately
blank_title   <- function(gg) gg + labs(title = " ") + theme(plot.title = element_text(size = 18))

# Change elements of the plots to make it easier to read
increase_font <- function(gg) gg + theme(axis.title = element_text(size = 17), axis.text = element_text(size = 15))
strip_size <- 17

risk_plot <- (risk_plot + theme(strip.text = element_text(size = strip_size))) %>% increase_font
estimates_plot <- estimates_plot %>% blank_title %>% increase_font

ggsave(
  ggarrange(risk_plot, estimates_plot,
            common.legend = TRUE, legend = "right", widths = c(2, 1)),
  file = str_interp("risk_and_scatterplot.pdf"), device = "pdf",
  width = 12, height = 4, path = img_path
)

