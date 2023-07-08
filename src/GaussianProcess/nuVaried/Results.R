source("src/PlottingFunctions.R")

model <- "GaussianProcess/nuVaried"
intermediates_path  <- paste0("intermediates/", model, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")
img_path            <- paste0("img/", model)

# NB: Order of this list gives the order of the facets in the MAD plot
param_labels <- c(
  "σ"  = expression(sigma[epsilon]),
  "ρ"  = expression(rho),
  "ν"  = expression(nu)
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
    gg, width = 12, height = 4, device = "pdf", path = img_path,
    file = paste0(loss, "_vs_m.pdf"),
  )
})

# ---- Joint distribution ----

# Load in estimates + true parameters
df     <- estimates_path %>% paste0("merged_scenarios.csv") %>% read.csv
likelihood_path <- paste0(estimates_path, "merged_likelihood_scenarios.csv")
if (file.exists(likelihood_path)) df <- rbind(df, read.csv(likelihood_path))

# Load realisations from the model
fields_scenarios <- paste0(intermediates_path, "fields.csv") %>% read.csv

# filter the estimates to only a subset of estimators
all_m  <- unique(df$m)
df <- df %>% filter(estimator %in% estimators, m %in% all_m)

# Single scenario:
plotlist <- scatterplotlong(df %>% filter(m == 150, k == df$k[1]), param_labels)
plotlist <- lapply(plotlist, function(gg) gg + theme(legend.text.align = 0))
figure   <- ggarrange(plotlist = plotlist, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(figure, file = "Scatterplot_m150_singleScenario.pdf",
       width = 9, height = 2.8, path = img_path, device = "pdf")

# Several scenarios and with a field realisation:
plot_scenario <- function(df, field) {

  fieldplot    <- field_plot(field) + theme(legend.position = "top")
  scatterplots <- scatterplotlong(df, param_labels)

  legend.grob <<- get_legend(scatterplots[[1]]) # cheeky super-assignment

  plotlist <- c(list(fieldplot), scatterplots)

  return(plotlist)
}

plotlist <- lapply(all_m, function(j) {

  plotlist <- lapply(unique(df$k), function(scen) {

    df     <- df %>% filter(k == scen, m == j)
    field  <- fields_scenarios %>% filter(scenario == scen, replicate == 1)

    ggarrange(plotlist = plot_scenario(df, field), legend = "none", nrow = 1, align = "hv")
  })

  fig <- ggarrange(plotlist = plotlist, legend.grob = legend.grob, ncol = 1)

  ggsave(fig, file = str_interp("Scatterplot_m${j}.pdf"),
         width = 8.3, height = 1.75 * length(unique(df$k)), path = img_path, device = "pdf")

  return(fig)
})
