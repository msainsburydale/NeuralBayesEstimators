source("src/PlottingFunctions.R")

model <- "ConditionalExtremes"
intermediates_path  <- paste0("intermediates/", model, "/")
estimates_path      <- paste0(intermediates_path, "Estimates/")
img_path            <- paste0("img/", model)

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

## Subset of estimators
estimators <- c("N_m1",  "ND_piecewise", "ML")



# ---- Variable sample size ----

# Load in estimates + true parameters
df     <- estimates_path %>% paste0("merged_test.csv") %>% read.csv
likelihood_path <- paste0(estimates_path, "merged_likelihood_test.csv")
if (file.exists(likelihood_path)) df <- rbind(df, read.csv(likelihood_path))
df$parameter <- gsub("δ₁", "δ1", df$parameter)

x <- lapply(c("MAE", "RMSE", "MSE", "MAD", "zeroone"), function(loss) {

  gg <- df %>%
    filter(estimator %in% estimators) %>%
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

# Load realisations from the model
fields <- paste0(intermediates_path, "fields.csv") %>% read.csv

# filter the estimates to only a subset of estimators
all_m  <- df$m %>% unique
df <- df %>% filter(estimator %in% estimators, m %in% all_m)

# Sanity check:
if (!all(names(param_labels) %in% unique(df$parameter)) ) {
  stop("Not all parameters have been given labels")
}

n_params <- length(param_labels)
fields <- fields %>% filter(replicate <= choose(n_params, 2))
num_scenarios <- length(unique(df$k))

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

# ---- Marginal distributions ----

densityplot <- function(df, param_labels, truth_colour = "red") {

  param_labeller <- label_parsed
  df <- mutate_at(df, .vars = "parameter", .funs = factor, levels = names(param_labels), labels = param_labels)

  gg <- ggplot(df) +
    geom_line(aes(x = estimate, group = estimator, colour = estimator), stat = "density") +
    geom_vline(aes(xintercept = truth), colour = truth_colour, linetype = "dashed") +
    facet_wrap(parameter ~ ., scales = "free", labeller = param_labeller) +
    labs(colour = "") +
    scale_colour_estimator(df) +
    theme_bw() +
    theme(legend.text.align = 0,
          panel.grid = element_blank(),
          strip.background = element_blank())

  return(gg)
}


plotlist <- lapply(all_m, function(j) {
  plotlist <- lapply(1:num_scenarios, function(i) {

    figure <- df %>% filter(k == i, m == j) %>% densityplot(param_labels = param_labels)

    figure$facet$params$nrow <- 2

    # Change elements of the plots to make it easier to read
    figure <- figure +
      theme(strip.text = element_text(size = 17),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 17),
            legend.key.height = unit(1, "cm"))

    ggsave(figure,
           file = str_interp("Densityplot_m${j}_scenario${i}.pdf"), device = "pdf",
           width = 12, height = 6, path = img_path)
  })
})
