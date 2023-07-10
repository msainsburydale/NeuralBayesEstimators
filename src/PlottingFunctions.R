suppressMessages({
  library("ggplot2")
  library("dplyr")
  library("tibble")
  library("viridis")
  library("reshape2")
  library("ggpubr")
  library("stringr")  # str_interp for string interpolation
  library("gridExtra")
  library("scales")   # to access break formatting functions
  library("purrr")    # map()
  library("combinat") # combn()
  library("tidyr")    # separate()
})

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# ---- Estimator labels and colours ----

estimator_colours <- c(
  # The "optimal" likelihood-based estimator used for comparison
  "ML"    = "#440154FF",
  "PL3"    = "#440154FF",
  # Neural estimators
  "N_m1"    = "#21908CFF",
  "ND_m30"  = "orange",
  "ND_m150"  = "green",
  "ND_m1to150"  = "pink",
  "ND_m1to150_expert"  = "#FDE725FF",
  "ND_m75"  = "#FDE725FF",
  "ND_m150" = "orange",
  "ND_piecewise" = "orange",
  # pairwise likelihood estimators used to check optimality of the PL3 estimator
  "PL2"    = "#21908CFF",
  "PL6" = "#FDE725FF",
  "PL9" = "red",
  "PLF" = "blue",
  # Estimators used in theoretical study
  "BayesEstimator" = "#21908CFF",
  "q = 5"    = "#FDE725FF",
  "q = 5_vanilla"    = "pink",
  "q = 5_vanilla_m5"    = "#FDE725FF",
  "q = 5_vanilla_m150"    = "orange",
  "q = 5_vanilla_m1to30"    = "red",
  "q = 5_vanilla_m1to150"    = "#440154FF",
  "q = 200"  = "orange",
  "q = 200_expert"  = "red",
  "q = 5_m5" = "#FDE725FF",
  "q = 5_m150" = "orange",
  "q = 5_m1to30" = "red",
  "NN_m10" = "red"
)

# NB Graphical devices are often OS-variable; bquote(bold(theta)) does not
# make theta bold, for some reason, but some computers don't like unicode, so 
# I'm just going with bquote(bold(theta)).
# boldtheta <- bquote(bold("\U03B8"))
boldtheta <- bquote(bold(theta))

estimator_labels <- c(
  "N_m1" = bquote(hat(.(boldtheta))[0]("·")),
  "ND_m1" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m==1),
  "ND_m10" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m==10),
  "ND_m30" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m==30),
  "ND60" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m==60),
  "ND_m75" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m==75),
  "ND_m130" = bquote(hat(.(boldtheta))[DS]("·")),
  "ND_m150" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m==150),
  "ND1to30" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m%in%1:30),
  "ND_m1to150" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m%in%1:150),
  "ND_m1to150_expert" = bquote(hat(.(boldtheta))[DS]("·") * " : " * m%in%1:150 * ", S(Z) = |Z|"),
  "ND_piecewise" = bquote(hat(.(boldtheta))[DS]("·")),
  "PL2" = expression(PL[2]),
  "PL3" =  bquote(hat(.(boldtheta))[PMAP]("·")),
  "PL6" = expression(PL[6]),
  "PL9" = expression(PL[9]),
  "PLF" = expression(PL[infinity]),
  "ML" = bquote(hat(.(boldtheta))[MAP]("·")),
  "BayesEstimator" = "Bayes estimator",
  "q = 5" = expression(hat(theta)[gamma*"*"]("·") * " : " * q[t] * " = 5"),
  "q = 5_vanilla" = expression("Vanilla Deep Set with " * q[t] == 5),
  "q = 200" = expression(hat(theta)[gamma*"*"]("·") * " : " * q[t] * " = 200"),
  "q = 200_expert" = expression(hat(theta)[gamma*"*"]("·") * " : " * q[t] * " = 200, max(·) in S(·)"),
  "q = 5_m5" = expression(hat(theta)[m==5]("·")),
  "q = 5_m150" = expression(hat(theta)[m==150]("·")),
  "q = 5_m1to30" = expression(hat(theta)[m%in%1:30]("·")),
  "q = 5_vanilla_m5" = expression(hat(theta)[m==5]("·")),
  "q = 5_vanilla_m150" = expression(hat(theta)[m==150]("·")),
  "q = 5_vanilla_m1to30" = expression(hat(theta)[m%in%1:30]("·")),
  "NN_m10" = expression(hat(theta)[gamma*"*"]("·")),
  "NN_m100" = expression(hat(theta)[gamma*"*"]("·"))
)

# Specifies the order that the estimators should appear in the plot legends.
estimator_order <- names(estimator_labels)

# Custom colour scale to enforce consistency for the estimators across the project
scale_colour_estimator <- function(df, ...) {
  estimators <- unique(df$estimator)
  ggplot2:::manual_scale(
    'colour',
    values = estimator_colours[estimators],
    labels = estimator_labels,
    breaks = estimator_order,
    ...
  )
}

scale_estimator <- function(df, scale = "colour", values = estimator_colours, ...) {
  estimators <- unique(df$estimator)
  ggplot2:::manual_scale(
    scale,
    values = values[estimators],
    labels = estimator_labels,
    breaks = estimator_order,
    ...
  )
}

# ---- Variable sample size plots ----

MAE <- function(x, y) mean(abs(x - y))
MSE <- function(x, y) mean((x - y)^2)
RMSE <- function(x, y) sqrt(mean((x - y)^2))
MAD <- mad
zeroone <- function(x, y, eps = y/10) mean(abs(x - y) > eps)

diagnosticplot <- function(df, param_labels, loss = MAE) {
  
  df <- df %>% mutate(residual = estimate - truth)
  
  all_m <- unique(df$m)
  breaks <- all_m[which(all_m %in% c(1, 30, 60, 90, 120, 150))]
  
  param_labeller <- label_parsed
  df <- mutate_at(df, .vars = "parameter", .funs = factor,  levels = names(param_labels), labels = param_labels)
  
  # Compute global risk for each combination of estimator and sample size m
  df <- df %>%
    group_by(estimator, parameter, m) %>%
    summarize(loss = loss(estimate, truth))
  
  # Plot risk vs. m
  gg <- ggplot(data = df, aes(x = m, y = loss, colour = estimator, group = estimator)) +
    geom_point() +
    geom_line() +
    facet_wrap(parameter ~ ., scales = "free", labeller = param_labeller) +
    labs(colour = "", x = expression(m), y = expression(r[Omega](hat(theta)("·")))) +
    scale_x_continuous(breaks = breaks) +
    scale_colour_estimator(df) +
    theme_bw() +
    theme(legend.text=element_text(size = 14),
          legend.text.align = 0,
          panel.grid = element_blank(),
          strip.background = element_blank())
  
  return(gg)
}


# Changes the limits of a ggplot legend without changing anything else of the scale
change_legend_limits <- function(gg, limits, aesthetic = "fill") {
  
  ## Find the scales associated with the specifed aesthetic
  sc <- as.list(gg$scales)$scales
  all_aesthetics <- sapply(sc, function(x) x[["aesthetics"]][1])
  idx <- which(aesthetic == all_aesthetics)
  
  ## Overwrite the breaks of the specifed aesthetic
  gg$scales$scales[[idx]][["limits"]] <- limits
  
  return(gg)
}


# ---- Joint distribution of estimators ----

scatterplotlong <- function(df, param_labels) {
  
  n_params    <- length(param_labels)
  param_names <- unique(df$parameter)
  
  if (n_params != length(param_names)) stop("The number of parameter labels differs to the number of parameters: Please ensure length(unique(df$parameter)) == length(param_labels)")
  if (!all(param_names %in% names(param_labels))) stop("Some parameters have not been given parameter labels: Please ensure all(unique(df$parameter) %in% names(param_labels))")
  
  # convert to wide form
  df <- df %>% pivot_wider(
    names_from = parameter, values_from = c("estimate", "truth")
  ) %>%
    as.data.frame
  
  combinations <- param_names %>% combn(2) %>% as.matrix
  
  # Generate the scatterplot estimation panels
  scatterplots <- apply(combinations, 2, function(p) {
    ggplot(data = df[sample(nrow(df)), ]) +
      geom_point(
        aes_string(
          x = paste("estimate", p[1], sep = "_"),
          y = paste("estimate", p[2], sep = "_"),
          colour = "estimator"
        ),
        alpha = 0.5, size = 1) +
      geom_point(
        aes_string(
          x = paste("truth", p[1], sep = "_"),
          y = paste("truth", p[2], sep = "_")
        ),
        shape = "+", size = 8, colour = "red"
      ) +
      labs(colour = "", x = param_labels[[p[1]]], y = param_labels[[p[2]]]) +
      scale_colour_estimator(df) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      scale_x_continuous(n.breaks = 4) +
      theme_bw() +
      theme(legend.position = "top")
  })
  
  if (length(scatterplots) == 1) {
    scatterplots <- scatterplots[[1]] +
      theme(
        legend.position = "right",
        legend.text.align = 0,
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  return(scatterplots)
}

# ---- Joint distribution of estimators as m varies (applicable only for two-parameter models) ----

# If we choose to use a common colour scale for the fields, the overall
# alignment is better since we're able to use facet_wrap.
distribution_vs_m <- function(df, fields, scenarios, param_labels, common_field_colour_scale = FALSE) {
  
  all_m  <- unique(df$m)
  
  df <- filter(df, k %in% scenarios)
  df$scenario <- df$k
  fields <- fields %>% filter(scenario %in% scenarios, replicate == 1)
  
  scatterplots <- lapply(all_m, function(i) {
    
    # filter the estimates to only a subset of m and estimators
    xi_hat <- xi_hat %>% filter(m == i, estimator %in% estimators)
    
    # Permute the rows to prevent one estimator dominating another
    xi_hat <- xi_hat %>% .[sample(nrow(.)), ]
    
    # NB: Can't use facet_grid() each panel needs free x and y scales
    figure <- scatterplotlong(df, param_labels) +
      facet_wrap(scenario ~ ., scales = "free")  +
      theme(legend.position = "top")
    
    legend.grob <<- get_legend(figure) # cheeky super-assignment
    
    if (i < max(all_m)) figure <- figure + xlab("")
    
    figure
  })
  
  if (!common_field_colour_scale) {
    
    fieldplot <- fields %>%
      group_split(scenario) %>%
      map(
        ~ggplot(.) +
          geom_tile(aes(x = x, y = y, fill = Z)) +
          labs(fill = "", x = " ", y = " ") +
          scale_fill_viridis_c(option = "magma") +
          facet_wrap(scenario ~ ., scales = "free") +
          theme_bw() +
          theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks   =  element_blank(),
            axis.text    = element_blank(),
            legend.position = "top",
            strip.background = element_blank(),
            strip.text.x = element_blank()
          )
      ) %>%
      ggarrange(plotlist = ., align = 'hv', nrow = 1)
    
    
    scatterplots <- ggarrange(plotlist = scatterplots, nrow = length(all_m),
                              legend.grob = legend.grob, align = "v")
    figure <- ggarrange(fieldplot, scatterplots, nrow = 2, legend = "none",
                        heights = c(1, 2))
    
  } else {
    
    fieldplot <- ggplot(fields) +
      geom_tile(aes(x = x, y = y, fill = Z)) +
      facet_wrap(scenario ~ ., scales = "free") +
      labs(fill = "", x = " ", y = " ") +
      scale_fill_viridis_c(option = "magma") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks   =  element_blank(),
        axis.text    = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
    
    plotlist <- c(list(fieldplot), scatterplots)
    figure <- ggarrange(plotlist = plotlist, nrow = length(plotlist),
                        legend.grob = legend.grob, align = "v")
    
  }
  
  return(figure)
}


# ---- Simulated fields ----


# The simulations may be highly varied in magnitude, so we need to
# use an independent colour scale. This means that we can't use facet_wrap().
field_plot <- function(field, regular = TRUE) {
  
  gg <- ggplot(field, aes(x = x, y = y))
  
  if (regular) {
    gg <- gg +
      geom_tile(aes(fill = Z)) +
      scale_fill_viridis_c(option = "magma")
  } else {
    gg <- gg +
      geom_point(aes(colour = Z)) +
      scale_colour_viridis_c(option = "magma")
  }
  
  gg <- gg +
    labs(fill = "", x = " ", y = " ") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks   =  element_blank(),
      axis.text    = element_blank(),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  return(gg)
}

plot_fields <- function(fields) {
  
  if(max(fields$replicate) < 8) stop("plot_fields() needs at least 8 field replicates for each parameter configuration")
  
  lapply(1:8, function(r) fields %>% filter(replicate == r) %>% field_plot) %>%
    ggarrange(plotlist = ., nrow = 2, ncol = 4)
  
}

# ---- Unused code that may be of use in the future ----

# cbrt <- function(x) sign(x) * abs(x)^(1/3)
# symmetriclog <- function(x) sign(x) * log(abs(x))

# Adds brackets if x starts with ln or log:
# "lna"  %>% add_brackets
# "loga" %>% add_brackets
# "la"   %>% add_brackets
# c("lna", "la") %>% add_brackets
# c("lna", "la") %>% sapply(add_brackets)
# add_brackets <- function(x) {
#   if (grepl("ln", x)) {
#     paste0(sub("ln", "ln(", x), ")")
#   } else if(grepl("log", x)) {
#     paste0(sub("log", "log(", x), ")")
#   } else {
#     x
#   }
# }

# .param_labeller <- as_labeller(function(x) sapply(x, add_brackets))

# ---- Removing points ----

# I thought I needed this for the likelihood estimators, since some of their
# estimates were outliers (vanilla ML or PL estimators) or right on the boundary
# of the support (MAP or PMAP estimators). However, I didn't end up needing
# this, but I'll keep it for future reference.

# Ω = list(
#   σ = c(0.1, 1),
#   ρ = c(2.0, 10.0),
#   ν = c(0.5, 2.0)
# )

# removeoutliers <- function(df, Ω, scaling_factor = 2) {
#
#   # Expand the support by scaling_factor
#   if (scaling_factor < 1) stop("scaling_factor should be greater than 1")
#   scaling_factor = c(1/scaling_factor, scaling_factor)
#   Ω = lapply(Ω, function(x) x * scaling_factor)
#
#   # First, add a two-column matrix to the data frame; the ith row of this matrix
#   # specifies the bounds for the parameter in the ith row of df.
#   bounds <- do.call(rbind, Ω[df$parameter])
#   bounds <- as.data.frame(bounds)
#   names(bounds) <- c("lower", "upper")
#   rownames(bounds) <- NULL
#   df <- cbind(df, bounds)
#
#   # Now find which estimates are outside of the bounds (i.e., which are outliers):
#   outliers <- which(with(df, lower > estimate | estimate > upper))
#   cat("Number of likelihood-based estimates that are outliers:", length(outliers), "\n")
#
#   # Find and remove parameter configurations associated with these outliers.
#   # We remove entire pairs of parameter configuration and sample size even if
#   # only one of the parameters is considered an outlier. This is because the
#   # failure of one parameter would likely compromise the estimation of the
#   # others.
#   bad_mk <- df[outliers, c("m", "k")]
#   df <- anti_join(df, bad_mk, by = c("m", "k"))
#
#   # remove lower and upper
#   df$lower <- NULL
#   df$upper <- NULL
#
#   df
# }
#

# scaling_factor = list("σ" = c(2, 1))

# removeboundarypoints <- function(df, Ω, scaling_factor = c(2, 1)) {
#
#   # # Shrink the support by scaling_factor
#   # if (length(scaling_factor) > 1) stop("length(scaling_factor) must equal 1")
#   # if (!(names(scaling_factor) %in% names(Ω))) stop("the element in scaling_factor should be named, and that name should appear in Ω")
#   # Ω[[names(scaling_factor)]] <- Ω[names(scaling_factor)][[1]] * scaling_factor[[1]]
#   #
#   # # First, add a two-column matrix to the data frame; the ith row of this matrix
#   # # specifies the bounds for the parameter in the ith row of df.
#   # bounds <- do.call(rbind, Ω[df$parameter])
#   # bounds <- as.data.frame(bounds)
#   # names(bounds) <- c("lower", "upper")
#   # rownames(bounds) <- NULL
#   # df <- cbind(df, bounds)
#
#   # Shrink the support by scaling_factor
#   if (scaling_factor[1] < 1) stop("scaling_factor[1] should be >= 1")
#   if (scaling_factor[2] > 1) stop("scaling_factor[2] should be <= 1")
#   Ω = lapply(Ω, function(x) x * scaling_factor)
#
#   # First, add a two-column matrix to the data frame; the ith row of this matrix
#   # specifies the bounds for the parameter in the ith row of df.
#   bounds <- do.call(rbind, Ω[df$parameter])
#   bounds <- as.data.frame(bounds)
#   names(bounds) <- c("lower", "upper")
#   rownames(bounds) <- NULL
#   df <- cbind(df, bounds)
#
#   # Now find which parameters are outside the support:
#   outside_idx <- which(with(df, truth < lower | truth > upper))
#
#   # Find and remove parameter configurations associated with these outliers.
#   # We remove entire pairs of parameter configuration and sample size even if
#   # only one of the parameters is considered an outlier. This is because the
#   # failure of one parameter would likely compromise the estimation of the
#   # others.
#   bad_mk <- df[outside_idx, c("m", "k")]
#   df <- anti_join(df, bad_mk, by = c("m", "k"))
#
#   # remove lower and upper
#   df$lower <- NULL
#   df$upper <- NULL
#
#   df
# }



# ---- Parallel coordinates ----

#' Parallel coordinate plot with individual y-axis for each variable
#'
#' Similar to `GGally::ggparcoord`, but creates an individual y-axis for each variable.
#'
#' @param truth a `data.frame` containing points to plot for each variable
#'
#' @examples
#' ggparcoord_ind_yaxis(iris, columns = 4:1)
#' ggparcoord_ind_yaxis(iris, columns = 4:1, groupColumn = "Species", alphaLines = 0.5)
#'
#' # Add median of each variable:
#' truth <- iris %>% select(4:1) %>% apply(2, median, simplify = FALSE) %>% data.frame
#' ggparcoord_ind_yaxis(iris, truth = truth, columns = 4:1, groupColumn = "Species", alphaLines = 0.5)
# ggparcoord_ind_yaxis <- function(
#   data,
#   columns,
#   facetColumn = NULL, # I don't actually use this argument in this project
#   groupColumn = NULL,
#   truth = NULL,
#   truthPointSize = 2,
#   alphaLines = 1,
#   nbreaks = 4,
#   axis_font_size = 3
# ) {
# 
# 
#   if (is.null(facetColumn)) {
#     facetColumn <- "panel"
#     data$panel <- "dummy"
#     if (!is.null(truth)) {
#       truth$panel <- "dummy"
#     }
# 
#   } else {
#     # Check that facetColumn is in both the data and truth data frames
#     if (!(facetColumn %in% names(data))) {
#       stop("facetColumn must be in data")
#     }
# 
#     if (!is.null(truth) && !(facetColumn %in% names(truth))) {
#       stop("facetColumn must be in truth dataframe")
#     }
# 
#     # Rename facet column to "panel" and group column to "network"
#     data <- rename(data, panel = facetColumn, network = groupColumn)
#     truth <- rename(truth, panel = facetColumn)
#   }
# 
#   # Rename group column to "network"
#   if (is.null(groupColumn)) {
#     data$network = "dummy"
#   } else {
#     data <- rename(data, panel = facetColumn, network = groupColumn)
#   }
# 
# 
#   # select the variables to plot and facet with
#   data_subset <- data %>% select(columns, facetColumn)
# 
#   # re-order truth to match columns
#   col_names <- data_subset %>% names
#   if (!is.null(truth)) {
#     truth <- truth %>% select(col_names)
#     data_subset <- data_subset %>% rbind(truth)
#   }
# 
#   # Calculate the axis breaks for each variable on the *original* scale.
#   # Note that the breaks computed by pretty() are guaranteed to contain all of
#   # the data. We include truth in these breaks, just in case one of the true
#   # points falls outside the range of the data (can easily happen in the context
#   # of comparing parameter estimates to the true values).
#   breaks_df <- data_subset %>%
#     reshape2::melt(id.vars = "panel", variable.name = "ind", value.name = "values") %>%
#     group_by(ind, panel) %>%   # group by the plotting variables AND the facet variable
#     summarize(breaks = pretty(values, n = nbreaks))
# 
#   # Normalise the breaks to be between 0 and 1, and set the coordinates of the
#   # tick marks. Importantly, if we want the axis heights to be the same, the
#   # breaks need to be normalised to be between exactly 0 and 1.
#   axis_df <- breaks_df %>%
#     mutate(yval = (breaks - min(breaks))/(max(breaks) - min(breaks))) %>%
#     mutate(xmin = as.numeric(ind) - 0.05,
#            xmax = as.numeric(ind),
#            x_text = as.numeric(ind) - 0.2)
# 
#   # Calculate the co-ordinates for our axis lines:
#   axis_line_df <- axis_df %>%
#     group_by(ind, panel) %>%
#     summarize(min = min(yval), max = max(yval))
# 
#   # Getting the minimum/maximum breaks on the original scale, to scale the
#   # data in the same manner that we scaled the breaks
#   minmax_breaks <- breaks_df %>%
#     summarize(min_break = min(breaks), max_break = max(breaks)) %>%
#     group_by(ind, panel)
# 
# 
#   # Normalise the original data in the same way that the breaks were normalised.
#   # This ensures that the scaling is correct.
#   # Do the same for the truth points, if they exist.
#   lines_df <- data %>%
#     select(columns, panel, network) %>%
#     mutate(row = row_number()) %>% # need row information to group individual rows in the plots
#     reshape2::melt(id.vars = c("panel", "network", "row"),
#                    variable.name = "ind", value.name = "values") %>%
#     group_by(ind, panel)
# 
#   if (!is.null(truth)) {
#     truth <- truth %>%
#       select(columns, panel) %>%
#       reshape2::melt(id.vars = "panel", variable.name = "ind", value.name = "values") %>%
#       group_by(ind, panel)
#   }
# 
#   groupRows <- lines_df %>% group_rows()
#   for (group in 1:nrow(minmax_breaks)) { # For every group
#     idx <- groupRows[[group]]
#     minb = minmax_breaks[group, "min_break"] %>% as.numeric
#     maxb = minmax_breaks[group, "max_break"] %>% as.numeric
#     lines_df[idx, "values"] <- (lines_df[idx, "values"] - minb) / (maxb -  minb)
#     if (!is.null(truth)) {
#       truth[group, "values"] <- (truth[group, "values"] - minb) / (maxb -  minb)
#     }
#   }
# 
# 
#   # Now plot:
#   gg <- ggplot() +
#     geom_line(data = lines_df %>% data.frame,
#               aes_string(x = "ind", y = "values", group = "row", colour = "network"),
#               alpha = alphaLines) +
#     geom_segment(data = axis_line_df, aes(x = ind, xend = ind, y = min, yend = max),
#                  inherit.aes = FALSE) +
#     geom_segment(data = axis_df, aes(x = xmin, xend = xmax, y = yval, yend = yval),
#                  inherit.aes = FALSE) +
#     geom_text(data = axis_df, aes(x = x_text, y = yval, label = breaks),
#               inherit.aes = FALSE, size = axis_font_size)
# 
#   if (!is.null(truth)) {
#     gg <- gg + geom_point(data = truth, aes(x = ind, y = values),
#                           inherit.aes = FALSE, colour = "red", size = truthPointSize)
#   }
# 
#   gg <- gg + facet_wrap(panel ~ ., labeller = param_labeller)
# 
#   gg <- gg + labs(colour = "")
# 
#   gg <- gg + theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       panel.border = element_blank(),
#       axis.title = element_blank(),
#       axis.ticks =  element_blank(),
#       axis.text.y = element_blank(),
#       strip.background = element_blank(),
#       strip.text.x = element_blank()
#     )
# 
#   # Remove the legend if there's only one group
#   n_network <- lines_df %>% data.frame() %>% select(network) %>% unique %>% length
#   if (n_network == 1) {
#     gg <- gg + theme(legend.position = "none")
#   }
# 
#   return(gg)
# }



# ---- time vs m ----

# time_vs_m <- function(df) {
#   
#   all_m <- unique(df$m)
#   breaks <- all_m[which(all_m %in% c(1, 30, 60, 90, 120, 150))]
#   
#   ggplot(df, aes(x = m, y = time, colour = estimator)) +
#     geom_line() +
#     geom_point() +
#     scale_colour_estimator(df) +
#     labs(colour = "", x = expression(m[e]), y = "time (s)") +
#     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                   labels = trans_format("log10", math_format(10^.x))) +
#     scale_x_continuous(breaks = breaks) +
#     theme_bw() +
#     theme(legend.text=element_text(size = 14),
#           legend.text.align = 0,
#           panel.grid = element_blank(),
#           strip.background = element_blank())
# }