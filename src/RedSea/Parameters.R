# NB: This is similar to ConditionalExtremes/ParameterConfigurations.R,
# except for the fact that we use an irregular grid of spatial locations defined
# by the Red Sea data set.

quick <- identical(commandArgs(trailingOnly=TRUE)[1], "--quick")

source("./src/GaussianProcess/nuFixed/ParameterFunctions.R", chdir = TRUE)

suppressMessages({
library("dplyr")
library("abind")
library("ggplot2")
library("viridis")
library("KScorrect")
})

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


# ---- LHS design for the parameter space ----

# NB: theta corresponds to rho in the manuscript, the range parameter.
# Here, lambda corresponds to the parameter in the function a(.), NOT the
# measurement error variance, which is denoted by lambda in the other .R scripts.
# The measurement error variance is sigma2e, as per the manuscript.
# Confusingly, the functions write Lambda as the measurement error variance.
# Note also that these parameter ranges are given on the original scale.

lower <- c(
  # parameters associated with a(.) and b(.)
  kappa  = 1,
  lambda = 1,
  beta   = 0.05,
  # Covariance parameters associated with the Gaussian process
  rho = 0.2,
  nu  = 0.4,
  # Location and scale parameters for the Subbotin distribution
  mu  = -0.5,
  tau = 0.2,
  # Parameters used to construct the shape parameter delta
  delta1 = 1
)

upper <- c(
  # parameters associated with a(.) and b(.)
  kappa  = 2,
  lambda = 5,
  beta   = 1,
  # Covariance parameters associated with the Gaussian process
  rho    = 5,
  nu     = 2.0,
  # Location and scale parameters for the Subbotin distribution
  mu     = 0.5,
  tau    = 1,
  # Parameters used to construct the shape parameter delta
  delta1 = 3
)

ranges <- data.frame(lower = lower, upper = upper)

## Number of unique parameter configurations to generate:
if (quick) {
  n <- c(train = 500, val = 60, test = 60, scenarios = 6)
} else {
  n <- c(train = 10000, val = 2000, test = 1000, scenarios = 6)
}

all_sets <- names(n)

set.seed(1)
# param_grid <- lapply(all_sets, function(set) lhs_design(n[set], t(ranges)) )
param_grid <- lapply(all_sets, function(set) {
  sapply(names(upper), function(x) runif(n[set], min = lower[x], max = upper[x]))
})
names(param_grid) <- all_sets


generate_parameter_configurations <- function(type) {

  cat("\nGenerating parameter configurations for the", type, "Red Sea data set...\n")

  # Directories
  data_path <- paste0("./data/RedSea/", type)
  path      <- paste0("./intermediates/RedSea/", type,"/parameter_configurations/")
  img_path  <- paste0("img/", type, "/RedSea")
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  dir.create(img_path, showWarnings = FALSE, recursive = TRUE)

  # ---- Generate Cholesky Factors ----

  load(file = paste0(data_path, "/S.rda"))

  params <- lapply(all_sets, function(set) {
    prepare_wrapper(
      S = S,
      path = paste0(path, set),
      nuGrid = param_grid[[set]][, "nu"],
      thetaGrid = param_grid[[set]][, "rho"],
      completeGrid = FALSE, saveLambda = FALSE, lambdaGrid = 0
    )
  })

  # ---- Save the other parameters ----

  names(param_grid) <- all_sets

  xi <- lapply(all_sets, function(set) {

    xi <- param_grid[[set]]
    save(file = paste0(path, set, "_xi.rda"), xi)

    xi
  })

  return(xi)
}

xi <- lapply(c("irregular", "regular"), generate_parameter_configurations)
