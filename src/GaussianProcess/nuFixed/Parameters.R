library("optparse")
option_list = list(
  make_option(c("-q", "--quick"), action="store_true", default = FALSE,
              help = "Generate small number of parameter", type = "logical")
)
opt_parser <- OptionParser(option_list = option_list)
opt   <- parse_args(opt_parser)
quick <- opt$quick


model    <- file.path("GaussianProcess", "nuFixed")
int_path <- file.path("intermediates", model, "parameter_configurations")
img_path <- file.path("img", model)
dir.create(img_path, recursive = TRUE, showWarnings = FALSE)

source(file.path("src", model, "ParameterFunctions.R"), chdir = TRUE)

all_sets <- c("train", "val", "test")

rangeTheta <- rbind(train = c(min = 0.2, max = 50),
                    val   = c(min = 1.8, max = 27),
                    test  = c(min = 2.0, max = 25))

dfRange <- rbind(train = c(min = 1, max = 255),
                 val   = c(min = 35, max = 224),
                 test  = c(min = 40, max = 216))


## Number of combinations of theta and nu to use, and the number of values of
## lambda to use with each theta-nu pair:
if (quick) {
  nTheta <- c(train = 100, val = 55, test = 55)
  nLambda  <- c(train = 10, val = 5, test = 5)
} else {
  nTheta <- c(train = 200, val = 40, test = 30)
  nLambda  <- c(train = 201, val = 50, test = 30)
}

params <- lapply(all_sets,  function(set, ...) {
  prepare_wrapper(
    path = file.path(int_path, set),
    nuGrid = 1,
    nLambda = nLambda[set],
    nTheta = nTheta[set],
    rangeTheta = rangeTheta[set, ],
    dfRange = dfRange[set, ],
    ...
  )
})
names(params) <- all_sets


# ---- Scenario parameters ----

# Here, we select a subset of tests parameters that will be used for the joint
# distribution plot. We use a complete grid of parameter values to facilitate
# a facetted plot.
params[["scenarios"]] <- prepare_wrapper(
  path = file.path(int_path, "scenarios"),
  nLambda = 3,  nTheta = 3, nuGrid = 1,
  dfRange = dfRange["test", ],
  rangeTheta = rangeTheta["test", ],
  completeGrid = TRUE
)

