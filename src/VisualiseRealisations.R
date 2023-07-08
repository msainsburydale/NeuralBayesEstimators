library("optparse")

option_list = list(
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="Model for this application.", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

model  <- opt$model
intermediates_path <- paste0("./intermediates/", model)
img_path <- paste0("./img/", model)
dir.create(img_path, recursive = TRUE, showWarnings = FALSE)

suppressMessages({
  library("ggplot2")
  library("dplyr")
  library("tibble")
  library("viridis") # viridis colour scale
  library("reshape2")
  library("ggpubr")
  library("gridExtra")
  library("stringr") # str_interp for string interpolation
  options(dplyr.summarise.inform = FALSE) # Suppress summarise info
  source("src/PlottingFunctions.R")
})


## load the parameters
if (model == "GaussianProcess/nuFixed") {
  params_path <- paste0(intermediates_path, "/parameter_configurations/")
  load(paste0(params_path, "scenarios_xi.rda"))
  theta <- xi # deprecation coercion
} else {
  theta <- read.csv(paste0(intermediates_path, "/scenario_parameters.csv"))
}

## Visualise the fields
fields <- paste0(intermediates_path, "/fields.csv") %>% read.csv

# Use only some of the fields for this plot
fields <- fields %>% filter(replicate <= 4)

# Some models produce fields that are highly varied. Hence, we use a separate
# colour scale for each field... This means that we can't use facet_wrap().
fieldplots <- lapply(unique(fields$scenario), function(j) {

  fields <- fields %>% filter(scenario == j)

  fieldplots <- lapply(unique(fields$replicate), function(i) {

    fields <- fields %>% filter(replicate == i)

    return(field_plot(fields))
  })

  fieldplots <- ggarrange(plotlist = fieldplots, nrow = 1, ncol = 4)

  # Draw the text that describes the parameters 
  text <- theta[j, ] %>%
    round(2) %>%
    paste(names(.), ., sep = " = ", collapse = ", ")

  # Unicode doesn't work
  text <- text %>%
    str_replace("σ", "sigma[epsilon]") %>%
    str_replace("ρ", "rho") %>%
    str_replace("ν", "nu") %>%
    str_replace("κ", "kappa") %>%
    str_replace("λ", "lambda") %>%
    str_replace("β", "beta") %>%
    str_replace("μ", "mu") %>%
    str_replace("τ", "tau" ) %>%
    str_replace("δ1", "delta1")

  # A bit more wrangling so that parse works (see ?plotmath)
  text <- paste0("list(", text, ")")
  text <- gsub("=", "==", text)

  textplot <- ggplot() +
    theme_void() +
    geom_text(aes(x = 0, y = 0, label = text), size = 10, parse = TRUE)

  # Combine the title plot and the fieldplots
  return(ggarrange(textplot, fieldplots, ncol = 1, heights = c(1,5)))
})

figure <- ggarrange(plotlist = fieldplots, ncol = 1)

suppressWarnings(
  ggsave(figure, width = 12, height = 3 * nrow(theta),
         file = "fields.pdf", device = "pdf", path = img_path)

)

