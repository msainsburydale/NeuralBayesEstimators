library("optparse")
option_list <- list(
  make_option("--path", type="character", default=NULL,
              help="Relative path for this run.", metavar="character"),
  make_option("--excluded_m", type="integer", default=NULL,
              help="Networks that we wish to exclude characterised by the number of replicates available during training, m.")
)
opt_parser <- OptionParser(option_list=option_list)
path       <- parse_args(opt_parser)$path
path       <- gsub("/", .Platform$file.sep, path)
excluded_m <- parse_args(opt_parser)$excluded_m

intermediates_path <- file.path("intermediates", path)
img_path <- file.path("img", path)
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)

source(file.path("src", "PlottingFunctions.R"))

# Find the runs folders:
all_dirs <- list.dirs(path = intermediates_path, recursive = TRUE)
runs_dirs <- all_dirs[which(grepl("runs_", all_dirs))]
if (!is.null(excluded_m)) runs_dirs <- runs_dirs[!grepl(excluded_m, runs_dirs)]

# Load the loss function per epoch files:
loss_per_epoch_list <- lapply(runs_dirs, function(x) {
  loss_per_epoch <- read.csv(file.path(x, "loss_per_epoch.csv"), header = FALSE)
  colnames(loss_per_epoch) <- c("training", "validation")
  loss_per_epoch$epoch <- 0:(nrow(loss_per_epoch) - 1)
  return(loss_per_epoch)
})

# Extract the title of each network:
network <- sub(".*runs_", "", runs_dirs)

# Extract the number of replicates used during training for each estimator:
m <- regmatches(network, gregexpr("[[:digit:]]+", network)) %>% as.numeric

# Sanity check:
if (!all(network %in% names(estimator_labels))) stop("Not all neural estimators have been given a label")

# Combine the loss matrices into a single data frame:
df <- do.call("rbind", loss_per_epoch_list)
df$network <- rep(network, sapply(loss_per_epoch_list, nrow))
df$m       <- rep(m, sapply(loss_per_epoch_list, nrow))
df <- df %>%
  melt(c("epoch", "network", "m"), variable.name = "set", value.name = "loss") %>%
  mutate_at(.vars = "network", .funs = factor, labels = estimator_labels[network])

# Create y limits using the minimum loss over all sets for a given m, so that
# the panels are directly comparable when appropriate.
df <- df %>%
  group_by(m) %>%
  mutate(ymin = min(loss), ymax = max(loss))

# Compute the minimum validation loss for a given m.
min_df <- df %>% filter(set == "validation") %>% summarise(min_val_loss = min(loss))
min_val_loss <- setNames(min_df$min_val_loss, min_df$m)
df$min_val_loss <- min_val_loss[as.character(df$m)]

# Plot the loss functions:
figure <- ggplot(df) +
  geom_line(aes(x = epoch, y = loss, colour = set)) +
  scale_color_manual(values = c("blue", "red")) +
  facet_wrap(~network, scales = "free", labeller = label_parsed, nrow = 1) +
  labs(colour = "", y = "loss") +
  geom_hline(aes(yintercept = min_val_loss), colour = "red", alpha = 0.3, linetype = "dashed") +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_bw()  +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

n_networks <- length(network)
width <- 3.3 * n_networks

ggsave(
  figure, width = width, height = 3.3,
  file = paste0("loss.pdf"), device = "pdf", path = img_path
)
