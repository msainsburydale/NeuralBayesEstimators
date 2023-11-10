path <- "Pretraining"

suppressMessages({
  library("dplyr")
  library("ggplot2")
  library("reshape2")
  library("stringr")
  options(dplyr.summarise.inform = FALSE) # Suppress summarise info
})

# Load the loss function per epoch files:
all_dirs  <- list.dirs(path = file.path("intermediates", path), recursive = TRUE)
runs_dirs <- all_dirs[which(grepl("runs_", all_dirs))]
# determine models automatically from runs_dirs (every unique value between "Pretraining/' and "/NotPretrained")
schemes   <- c("Pretrained", "NotPretrained")
loss_per_epoch_list <- lapply(schemes, function(scheme) {
  path <- file.path("intermediates", path, scheme, "runs_N30", "loss_per_epoch.csv")
  loss_per_epoch <- read.csv(path, header = FALSE)
  colnames(loss_per_epoch) <- c("training", "validation")
  loss_per_epoch$epoch <- 0:(nrow(loss_per_epoch) - 1)
  loss_per_epoch$scheme <- scheme
  loss_per_epoch
})

df <- do.call("rbind", loss_per_epoch_list)
df <- df %>% melt(c("epoch", "scheme"), variable.name = "set", value.name = "loss")

# Edit the factor levels to control the facet order and display
df$scheme <- factor(
  df$scheme,
  levels = c("Pretrained", "NotPretrained"),
  labels = c("Pre-trained", "Not pre-trained")
)

## Remove training set here since we're not testing overfitting.
df <- df %>% filter(set == "validation")

## Compute the global minimum validation risk
min_val_loss <- min(df$loss)

fig <- ggplot(df) +
  geom_line(aes(x = epoch, y = loss, linetype = scheme), colour = "red") +
  scale_linetype_manual(values = c("solid", "twodash")) +
  labs(colour = "", y = expression(r[Omega](hat(theta))), linetype = "") +
  geom_hline(aes(yintercept = min_val_loss), colour = "red", alpha = 0.3, linetype = "dashed") +
  theme_bw() +
  theme(strip.background = element_blank())

fig <- fig + coord_cartesian(ylim = c(0.2, 1))

img_path <- file.path("img", path)
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)
ggsave(
  fig, width = 6, height = 3.3,
  file = "risk.pdf", device = "pdf", path = img_path
)
