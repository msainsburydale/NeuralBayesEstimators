path <- "SimulationOnTheFly"

# ---- Results ----

suppressMessages({
  library("dplyr")
  library("ggplot2")
  library("reshape2")
  options(dplyr.summarise.inform = FALSE) # Suppress summarise info
})

# Load the loss function per epoch files:
all_dirs  <- list.dirs(path = paste0("intermediates/", path), recursive = TRUE)
runs_dirs <- all_dirs[which(grepl("runs_", all_dirs))]
networks  <- c("smallNetwork", "largeNetwork")
schemes   <- c("everyEpoch", "someEpochs", "noEpochs", "recycled")
loss_per_epoch_list <- expand.grid(networks, schemes) %>% 
  apply(1, function(study) {
    network = study[1]
    scheme = study[2]
    path <- paste("intermediates", path, network, scheme, "runs_N1/loss_per_epoch.csv", sep = "/")
    loss_per_epoch <- read.csv(path, header = FALSE)
    colnames(loss_per_epoch) <- c("training", "validation")
    loss_per_epoch$epoch <- 0:(nrow(loss_per_epoch) - 1)
    loss_per_epoch$network <- network
    loss_per_epoch$scheme <- scheme
    loss_per_epoch
  })

df <- do.call("rbind", loss_per_epoch_list)
df <- df %>% melt(c("epoch", "scheme", "network"), variable.name = "set", value.name = "loss")

# Edit the factor levels to control the facet order and display
refresh_epoch <- 30
df$scheme <- factor(
  df$scheme, 
  levels = c("everyEpoch", "someEpochs", "noEpochs", "recycled"), 
  labels = c("Simulate training data every epoch", 
             paste("Simulate training data every", refresh_epoch, "epochs"), 
             "Simulate training data once only", 
             "Cycle over large fixed training set")
)

library("stringr")
df <- df %>% 
  mutate(network = str_to_title(network)) %>% 
  mutate(set = str_to_title(set))
df$network <- gsub("network", " network", df$network)

df <- df %>% 
  mutate(study = paste(network, set, sep = ": "))

minimum_loss <- df %>% filter(set == "Validation") %>% summarise(min(loss)) %>% unlist

size <- 0.3 

img_path <- paste0("img/", path)
dir.create(img_path, showWarnings = FALSE, recursive = TRUE)

# Grid of plots
fig <- ggplot(df) + 
  geom_line(aes(x = epoch, y = pmin(loss, 0.9), colour = set), size = size) + 
  scale_color_manual(values = c("blue", "red")) +
  geom_hline(aes(yintercept = minimum_loss), colour = "red", alpha = 0.3, linetype = "dashed") +
  facet_grid(network ~ scheme) + 
  labs(colour = "", y = expression(r[Omega](hat(theta))), linetype = "") + 
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position="none") + 
  scale_y_continuous(expand = c(0, 0.01)) 


width = 3.3 * length(schemes)
ggsave(
  fig, width = width, height = 5.3,
  file = "loss_grid.pdf", device = "pdf", path = img_path
)

# Save the minimum loss function score:
# minima <- df %>% 
#   group_by(set, scheme, network) %>% 
#   summarise(min = min(loss)) %>% 
#   as.data.frame()
# 
# results_path <- paste0("results/", path)
# dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
# write.csv(minima, file = paste(results_path, "loss_minima.csv", sep = "/"))
