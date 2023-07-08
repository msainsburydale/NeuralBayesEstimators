suppressMessages({
library("dplyr")
library("ggplot2")
library("viridis")
})

model <- "GaussianProcess/nuFixed"
params_path <- paste0("intermediates/", model, "/parameter_configurations/")
img_path    <- paste0("img/", model)
dir.create(img_path, recursive = TRUE, showWarnings = FALSE)

all_sets <- c("train", "val", "test", "scenarios")

# ---- Load the parameter configurations ----

xi <- lapply(all_sets, function(set) {
  load(file = paste0(params_path, set, "_xi.rda"))
  xi <- data.frame(xi, set = set)
  return(xi)
})
df <- do.call(rbind, xi)

# Drop any parameters that contain only a single value (e.g., fixed smoothness):
df <- df[, which(sapply(df, function(x) length(unique(x))) > 1)]

# Convert sigma to log(sigma^2), for consistency with Gerber and Nychka (2021):
df$logsigma2 <- log(df$sigma^2)
df$sigma <- NULL


# ---- Plot the parameter configurations ----

# aesthetics for the plot. See, e.g., show_col(viridis(6))
set_colours <- c("#440154FF", "#21908CFF", "#FDE725FF", "red")
set_labels <- c(expression(theta[train]), expression(theta[val]), expression(theta[test]), expression(theta[scenarios]))
set_sizes <- c(1, 1, 1, 2) / 2
names(set_labels) <- names(set_sizes) <- names(set_colours) <- all_sets

figure <- ggplot(data = df, aes(x = logsigma2, y = rho, colour = set, size = set)) +
  geom_point() +
  labs(colour = "", x = expression(log(sigma[epsilon]^2)), y = expression(rho)) +
  scale_color_manual(values = set_colours,
                     labels = set_labels,
                     breaks = all_sets[which(all_sets != "scenarios")]) +
                     # breaks = all_sets) +
  scale_size_manual(values = set_sizes) +
  guides(size = "none") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.text = element_text(size = 13))

ggsave(
  figure, filename = "ParameterConfigurations.png",
  path = img_path, device = "png", width = 5, height = 3
)

ggsave(
  figure, filename = "ParameterConfigurations.pdf",
  path = img_path, device = "pdf", width = 5, height = 3
)
