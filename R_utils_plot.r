suppressPackageStartupMessages({
  library(ggplot2)
})

theme_clean <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

save_plot <- function(p, out_path, width = 8, height = 5, dpi = 300) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, plot = p, width = width, height = height, dpi = dpi)
  message("Saved plot: ", out_path)
}