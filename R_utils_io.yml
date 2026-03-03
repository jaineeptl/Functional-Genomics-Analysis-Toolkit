suppressPackageStartupMessages({
  library(yaml)
})

read_config <- function(path = "config.yml") {
  if (!file.exists(path)) stop("config.yml not found at repo root.")
  cfg <- yaml::read_yaml(path)
  cfg
}

ensure_dir <- function(x) {
  if (is.null(x) || is.na(x) || x == "") stop("Directory path is empty.")
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(x)
}

stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

save_table <- function(df, out_path) {
  ensure_dir(dirname(out_path))
  utils::write.csv(df, out_path, row.names = FALSE)
  message("Wrote: ", out_path)
}

save_rds <- function(obj, out_path) {
  ensure_dir(dirname(out_path))
  saveRDS(obj, out_path)
  message("Wrote: ", out_path)
}