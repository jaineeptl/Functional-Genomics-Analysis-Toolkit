load_seurat_object <- function(cfg) {
  in_path <- cfg$input_rds
  if (is.null(in_path) || is.na(in_path)) {
    message("Choose a Seurat .rds file...")
    in_path <- file.choose()
  }
  if (!file.exists(in_path)) stop("Input .rds does not exist: ", in_path)

  obj <- readRDS(in_path)
  if (!inherits(obj, "Seurat")) stop("Loaded object is not a Seurat object.")

  message("Loaded: ", in_path)
  message("Assays: ", paste(Assays(obj), collapse = ", "))
  message("Default assay: ", DefaultAssay(obj))
  message("Metadata columns: ", paste(colnames(obj@meta.data), collapse = ", "))

  obj
}