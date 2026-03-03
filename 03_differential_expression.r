run_de <- function(obj, cfg) {
  DefaultAssay(obj) <- cfg$default_assay %||% DefaultAssay(obj)

  ident_col <- cfg$ident_col
  if (!ident_col %in% colnames(obj@meta.data)) stop("ident_col not in metadata: ", ident_col)

  obj <- SetIdent(obj, value = obj@meta.data[[ident_col]])

  a <- as.character(cfg$group_a)
  b <- as.character(cfg$group_b)

  if (!a %in% levels(Idents(obj))) stop("group_a not found in identities: ", a)
  if (!b %in% levels(Idents(obj))) stop("group_b not found in identities: ", b)

  message("DE: ", ident_col, " | ", a, " vs ", b, " | assay=", DefaultAssay(obj))

  de <- FindMarkers(
    object = obj,
    ident.1 = a,
    ident.2 = b,
    test.use = cfg$de$test_use %||% "wilcox",
    min.pct = as.numeric(cfg$de$min_pct %||% 0.1),
    logfc.threshold = as.numeric(cfg$de$logfc_threshold %||% 0),
    only.pos = isTRUE(cfg$de$only_pos %||% FALSE)
  )

  de <- de %>%
    tibble::rownames_to_column("gene")

  de
}