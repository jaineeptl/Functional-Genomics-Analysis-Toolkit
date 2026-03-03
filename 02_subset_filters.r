apply_subset_rules <- function(obj, cfg) {
  if (is.null(cfg$subset_rules$enable) || !isTRUE(cfg$subset_rules$enable)) return(obj)

  rules <- cfg$subset_rules$rules
  if (is.null(rules) || length(rules) == 0) return(obj)

  md <- obj@meta.data

  keep <- rep(TRUE, nrow(md))
  for (r in rules) {
    col <- r$col
    if (!col %in% colnames(md)) stop("subset_rules col not found in metadata: ", col)

    if (!is.null(r$equals)) {
      keep <- keep & (as.character(md[[col]]) == as.character(r$equals))
    } else if (!is.null(r$regex)) {
      keep <- keep & grepl(r$regex, as.character(md[[col]]))
    } else {
      stop("Each subset rule must have 'equals' or 'regex'.")
    }
  }

  obj2 <- obj[, keep]
  message("Subset kept cells: ", ncol(obj2), " / ", ncol(obj))
  obj2
}

apply_gene_filter <- function(obj, cfg) {
  gf <- cfg$gene_filter
  if (is.null(gf$enable) || !isTRUE(gf$enable)) return(obj)

  gene <- as.character(gf$gene)
  min_expr <- as.numeric(gf$min_expr %||% 0)

  if (!gene %in% rownames(obj)) stop("gene_filter gene not found in object: ", gene)

  expr <- FetchData(obj, vars = gene)[, 1]
  obj2 <- obj[, expr > min_expr]
  message("Gene filter (", gene, " > ", min_expr, ") kept cells: ", ncol(obj2), " / ", ncol(obj))
  obj2
}

`%||%` <- function(a, b) if (!is.null(a)) a else b