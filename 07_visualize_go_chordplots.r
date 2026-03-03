make_go_chord <- function(go_df, de_tbl, title, out_path = NULL,
                          top_terms = 5, max_genes_per_term = 12,
                          min_abs_logfc = 0.25, label_top_genes_total = 40,
                          gene_scale = c("white", "red3")) {

  if (is.null(go_df) || nrow(go_df) == 0) {
    warning("GO dataframe is empty; skipping chordplot.")
    return(NULL)
  }
  stopifnot(all(c("Description", "p.adjust", "geneID") %in% colnames(go_df)))

  top_terms_df <- go_df %>%
    dplyr::transmute(
      term = as.character(Description),
      padj = as.numeric(p.adjust),
      geneID = as.character(geneID)
    ) %>%
    dplyr::filter(!is.na(term), term != "", is.finite(padj), !is.na(geneID), geneID != "") %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(term, .keep_all = TRUE) %>%
    dplyr::slice_head(n = top_terms)

  edges <- top_terms_df %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(gene = geneID) %>%
    dplyr::inner_join(de_tbl, by = "gene") %>%
    dplyr::filter(is.finite(logFC), abs(logFC) >= min_abs_logfc) %>%
    dplyr::group_by(term) %>%
    dplyr::slice_max(order_by = abs(logFC), n = max_genes_per_term, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(term, gene, .keep_all = TRUE)

  if (nrow(edges) < 2) {
    warning("Not enough edges after filtering for chord plot.")
    return(NULL)
  }

  # label only top genes by max abs(logFC)
  genes_keep <- edges %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(max_abs = max(abs(logFC)), .groups = "drop") %>%
    dplyr::arrange(desc(max_abs)) %>%
    dplyr::slice_head(n = min(label_top_genes_total, nrow(.))) %>%
    dplyr::pull(gene)

  mat <- edges %>% dplyr::transmute(term = term, gene = gene, value = 1)

  term_levels <- unique(mat$term)
  gene_levels <- unique(mat$gene)

  term_pal <- brewer.pal(min(12, length(term_levels)), "Set3")
  if (length(term_pal) < length(term_levels)) {
    term_pal <- c(term_pal, circlize::rand_color(length(term_levels) - length(term_pal)))
  }
  term_cols <- setNames(term_pal[seq_along(term_levels)], term_levels)

  gene_logfc <- edges %>% dplyr::distinct(gene, logFC)
  lim <- max(abs(gene_logfc$logFC), na.rm = TRUE)
  lim <- max(lim, 0.25)

  col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("blue3", "white", gene_scale[2]))
  gene_cols <- setNames(col_fun(gene_logfc$logFC), gene_logfc$gene)

  grid.col <- c(term_cols, gene_cols)
  link_cols <- adjustcolor(term_cols[mat$term], alpha.f = 0.75)

  if (!is.null(out_path)) {
    ensure_dir(dirname(out_path))
    grDevices::png(out_path, width = 1800, height = 1400, res = 200)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4, 4, 4, 14))

  circos.clear()
  circos.par(start.degree = 90, track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0))

  chordDiagram(
    x = mat,
    grid.col = grid.col,
    col = link_cols,
    transparency = 0.4,
    link.border = "grey25",
    link.lwd = 1,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )

  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")

      is_gene <- !sector %in% term_levels
      if (is_gene && sector %in% genes_keep) {
        circos.text(mean(xlim), ylim[1] + 0.12, sector,
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 0.45)
      }
    },
    bg.border = NA
  )

  title(title)
  invisible(TRUE)
}