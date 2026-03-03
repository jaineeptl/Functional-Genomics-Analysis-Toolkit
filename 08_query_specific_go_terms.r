extract_genes_from_go_term <- function(ego_result_df, go_id) {
  if (is.null(ego_result_df) || nrow(ego_result_df) == 0) return(character(0))
  stopifnot("ID" %in% colnames(ego_result_df))
  row <- ego_result_df %>% dplyr::filter(.data$ID == go_id)
  if (nrow(row) == 0) return(character(0))
  genes <- row$geneID[[1]]
  unique(unlist(strsplit(genes, "/")))
}

plot_gene_logfc <- function(de_tbl, genes, title, out_path = NULL) {
  df <- de_tbl %>%
    dplyr::filter(gene %in% genes) %>%
    dplyr::arrange(logFC)

  if (nrow(df) < 2) {
    warning("Too few genes found in DE table for plot.")
    return(NULL)
  }

  p <- ggplot(df, aes(x = logFC, y = reorder(gene, logFC), color = logFC)) +
    geom_segment(aes(x = 0, xend = logFC, yend = gene), linewidth = 0.8, color = "grey65") +
    geom_point(size = 3) +
    scale_color_gradient2(midpoint = 0, name = "logFC") +
    labs(title = title, x = "log fold change", y = NULL) +
    theme_clean(10)

  if (!is.null(out_path)) save_plot(p, out_path, width = 8, height = 6)
  p
}