make_go_barplot <- function(go_df, de_tbl, title, out_path = NULL, max_terms = 30) {
  if (is.null(go_df) || nrow(go_df) == 0) {
    warning("GO dataframe is empty; skipping barplot.")
    return(NULL)
  }
  stopifnot(all(c("Description", "Count", "p.adjust", "geneID") %in% colnames(go_df)))

  go_df2 <- go_df %>%
    dplyr::transmute(
      term = as.character(Description),
      count = as.numeric(Count),
      padj = as.numeric(p.adjust),
      geneID = as.character(geneID)
    ) %>%
    dplyr::filter(is.finite(count), count > 0, is.finite(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = max_terms)

  # Expand term->genes and compute mean logFC for genes in term (if overlap exists)
  go_long <- go_df2 %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(gene = geneID)

  term_fc <- go_long %>%
    dplyr::inner_join(de_tbl, by = "gene") %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(term_logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

  go_plot <- go_df2 %>%
    dplyr::left_join(term_fc, by = "term") %>%
    dplyr::mutate(term_logFC = ifelse(is.na(term_logFC), 0, term_logFC))

  p <- ggplot(go_plot, aes(x = reorder(term, count), y = count, fill = term_logFC)) +
    geom_col(width = 0.85) +
    coord_flip() +
    scale_fill_gradient2(midpoint = 0, name = "Mean logFC\n(genes in term)") +
    labs(title = title, x = NULL, y = "Gene count") +
    theme_clean(10)

  if (!is.null(out_path)) save_plot(p, out_path, width = 9, height = 6)
  p
}