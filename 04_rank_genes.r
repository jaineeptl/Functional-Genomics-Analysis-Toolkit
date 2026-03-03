rank_de_genes <- function(de) {
  if (!("p_val" %in% colnames(de))) stop("DE table must include p_val.")
  fc_col <- if ("avg_log2FC" %in% colnames(de)) "avg_log2FC" else if ("avg_logFC" %in% colnames(de)) "avg_logFC" else NULL
  if (is.null(fc_col)) stop("DE table must include avg_log2FC or avg_logFC.")

  de2 <- de %>%
    dplyr::mutate(
      logFC = as.numeric(.data[[fc_col]]),
      p_val = as.numeric(.data$p_val),
      rank_score = abs(logFC) * -log10(p_val + 1e-300)
    ) %>%
    dplyr::filter(is.finite(rank_score)) %>%
    dplyr::arrange(desc(rank_score))

  list(de_ranked = de2, fc_col = fc_col)
}

split_up_down <- function(de_ranked, fc_col, top_n = 1000) {
  top <- head(de_ranked, top_n)

  up_a <- top %>% dplyr::filter(.data[[fc_col]] > 0)
  up_b <- top %>% dplyr::filter(.data[[fc_col]] < 0)

  list(top = top, up_a = up_a, up_b = up_b)
}