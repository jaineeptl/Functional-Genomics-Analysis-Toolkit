run_go <- function(genes_symbol, cfg) {
  genes_symbol <- unique(genes_symbol)
  genes_symbol <- genes_symbol[!is.na(genes_symbol) & genes_symbol != ""]
  if (length(genes_symbol) < 5) {
    warning("Too few genes for GO enrichment (n < 5). Returning NULL.")
    return(NULL)
  }

  conv <- bitr(
    genes_symbol,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (is.null(conv) || nrow(conv) == 0) {
    warning("No SYMBOL -> ENTREZID mappings found.")
    return(NULL)
  }

  ego <- enrichGO(
    gene = conv$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = cfg$go$ont %||% "BP",
    pAdjustMethod = cfg$go$p_adjust_method %||% "BH",
    pvalueCutoff = as.numeric(cfg$go$pvalue_cutoff %||% 0.05),
    qvalueCutoff = as.numeric(cfg$go$qvalue_cutoff %||% 0.2),
    readable = isTRUE(cfg$go$readable %||% TRUE)
  )

  ego
}

enrich_to_df <- function(ego, cfg) {
  if (is.null(ego)) return(NULL)
  df <- ego@result
  if (!is.null(cfg$go$subtree_root)) {
    df <- filter_enrich_result_to_subtree(df, cfg$go$subtree_root)
  }
  df
}