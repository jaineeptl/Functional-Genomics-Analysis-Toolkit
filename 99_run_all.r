source("00_setup.R")
source("01_load_seurat.R")
source("02_subset_filters.R")
source("03_differential_expression.R")
source("04_rank_genes.R")
source("05_go_enrichment.R")
source("06_visualize_go_barplots.R")
source("07_visualize_go_chordplots.R")
source("08_query_specific_go_terms.R")

cfg <- read_config("config.yml")
ensure_dir(cfg$out_dir)
ensure_dir(cfg$fig_dir)

# 1) Load
obj <- load_seurat_object(cfg)

# 2) Set assay
if (!is.null(cfg$default_assay)) DefaultAssay(obj) <- cfg$default_assay

# 3) Filter
obj_f <- obj
obj_f <- apply_subset_rules(obj_f, cfg)
obj_f <- apply_gene_filter(obj_f, cfg)

# 4) DE
de <- run_de(obj_f, cfg)

# Standardize DE for downstream
ranked <- rank_de_genes(de)
de_ranked <- ranked$de_ranked
fc_col <- ranked$fc_col

de_tbl <- de_ranked %>% dplyr::transmute(gene = gene, logFC = as.numeric(.data[[fc_col]]), p_val = p_val, p_adj = p_val_adj)

save_table(de, file.path(cfg$out_dir, paste0("de_", cfg$group_a, "_vs_", cfg$group_b, ".csv")))
save_table(de_ranked, file.path(cfg$out_dir, paste0("de_ranked_", cfg$group_a, "_vs_", cfg$group_b, ".csv")))

# 5) Pick genes for GO
sig <- de_ranked %>% dplyr::filter(is.finite(p_val), p_val <= as.numeric(cfg$rank$min_p %||% 0.05))
use_fallback <- nrow(sig) == 0

spl <- split_up_down(de_ranked, fc_col, top_n = as.numeric(cfg$rank$top_n %||% 1000))
genes_up_a <- spl$up_a$gene
genes_up_b <- spl$up_b$gene

# 6) GO enrichment (up in A and up in B)
ego_a <- run_go(genes_up_a, cfg)
ego_b <- run_go(genes_up_b, cfg)

go_a <- enrich_to_df(ego_a, cfg)
go_b <- enrich_to_df(ego_b, cfg)

if (!is.null(go_a)) save_table(go_a, file.path(cfg$out_dir, paste0("go_up_", cfg$group_a, ".csv")))
if (!is.null(go_b)) save_table(go_b, file.path(cfg$out_dir, paste0("go_up_", cfg$group_b, ".csv")))

# 7) Barplots
if (!is.null(go_a)) {
  make_go_barplot(
    go_df = go_a,
    de_tbl = de_tbl,
    title = paste0("GO (", cfg$go$ont, ") | Up in ", cfg$group_a),
    out_path = file.path(cfg$fig_dir, paste0("go_bar_up_", cfg$group_a, ".png")),
    max_terms = as.numeric(cfg$plots$barplot_max_terms %||% 30)
  )
}

if (!is.null(go_b)) {
  make_go_barplot(
    go_df = go_b,
    de_tbl = de_tbl,
    title = paste0("GO (", cfg$go$ont, ") | Up in ", cfg$group_b),
    out_path = file.path(cfg$fig_dir, paste0("go_bar_up_", cfg$group_b, ".png")),
    max_terms = as.numeric(cfg$plots$barplot_max_terms %||% 30)
  )
}

# 8) Chord plots
if (!is.null(go_a)) {
  make_go_chord(
    go_df = go_a,
    de_tbl = de_tbl,
    title = paste0("Top GO terms (Up in ", cfg$group_a, ") ↔ genes (logFC)"),
    out_path = file.path(cfg$fig_dir, paste0("go_chord_up_", cfg$group_a, ".png")),
    top_terms = as.numeric(cfg$plots$chord_top_terms %||% 5),
    max_genes_per_term = as.numeric(cfg$plots$chord_max_genes_per_term %||% 12),
    min_abs_logfc = as.numeric(cfg$plots$chord_min_abs_logfc %||% 0.25),
    label_top_genes_total = as.numeric(cfg$plots$chord_label_top_genes_total %||% 40)
  )
}

if (!is.null(go_b)) {
  make_go_chord(
    go_df = go_b,
    de_tbl = de_tbl,
    title = paste0("Top GO terms (Up in ", cfg$group_b, ") ↔ genes (logFC)"),
    out_path = file.path(cfg$fig_dir, paste0("go_chord_up_", cfg$group_b, ".png")),
    top_terms = as.numeric(cfg$plots$chord_top_terms %||% 5),
    max_genes_per_term = as.numeric(cfg$plots$chord_max_genes_per_term %||% 12),
    min_abs_logfc = as.numeric(cfg$plots$chord_min_abs_logfc %||% 0.25),
    label_top_genes_total = as.numeric(cfg$plots$chord_label_top_genes_total %||% 40)
  )
}

# 9) Example: extract genes from one GO term and plot their logFC
# (Only runs if subtree root is set and GO results exist)
example_go <- cfg$go$subtree_root
if (!is.null(go_a) && !is.null(example_go)) {
  genes_term <- extract_genes_from_go_term(go_a, example_go)
  if (length(genes_term) > 0) {
    plot_gene_logfc(
      de_tbl = de_tbl,
      genes = genes_term,
      title = paste0("Genes in ", example_go, " (from GO results)"),
      out_path = file.path(cfg$fig_dir, paste0("genes_", gsub(":", "_", example_go), ".png"))
    )
  }
}

message("Done.")