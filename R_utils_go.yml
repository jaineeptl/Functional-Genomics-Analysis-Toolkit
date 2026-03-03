suppressPackageStartupMessages({
  library(dplyr)
  library(GO.db)
})

get_go_subtree <- function(root_go_id) {
  if (is.null(root_go_id) || is.na(root_go_id) || root_go_id == "") return(NULL)
  subtree <- unique(c(root_go_id, as.character(GOBPOFFSPRING[[root_go_id]])))
  subtree
}

filter_enrich_result_to_subtree <- function(enrich_df, root_go_id) {
  subtree <- get_go_subtree(root_go_id)
  if (is.null(subtree)) return(enrich_df)
  if (!("ID" %in% colnames(enrich_df))) stop("enrich_df must include an ID column (GO IDs).")
  enrich_df %>% dplyr::filter(.data$ID %in% subtree)
}