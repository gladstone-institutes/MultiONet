#' compare_gene_sets
#'
#' Visualize overlap of significantly enriched terms across multiple gene sets using an UpSet plot.
#'
#' @param gene_sets A named list of ORA result data frames (e.g., from clusterProfiler::enrichGO), each containing at least `ID`, `pvalue`, and `qvalue` columns.
#' @param pvalueCutoff Numeric. P-value threshold to consider a term significant (default 0.05).
#' @param qvalueCutoff Numeric. Q-value threshold to consider a term significant (default 0.05).
#'
#' @return A binary matrix used to generate an UpSet plot. The plot is shown automatically.
#' @export
#' 
compare_gene_sets <- function(gene_sets, pvalueCutoff = 0.05, qvalueCutoff = 0.05) {
  require(UpSetR)
  
  if (length(gene_sets) < 2) {
    stop("Need at least two gene sets.")
  }
  
  # Filter each ORA result based on pvalue and qvalue thresholds
  filtered_sets <- lapply(gene_sets, function(df) {
    df=as.data.frame(df)
    df <- df[df$pvalue < pvalueCutoff & df$qvalue < qvalueCutoff, ]
    unique(df$ID)
  })
  
  # Collect all unique terms across sets
  all_terms <- unique(unlist(filtered_sets))
  
  # Create a binary matrix where rows are terms and columns are gene sets
  bin_mat <- sapply(filtered_sets, function(set) {
    as.integer(all_terms %in% set)
  })
  
  rownames(bin_mat) <- all_terms
  bin_df <- as.data.frame(bin_mat)
  
  # Generate the UpSet plot
  upset(bin_df,
        sets = names(filtered_sets),
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Shared Gene Set Count")
  
  return(invisible(bin_df))
}
