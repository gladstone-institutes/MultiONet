#' run_ora
#'
#' Perform Over-Representation Analysis (ORA) using Gene Ontology (GO) and optionally plot cnetplot
#'
#' @param genes A character vector of gene symbols
#' @param ontology Character. Ontology to use for GO analysis ("BP", "MF", or "CC")
#' @param pvalueCutoff Numeric. P-value cutoff for enrichment (default = 0.05)
#' @param qvalueCutoff Numeric. Q-value (adjusted p-value) cutoff for enrichment (default = 0.05)
#' @param showCategory Integer. Number of top enriched terms to show in the plot (default = 10)
#' @param plot Logical. Whether to display a cnetplot of enriched terms (default = FALSE)
#'
#' @return A data frame of enriched GO terms. Returns NULL if fewer than 5 genes are mapped or no terms are enriched.
#' @export
run_ora <- function(genes, ontology = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, showCategory = 10, plot = FALSE){
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(enrichplot)
  require(igraph)
  require(ggraph)
  
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_ids <- gene_df$ENTREZID
  
  message("Number of matched gene IDs: ", length(gene_ids))
  
  if (length(gene_ids) < 5) {
    warning("Too few mapped genes for enrichment.")
    return(NULL)
  }
  
  ego <- enrichGO(
    gene = gene_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(ego) == 0) {
    message("No enriched terms found.")
    return(NULL)
  }
  
  if (plot) {
    p <- cnetplot(ego, showCategory = showCategory, circular = FALSE, colorEdge = TRUE)
    print(p)
  }
  
  return(as.data.frame(ego))
}