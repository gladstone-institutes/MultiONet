#' Reduce GO Enrichment Result Using Semantic Similarity
#'
#' This function reduces redundancy in GO enrichment results by grouping similar GO terms
#' based on semantic similarity using the `rrvgo` and `GOSemSim` packages.
#'
#' @param enrichment_result A data.frame or enrichResult object with columns \code{ID} and \code{p.adjust}.
#'
#' @return A list with:
#' \describe{
#'   \item{processed_result}{All terms mapped to parents with semantic similarity}
#'   \item{reduced_result}{Subset with unique parent terms only}
#' }
#' @import rrvgo org.Hs.eg.db GO.db GOSemSim
#' @export
#'
#' @examples
#' \dontrun{
#'   library(clusterProfiler)
#'   ego <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "BP")
#'   reduced <- reduce_enrichment_result(ego@result)
#' }
#' 

reduce_enrichment_result <- function(enrichment_result) {
  if (!requireNamespace("rrvgo", quietly = TRUE)) stop("Package 'rrvgo' is required.")
  if (!requireNamespace("GOSemSim", quietly = TRUE)) stop("Package 'GOSemSim' is required.")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Package 'org.Hs.eg.db' is required.")
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("Package 'GO.db' is required.")
  
  terms <- enrichment_result$ID
  scores <- -log10(enrichment_result$p.adjust)
  
  semData <- GOSemSim::godata('org.Hs.eg.db', ont = "BP")
  flag <- terms %in% names(semData@IC)
  
  terms <- terms[flag]
  scores <- scores[flag]
  
  simMatrix <- rrvgo::calculateSimMatrix(terms, orgdb = org.Hs.eg.db, method = "Wang")
  names(scores) <- rownames(simMatrix)
  
  reduced <- rrvgo::reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = org.Hs.eg.db)
  reduced_unique <- reduced[reduced$go == reduced$parent, ]
  
  return(list(processed_result = reduced, reduced_result = reduced_unique))
}
