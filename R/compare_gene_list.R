#' compare_gene_list
#'
#' Visualize overlaps among gene lists using a Venn diagram.
#'
#' @param gene_sets A named list of gene vectors. Each element should be a character vector of gene names.
#'
#' @return A ggplot2 Venn diagram object is printed to the console.
#' @export
compare_gene_list <- function(gene_sets) {
  require(ggVennDiagram)
  require(ggplot2)
  
  if (length(gene_sets) > 5) {
    warning("Number of sets should be less than or equal to 5 for Venn diagram visualization.")
  }
  
  p <- ggVennDiagram(gene_sets) +
    scale_fill_gradient()
  print(p)
  
  invisible(p)
}