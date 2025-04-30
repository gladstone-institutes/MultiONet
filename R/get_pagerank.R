#' get_pagerank
#'
#' Compute PageRank scores for a network and visualize using ggnet2.
#'
#' @param network An igraph object representing a gene or protein interaction network.
#'
#' @return A named numeric vector of PageRank scores for each node. A ggnet2 plot is displayed.
#' @export
get_pagerank <- function(network, label_size=2) {
  require(GGally)
  require(igraph)
  
  weights <- page_rank(network)$vector
  
  gp <- ggnet2(network,
               label = TRUE,
               color = "steelblue",
               size = weights,
               label.size = label_size,
               legend.size = FALSE,
               legend.color = FALSE) +
    theme(legend.position = "none")
  
  print(gp)
  return(weights)
}
