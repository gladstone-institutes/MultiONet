#' Create Directed Graph with Omnipath Directionality
#'
#' This function takes an igraph network, retrieves directionality from Omnipath, adjusts edge directions accordingly, and plots a directed network.
#'
#' @param network igraph object representing the initial network
#' @param node_size Numeric, size of the network nodes
#' @param label_genes Logical, whether to label genes (nodes) or not
#' @param label_size Numeric, font size for node labels
#'
#' @return A directed igraph network visualized using ggnet2
#'
#' @examples
#' make_directed_graph(network = your_network_object)
#'
#' @export
make_directed_graph <- function(network,
                                node_size = 5,
                                label_genes = TRUE,
                                label_size = 2) {
  
  # Load required packages
  require(OmnipathR)
  require(igraph)
  require(GGally)
  require(dplyr)
  require(magrittr)
  
  # Import directional interactions from Omnipath
  interactions <- OmnipathR::import_omnipath_interactions()
  directions <- interactions[, c("source_genesymbol", "target_genesymbol")]
  
  # Get edge list from the provided network
  network_edges <- as.data.frame(get.edgelist(network), stringsAsFactors = FALSE)
  colnames(network_edges) <- c("from", "to")
  
  # Identify edges that match Omnipath interactions (in correct or reverse direction)
  edges_forward <- paste0(network_edges$from, "_", network_edges$to)
  edges_backward <- paste0(network_edges$to, "_", network_edges$from)
  omni_edges <- paste0(directions$source_genesymbol, "_", directions$target_genesymbol)
  
  # Identify edges to keep direction
  directed_forward <- network_edges[edges_forward %in% omni_edges, ]
  directed_backward <- network_edges[edges_backward %in% omni_edges, ]
  
  # Reverse backward edges
  if (nrow(directed_backward) > 0) {
    for (i in seq_len(nrow(directed_backward))) {
      edge_id <- get.edge.ids(network, c(directed_backward$from[i], directed_backward$to[i]))
      network <- delete.edges(network, edge_id)
      network <- add.edges(network, c(directed_backward$to[i], directed_backward$from[i]))
    }
  }
  
  # Finalize directed edges
  final_edges <- rbind(directed_forward, directed_backward[, c("to", "from")])
  
  # Create final directed graph
  directed_network <- graph_from_data_frame(final_edges, directed = TRUE)
  
  if(dim(get.edgelist(directed_network))[1]>0){
    directed_network_simple <- igraph::simplify(directed_network, remove.multiple = TRUE, remove.loops = TRUE)
    
    # Plot the directed network
    GGally::ggnet2(directed_network_simple,
                   label = label_genes,
                   size = node_size,
                   label.size = label_size,
                   color = "steelblue",
                   arrow.size = 10,
                   arrow.gap = 0.025)
  }else{
    print("There is no direction assgiend to edges")
  }
  

}
