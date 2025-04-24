#' find_mcl
#'
#' Detect MCL (Markov Cluster Algorithm) subnetworks in a graph and optionally plot the results.
#'
#' @param network An igraph object representing the input network.
#' @param min_node Integer. Minimum number of nodes to retain a cluster (default = 10).
#' @param expansion Integer. MCL expansion parameter (default = 2).
#' @param inflation Numeric. MCL inflation parameter (default = 1.5).
#' @param plot_network Logical. Whether to plot the resulting MCL subnetworks (default = FALSE).
#' @param label_genes Logical. Whether to display vertex labels in the plot (default = TRUE).
#'
#' @return A named list of igraph objects representing the largest filtered MCL subnetworks.
#' @export
find_mcl <- function(network, min_node = 10, expansion = 2, inflation = 1.5, plot_network = FALSE, label_genes = TRUE) {
  require(igraph)
  require(MCL)
  
  network_simplified <- simplify(network, remove.multiple = TRUE, remove.loops = FALSE)
  adj <- as.matrix(as_adjacency_matrix(network_simplified))
  mcl_result <- mcl(adj, addLoops = TRUE, expansion = expansion, inflation = inflation)
  clusters <- mcl_result$Cluster
  
  V(network)$cluster <- clusters
  
  subnetworks <- lapply(unique(clusters), function(i) {
    nodes_in_comp <- which(mcl_result$Cluster == i)
    if (length(nodes_in_comp) >= min_node) {
      subgraph <- induced_subgraph(network, nodes_in_comp)
      subgraph <- delete_vertices(subgraph, which(degree(subgraph) == 0))
      components <- clusters(subgraph, mode = "weak")
      biggest_cluster_id <- which.max(components$csize)
      vert_ids <- V(subgraph)[components$membership == biggest_cluster_id]
      largest <- induced_subgraph(subgraph, vert_ids)
      if (vcount(largest) > min_node) largest else NULL
    } else {
      NULL
    }
  })
  
  subnetworks <- Filter(Negate(is.null), subnetworks)
  names(subnetworks) <- paste0("mcl", seq_along(subnetworks))
  
  message("Number of MCL networks: ", length(subnetworks))
  
  largest_network_summary <- data.frame(
    Node_Count = sapply(subnetworks, vcount),
    Edge_Count = sapply(subnetworks, ecount)
  )
  print(largest_network_summary)
  
  if (plot_network && length(subnetworks) > 0) {
    g_combined <- Reduce(igraph::union, subnetworks)
    require(ggnetwork)
    require(GGally)
    require(ggplot2)
    
    gp <- ggnet2(g_combined,
                 label = label_genes,
                 color = "steelblue",
                 size = 5,
                 label.size = 3) +
      labs(title = "MCL networks")
    print(gp)
  }
  
  return(subnetworks)
}
