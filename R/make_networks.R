#' make_networks
#'
#' Construct and visualize gene interaction networks using STRINGdb data
#'
#' @param gene_sets A named list of gene vectors, each vector a gene set
#' @param indirect_neighbors Integer: how many network steps to expand (default 3)
#' @param score_threshold Minimum STRING score to include edge (default 700)
#' @param network_file_prefix Prefix for CSV files written to disk (default "network")
#' @param string_link_file Path to STRING links file (gzipped)
#' @param string_alias_file Path to STRING aliases file (gzipped)
#' @param plot_network Logical: whether to plot each network (default FALSE)
#' @param id_source ID source for STRING aliases (e.g., "Ensembl_UniProt_GN")
#' @param node_size Size of nodes in the network plot
#' @param label_genes Logical: whether to label genes in plots
#'
#' @return A list containing two elements: networks and largest_networks, each a list of igraph objects
#' @export
make_networks <- function(gene_sets, indirect_neighbors=3, score_threshold=700,
                          network_file_prefix="network",
                          string_link_file="9606.protein.links.full.v11.5.txt.gz",
                          string_alias_file="9606.protein.aliases.v11.5.txt.gz",
                          plot_network=FALSE, plot_network_node_threshold = 100, id_source="Ensembl_UniProt_GN",
                          node_size=5, label_genes=TRUE) {
  
  require(ggnetwork)
  require(dplyr)
  require(igraph)
  require(GGally)
  
  aliases <- read.delim(gzfile(string_alias_file), stringsAsFactors = FALSE)
  aliases <- aliases[aliases$source == id_source, ]
  links <- read.delim(gzfile(string_link_file), sep=" ")
  filtered <- links[links$combined_score > score_threshold, ]
  
  unique_pairs <- filtered[, c("protein1", "protein2")] %>%
    mutate(
      col_min = pmin(protein1, protein2),
      col_max = pmax(protein1, protein2)
    ) %>%
    distinct(col_min, col_max, .keep_all = TRUE) %>%
    select(-col_min, -col_max)
  
  g <- graph_from_data_frame(unique_pairs, directed = FALSE)
  
  current_names <- V(g)$name
  name_map <- setNames(aliases$alias, aliases$X.string_protein_id)
  V(g)$name <- ifelse(current_names %in% names(name_map), name_map[current_names], current_names)
  
  duplicated_names <- duplicated(V(g)$name)
  if (any(duplicated_names)) V(g)$name <- make.unique(V(g)$name)
  
  edge_table <- igraph::as_data_frame(g, what = "edges")
  node_table <- igraph::as_data_frame(g, what = "vertices")
  
  intersection=Reduce(intersect, gene_sets)
  if(length(intersection)>0)  gene_sets <- c(gene_sets, list(intersection=intersection) )
    
  gene_sets <- c(gene_sets, list(combined = as.matrix(unique(unlist(gene_sets)))))
  

  
  networks <- list()
  largest_networks <- list()
  
  for (gsi in seq_along(gene_sets)) {
    cat("Processing", names(gene_sets)[gsi], "...\n")
    mapped <- merge(gene_sets[[gsi]], aliases, by.x = "V1", by.y = "alias", all.x = TRUE)
    start_nodes <- unique(mapped$V1)
    valid_start_nodes <- start_nodes[start_nodes %in% V(g)$name]
    
    if (length(valid_start_nodes) == 0) {
      cat("No mapped genes in", names(gene_sets)[gsi], "\n")
      next
    }
    
    hop_nodes <- unique(unlist(ego(g, order = indirect_neighbors, nodes = valid_start_nodes, mode = "all")))
    sub_g <- induced_subgraph(g, hop_nodes)
    sub_g <- delete_vertices(sub_g, which(degree(sub_g) == 0))
    
    if (plot_network && vcount(sub_g) <= plot_network_node_threshold) {
      net <- intergraph::asNetwork(sub_g)
      gp <- ggnet2(net, label = label_genes, color = "steelblue", size = node_size, label.size = 3) +
        labs(title = paste0("Subnetwork for ", names(gene_sets)[gsi]))
      print(gp)
    }
    
    comps <- clusters(sub_g, mode = "weak")
    largest <- induced_subgraph(sub_g, V(sub_g)[comps$membership == which.max(comps$csize)])
    
    if (plot_network && vcount(largest) <= plot_network_node_threshold) {
      net <- intergraph::asNetwork(largest)
      gp <- ggnet2(net, label = label_genes, color = "steelblue", size = node_size, label.size = 3) +
        labs(title = paste0("Largest subnetwork for ", names(gene_sets)[gsi]))
      print(gp)
    }
    
    networks[[gsi]] <- sub_g
    largest_networks[[gsi]] <- largest
    
    write.csv(edge_table, paste0(network_file_prefix, "_", names(gene_sets)[gsi], "_edges.csv"), row.names = FALSE)
    write.csv(cbind(node_table, GeneSet = names(gene_sets)[gsi]),
              paste0(network_file_prefix, "_", names(gene_sets)[gsi], "_nodes.csv"), row.names = FALSE)
  }
  
  # network statstics
  network_summary <- data.frame(
    Network = paste0("Network from ", names(gene_sets), seq_along(networks)),
    Node_Count = sapply(networks, vcount),
    Edge_Count = sapply(networks, ecount)
  )
  
  print(network_summary)
  
  largest_network_summary <- data.frame(
    Network = paste0("Largest sub network from ", names(gene_sets), seq_along(largest_networks)),
    Node_Count = sapply(largest_networks, vcount),
    Edge_Count = sapply(largest_networks, ecount)
  )
  
  print(largest_network_summary)
  
  # pairwise network statistics
  n <- length(gene_sets)
  node_intersections <- matrix(0, n, n)
  edge_intersections <- matrix(0, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      nodes_i <- V(networks[[i]])$name
      nodes_j <- V(networks[[j]])$name
      node_intersections[i, j] <- length(intersect(nodes_i, nodes_j))
      
      # Get edges as sets of strings (unordered pair representation)
      edges_i <- apply(as_edgelist(networks[[i]]), 1, function(x) paste(sort(x), collapse = "_"))
      edges_j <- apply(as_edgelist(networks[[j]]), 1, function(x) paste(sort(x), collapse = "_"))
      edge_intersections[i, j] <- length(intersect(edges_i, edges_j))
    }
  }
  
  # Convert to data frames
  colnames(node_intersections) <- rownames(node_intersections) <- names(gene_sets)
  colnames(edge_intersections) <- rownames(edge_intersections) <- names(gene_sets)
  
  node_df <- as.data.frame(node_intersections)
  edge_df <- as.data.frame(edge_intersections)
  
  cat("Node Intersections:\n")
  print(node_df)
  
  cat("\nEdge Intersections:\n")
  print(edge_df)
  
  
  return(list(networks = networks, largest_networks = largest_networks))
}
