% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_networks.R
\name{make_networks}
\alias{make_networks}
\title{make_networks}
\usage{
make_networks(
  gene_sets,
  indirect_neighbors = 3,
  score_threshold = 700,
  number_of_largest_networks = 1,
  network_file_prefix = "network",
  string_link_file = "9606.protein.links.full.v11.5.txt.gz",
  string_alias_file = "9606.protein.aliases.v11.5.txt.gz",
  plot_network = FALSE,
  plot_network_node_threshold = 100,
  id_source = "Ensembl_UniProt_GN",
  node_size = 5,
  label_genes = TRUE,
  label_size = 2
)
}
\arguments{
\item{gene_sets}{A named list of gene vectors, each vector a gene set}

\item{indirect_neighbors}{Integer: how many network steps to expand (default 3)}

\item{score_threshold}{Minimum STRING score to include edge (default 700)}

\item{number_of_largest_networks}{Integer: how many largest sub-networks to make (default 1)}

\item{network_file_prefix}{Prefix for CSV files written to disk (default "network")}

\item{string_link_file}{Path to STRING links file (gzipped)}

\item{string_alias_file}{Path to STRING aliases file (gzipped)}

\item{plot_network}{Logical: whether to plot each network (default FALSE)}

\item{id_source}{ID source for STRING aliases (e.g., "Ensembl_UniProt_GN")}

\item{node_size}{Size of nodes in the network plot}

\item{label_genes}{Logical: whether to label genes in plots}
}
\value{
A list containing two elements: networks and largest_networks, each a list of igraph objects
}
\description{
Construct and visualize gene interaction networks using STRINGdb data
}
