% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_gene_sets.R
\name{compare_gene_sets}
\alias{compare_gene_sets}
\title{compare_gene_sets}
\usage{
compare_gene_sets(gene_sets, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
}
\arguments{
\item{gene_sets}{A named list of ORA result data frames (e.g., from clusterProfiler::enrichGO), each containing at least \code{ID}, \code{pvalue}, and \code{qvalue} columns.}

\item{pvalueCutoff}{Numeric. P-value threshold to consider a term significant (default 0.05).}

\item{qvalueCutoff}{Numeric. Q-value threshold to consider a term significant (default 0.05).}
}
\value{
A binary matrix used to generate an UpSet plot. The plot is shown automatically.
}
\description{
Visualize overlap of significantly enriched terms across multiple gene sets using an UpSet plot.
}
