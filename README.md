# MultiONet

**MultiONet** is an R package designed for **multi-omics network analysis**. It integrates differential expression and chromatin accessibility (or other omics) data through gene networks to support functional enrichment, graph analysis, and hypothesis generation.

This package enables researchers to generate interaction networks, reduce redundancy in GO terms, compare enriched pathways across datasets, and visualize subnetworks with biological relevance.

---

##  Features

- Construct STRING-based gene interaction networks for multiple omics inputs.
- Perform GO term over-representation analysis (ORA).
- Reduce redundancy in GO enrichment results using semantic similarity.
- Identify shared and unique GO terms across omics data.
- Visualize and compare subnetworks.
- Compute node centrality metrics (e.g., PageRank) for network prioritization.

---

##  Example: Kidney Injury Multi-omics Data

This example walks through using **MultiONet** on public kidney injury data from [GSE254185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185) (Maximilian Reck et al., *Nature Communications*, 2025). The goal is to compare transcriptomic and chromatin accessibility signatures in **proximal tubule (PT) injured cells**.

```r
# Load differential expression and peak data
data_a <- deg %>% filter(Cluster == 1, adj..p.value < 0.05)
data_b <- peak %>% filter(Cluster == 1, adj..p.value < 0.05)

# Compare DEG and peak gene lists to identify overlaps and differences
input <- list(deg = data_a$Gene, peak = data_b$Closest.Gene)
result <- compare_gene_list(input)

# ORA: individual and combined gene sets
ora_deg <- run_ora(data_a$Gene)
ora_peak <- run_ora(data_b$Closest.Gene)
ora_intersect <- run_ora(intersect(data_a$Gene, data_b$Closest.Gene))
ora_union <- run_ora(union(data_a$Gene, data_b$Closest.Gene))

# Network-based integration
input <- list(deg = data_a$Gene, peak = data_b$Closest.Gene)
result <- make_networks(input,
                        string_link_file = "9606.protein.links.full.v11.5.txt.gz",
                        string_alias_file = "9606.protein.aliases.v11.5.txt.gz",
                        score_threshold = 997,
                        plot_network = TRUE)

# Extract ORA from largest networks
ora_net_deg <- run_ora(V(result$largest_networks$Largest_network_of_deg_1)$name)
ora_net_peak <- run_ora(V(result$largest_networks$Largest_network_of_peak_1)$name)

# Perform ORA for genes in the largest networks then plot a upset plot to examine the similarity between the two results
compare_gene_sets(list(set1 = ora_net_deg, set2 = ora_net_peak))

# Reduce redundancy in GO terms
combined_raw <- bind_rows(ora_net_deg %>% mutate(source = "deg"),
                          ora_net_peak %>% mutate(source = "peak"))
reduced <- reduce_enrichment_result(combined_raw)

# Heatmap of top GO terms
reduced_heatmap(reduced$processed_result, combined_raw,
                set1_processed_ora = reduce_enrichment_result(ora_net_deg)$processed_result,
                set2_processed_ora = reduce_enrichment_result(ora_net_peak)$processed_result)

# Visualize directed networks with Omnipath, a database built from above 100 resources:
make_directed_graph(result[['largest_networks']][["Largest_network_of_deg_1"]])

# PageRank on combined network
get_pagerank(result$largest_networks$Largest_network_of_combined_1)
