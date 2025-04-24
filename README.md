# MultiONet

## download STRING networks
download.file(
  "https://stringdb-static.org/download/protein.links.full.v11.5/9606.protein.links.full.v11.5.txt.gz",
  destfile = "9606.protein.links.full.v11.5.txt.gz"
)

download.file(
  "https://stringdb-static.org/download/protein.alias.v11.5/9606.protein.aliases.v11.5.txt.gz",
  destfile = "9606.protein.aliases.v11.5.txt.gz"
)


## make gene sets
gene_sets=list(Set1=as.matrix(set1),Set2=as.matrix(set2))

## Venn diagram
compare_gene_sets(gene_sets)

## Generate networks
network_result=make_networks(gene_sets, indirect_neighbors=1, score_threshold=900, network_file_prefix="network", string_link_file="9606.protein.links.full.v11.5.txt.gz", string_alias_file="9606.protein.aliases.v11.5.txt.gz", 
                         plot_network=TRUE, label_genes=TRUE)

## Find mcl subnetworks
mcl_result=find_mcl(network_result[["largest_networks"]][[2]], min_node=50, expansion = 5, inflation = 10, plot_network = TRUE, label_genes = TRUE)

## Perform ORA
ora_result1=run_ora(V(mcl_result[[1]])$name, pvalueCutoff = 0.05, qvalueCutoff = 0.05, showCategory = 10, plot = FALSE)

ora_result2=run_ora(V(mcl_result[[2]])$name, pvalueCutoff = 0.05, qvalueCutoff = 0.05, showCategory = 10, plot = FALSE)

## Make an upset plot
compare_gene_sets(list(set1=ora_result1, set2=ora_result2), pvalueCutoff = 0.05, qvalueCutoff = 0.05)

## Find page ranks
pageranks=get_pagerank(mcl_result[[1]])
