#' Reduced Heatmap Plot for GO Term Enrichment
#'
#' This function generates a clustered heatmap for GO enrichment results across two sets,
#' optionally using GO IDs or descriptions as column labels.
#'
#' @param combined_processed_ora A data.frame of processed combined ORA results.
#' @param combined_original_ora A data.frame of original combined ORA results.
#' @param set1_processed_ora A data.frame of processed ORA results for set1.
#' @param set2_processed_ora A data.frame of processed ORA results for set2.
#' @param file_name Filename (without extension) for output PDF.
#' @param max_per_set Number of top terms per set to retain (default: 20).
#' @param use_ID Logical; if TRUE, uses GO ID instead of Description (default: FALSE).
#' 
#' @return A matrix of padj_log scores used in the heatmap.
#' @export
#' 
#' 


reduced_heatmap <- function(combined_processed_ora, combined_original_ora, set1_processed_ora, set2_processed_ora, file_name, max_per_set = 20, use_ID = FALSE) {
  require(gplots)
  require(reshape2)
  require(GOSemSim)
  require(tidyr)
  require( rrvgo)
  
  # Join original and parent info
  joined <- right_join(combined_processed_ora%>%mutate(ID=go), combined_original_ora[, c("ID", "source")], by = "ID")
  
  # Get shared parents between two sets
  shared_parents <- joined %>%
    group_by(parent) %>%
    summarise(n = n_distinct(source)) %>%
    filter(n == 2) %>%         # Keep only those that have both 'a' and 'b'
    pull(parent)
  
  # Get 1. shared parents 2. unique parents from each set
  Largest_network_of_set1_set2_ora_common_results=combined_original_ora[combined_original_ora$ID %in% shared_parents, ]
  Largest_network_of_set2_ora_unique_results=combined_original_ora[!combined_original_ora$ID %in% shared_parents 
                                                                   & combined_original_ora$ID %in% set2_processed_ora$parent, ]
  
  Largest_network_of_set1_ora_unique_results=combined_original_ora[!combined_original_ora$ID %in% shared_parents 
                                                                   & combined_original_ora$ID %in% set1_processed_ora$parent, ]
  
  dim(Largest_network_of_set1_set2_ora_common_results)
  dim(Largest_network_of_set2_ora_unique_results)
  dim(Largest_network_of_set1_ora_unique_results)
  
  
  derived_result=rbind(Largest_network_of_set1_set2_ora_common_results, Largest_network_of_set2_ora_unique_results,
                       Largest_network_of_set1_ora_unique_results)
  
  
  heat_data <- derived_result %>%
    mutate(padj_log = -log10(p.adjust)) %>%
    dplyr::select(ID, source, padj_log) %>%
    group_by(source, ID) %>%
    summarise(padj_log = max(padj_log, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = ID, values_from = padj_log)
  
  #heat_data <- derived_result %>%
  #  mutate(padj_log = -log10(p.adjust)) %>%
  #  select(ID, source, padj_log) %>%
  #  pivot_wider(names_from = ID, values_from = padj_log)
  
  #  Convert back to long format for ggplot heatmap
  heat_long <- heat_data %>%
    pivot_longer(-source, names_to = "GO_ID", values_to = "padj_log")
  
  # get top N terms from each result
  top_terms <- heat_long %>%
    group_by(source) %>%
    slice_max(order_by = padj_log, n = max_per_set, with_ties = FALSE) %>%
    ungroup()
  
  
  mat <- reshape2::dcast(top_terms, source ~ GO_ID, value.var = "padj_log")
  
  if(use_ID==FALSE){
    
    
    # Set rownames and remove first column (which is now rownames)
    colnames(mat) <- derived_result$Description[match(colnames(mat), derived_result$ID)]
    rownames(mat) = mat[,1]
    mat <- mat[, -1]
    
    # Convert to numeric matrix (important!)
    mat_new=matrix(as.numeric(unlist(mat)), nrow = nrow(mat), ncol = ncol(mat))
    colnames(mat_new)=colnames(mat)
    rownames(mat_new)=rownames(mat)
    storage.mode(mat_new) <- "numeric"
    mat=mat_new
    
    terms=top_terms$GO_ID
    semData <- GOSemSim::godata('org.Hs.eg.db', ont = "BP")
    flag=terms %in% names(semData@IC)
    sum(flag)
    length(terms)
    terms=terms[flag]
    
    simMatrix <-  rrvgo::calculateSimMatrix(terms, orgdb = org.Hs.eg.db, method = "Wang")
    
    
    # Hierarchical clustering ===
    dist_mat <- as.dist(1 - simMatrix)  # Convert similarity to distance
    
    #dist_mat=as.matrix(dist_mat)
    #colnames(dist_mat)=derived_result$Description[match(colnames(dist_mat), derived_result$ID)]
    #rownames(dist_mat)=colnames(dist_mat)
    
    hc_reps <- hclust(dist_mat, method = "average")
    hc_reps[['labels']]=paste0(derived_result$Description[match(hc_reps[['labels']], derived_result$ID)])
    
    minv=min(mat, na.rm=T)
    mat[is.na(mat)]=-99
    real_range <- range(mat, na.rm = TRUE)
    breaks <- c(minv-1, seq(minv, real_range[2], length.out = 21))
    
    
    pdf(paste0(file_name,".pdf"))
    gp=gplots::heatmap.2(mat, labRow=rownames(mat), cexRow = 0.3, cexCol = 0.5, breaks=breaks, Rowv=NA, Colv=as.dendrogram(hc_reps), col = c("grey",colorRampPalette(c("white", "red"))(20)),
                         density.info="none", trace="none",labCol=hc_reps$labels, margins = c(30,5), srtCol=60)
    dev.off()
    
    print(gp)
    
  }else{# use ID
    
    
    # Set rownames and remove first column (which is now rownames)
    rownames(mat) = mat[,1]
    mat <- mat[, -1]
    
    # Convert to numeric matrix (important!)
    mat_new=matrix(as.numeric(unlist(mat)), nrow = nrow(mat), ncol = ncol(mat))
    colnames(mat_new)=colnames(mat)
    rownames(mat_new)=rownames(mat)
    storage.mode(mat_new) <- "numeric"
    mat=mat_new
    
    terms=top_terms$GO_ID
    semData <- GOSemSim::godata('org.Hs.eg.db', ont = "BP")
    flag=terms %in% names(semData@IC)
    sum(flag)
    length(terms)
    terms=terms[flag]
    
    simMatrix <-  rrvgo::calculateSimMatrix(terms, orgdb = org.Hs.eg.db, method = "Wang")
    
    
    # Hierarchical clustering ===
    dist_mat <- as.dist(1 - simMatrix)  # Convert similarity to distance
    
    #dist_mat=as.matrix(dist_mat)
    #colnames(dist_mat)=derived_result$Description[match(colnames(dist_mat), derived_result$ID)]
    #rownames(dist_mat)=colnames(dist_mat)
    
    hc_reps <- hclust(dist_mat, method = "average")
    hc_reps[['labels']]=paste0(derived_result$ID[match(hc_reps[['labels']], derived_result$ID)])
    
    minv=min(mat, na.rm=T)
    mat[is.na(mat)]=-99
    real_range <- range(mat, na.rm = TRUE)
    breaks <- c(minv-1, seq(minv, real_range[2], length.out = 21))
    
    
    pdf(paste0(file_name,".pdf"))
    gp=gplots::heatmap.2(mat, labRow=rownames(mat), cexRow = 0.3, cexCol = 0.5, breaks=breaks, Rowv=NA, Colv=as.dendrogram(hc_reps), col = c("grey",colorRampPalette(c("white", "red"))(20)),
                         density.info="none", trace="none",labCol=hc_reps$labels, margins = c(30,5), srtCol=60)
    dev.off()
    
    print(gp)
    
  }
  return(mat)
}
