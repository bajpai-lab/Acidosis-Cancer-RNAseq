# Helper function to merge only two datasets with specified suffixes
merge_intersecting_genes_two_sets <- function(set1, set2, name1, name2){
  # Handling suffixes explicitly since I can't get merge() to do it consistently
  names(set1) <- ifelse(names(set1)=="gene_name", names(set1), str_c(names(set1),name1))
  names(set2) <- ifelse(names(set2)=="gene_name", names(set2), str_c(names(set2),name2))
  
  # Merging the two datasets by gene_name entries
  # Default behavior for merge() is filtering out non-matches
  merged <- merge(set1,
                  set2,
                  by = "gene_name")
  return(merged)
}

# Main function to merge arbitrary numbers of datasets from a list
# Using list entry names for suffixes.
merge_intersecting_genes <- function(results_list){
  if (length(results_list) < 2){
    return (FALSE)
  }
  
  outcome <- results_list[[1]]
  names(outcome) <- ifelse(names(outcome)=="gene_name", names(outcome), str_c(names(outcome),paste0(".",names(results_list)[1])))
  
  for (name in names(results_list[-c(1)])){
    outcome <- merge_intersecting_genes_two_sets(outcome, results_list[[name]], "", paste0(".",name))
  }
  
  return (outcome)
}

fuzzy_evaluate <- function(data, common_prefix, cutoff, operator=c("<",">")){
  operator <- match.arg(operator)
  
  if (operator == "<"){
    data[apply(data[startsWith(colnames(data),common_prefix)], 1, max) < cutoff,]
  } else if (operator == ">"){
    data[apply(data[startsWith(colnames(data),common_prefix)], 1, min) > cutoff,]
  }
}

directional_intersect <- function(data, cutoff_FDR = 0.1){
  inter <- merge_intersecting_genes(data)
  inter_UP <- fuzzy_evaluate(data = inter, common_prefix = "log2FoldChange", operator = ">", cutoff = 0)
  inter_DOWN <- fuzzy_evaluate(data = inter, common_prefix = "log2FoldChange", operator = "<", cutoff = 0)
  sig_inter <- fuzzy_evaluate(data = inter, common_prefix = "padj", operator = "<", cutoff = cutoff_FDR)
  sig_inter_UP <- fuzzy_evaluate(data = sig_inter, common_prefix = "log2FoldChange", operator = ">", cutoff = 0)
  sig_inter_DOWN <- fuzzy_evaluate(data = sig_inter, common_prefix = "log2FoldChange", operator = "<", cutoff = 0)
  sig_inter_BOTH <- rbind(sig_inter_UP,sig_inter_DOWN)
  
  return (list(
    intersect = inter,
    intersect_UP = inter_UP,
    intersect_DOWN = inter_DOWN,
    sig = sig_inter,
    sig_UP = sig_inter_UP,
    sig_DOWN = sig_inter_DOWN,
    sig_BOTH = sig_inter_BOTH
  ))
}