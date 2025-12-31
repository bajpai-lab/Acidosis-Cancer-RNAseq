# source('intersect.R')
# Creates l2fc vs l2fc plot with only the subset colored
create_compare_plot_inverse<- function(set1, set2, name1, name2, padj_cutoff, out_path, subset = c(), subset_name = ""){
  print(paste("Creating l2fc plot for",name1,"vs",name2,"(",out_path,")"))
  
  # Create subset of all shared entries in sets 1 and 2
  merged_diff <- merge(set1,
                       set2,
                       by = "gene_name")
  
  fishertest <- with(merged_diff[merged_diff$gene_name %in% subset,],
                     fisher.test(table(padj.x<padj_cutoff,padj.y<padj_cutoff)))
  
  # Silly implicitly-assigned way to do colors
  merged_diff <- merged_diff %>% mutate(hack=case_when(
    ((gene_name %in% subset) & !((padj.y < padj_cutoff) & (padj.x < padj_cutoff))) ~
      paste0("Subset",(padj.y < padj_cutoff),(padj.x < padj_cutoff)),
    ((gene_name %in% subset) & (padj.y < padj_cutoff) & (padj.x < padj_cutoff)) ~
      paste0("Subset",(padj.y < padj_cutoff),(padj.x < padj_cutoff),(sign(log2FoldChange.x) == sign(log2FoldChange.y))),
    !(gene_name %in% subset) ~
      paste0("signotSubset"),
    TRUE ~
      ""
  ))
  
  # Exclude genes insignificant in both sets from plot
  merged_diff_source <- merged_diff[!merged_diff$hack == "SubsetFALSEFALSE",]
  
  
  permutations <- list(
    "5-colors" = c(1:5),
    "both-sig" = c(1,4,5),
    "both-sig-same-dir" = c(1,5)
  )
  
  for (permutation in names(permutations)){
    # Create text labels for gene groups
    allowed_groups <- levels(factor(merged_diff_source$hack))[permutations[[permutation]]]
    merged_diff <- merged_diff_source[merged_diff_source$hack %in% allowed_groups,]
    related_labels = c(
      "signotSubset"=paste0("padj < 0.1 either set. NON-", subset_name),
      "SubsetFALSETRUE"=paste0("padj < 0.1 only (",name1,") and ",subset_name," [", length(merged_diff$hack[merged_diff$hack == "SubsetFALSETRUE"]), " genes]"),
      "SubsetTRUEFALSE"=paste0("padj < 0.1 only (",name2,") and ",subset_name," [", length(merged_diff$hack[merged_diff$hack == "SubsetTRUEFALSE"]), " genes]"),
      "SubsetTRUETRUEFALSE"=paste0("padj < 0.1 both sets. DIFFERENT DIRECTION and ",subset_name," [", length(merged_diff$hack[merged_diff$hack == "SubsetTRUETRUEFALSE"]), " genes]"),
      "SubsetTRUETRUETRUE"=paste0("padj < 0.1 both sets. SAME DIRECTION and ",subset_name," [", length(merged_diff$hack[merged_diff$hack == "SubsetTRUETRUETRUE"]), " genes]")
    )
    related_labels = related_labels[permutations[[permutation]]]
    
    # Create colors for gene groups (ordered with labels set)
    related_outline_colors = c(
      "grey",
      "#d25bf5",
      "lightgreen",
      "red",
      "#4ec5f7")
    related_outline_colors = related_outline_colors[permutations[[permutation]]]
    related_fill_colors = c(
      "grey",
      "#d25bf5",
      "lightgreen",
      "red",
      "#4ec5f7")
    related_fill_colors = related_fill_colors[permutations[[permutation]]]
    
    # Create plot
    merged_diff <- merged_diff[order(merged_diff$hack), ]
    p <- ggplot(merged_diff,aes(x=log2FoldChange.x, y=log2FoldChange.y, color = hack, fill = hack)) +
      
      # Plot labels
      labs(title = paste0('log2 Fold Change comparison (',permutation,')'),
           subtitle = paste0('"', name1, '" vs "', name2, '"')) +
      xlab(paste0('log2FC "', name1, '"')) +
      ylab(paste0('log2FC "', name2, '"')) +
      
      geom_point(shape=21, size=1.5, alpha =0.4) +
      
      # (Pink, green, orange, red) color scheme  
      scale_fill_manual(values = related_fill_colors,
                        labels = related_labels)+  
      scale_color_manual(values = related_outline_colors,
                         labels = related_labels)+
      theme_classic()+
      theme(legend.title=element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            legend.title.position = "top",
            legend.text.position = "right",
            legend.text = element_text(size=6, hjust = 0, vjust = 0.5, angle = 0)
      )
    
    # Write out plot
    ggsave(plot=p, filename=(out_path %>% str_replace(".pdf",paste0("_",permutation,".pdf"))), width=6, height=6)
    print(paste("Exported",name1,"vs",name2,"(",(out_path %>% str_replace(".pdf",paste0("_",permutation,".pdf"))),")"))
  }

  
  fisher_report <- paste(fishertest$method,"\n",
        "odds ratio: ",format(fishertest$estimate, digits=4, scientific=FALSE),"\n",
        "reported P value: ",format(fishertest$p.value, digits=4, scientific=TRUE))
  
  write_file(fisher_report,paste0(out_path %>% str_replace(".pdf",""),"_fisher_report.txt"))
  
  
  print(paste0("Significant genes: ", dim(merged_diff[merged_diff$padj.x < padj_cutoff & merged_diff$padj.y < padj_cutoff,])[1]))
  print(paste0("Sig both pos:", dim(merged_diff[merged_diff$padj.x < padj_cutoff & merged_diff$padj.y < padj_cutoff & merged_diff$log2FoldChange.x > 0 & merged_diff$log2FoldChange.y > 0,])[1]))
  print(paste0("Sig both neg:", dim(merged_diff[merged_diff$padj.x < padj_cutoff & merged_diff$padj.y < padj_cutoff & merged_diff$log2FoldChange.x < 0 & merged_diff$log2FoldChange.y < 0,])[1]))
}