# source('intersect.R')
create_compare_plot<- function(set1, set2, name1, name2, padj_cutoff, out_path){
  print(paste("Creating l2fc plot for",name1,"vs",name2,"(",out_path,")"))
  
  # Create subset of all shared entries in sets 1 and 2
  merged_diff <- merge(set1,
                       set2,
                       by = "gene_name")
  
  padjtable <- with(merge_intersecting_genes_two_sets(set1,set2, ".1", ".2"),
                 table(padj.1<padj_cutoff,padj.2<padj_cutoff))
  
  fishertest <- with(merge_intersecting_genes_two_sets(set1,set2, ".1", ".2"),
                     fisher.test(table(padj.1<padj_cutoff,padj.2<padj_cutoff)))
  
  # Silly implicitly-assigned way to do colors
  merged_diff <- merged_diff %>% mutate(hack=paste0((padj.x < padj_cutoff),
                                                    (padj.y < padj_cutoff),
                                                    case_when(
                                                      ((padj.x < padj_cutoff) & (padj.y < padj_cutoff) & (sign(log2FoldChange.x) == sign(log2FoldChange.y))) ~ "both",
                                                      ((padj.x < padj_cutoff) & (padj.y < padj_cutoff) & (sign(log2FoldChange.x) != sign(log2FoldChange.y))) ~ "diff",
                                                      TRUE ~ ""
                                                    )))
  merged_diff <- merged_diff[!merged_diff$hack == "FALSEFALSE",]
  related_labels = c(#"FALSEFALSE"=paste0("padj > 0.1 (",name1,") and (",name2,") [", length(merged_diff$hack[merged_diff$hack == "FALSEFALSE"]), " genes]"),
    #"FALSETRUE"=paste0("padj < 0.1 only (",name2,") [", length(merged_diff$hack[merged_diff$hack == "FALSETRUE"]), " genes]"),
    #"TRUEFALSE"=paste0("padj < 0.1 only (",name1,") [", length(merged_diff$hack[merged_diff$hack == "TRUEFALSE"]), " genes]"),
    "FALSETRUE"=paste0("padj < 0.1 only (",name2,")"),
    "TRUEFALSE"=paste0("padj < 0.1 only (",name1,")"),
    "TRUETRUEboth"=paste0("padj < 0.1 (",name1,") and (",name2,") SAME DIRECTION [", length(merged_diff$hack[merged_diff$hack == "TRUETRUEboth"]), " genes]"),
    "TRUETRUEdiff"=paste0("padj < 0.1 (",name1,") and (",name2,") [", length(merged_diff$hack[merged_diff$hack == "TRUETRUEdiff"]), " genes]"))
  related_colors = c(#"lightpink",
    "lightgreen",
    "#d25bf5",
    "#4ec5f7",
    "red")
  
  # Create plot
  p <- ggplot(merged_diff,aes(x=log2FoldChange.x, y=log2FoldChange.y, color = hack, fill = hack)) +
    
    # Plot labels
    labs(title = paste0('log2 Fold Change comparison'),
         subtitle = paste0('"', name1, '" vs "', name2, '"')) +
    xlab(paste0('log2FC "', name1, '"')) +
    ylab(paste0('log2FC "', name2, '"')) +
  
    geom_point(shape=21, size=1.5, alpha =0.4) +
  
    # (Pink, green, orange, red) color scheme  
    scale_fill_manual(values = related_colors,
                      labels = related_labels)+  
    scale_color_manual(values = related_colors,
                      labels = related_labels)+
    theme_classic()+
    theme(legend.title=element_blank(),
          legend.position = "top",
          legend.direction = "vertical",
          legend.title.position = "top",
          legend.text.position = "right",
          legend.text = element_text(size=7, hjust = 0, vjust = 0.5, angle = 0)
    )

  
  fisher_report <- paste(fishertest$method,"\n",
        "odds ratio: ",format(fishertest$estimate, digits=4, scientific=FALSE),"\n",
        "reported P value: ",format(fishertest$p.value, digits=4, scientific=TRUE))
  
  write_file(fisher_report,paste0(out_path %>% str_replace(".pdf",""),"_fisher_report.txt"))
  
  
  print(paste0("Significant genes: ", dim(merged_diff[merged_diff$padj.x < padj_cutoff & merged_diff$padj.y < padj_cutoff,])[1]))
  print(paste0("Sig both pos:", dim(merged_diff[merged_diff$padj.x < padj_cutoff & merged_diff$padj.y < padj_cutoff & merged_diff$log2FoldChange.x > 0 & merged_diff$log2FoldChange.y > 0,])[1]))
  print(paste0("Sig both neg:", dim(merged_diff[merged_diff$padj.x < padj_cutoff & merged_diff$padj.y < padj_cutoff & merged_diff$log2FoldChange.x < 0 & merged_diff$log2FoldChange.y < 0,])[1]))
  
  # Write out and return plot
  ggsave(plot=p, filename=out_path, width=6, height=6)
  print(paste("Exported",name1,"vs",name2,"(",out_path,")"))
  return(p)
}