# source('intersect.R')
create_compare_plot<- function(set1, set2, name1, name2, padj_cutoff, out_path, subset_up = c(), subset_down = c()){
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
                                                      ((padj.x < padj_cutoff) & (padj.y < padj_cutoff) & (sign(log2FoldChange.x) == sign(log2FoldChange.y)) & (gene_name %in% subset_down)) ~ "bothSUBDOWN",
                                                      ((padj.x < padj_cutoff) & (padj.y < padj_cutoff) & (sign(log2FoldChange.x) == sign(log2FoldChange.y)) & (gene_name %in% subset_up)) ~ "bothSUBUP",
                                                      ((padj.x < padj_cutoff) & (padj.y < padj_cutoff) & (sign(log2FoldChange.x) == sign(log2FoldChange.y))) ~ "both",
                                                      ((padj.x < padj_cutoff) & (padj.y < padj_cutoff) & (sign(log2FoldChange.x) != sign(log2FoldChange.y))) ~ "diff",
                                                      TRUE ~ ""
                                                    )))
  
  # Exclude genes insignificant in both sets from plot
  merged_diff <- merged_diff[!merged_diff$hack == "FALSEFALSE",]
  
  # Create text labels for gene groups
  related_labels = c(
    "FALSETRUE"=paste0("padj < 0.1 only (",name2,")"),
    "TRUEFALSE"=paste0("padj < 0.1 only (",name1,")"),
    "TRUETRUEboth"=paste0("padj < 0.1 (",name1,") and (",name2,") SAME DIRECTION [", length(merged_diff$hack[merged_diff$hack == "TRUETRUEboth" | merged_diff$hack == "TRUETRUEbothSUBUP" | merged_diff$hack == "TRUETRUEbothSUBDOWN"]), " genes]"),
    "TRUETRUEbothSUBDOWN"=paste0("padj < 0.1 (",name1,") and (",name2,") SAME DIRECTION and WORSE OS [", length(merged_diff$hack[merged_diff$hack == "TRUETRUEbothSUBDOWN"]), " genes]"),
    "TRUETRUEbothSUBUP"=paste0("padj < 0.1 (",name1,") and (",name2,") SAME DIRECTION and BETTER OS [", length(merged_diff$hack[merged_diff$hack == "TRUETRUEbothSUBUP"]), " genes]"),
    "TRUETRUEdiff"=paste0("padj < 0.1 (",name1,") and (",name2,") [", length(merged_diff$hack[merged_diff$hack == "TRUETRUEdiff"]), " genes]")
    )
  
  # Create colors for gene groups (ordered with labels set)
  related_outline_colors = c(
    "lightgreen",
    "#d25bf5",
    "#4ec5f7",
    "#4ec5f7",
    "#4ec5f7",
    "red")
  related_fill_colors = c(
    "lightgreen",
    "#d25bf5",
    "#4ec5f7",
    "black",
    "black",
    "red")
  
  # Create plot
  merged_diff <- merged_diff[order(merged_diff$hack), ]
  p <- ggplot(merged_diff,aes(x=log2FoldChange.x, y=log2FoldChange.y, color = hack, fill = hack)) +
    
    # Plot labels
    labs(title = paste0('log2 Fold Change comparison'),
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