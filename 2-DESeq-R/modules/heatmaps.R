plot_heatmap <- function(INPUT,
                         TAG,
                         THE_PALETTE="DEFAULT",
                         THE_TITLE,
                         GLOBAL_MIN,
                         GLOBAL_MAX,
                         FONT_SIZE=7,
                         TOP_GENES=0,
                         TOP_GENES_LOW=0,
                         TOP_GENES_HIGH=0,
                         KEEP_ENSG=TRUE,
                         NUMBERS=FALSE,
                         FIXED_ROWS=c(),
                         CELL_VALUE_UNIT="log2FoldChange",
                         CELL_HEIGHT=NA,
                         CELL_WIDTH=NA,
                         TRIM_LEGEND=FALSE,
                         SHOW_LEGEND=TRUE){
  print(paste("Evaulating (",TOP_GENES_LOW,",",TOP_GENES,",",TOP_GENES_HIGH,") heat map for",TAG))
  print(dim(INPUT))
  write.csv(INPUT, file=paste0("results/heatmap_hits_",TAG,".csv"))
  
  # Input should be a matrix of gene_name and arbitrary number of l2fcs from results files
  if (!dim(INPUT)[1]>1){
    print(paste("Cannot create heat map for",TAG,"due to insufficient entries"))
    return(FALSE)
  }
  
  # Remove ENSG ID "gene names"
  if (!KEEP_ENSG){
    INPUT <- INPUT[!startsWith(INPUT$gene_name,"ENSG"),]
  }
  
  # Set row names for heatmap
  rownames(INPUT) <- INPUT$gene_name
  
  # Create an ordered list of entries based on average L2FC
  SELECT <- order(rowMeans(INPUT[!colnames(INPUT) %in% c("gene_name")]),
                  decreasing=TRUE)
  
  # Top # of genes
  if (TOP_GENES > 0){
    TOP_GENES_LOW <- TOP_GENES
    TOP_GENES_HIGH <- TOP_GENES
  }
  if (dim(INPUT)[1] >= TOP_GENES_LOW + TOP_GENES_HIGH){
    SELECT_LOW <- c()
    SELECT_HIGH <- c()
    if (TOP_GENES_LOW > 0){
      SELECT_LOW <- SELECT[(length(SELECT)-TOP_GENES_LOW+1):length(SELECT)]
    }
    if (TOP_GENES_HIGH > 0){
      SELECT_HIGH <- SELECT[1:TOP_GENES_HIGH]
    }
    SELECT <- c(SELECT_HIGH, SELECT_LOW)
  }
  
  # Finding range of current dataset
  INPUT_MIN <- min(INPUT[!colnames(INPUT) %in% c("gene_name")])
  INPUT_MAX <- max(INPUT[!colnames(INPUT) %in% c("gene_name")])
  
  print(paste0("INPUT_MIN: ", INPUT_MIN, "  INPUT_MAX: ",INPUT_MAX))

  # Setting undefined variables
  if (missing(GLOBAL_MIN)){
    print(paste0("No global min set, using ", INPUT_MIN))
    GLOBAL_MIN <- INPUT_MIN
  }
  if (missing(GLOBAL_MAX)){
    print(paste0("No global max set, using ", INPUT_MAX))
    GLOBAL_MAX <- INPUT_MAX
  }
  if (missing(FONT_SIZE)){
    print(paste0("Font size set, using 1"))
    FONT_SIZE <- 1
  }
  
  # Used to index the colors before providing to pheatmap
  #color_range <- seq(INPUT_MIN,INPUT_MAX,length.out=100)#dim(INPUT)[1])
  color_range <- seq(GLOBAL_MIN,GLOBAL_MAX,length.out=100)
  # Defining a (min=RED, 0=WHITE, max=BLUE) palette
  if (THE_PALETTE == "DEFAULT"){
    print("----PALETTE IS DEFAULT")
    THE_PALETTE<-c("#721d1a","orange","#3ada35","#5ebfe0","#404597")
    if (GLOBAL_MIN < 0 & GLOBAL_MAX > 0){
      print("----PALETTE CAN 0 SPLIT AS:")
      breaks_range <- c(seq(GLOBAL_MIN,
                            0,
                            length.out=3),
                        seq(0,
                            GLOBAL_MAX,
                            length.out=3)[-1])
      print(breaks_range)
    }else if (GLOBAL_MAX <= 0){
      breaks_range <- c(seq(GLOBAL_MIN,
                            0,
                            length.out=5))
    }else if (GLOBAL_MIN >= 0){
      breaks_range <- c(seq(0,
                            GLOBAL_MAX,
                            length.out=5))
    }
  }else{
    breaks_range <- seq(INPUT_MIN,INPUT_MAX,length.out=length(THE_PALETTE))
  }
  # Initialize using SELECT of INPUT
  final_dataset <- INPUT[SELECT,]
  # Append any fixed rows without duplicating existing entries
  final_dataset <- rbind(final_dataset,INPUT[INPUT$gene_name %in% c(FIXED_ROWS) & !(INPUT$gene_name %in% final_dataset$gene_name),])
  # Remove "gene_name" column regardless of its index in the matrix
  final_dataset <- final_dataset[,!colnames(INPUT) %in% c("gene_name")]
  
  SUBTITLE <- paste0("Stats: ",dim(final_dataset)[1], " rows, ",
                     dim(final_dataset)[2]," cols, ",
                     format(INPUT_MIN, digits=3, scientific=FALSE), "min, ",
                     format(INPUT_MAX, digits=3, scientific=FALSE), "max")
  
  
  p <- pheatmap(final_dataset, cluster_rows=FALSE, show_rownames=TRUE,
                labels_row = rownames(final_dataset),
                fontsize_row = FONT_SIZE,
                breaks = seq(GLOBAL_MIN,GLOBAL_MAX,length.out=101),
                color=circlize::colorRamp2(breaks_range, THE_PALETTE)(color_range),
                legend_breaks = round(seq(GLOBAL_MIN,GLOBAL_MAX,length.out=5),digits=2),
                border_color = NA,
                main = "PLACEHOLDER_TITLE",
                angle_col = "0",
                cluster_cols=FALSE,
                display_numbers=NUMBERS,
                fontsize_number=FONT_SIZE,
                number_color="white",
                cellheight = CELL_HEIGHT,
                cellwidth = CELL_WIDTH,
                drop_levels = TRIM_LEGEND,
                legend = SHOW_LEGEND)
  
  p$gtable$grobs[[1]] <- grid::textGrob(paste0(
    THE_TITLE,"\n",SUBTITLE
  ), gp = gpar(fontface="bold"))
  
  ggsave(filename=paste0("results/heatmap_",TAG,".pdf"),
         p$gtable,
         width=(2+dim(final_dataset)[2]*CELL_WIDTH*0.015),
         height=(1+dim(final_dataset)[1]*CELL_HEIGHT*0.01375))
  #ggsave(filename=paste0("results/heatmap_",TAG,".pdf"), p, width=6, height=4)
  return(p)
}