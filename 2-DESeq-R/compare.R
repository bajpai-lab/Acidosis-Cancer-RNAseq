library(tximport)
library(stringr)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(eulerr)
library(pheatmap)
library(grid)
library(gridExtra)
library(zeallot)
library(PCAtools)
library(readr)

##### Configurations
##### Configurations
##### Configurations
##### Configurations
{
  # FDR cutoff for significant results
  cutoff_FDR=0.1
}

##### DDS creation
##### DDS creation
##### DDS creation
##### DDS creation
{
  source ('modules/run_DESeq.R')
  
  # MDA-MB-231 24h + 10w
  {
    dds_mda_24h10w_samples <- get_dds(group_name = "mda_24h10w_samples",
                                      combined_categories = list(
                                        c("duration","condition")
                                      ),
                                   design = ~ combined_1,
                                   samples_exclude = c(7:12,21:44), # MDA 23d, Corbet 6w, and SUM-159
                                   test_name = "test_mda_24h",
                                   test_contrast = c("combined_1", "24h_acidosis", "24h_normal"))
  }
}

##### Calling results and significance
##### Calling results and significance
##### Calling results and significance
##### Calling results and significance
{
  # Requires sourcing run_DESeq.R for get_results()
  
  # Calling results from dds_mda_samples
  {
    res <- list(
      # dds_mda_24h10w_samples: ~ (duration)_(condition)
      mda24h10w_24hAvN = get_results(dds = dds_mda_24h10w_samples, contrast = c("combined_1", "24h_acidosis", "24h_normal"), group_name = "mda_24h10w_samples", results_name = "mda24h10w_24h_AvN"),
      mda24h10w_10wAvN = get_results(dds = dds_mda_24h10w_samples, contrast = c("combined_1", "10w_acidosis", "10w_normal"), group_name = "mda_24h10w_samples", results_name = "mda24h10w_10w_AvN")
    )
    
    # Building FDR-cutoff-filters of results
    {
      sig_res <- list()
      for (entry_name in names(res)){
        print(paste0("Dims pre- and post-FDR cutoff for ", entry_name,":"))
        sig_res[[entry_name]] <- res[[entry_name]] %>% filter(padj < cutoff_FDR)
        print(dim(res[[entry_name]]))
        print(dim(sig_res[[entry_name]]))
      }
    }
  }
}


##### Matching to prog genes
##### Matching to prog genes
##### Matching to prog genes
##### Matching to prog genes
{
  # Loading results from new TCGA analysis
  prog_BETTER = read.csv('../3-TCGA-R/PlanA/gene-list_planA_BOTHfcUP_progBETTER.csv')
  prog_WORSE = read.csv('../3-TCGA-R/PlanA/gene-list_planA_BOTHfcDOWN_progWORSE.csv')
}

##### Generating intersects
##### Generating intersects
##### Generating intersects
##### Generating intersects
{
  source ('modules/intersect.R')
  # MDA-only comparisons
  {
    # "(MDA24h10w 24h) U (MDA24h10w 10w)"
    {
      planA <- directional_intersect(list(
        MDA24h10w_24h=res$mda24h10w_24hAvN,
        MDA24h10w_10w=res$mda24h10w_10wAvN
      ))
      write.csv(planA$sig_BOTH, file = "results/intersect_PlanA_(10w)U(24h)_BOTH-FC-by-acidosis_0.1FDR.csv")
      write.csv(planA$sig, file = "results/intersect_PlanA_(10w)U(24h)_ANY-FC-by-acidosis_0.1FDR.csv")
    }
  }
}

##### Venn/Euler diagrams
##### Venn/Euler diagrams
##### Venn/Euler diagrams
##### Venn/Euler diagrams
{
  source ('modules/euler_diagram.R')
  
  # Plan A SUM MDA
  {
    euler_plot <- create_eulerr(
      sig_res$mda24h10w_24hAvN$gene_name, "Short (MDA 24h)",
      sig_res$mda24h10w_10wAvN$gene_name, "Long (MDA 10w)",
      planA$sig_BOTH$gene_name, "Common BOTH"
    )
    euler_plot <- grid.arrange(grobs = list(euler_plot), top = "(STAR) Plan A - Overlap of genes altered by acidosis (FDR < 0.1)")
    ggsave(euler_plot, filename="results/euler_diagram_planA.pdf", width=6, height=5)
  }
}

#### Comparing results with l2fc vs l2fc plots
#### Comparing results with l2fc vs l2fc plots
#### Comparing results with l2fc vs l2fc plots
#### Comparing results with l2fc vs l2fc plots
{
  source ('modules/l2fc.R')
  
  create_compare_plot(set1 = res$mda24h10w_24hAvN,
                      set2 = res$mda24h10w_10wAvN,
                      name1 = "(pHe 6.4 / 7.4) MDA-MB-231 24 h",
                      name2 = "(pHe 6.4 / 7.4) MDA-MB-231 10 w",
                      padj_cutoff = cutoff_FDR,
                      out_path = "results/l2fc_planA_MDA24h10w.pdf")
}

#### Selective Heatmaps
#### Selective Heatmaps
#### Selective Heatmaps
#### Selective Heatmaps
{
  source ('modules/heatmaps.R')
  
  ## Larger aggregation of genes from different pathway databases
  {
    LA_cc = read.csv('path_groups/sort-cc.csv',header = FALSE)$V1 # Cell cycle genes
    LA_hip = read.csv('path_groups/sort-HIPPO.csv',header = FALSE)$V1 # HIPPO genes
    LA_ILK = read.csv('path_groups/sort-ILK.csv',header = FALSE)$V1 # ILK signalling genes
    LA_int = read.csv('path_groups/sort-integrin.csv',header = FALSE)$V1 # Integrin-cytoskeleton signaling genes
    LA_PI3K = read.csv('path_groups/sort-PI3K.csv',header = FALSE)$V1 # PI3K signalling genes
    LA_glyc = read.csv('path_groups/sort-glycolysis.csv',header = FALSE)$V1 # Glycolysis-related genes
  }
  
  # Function for plotting HIPPO, HIPPO UP, ILK+, ILK+ DOWN
  # assumes hmsource, fullsource, progs, pathways, cell_height/width, global min/max are defined
  get_hm = function(pathways=NA, prog=NA, tag, title=NA, TOP_GENES=30, SHOW_LEGEND=TRUE, FIXED_ROWS=c()){
    path_filter <- fullsource$gene_name #placeholder that returns full set
    prog_filter <- fullsource$gene_name #placeholder that returns full set
    TOP_GENES_HIGH <- 0
    TOP_GENES_LOW <- 0
    
    if (!is.na(pathways)){
      switch(pathways,
             HIPPO={
               path_filter <- c(LA_hip)
               TOP_GENES_HIGH <- TOP_GENES
               print("Added HIPPO filter")
             },
             ILKmany={
               path_filter <- c(LA_ILK,LA_cc,LA_int,LA_PI3K,LA_glyc)
               TOP_GENES_LOW <- TOP_GENES
               print("Added ILK filter")},
             {stop("No matching path_filter preset")}) #fails aggressively to encourage parameter checking
    }
    if (!is.na(prog)){
      switch(prog,
             UP={
               prog_filter <- prog_BETTER$gene_name
               print("Added TCGA BETTER filter")
             },
             DOWN={
               prog_filter <- prog_WORSE$gene_name
               print("Added TCGA WORSE filter")
             },
             {stop("No matching prog_filter preset")}) #fails aggressively to encourage parameter checking
    }
    if (is.na(title)){
      title <- paste0(tag,
                      (if (!is.na(pathways)) paste0("-U-(",pathways,")") else ""),
                      (if (!is.na(prog)) paste0("-U-(TCGA ",prog,")") else ""))
    }
    if (TOP_GENES > 0 && is.na(TOP_GENES_HIGH) && is.na(TOP_GENES_LOW)){
      TOP_GENES_HIGH <- TOP_GENES
      TOP_GENES_LOW <- TOP_GENES
    }
    
    print(paste0("For ",tag," we have (", dim(hmsource)[1], ",", dim(hmsource)[2],")"))
    hm <- hmsource[hmsource$gene_name %in% path_filter
                   & hmsource$gene_name %in% prog_filter,]
    
    print(paste0("Reduced to (", dim(hm)[1], ",", dim(hm)[2],")"))
    print(head(hm))
    full <- fullsource[fullsource$gene_name %in% path_filter
                       & fullsource$gene_name %in% prog_filter,]
    
    write.csv(full, file = paste0("results/heatmapEXTENDED_",
                                  tag,
                                  (if (!is.na(pathways)) paste0("_",pathways) else ""),
                                  (if (!is.na(prog)) paste0("_TCGA-",prog) else ""),
                                  "_anyq_anyFC.csv"))
    print(paste0("  Plotting top ", TOP_GENES))
    pTOP <- plot_heatmap(INPUT = hm,
                 TAG = paste0(tag,
                              (if (!is.na(pathways)) paste0("_",pathways) else ""),
                              (if (!is.na(prog)) paste0("_TCGA-",prog) else "")),
                 THE_TITLE = title,
                 FONT_SIZE = FONT_SIZE,
                 TOP_GENES_HIGH = TOP_GENES_HIGH,
                 TOP_GENES_LOW = TOP_GENES_LOW,
                 FIXED_ROWS=FIXED_ROWS,
                 GLOBAL_MIN = hmmin,
                 GLOBAL_MAX = hmmax,
                 CELL_HEIGHT = cell_height,
                 CELL_WIDTH = cell_width,
                 SHOW_LEGEND = SHOW_LEGEND)
    print(paste0("  Plotting all genes"))
    pALL <- plot_heatmap(INPUT = hm,
                 TAG = paste0(tag,
                              (if (!is.na(pathways)) paste0("_",pathways) else ""),
                              (if (!is.na(prog)) paste0("_TCGA-",prog) else ""),
                              "_ALL-GENES"),
                 THE_TITLE = paste0(title," ALL GENES"),
                 FONT_SIZE = FONT_SIZE,
                 TOP_GENES_HIGH = TOP_GENES_HIGH*10000, # TODO: Fix plot_heatmap() failure on default TOP_GENE values
                 TOP_GENES_LOW = TOP_GENES_LOW*10000,
                 GLOBAL_MIN = hmmin,
                 GLOBAL_MAX = hmmax,
                 CELL_HEIGHT = cell_height,
                 CELL_WIDTH = cell_width,
                 SHOW_LEGEND = SHOW_LEGEND)
    # return (list(pTOP=pTOP,pALL=pALL))
  }
  
  hmmin = -1.1
  hmmax = 1.1
  cell_height = 4
  cell_width = 70
  FONT_SIZE = 4
  FIXED_ROWS_HIPPO = c()
  FIXED_ROWS_ILK = c("SLC9A1")
  
  # (24h)U(10w)
  {
    fullsource <- planA$intersect
    tag <- "24h-U-10w"
    
    hmsource <- planA$sig_UP[c(1,3,9)]
    colnames(hmsource) <- c("gene_name", "MDA 24h", "MDA 10w")
    get_hm(pathways="HIPPO", tag=tag, FIXED_ROWS=FIXED_ROWS_HIPPO)
    get_hm(pathways="HIPPO",prog="UP", tag=tag, FIXED_ROWS=FIXED_ROWS_HIPPO)
    
    hmsource <- planA$sig_DOWN[c(1,3,9)]
    colnames(hmsource) <- c("gene_name", "MDA 24h", "MDA 10w")
    get_hm(pathways="ILKmany", tag=tag, FIXED_ROWS=FIXED_ROWS_ILK)
    get_hm(pathways="ILKmany",prog="DOWN", tag=tag, FIXED_ROWS=FIXED_ROWS_ILK)
  }
}

#### Pathway Analysis
#### Pathway Analysis
#### Pathway Analysis
#### Pathway Analysis

{
  source("modules/pathways.R")
  
  # "(MDA 24h) U (MDA 10w)"
  # IPA_PATHWAYS_STAR_MDA24h_vs10W.csv
  toppathwaysIPA(input = "GO/IPA_PATHWAYS_STAR_MDA24h_vs10W.csv", tag = "((MDA 24h) U (MDA 10w)) BOTH Common", limit = 50, exclude_from_top_by_index = c())
  toppathwaysIPA(input = "GO/IPA_PATHWAYS_STAR_MDA24h_vs10W.csv", tag = "((MDA 24h) U (MDA 10w)) BOTH Common REDUCED", limit = 33, exclude_from_top_by_index = c(3,5,6,11,12,19,20,23:26,28:30,42,46,48),height=7)
  toppathwaysIPA(input = "GO/IPA_PATHWAYS_STAR_MDA24h_vs10W.csv", tag = "((MDA 24h) U (MDA 10w)) BOTH Common TOP20plus", limit = 26, exclude_from_top_by_index = c(21:26,28:32,34:44,46,48),height=6)
  toppathwaysIPA(input = "GO/IPA_PATHWAYS_STAR_MDA24h_vs10W.csv", tag = "((MDA 24h) U (MDA 10w)) BOTH Common TOP10plus", limit = 17, exclude_from_top_by_index = c(11:13,15:26,28:32,34:44,46,48),height=6)
}


#### Trick for nice rasters of pdfs:

# ls *.pdf | parallel "convert -verbose -density 300 -trim {} -quality 100 -flatten -sharpen 0x1.0 {.}_render.png"
