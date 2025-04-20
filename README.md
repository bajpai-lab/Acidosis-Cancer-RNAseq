# Acidosis-Cancer-RNAseq

## Summary of repository
This analysis is broken into three major parts each with its own README.md within with more details:

1. `1-RNAseq-STAR`
    - Processing, aligning, and quantifying of RNA-seq data using STAR and Subread. This stage is best replicated in a Unix environment.
2. `2-DESeq-R`
    - Analysis of RNA-seq count data using the DESeq2 package in R. Includes code for visualizations. Some require the completion of TCGA-R.
3. `3-TCGA-R`
    - Identifying of genes associated with acidosis and either better or worse prognostic outcomes for patients recorded in clinical data available between TCGA and Liu, J., et al., Cell 2018.

## Lists of package versions used
#### R (for 2-DESeq-R and 3-TCGA-R)
```R
> sessionInfo()
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Manjaro Linux

Matrix products: default
BLAS:   /usr/lib/libblas.so.3.12.0 
LAPACK: /usr/lib/liblapack.so.3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Chicago
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gtsummary_2.1.0             readxl_1.4.5                survminer_0.5.0             ggpubr_0.6.0                survival_3.8-3              org.Hs.eg.db_3.20.0        
 [7] AnnotationDbi_1.68.0        lubridate_1.9.4             forcats_1.0.0               purrr_1.0.4                 tidyr_1.3.1                 tibble_3.2.1               
[13] tidyverse_2.0.0             TCGAbiolinks_2.34.1         DT_0.33                     readr_2.1.5                 PCAtools_2.18.0             zeallot_0.1.0              
[19] gridExtra_2.3               pheatmap_1.0.12             eulerr_7.0.2                ggrepel_0.9.6               ggplot2_3.5.2               dplyr_1.1.4                
[25] DESeq2_1.46.0               SummarizedExperiment_1.36.0 Biobase_2.66.0              MatrixGenerics_1.18.1       matrixStats_1.5.0           GenomicRanges_1.58.0       
[31] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0         stringr_1.5.1               tximport_1.34.0            
```
#### Conda (for 1-RNAseq-STAR)
For this process, we use three conda environments for our packages. As an example for determining environment, `/home/lexic/.conda/envs/samtools/bin/samtools` has `envs/samtools` which means it was the `samtools` environment.
- For "pre" see: `conda_env_pre.txt`
- For "star" see: `conda_env_star.txt`
- For "samtools" see: `conda_env_samtools.txt`