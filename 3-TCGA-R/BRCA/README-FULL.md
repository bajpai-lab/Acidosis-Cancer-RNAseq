# TNBC acidosis RNA-seq and patients' survival analysis
The archived patients' sequencing data and their corresponding clinical outcomes enable us to examine the in vivo relevance of acidosis-activated and repressed genes defined in the MDA-MB-231 cells (triple-negative breast cancer, or TNBC, cell line). We further isolated a group of acidosis-induced genes with greater prognostic significance based on their predicted hazardous ratio (HR) on TNBC patients' overall survival (OS). Patients' sequence data and clinical records are extracted from The Cancer Genome Atlas (TCGA) via the TCGAbiolinks (version 2.22.4) package. The patients’ clinical outcome records are updated in Liu, J. et al., Cell 2018. Breast cancer patients (BRCA) were selected. We select patients with the BC Basal subtype for TNBC. </br>
## RNA-seq analysis
First, we will generate significant differentially expressed genes (DEGs) for the cells cultured under acidic condition compared to the normal PH environment. We normalize the raw RNA sequencing reads, generated from nf-co pipeline, and perform principal component analysis (PCA) using the DSEq2 package. (Code: RNA_seq_analysis.R)
```
#Iinstall/load necessary packages, including DESeq2:
rm(list=ls()) 

#BiocManager::install("DESeq2")

library( "DESeq2" )
library( "dplyr" )
library( "tibble" )
library( "stringr" )
library("readr")
library( "ggplot2" )
library( "ggrepel" )
library("org.Hs.eg.db")
setwd("~/TNBC_acidosis_RNA_seq_TCGA_analysis")

```
Then, we read in the raw read file, located in the parent directory, and get rid of reads related to PH = 6.7: 
```
# Read in file, define rownames, columnnames and conditions
setwd("~/TNBC_acidosis_RNA_seq_TCGA_analysis")
df2 <- read.delim( "raw_reads_24hrs.txt")
colnames( df2 ) <- str_replace( colnames(df2), ".sorted.bam", "" )
rownames( df2 ) <- df2$Geneid
df2 <- df2[,-c(1,2)]
## Identify the conditions
df_samples2 <- data.frame(
  condition = str_split( colnames(df2), "[_]", simplify=TRUE )[,1],
  row.names = colnames(df2)
)
## Getting rid of the outlines based on PCA plots
df_samples <- df_samples2[-c(2,5,8,11),,drop=FALSE]
df <- df2[,rownames(df_samples)]
rm(df_samples2)
rm(df2)


df_samples <- data.frame(
  condition = str_split( colnames(df), "[_]", simplify=TRUE )[,1],
  row.names = colnames(df)
)
```
Then, we generate the PCA plots using DESeq2 package:
```

# DSEQ2 package normalization, pca plot
dds <- DESeqDataSetFromMatrix( df, df_samples, design=~condition )
dds <- DESeq( dds )
vsd <- vst( dds )
plotPCA( vsd )
```
Gene expressions of cells under normal PH conditions and acidosis environment are well-seperated with PC1 score, which takes over 70% of the differences:</br>
![Untitled-10](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/8e7f7615-f856-4c1d-921f-5d79b1d0d11a)
 </br>

Then, we sort out conditions: normal PH vs. acidosis (labeled as "Control" and "High" in the read file, respectively) for gene expression comparison. 
```

##Sorting out conditions and output the results##

res_Control_vs_HighDosage <- results( dds, contrast=c( "condition", "Control", "High" ) )
dfk <- data.frame( res_Control_vs_HighDosage ) %>%
  rownames_to_column() %>%
  arrange( log2FoldChange )
```
We then find significant DEGs based on FDR (padj). We will record genes at any Fold change (FC) with FDR≤0.1 and then, a more stringent FDR≤0.05 at FC≥1.25 either up and down. We use org.Hs.eg.db package to match genes' Ensemble IDs to gene names. 
```

dfk<-dfk[!dfk$padj %in% NA,]
dfk<-dfk[dfk$padj<=0.1,]
annots <- select(org.Hs.eg.db, keys=dfk$rowname, columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[annots$ENSEMBL %in% dfk$rowname,]
ds <- annots[duplicated(annots$ENSEMBL),]
dfk <- dfk[!dfk$rowname %in% ds$ENSEMBL,]
annots <- annots[!annots$ENSEMBL %in% ds$ENSEMBL,]
annots <- annots[match(dfk$rowname, annots$ENSEMBL),]
ENSEMBL <- dfk$rowname
GeneSymbol <- annots$SYMBOL
log2FC_PH_7_vs_6 <- dfk$log2FoldChange
q_value <- dfk$padj
rm(res, dfk, annots, ds)
dfba <- data.frame(ENSEMBL,GeneSymbol,log2Fold,q_value)
dfba <- dfba[!dfba$GeneSymbol %in% NA,]
write_csv(dfba,'Control_vs_acidoses_q0.1_anyFC.csv')
dfba <- dfba[abs(dfba$log2FC_PH_7_vs_6)>=0.322,]
dfba <- dfba[dfba$q_value<=0.05,]
write_csv(dfba,'Control_vs_acidoses_q0.05_1.25_FC.csv')
```
To understand the upstream activation/inhibition and possible biological pathways enrichment/inhibition by acidosis, we conducted the core anlaysis using IPA. The outputs for the pathway analysis and upstream analysis are stored in the directory: "IPA_output". Running the MATLAB file: plotting_top_pathways.m in the same folder, we can visually observe the 12 inhibited pathways by acidosis with p≤0.05 and combined score (-log(p)*z) ≥5. These include cell cycle regulation, PTEN (related to PI3K/AKT), mTOR, and ILK signaling pathways. The corresponding top 10 inhibited upstreams by acidosis incluce PI3K family proteins and MAPK-related genes. On the countary, HIPPO pathway is enriched by acidosis. </br>
![Untitled-13](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/c9c9ded0-e8b3-4821-9b2b-e146b629fff8)

</br>

## Patients' survival analysis using patients' sequencing data and clinical outcomes from TCGA
Now, we will exaimine the overall survival outcomes of the genes upregualted by acidosis. General survival analysis (Kaplan-Meier plots and Cox regression analysis) is performed using survminer package (version: 0.4.9) ggsurvplot function (contained within survminer v.0.4.9) is used to generate Kaplan-Meier plot. Cox regressional analysis gives the HR and 95% confidence interval (CI). (TCGA_survival_curve_plotting_patient_TPM_recording.R) First, install and load the necessary libraries.
```
rm(list=ls()) 
# Install BiocManager package and TCGAbiolink
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")
#TCGAbiolinks version 2.22.4
#BiocManager::install("EDASeq")
#BiocManager::install("org.Hs.eg.db")
#install.packages("survminer")
#install.packages('DT')
#install.packages("tidyverse")
#BiocManager::install("SummarizedExperiment")
#install.packages("gtsummary")
library(SummarizedExperiment)
library(DT)
library(TCGAbiolinks)
library(tidyverse)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(readxl)
library(gtsummary)
#  # Query platform Illumina HiSeq 
# prep data from TCGA
setwd("~/TNBC_acidosis_RNA_seq_TCGA_analysis")
```
Then, prepare data from TCGA and extract RNA seq data. Download the TCGA data set to a local directory named "survivalanalysis". (Please create the directory "survivalanalysis" prior to performing this step.) The project name is "TCGA-BRCA"

```
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  experimental.strategy = "RNA-Seq", 
  access = "open", 
  legacy = FALSE
)
GDCdownload(
  query, 
  directory = "survivalanalysis"
)
#  # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
data <- GDCprepare(
  query, 
  directory = "survivalanalysis", 
  save = TRUE, 
  save.filename = "survivalanalysis/tcga_brca_data_2.15.RDS"
)
```
We then create expression matrix with geneID as rows and patientID (barcode) as columns. Using org.Hs.eg.db package, we can convert the gene ensemble ID to gene symbol. We normally use transcripts per million reads (TPM) normalization for the patients' gene expression, although Fragments Per Kilobase of transcripts per million mapped reads (FPKM) normalization is also acceptable. We also need to delete genes with very low expression (~0) and undetermined gene names. 
```
# create expression matrix with geneID in rows and samples (barcode) as col
data.mtx <- assay(data, "tpm_unstrand")
colnames(data.mtx) <- gsub("\\-\\w*\\-\\w*\\-\\w*\\-\\w*$", "", colnames(data.mtx))
rows <- mapIds(org.Hs.eg.db, gsub('\\.\\d*', '',rownames(data.mtx)), keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
rows <- unname(rows)
rownames(data.mtx) <- rows
## getting rid of low-expression rows and rows with name "NA"
data.mtx <- data.mtx[!is.na(rownames(data.mtx)),]
data.mtx <- data.mtx[!rowMedians(data.mtx)<=5,]
```
We then extract the patients' clinical outcomes information from the TCGABioLink package, and select TNBC (Subtype: Basal). Then, we match the TNBC patients with the updated pateints' information (Liu, J. et al., Cell 2018. The data, named as "patients.xlsx" is in the parent's directory). 
```

# get subtype information
data.subtype <- TCGAquery_subtype(tumor = "BRCA")
data.subtype <- data.subtype[data.subtype$BRCA_Subtype_PAM50 %in% c('Basal'),] #Basal: TNBC
## Upload updated clinical data. Only use the data.subtype if the information is not available 
## patient data is found in: Liu, J., et al., Cell 2018
Clinical_data <- read_excel("patient_data.xlsx")
Clinical_data <- Clinical_data[Clinical_data$bcr_patient_barcode %in% (data.subtype$patient),]
Clinical_data <- Clinical_data %>%
  mutate(OS.time = abs(as.numeric(Clinical_data$OS.time)/365.25))
Clinical_data <- Clinical_data %>%
  mutate(PFI = ifelse(is.na(Clinical_data$PFI), Clinical_data$OS, Clinical_data$PFI))
Clinical_data <- Clinical_data %>%
  mutate(PFI.time = abs(as.numeric(Clinical_data$PFI.time)/365.25))
rm(data.subtype)
```
We then construct the survival matrix for the subsequent overall survival Kaplan-Meier analysis. If you want to conduct progression free analysis, you can swithc OS to PFI, OS.time to PFI.time.
```
survivaldf <- data.frame(patient = Clinical_data$bcr_patient_barcode, 
                         status = Clinical_data$OS,
                         time = Clinical_data$OS.time)
survivaldf <- survivaldf[survivaldf$patient %in% colnames(data.mtx),]
Gene <- rownames(data.mtx)
```
We set the optimal cut-point for the gene expression in each patient via surv_cutpoint function (within survminer package) with the package's default minimal proportion of observation (10%). This converts each gene expression into a binary output ("high" versus "low") based on the gene expression cutpoint. We are specifically intrested in acidosis-induced genes that will bring a reasonable seperation of survival curve (p≤0.1). However, genes whose high expression causes significant survival curve difference (p≤0.05, evaluated by surv_pvalue function within survminer package) and bring a HR of more than 1.2 (evaluated by Cox regression) are selected as candidate genes for bringing worst OS outcomes. On the countary, genes whose high expression casues OS HR ≤ 1/1.2 are considered candidate genes for bringing better OS outcomes. 
```
survival_genes <- character()
survival_genes_better <- character()
for (i in 1:length(Gene)) {
  ## Extracting genes and tpms
  gene <- Gene[i]
  tpms <- data.mtx[i,]
  df <- survivaldf 
  df$tpm <- tpms[match(df$patient,names(tpms))]
  ## Determine the optimal cut-point for each gene. Minimum portion of observation group: 10% (Default). 
  cut <- surv_cutpoint(df, time = "time", event = "status", variables = c('tpm'), minprop = 0.1)
  ## Divide each variable values based on the cut-point 
  data.cut <- surv_categorize(cut)
  ## Computes an estimate of a survival curve for censored data using either the Kaplan-Meier
  fit <- survfit(Surv(time, status) ~ tpm, data = data.cut)
  #Calculating HR for each gene: 
  df2s <- coxph(Surv(time, status)~tpm=='high', data = data.cut)
  df2s <- df2s$coefficients
  df2s <- exp(df2s)
  ## If the survival pvalue is smaller than 0.05, the gene is considered significant.
  if (surv_pvalue(fit)$pval <= 0.1) {
    survivaldf[gene] <- data.cut$tpm
    #Prognostic significant genes: HR>=1.2 p<=0.05
    if (surv_pvalue(fit)$pval <= 0.05) {
      if (df2s>=1.2) {
        survival_genes <- append(survival_genes, gene)
      }
      if (df2s<=1/1.2){
        survival_genes_better <- append(survival_genes_better,gene)
      }
    }
  }
  
  rm(gene, tpms, df, cut, data.cut, fit, df2s)
}
```
We now evaluate the clinical outcomes for the high expression of acidosis-activated genes with FC≥1.25. If a patient have high expression of over 50% of any acidosis-activated genes, she is considered to have high expression of the overall acidosis-activated genes. We also extract the list of genes that bring worse OS outcome and are repressed by acidosis and list of genes that bring better OS outcome and are upregulated by acidosis. 
```
Gene <- read_csv("Control_vs_acidoses_q0.05_1.25_FC_up_down.csv")
Gene_repressed_by_acidosis <- Gene[Gene$log2FC_PH_7_vs_6>0,]
Gene_up <- Gene[Gene$log2FC_PH_7_vs_6<0,]


## Survival curve plotting
surdf_overlap <- survivaldf[,colnames(survivaldf) %in% c(Gene_up$GeneSymbol,"patient","status","time")]
surdf_overlap$Row_mode<-apply(surdf_overlap[], 1, function(x) {names(which.max(table(factor(x,unique(x)))))})

fit <- survfit(Surv(time, status) ~ Row_mode, data = surdf_overlap) 

coxph(Surv(time, status) ~ Row_mode=='high', data = surdf_overlap) %>%
  tbl_regression(exp=TRUE)

ggsurvplot(fit, pval = T, risk.table = T, 
           xlim = c(0, 10), break.time.by = 1,
           tables.y.text = F)
```
Out of the roughly 2,300 genes upregualted by acidosis, over 50% (~1,300) genes are reported to seperate survival curves by p≤0.1. Therefore, the survival curves for these genes is meaningful. High expression of these genes bring better OS outcome. A copy of this plot is stored in the directory: "TCGA_survival_curves": </br>
![Untitled-9](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/b0124e6a-9b08-4304-a5c9-61135874eac9) </br>
Then, we isolate the 270+ genes that bring worse OS outcomes and are respressed by acidosis, and 300+ genes that bring better OS outcomes and are upregulated by acidosis. The gene lists are in the directory: TCGA_prog_sig_genes. 
```
Gene_repressed_by_acidosis <- Gene_repressed_by_acidosis[Gene_repressed_by_acidosis$GeneSymbol %in% survival_genes, ]
Gene_up <- Gene_up[Gene_up$GeneSymbol %in% survival_genes_better,]
write_csv(Gene_repressed_by_acidosis, "TCGA_prog_sig_genes/Prognostic_significant_genes_repressed_by_acidosis.csv")
write_csv(Gene_up,"/TCGA_prog_sig_genes/Prognostic_better_outcome_genes_activated_by_acidosis.csv")
```
Indeed, the 270+ genes bring far worse OS outcomes as shown below (file is in the directory: "TCGA_survival_curves"). </br>
![Untitled-9](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/5aa40416-be6f-4af3-8c1a-652837a8d74f) </br>
Similarly, the 500+ genes bring better OS outcomes as shown below (file is in the directory: "TCGA_survival_curves") </br>
![500_genes_better_outomce_up_by_acidosis](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/68b4c386-cc31-487c-8d9b-27a0ea185470)

 </br>

## Functional characterization of the acidosis-induced DEGs
We now extract the downstream genes related to cell proliferation and actin cytoskeleton signaling and conclude 19 of these downstream genes are repressed by acidosis and bring worse OS outcomes. On the countary, 10 Hippo downstream genes are upregulated by acidosis and bring better OS outcome. The gene list is also in the folder: "IPA_output". </br>
![Untitled-11](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/e6742475-83aa-4781-8e83-db8992dda8a6) </br>
Among these, we saw acidosis represses NHE-related exchanger genes such asSLC9A1, which is also related to RHO-GTPase signaling and subsequent ILK signaling. ILK downstream MYL (Myosin Light chain) genes are also among the acidosis-repressed genes with worse OS outcomes. On the countary, HIPPO downstream LATS1 is upregulated by acidosis, meaning in the acidic environment, YAP/TAZ is more likely to be trapped in the cytoplasm and proliferation-related genes being repressed. Additional genes related to PI3K/AKT pathways, such as INPP5K and LAMTOR1 (Late Endosomal/Lysosomal Adaptor, MAPK And MTOR Activator 1), a protein coding gene that enables GTPase binding activity and contribute to the PI3K/AKT activity. We can plot the OS curves for each of the aferomentioned genes: </br>
```
# Example of plotting KM curve for INPP5K;
# Switch INPP5K to other genes for all KM plots
fit <- survfit(Surv(time, status) ~ INPP5K, data = surdf_overlap) 

coxph(Surv(time, status) ~ INPP5K=='high', data = surdf_overlap) %>%
  tbl_regression(exp=TRUE)

ggsurvplot(fit, pval = T, risk.table = T, 
           xlim = c(0, 10), break.time.by = 1,
           tables.y.text = F)
```
![306356869-2f227782-0379-451d-a3a6-774abec2f732](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/92a1939c-e82a-4827-8806-0c9150c0d53b)
</br>
![Untitled-15](https://github.com/taojiaxiang1991/TNBC_acidosis_RNA_seq_TCGA_analysis/assets/111034240/934a0a0b-ed28-46c6-8d02-10a0ff36e428) </br>



