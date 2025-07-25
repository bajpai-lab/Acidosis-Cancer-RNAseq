# TCGA analysis and generation of prognosis genes
## Background
The main script of this subsection `./TCGA.R` uses [The Cancer Genome Atlas](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) to identify genes associated with better and worse survival outcomes for patients with Basal BCRA-associated tumors.
1. To do this, we first download "TCGA-BRCA" RNA-Seq transcriptome data using the [`TCGAbiolinks`](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) package.
    - _See `# Acquiring TCGA datasets` in the script for details._
2. Then, we create a dataset associating patient TPM values from TCGA transcriptome results with the patient's clinical metadata.
    - It is also here that we filter for patients with "Basal" subtype of BRCA-associated tumors.
    - _See `# Creating "data.mtx" and "survivaldf"` in the script for details_
3. With that dataset, we analyze Kaplan-Meier survival fits for each gene in `data.mtx` across all patients in `survivaldf`. If this gene's curve meets our criteria for a "better" or "worse" outcome, we append it to `survival_genes_better` or `survival_genes_worse`.
    - Currently, those criteria are in terms of curve fit statistics of `p <= 0.05` and Hazard Ratio (variable called `df2s`) of `>= 1.2` or `<= 1/1.2`.
    - _See `# Populate "survival_genes_worse", "survival_genes_better"` in the script for details_
4. With a list of prognosis-relevant genes from TCGA, we are now able to filter the genes which returned as hits in our DESeq analysis.
    - _Note: If you are replicating our work, you only need to run a portion of the DESeq analysis to obtain hits to pass to TCGA analysis. Skip lines which include "prog genes" if you do not have them yet. Otherwise, our provided list of prognosis genes should work fine for a full run of the DESeq script._
    - The function to generate these plots is stored under `# Shared function to create all survival plots` in the TCGA script. Run this to create the function definition first.
        - The function will automatically generate a plot, risk table, gene list, and "high" curve stats.
        - To plot a gene list, we are aggregating all gene "high" and "low" expresssions in a patient based on the mode of the two.
            - By itself, this does not account for the statistical relevance of a gene's "high" and "low" expression in each patient!
            - We filter by `survival_genes_worse` or `survival_genes_better` to obtain statistically-relevant genes.
    - Calls to this function are made under `# Creating Kaplan-Meier curves` in the script. Here, we load our DESeq gene hits and:
        1. Sort genes by those "activated" or "repressed" by acidosis using their log2 fold changes.
        2. Filter genes by those prognostically relevant from our TCGA analysis.
        3. Pass our list of genes to the plotting function.

# Summary
## Intake:
- [The Cancer Genome Atlas](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) or `/3-TCGA-R/survivalanalysis/tcga_brca_data_2.15.RDS`
    - Clinical data for determining prognostically-significant genes in a tumor type.
- `/3-DESeq-R/results/intersect_PlanA_(10w)U(24h)_BOTH-FC-by-acidosis_0.1FDR.csv`
    - Gene hits from DESeq2 analysis of acidosis-treated cell lines.

## Output:

- `/3-TCGA-R/survivalanalysis/Prognostic_better_outcome_genes_0.05q_1.2HR.csv`
    - Human genes significantly associated with a better outcome in cancer patients.
- `/3-TCGA-R/survivalanalysis/Prognostic_worse_outcome_genes_0.05q_1.2HR.csv`
    - Human genes significantly associated with a worse outcome in cancer patients.
- `/3-TCGA-R/PlanA/`
    - Kaplan-Meier survival analysis and [coxph](https://rdrr.io/cran/survival/man/coxph.html) statistics for genes from "Plan A" RNA-seq analysis in DESeq2.
    - **PlanA** = Intersection of 0.1 FDR gene hits found in both MDA-MB-231 10 weeks and 24 hours whose log2FoldChange act in the same direction.
        - `BOTH-FC(q0.1(MDA 10w) ∪ q0.1(MDA 24h))`

For more details on this analysis, see [**the full (partially outdated) README-FULL.md**](README-FULL.md). 