# RNA-seq raw-read processing, aligning, and quantifying
## Background
This section is broken into sequential steps similar to the root of this repository. Under each step, we include a `commands` file which stores the lines run. Each step was completed by running `bash commands` within the root of a step's directory (e.g. `cd ./2-trim; bash commands`).

**NOTE 1:** _This repository excludes files which are in excess of 10MB in an effort to keep this work light enough for Github. Raw sequence files are available under the usual channels (see paper)._

**NOTE 2:** _Most scripts in this section feature [absolute paths](https://www.redhat.com/en/blog/linux-path-absolute-relative) to programs or files. Please make note of these paths and **update them according to your particular environment**._

**NOTE 3:** _We use **job counts which may not match your system's capabilities**. Please check our scripts' `-j`, `-t`, and other multithreading parameters BEFORE running them!_

# Pipeline summary
## Peculiarities

Since in-house bioinformatics pipelines tend to have unusual variances in approach between labs, we will note some of ours. These are two potentially surprising design choices we made for ours:
- The use of [GNU parallel](https://en.wikipedia.org/wiki/GNU_parallel) to parallelize jobs and manage system load
    - This is an easy way to parallelize jobs within one system, but can result in surprising syntax where we must manipulate file strings prior to running a "main" program.
- The use of absolute paths to [Miniforge (Conda)](https://github.com/conda-forge/miniforge) environments' scripts.
    - This is a workaround to allow `parallel` to access scripts without the long delay of activating a conda environment within a job.
    - This approach can result in complications with nested dependencies. In this work, instances of this complication are resolved by the `--path_to_cutadapt` parameter in Trim Galore.

## 2-trim
Here, we use [Trim Galore](https://github.com/FelixKrueger/TrimGalore) to trim our fastq data by a quality cutoff of `30` and against adapter sequences `CTGTCTCTTATACACATCT` and `AGATGTGTATAAGAGACAG`.
We take raw reads located at `../0-fastq/*.fastq.gz` and remove `_S#...` endings from the file names to shorten them.

```bash
# Content of /1-RNAseq-STAR/2-trim/commands
find ../0-fastq/*.fastq.gz | sed -r 's/^(.*_S[0-9]).*$/\1/g' | sort -u | parallel -j 28 \
"/home/lexic/.conda/envs/pre/bin/trim_galore \
--path_to_cutadapt /home/lexic/.conda/envs/pre/bin/cutadapt \
-a 'CTGTCTCTTATACACATCT' \
-a 'AGATGTGTATAAGAGACAG'  \
-j 2 -q 30 {}*.gz \
--gzip -o . > {/}.log 2>&1"
```

## 3-align
We align our trimmed reads from `../2-trim/*.gz` (created the step prior) using [STAR](https://github.com/alexdobin/STAR) against the [Human - Release 47 GRCh38.p14](https://www.gencodegenes.org/human/) genome. Once again, we make some string manipulations with `sed` to prepare our fastq files for paired or unpaired alignment.

Our datasets were a mix of paired and unpaired sequence files, so we wrote slightly different alignments for each. When the unpaired reads were processed, it was without the presence of paired-end reads in `2-trim`. This may require more creative separating to now run together.

```bash
# Content of /1-RNAseq-STAR/2-trim/commands
# Unpaired
ls ../2-trim/*.gz | grep -Pv "/SUM" | sed -r 's/^(.*\/[^_]+_.+)_S[0-9].*gz$/\1/g' | sort -u | parallel --progress --joblog unpaired_parallel.log -j 4 "(time (/opt/miniforge/envs/star/bin/STAR --runThreadN 15 --genomeDir /blazer/references/grch38 --readFilesIn \$(ls -1 {}*.fq.gz | paste -sd ',' -) --outFileNamePrefix {/} --readFilesCommand zcat)) 2>&1 | tee {/.}.log"

# Paired end
ls ../2-trim/*.gz | grep -P "/SUM" | sed -r 's/^(.*)_R[12].*.gz/\1/g' | sort -u | parallel --progress --joblog single_parallel.log -j 4 "(time (/opt/miniforge/envs/star/bin/STAR --runThreadN 15 --genomeDir /blazer/references/grch38 --readFilesIn {}*.fq.gz --outFileNamePrefix {/} --readFilesCommand zcat)) 2>&1 | tee {/.}.log"
```

## 4-sambam
Next we compress the [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) files generated as `../3-align/*.sam` into smaller [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) files using [samtools](https://github.com/samtools/samtools).
```bash
# Content of /1-RNAseq-STAR/4-sambam/commands
ls ../3-align/*.sam | sed -r 's/^(.*\/[^_]+_.+)Aligned.*\.sam$/\1/g' | sort -u | parallel --progress --joblog parallel.log -j 50 "(time (/home/lexic/.conda/envs/samtools/bin/samtools view -S -b {}*.sam > {/}.bam)) 2>&1 | tee {/.}.log"
```

## 5-sort
Now we sort the BAM files created in `../4-sambam/*.bam` so our quantification process can more quickly parse them. This again uses samtools.

```bash
# Content of /1-RNAseq-STAR/5-sort/commands
ls ../4-sambam/*.bam | sed -r 's/^(.*\/[^_]+_.+)\.bam$/\1/g' | sort -u | parallel --progress --joblog parallel.log -j 8 "(time (/home/lexic/.conda/envs/samtools/bin/samtools sort -@ 7 -o {/}_sorted.bam {}.bam)) 2>&1 | tee {/.}.log"
```

## 6-count
Using the sorted BAM files in `../5-sort/*.bam`, we can now quantify the number of reads present on each gene as counts. Quantification is done using [Subread](https://subread.sourceforge.net/)'s [featureCounts](https://subread.sourceforge.net/featureCounts.html) program.

This program creates a count matrix for the list of BAM files we provide for every gene. We exclude unnamed genes via `grep -Pv "^ENSG"` to obtain only named genes for DESeq. Our DESeq pipeline will handle the uniting of these two count matrices into one.
```bash
# Content of /1-RNAseq-STAR/6-count/commands
# Single-end reads
(time (/opt/miniforge/envs/star/bin/featureCounts -a /blazer/references/grch38/gencode.v44.chr_patch_hapl_scaff.annotation.gtf -T 60 -o count_matrix_single.tsv $(ls ../5-sort/*.bam | grep -v "/SUM") -g 'gene_name')) 2>&1 | tee single_output.log

# Paired-end reads
(time (/opt/miniforge/envs/star/bin/featureCounts -a /blazer/references/grch38/gencode.v44.chr_patch_hapl_scaff.annotation.gtf -T 60 -o count_matrix_paired.tsv -p $(ls ../5-sort/*.bam | grep "/SUM") -g 'gene_name')) 2>&1 | tee paired_output.log

cat count_matrix_single.tsv | grep -Pv "^ENSG" | grep -Pv "#" > count_matrix_noENSG_single.tsv
cat count_matrix_paired.tsv | grep -Pv "^ENSG" | grep -Pv "#" > count_matrix_noENSG_paired.tsv
```
# Summary
## Intake:
- FASTQ format reads from various RNA-seq experiments
- GRCh38 human genome for alignment

## Output:
- `count_matrix_noENSG_single.tsv`
    - A [count matrix](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input) for all single-end RNA-seq data.
- `count_matrix_noENSG_paired.tsv`
    - A [count matrix](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input) for all paired-end RNA-seq data.