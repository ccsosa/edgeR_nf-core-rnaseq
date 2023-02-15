# edgeR_nf-core-rnaseq
This repository consists of a function to implement an edgeR R package pipeline which obtains differential expressed genes automatically through CPU-parallelized computing for pairwise contrasts. 

A full list of libraries needed to run this code is included below.

**Dependencies:** `R (>= 4.0.0)`

**Imports:** `base, utils, methods, stats, ggplot2, tximport, parallel`

## example files loaded
- INDEX_BASH_1.sh (Salmon order to get a salmon index folder with the files to run nf-core/rnaseq) 
    >This must be done before run nf-core/rnaseq pipeline!
- run_5.sh (Bash script to download the SRR samples from NCBI with fastq-dump)
- run5.csv (nf-core/rnaseq sample file to be used to get the salmon counts)
- run5_nf_test.sh (Bash script to run nextflow nf-core/rnaseq)
- nfcore_rna_seq.config (Config file to run nf-core/rnaseq. This file is in the folder nextflow_config)

# Description:

This code performs the following steps described graphically in the following schema :

![Code schema](https://raw.githubusercontent.com/ccsosa/edgeR_nf-core-rnaseq/main/images/edgeR_pipeline.drawio.png)

## Data preprocessing
- 1.) Create output-dir folder structure with two subfolders: csv and graphics
- 2.) Read Salmon outcomes (salmon_tx2gene.tsv, and tximport)
- 3.) Read the metadata file
- 4.) Create  groups for the samples based on the samples files provided for nf-core/rnaseq using combinatorics or using a group_vect object

## edgeR preprocessing
- 5.) Save raw counts obtained from salmon files
- 6.) Creates a sample count table and calculates normalizing factors
- 7.) Creates an exploratory plot (plotMDS) and save it (MDS.pdf)
- 8.) Creates the design object to be used for differential expressed genes step
- 9.) Estimate the dispersion  and save a plot of it (plotBCV.pdf)
- 10.) Creates the contrast object to be used for the differential expressed genes step (list)
## edgeR paralleling processing
- 11.) Run in parallel the next steps per contrast:
    - Use a Fit a quasi-likelihood negative binomial generalized log-linear model (glmQLFit) with the design and contrasts
    - Use glmQLFTest test to get the DEG outcome table
    - Run a Benjamini-Hochberg procedure to obtain false discovery rate values (FDR)
    - Split up and downregulated genes using the sign of the log2fold change
    - Save raw and filtered by a p-value results
    - Save plots of glmQLFit and biological coefficient variance
- 12.) Creates a summary file for all up and downregulated genes for all the contrasts provided
- 13.) Return a list with two slots: contrasts and the results of the contrasts (DEG)

# Parameters required:
- `folder_name` Folder name with the results obtained by nf-core/rnaseq results. 
    >If you have a folder named "2" this name will be used to call the metadata table and to name the folder with the outputs
- `data_dir` Folder name with the results obtained by nf-core/rnaseq for the Salmon counts
- `out_dir` Folder where the outcomes will be written
- `met_dir`  Folder where the metadata file will be loaded. Please name your metadata as folder name (e.g. run2.csv")
- `pval` p-value used to filter the results of differential expressed genes (default value = 0.05)
- `numCores` numeric, Number of cores to use for the process (default value numCores=2)
- `plot_MDS` This is a boolean value to indicate if the exploratory PCA plot should be displayed in the R session
- `group_vect` (default value numCores=2) This is a character object which represents the groups available in the metadata data.
    - If the value is NULL the script automatically will obtain groups using the sample names. 
    - For instance, if there are four samples named CONTROL_SRR1, CONTROL_SRR2,DROUGHT_SRR3, DROUGHT_SRR4, the groups will be control and drought respectively.
    - If values are provided, they must respect the order used in the metadata file. The group_vect object will be:
     ``` group_vect <- c("CONTROL","CONTROL","DROUGHT","DROUGHT") ```

## R code outcomes structure:
This code creates a folder with the name given in the `folder_name` parameter. The results have this structure:
![Directory structure](https://github.com/ccsosa/edgeR_nf-core-rnaseq/blob/main/images/edgeR_outcomes.jpg)

This code will create two subfolders:
- `csv`: CSV files obtained from the edgeR pipeline. Two files are created per contrast:
  - glmQLFTest_[CONTRAST]_pval_[pval]fulltable.csv (Outcome without p-value filter)
  - glmQLFTest_[CONTRAST]_pval_[pval]_filtered.csv (Outcome with a p-value threshold filter)
- `graphics`: subfolder to save the edgeR plots per contrast. Two files are created per contrast:
  - plotMD_glmQLFTest_[CONTRAST].pdf (Average log CPM vs LFC  value displaying differential expressed genes)
  - [CONTRAST]_BCV.pdf (Average Log2 CPM Vs Quarter root mean deviance)

In the main folder the next files are saved:
- Exploratory analysis:
  - abundance_drought_tpm_genelevel.txt (Raw tpm for the samples provided)
  - CPM_normalized.txt (Raw Counts per million per gene)
  - library_size.txt (library sizes per sample)
  - contrasts.csv (Contrast designed inferred by limma and edgeR)
  - plotBCV.pdf (Average log CPM Vs biological coefficient variation plot)
  - MDS.pdf (Multidimensional plot, each group is displayed in numbers)
- Summaries per contrast 
  - glmQLFTest_[CONTRAST]_pval_[pval]summary.csv (Summary of up and downregulated genes per contrast after applying a p-value threshold)
 
### Differential expressed genes outcome format:

"" | logFC | logCPM | F | PValue | FDR | status | status_name
------------ | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
TraesCS1A02G000200 | 0.116659367 | 1.677943702 | 0.118024516 | 0.739867646 | 0.820941735 | 1 | UP
TraesCS1A02G000400 | -0.2276499 | 8.025249611 | 2.010384676 | 0.193246796 | 0.313335743 | -1 | DOWN
TraesCS1A02G000900 | -0.146732262 | 9.310521499 | 0.737121095 | 0.415080596 | 0.545800092 | -1 | DOWN
TraesCS1A02G002500 | 0.211749279 | 5.076743544 | 4.355186469 | 0.06966184 | 0.147293307 | 1 | UP
TraesCS1A02G002700 | -0.517768365 | 5.575591155 | 11.08925656 | 0.010088987 | 0.035546276 | -1 | DOWN

- Gene (Gene id)
- logFC (Log2Fold change)
- logCPM (log of counts per million for each gene)
- F (F statistics value)
- Pvalue (p-value)
- FDR (Benjamini-Hochberg Procedure applied over the p-value)
- status (logFC sign: 1 and -1 for positive and negative values respectively)
- status_name : UP for positive logFC and DOWN for negative logFC)


# Steps to run nf-core/rnaseq stage 1 and 3 (This is mandatory to use this R code)

- Download the transcriptome and genome fasta files for your target species
- Prepare a Salmon index file according to your needs (A large genome probably will work better with no decoy aware option)  (see example file: INDEX_BASH_1.sh )
  > See https://salmon.readthedocs.io/en/latest/salmon.html
 - Download the SRR samples with the fastq files from NCBI using fastq-dump (see example file: run_5.sh)
- Prepare the samples file for running nf-core/rnaseq (please see the file added in the examples folder:`run5.csv`). The format is the following and contains the next columns:
![Sample_fig](https://github.com/ccsosa/edgeR_nf-core-rnaseq/blob/main/images/sample_file_rnaseq.jpg)
    -sample: This is the sample name used for nf-core/rnaseq and it represents a group and the SRR downloaded from NCBI. 
    >(PLEASE USE "_" as a separator always to use the automatic feature of the R code!)
    - fastq_1: Directory of the forward fastq file for the SRR sample downloaded
    - fastq_2: Directory of the reverse fastq file for the SRR sample downloaded
    - strandedness: Represent the strand of the RNA-Seq experiment. Leave as unstranded. Salmon will detect strandedness automatically
      > For relevant information see: https://nf-co.re/rnaseq
      
      >If you want that controls always appear in this way Treatment-Control, please add a letter after C (e.g. Aluminium would be Taluminium, this ensure that the contrast will be Taluminum-Control)
      
    - Configure resource to be used for the nf-core/rnaseq run (see in the nextflow_config file: nfcore_rna_seq.config)
- Submit to your HPC for processing using your resource management tool (An example is provided in the examples folder: run5_nf_test.sh).
```
    #!/bin/bash
module load  nextflow
nextflow run nf-core/rnaseq --input /scratch/bis_klpoe/chsos/data/sample_files/run5.csv --skip_alignment --outdir /scratch/bis_klpoe/chsos/analysis/ --pseudo_aligner 'salmon' -profile singularity -c /scratch/bis_klpoe/chsos/data/config_file/nfcore_rna_seq.config -w /scratch/bis_klpoe/chsos/analysis/work  --salmon_quant_libtype A
```
- Run parameters:
    - --skip_alignment (skip use alignment)
    -  --pseudo_aligner 'salmon' (Use Salmon for pseudo alignments and counts)
    -  --salmon_quant_libtype A (Allows salmon to find the samples strandedness)
    -  --input (sample file prepared previously
    -  --outdir (Folder  to save the nf-core/rnaseq results)
    -  -profile singularity (Use singularity docker)
    -  -c (Configuration file to be read)
    -  -w (workdir)


# How to run the R code:
-  Configure the config file to run nf-core/rnaseq pipeline stages 1 and 3 (https://nf-co.re/rnaseq) (An example is available in the folder examples)
-  Run nf-core/rnaseq pipeline
-  If the user has several salmon counts, please name the salmon file as salmon_[`folder_name`] (e.g. folder_name is "2")
-  Detect metadata file to use. This is the sample file used for nf-core/rnaseq pipeline
-  Include in `met_dir` the folder where your sample data is saved (This is used as metadata to run DEG in edgeR)
-  Select the number of cores to use `numCores`
  > The number of cores never would be equal to the number of cores available in your machine
- Define if you want to use a  group_vect object to assign samples into groups.
  > The order you provide the groups must match to the sample order used in the sample file to run nf-core/rnaseq pipeline! 
-  Define if you want to observe an exploratory plot in your R session with the parameter `plot_MDS`
- Load the code `edgeR_func.R`
- Run



```r


##############################################################################
#parameters
folder_name <- "5" #folder name
pval = 0.05 #p value for filtering
numCores <- 4 #number of cores to use in parallel
plot_MDS <- TRUE #if plot should be appears in R session
#groups if it not easy to get from the headers
group_vect <- c("LEAVES_D","LEAVES_D","LEAVES_D",
                "LEAVES_C","LEAVES_C","LEAVES_C") #there are six samples
#automatic option group_vect <- NULL
##############################################################################
#Directories
#where are the salmon files
data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
#defining output folder
out_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"
#path where the nf core metadata is available
met_dir <- "/scratch/bis_klpoe/chsos/data/sample_files/DONE"
##############################################################################
x <- DEG_edgeR_func(folder_name=folder_name,
                    data_dir=data_dir,
                    pval=pval,
                    out_dir=out_dir,
                    met_dir=met_dir,
                    plot_MDS=TRUE,
                    numCores=4,
                    group_vect=group_vect)



```

## Authors

Main:Chrystian C. Sosa

Other contributors: Indeewari Dissanayake, Gust Bilcke
