# edgeR_nf-core-rnaseq
This repository consists of a set of functions to use implement an edgeR R package pipeline which obtains differential expressed genes automatically throughtCPU-parallelized computing. 

A full list of libraries needed for run this code is included below.

**Dependencies:** `R (>= 4.0.0)`

**Imports:** `base, utils, methods, stats, ggplot2, tximport, parallel`



# Description:

This code perform the following steps described graphically in the following schema :

![Code schema](https://raw.githubusercontent.com/ccsosa/edgeR_nf-core-rnaseq/main/images/edgeR_pipeline.drawio.png)

## Data preprocessing
- 1.) Create output dir folder structure with two subfolders: csv and graphics
- 2.) Read Salmon outcomes (salmon_tx2gene.tsv, and tximport)
- 3.) Read the metadata file
- 4.) Create  groups for the samples based on the samples files provided for nf-core/rnaseq using combinatorics or using a group_vect object

## edgeR preprocessing
- 5.) Save raw counts obtained from salmon files
- 6.) Creates a sample count table and calculate normalizing factors
- 7.) Creates a exploratory plot (plotMDS) and save it (MDS.pdf)
- 8.) Creates the design object to be use for differential expressed genes step
- 9.) Estimate the dispersion  and save a plot of it (plotBCV.pdf)
- 10.) Creates the contrast object to be use for differential expressed genes step (list)
## edgeR paralleling processing
- 11.) Run in parallel the next steps per contrast:
    - Use a Fit a quasi-likelihood negative binomial generalized log-linear model (glmQLFit) with the design and contrasts
    - Use glmQLFTest test to get the DEG outcome table
    - Run a Benjamini-Hochberg procedure to obtain false discovery rate values (FDR)
    - Split up and downregulated genes using the sign of the log2fold change
    - Save raw and filtered by a p value results
    - Save plots of glmQLFit and biological coefficient variance
- 12.) Creates a summary file for all up and downregulated genes for all the contrasts provided
- 13.) Return a list with two slots: contrasts and the results of the contrasts (DEG)


# Nextflow (nf-core/rnaseq) steps
- Prepare the samples file for running nf-core/rnaseq (please see the file added in the examples folder 

# Outcomes structure:
This code creates a folder with the name given in the the `folder_name` parameter. The results have this structure:
![Directory structure](https://github.com/ccsosa/edgeR_nf-core-rnaseq/blob/main/images/edgeR_outcomes.jpg)

This code will create two subfolders:
- `csv`: CSV files obtained from the edgeR pipeline. Two files are created per contrast:
  - glmQLFTest_[CONTRAST]_pval_[pval]fulltable.csv (Outcome without p value filter)
  - glmQLFTest_[CONTRAST]_pval_[pval]_filtered.csv (Outcome with a p value threshold filter)
- `graphics`: subfolder to save the edgeR plots per contrast. Two files are created per contrast:
  - plotMD_glmQLFTest_[CONTRAST].pdf (Average log CPM vs LFC  value displaying differential expressed genes)
  - [CONTRAST]_BCV.pdf (Average Log2 CPM Vs Quarter root mean deviance)

In the main folder the next files are saved:
- Exploratory analysis:
  - abundance_drought_tpm_genelevel.txt (Raw tpm for the samples provided)
  - CPM_normalized.txt (Raw Counts per million per gene)
  - library_size.txt (library sizes per sample)
  - contrasts.csv (Contrast designed infered by limma and edgeR)
  - plotBCV.pdf (Average log CPM Vs biological coefficient variation plot)
  - MDS.pdf (Multidimensional plot, each group is displayed in numbers)
- Summaries per contrast 
  - glmQLFTest_[CONTRAST]_pval_[pval]summary.csv (Summary of up and downregulated genes per contrast after apply a p-value threshold)
 
## Differential expressed genes outcome format:

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




# How to run the code:
-  Configure the config file to run nf-core/rnaseq pipeline stages 1 and 3 (https://nf-co.re/rnaseq) (An example is available in the folder examples)
-  Run nf-core/rnaseq pipeline
-  If user have several salmon counts, please name the salmon file as salmon_[`folder_name`] (e.g. folder_name is "2")
-  Detect metadata file to use. This is the sample file used for nf-core/rnaseq pipeline
-  

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
