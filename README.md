# edgeR_nf-core-rnaseq
This repository consists of a set of functions to use implement an edgeR R package pipeline which obtains differential expressed genes automatically throughtCPU-parallelized computing. 

A full list of libraries needed for run this code is included below.

**Dependencies:** `R (>= 4.0.0)`

**Imports:** `base, utils, methods, stats, ggplot2, tximport, parallel`



## Description:
![Code schema](https://raw.githubusercontent.com/ccsosa/edgeR_nf-core-rnaseq/main/images/edgeR_pipeline.drawio.png)

### Steps
- Prepare the samples file for running nf-core/rnaseq (please see the file added in the examples folder

### Outcomes structure:

![Directory structure](https://github.com/ccsosa/edgeR_nf-core-rnaseq/blob/main/images/edgeR_outcomes.jpg)


### How to run the code:
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
