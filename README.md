# edgeR_nf-core-rnaseq
This repository consists of a set of functions to use implement an edgeR R package pipeline which obtains differential expressed genes automatically throughtCPU-parallelized computing. 

A full list of libraries needed for run this code is included below.

**Dependencies:** `R (>= 4.0.0)`

**Imports:** `base, utils, methods, stats, ggplot2, tximport, parallel`



## Description:
![Code schema](https://raw.githubusercontent.com/ccsosa/edgeR_nf-core-rnaseq/main/images/edgeR_pipeline.drawio.png)

### Steps
- Prepare the samples file for running nf-core/rnaseq (please see the file added in the examples folder
