#!/bin/bash
module load  nextflow
nextflow run nf-core/rnaseq --input /scratch/bis_klpoe/chsos/data/sample_files/run5.csv --skip_alignment --outdir /scratch/bis_klpoe/chsos/analysis/ --pseudo_aligner 'salmon' -profile singularity -c /scratch/bis_klpoe/chsos/data/config_file/nfcore_rna_seq.config -w /scratch/bis_klpoe/chsos/analysis/work  --salmon_quant_libtype A