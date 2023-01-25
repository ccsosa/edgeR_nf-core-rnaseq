#!/bin/bash

#salmon index -t gentrome.fa -d decoys.txt -p 8 -i salmon_index --gencode
#--decoys /scratch/bis_klpoe/chsos/data/salmon_index/decoys/decoys.txt
module load  salmon
salmon index -i transcripts_index --threads 6  --kmerLen 31 --tmpdir /scratch/bis_klpoe/chsos/data/salmon_index/tmp --keepFixedFasta --transcripts /scratch/bis_klpoe/chsos/data/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_transcripts.fasta