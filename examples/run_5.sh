#!/bin/bash
module load sra_toolkit
fastq-dump -I --split-files --gzip -M 36 SRR17729734 SRR17729735 SRR17729736 SRR17729737 SRR17729738 SRR17729739