# CRISPRt_gRNA Design for CN PAMs
[![DOI](https://zenodo.org/badge/567437408.svg)](https://zenodo.org/badge/latestdoi/567437408)

## Overview
This repository, a specialized adaptation of the `sgrna_design` suite, is tailored for designing sgRNAs with CN PAMs, 32 nucleotides in length, for CRISPRt Overexpression applications in bacterial systems. It extends John S. Hawkins's foundational work, which was modified by Ryan D. Ward.

## Quick-start Guide

### Cloning, Environment Setup, and Genomic Data Retrieval
- Clone the repository and set up the Conda environment.
- Retrieve genomic data from NCBI using the accession number.

### Executing the sgRNA Design Script
Run the design script to generate sgRNAs:
```
./CN32_design.py --input_genbank_genome_name ${ACC_NO}.gb --tsv_output_file ${ACC_NO}_sgrna.tsv
```

### Post-Processing with `awk` Functions
After sgRNA design, use `awk` commands to process the output:
```
awk 'NR == 1 {print; next} $2 <= -95 && $9 == "sense" && $10 == 39 {print $0 | "sort -k1,1V -k2,2nr"}' ${ACC_NO}_sgrna.tsv | 
awk 'NR == 1 {print; next} {guides[$1]++} guides[$1] <= 10 {print}' > top_ten_sgRNA_outputs.tsv #ranked
```
These commands filter and rank the sgRNAs based on specific criteria, prioritizing those with high specificity and limiting the selection to the top ten sgRNAs per gene.

## Modifications in This Approach
- Adapted for CN PAM compatibility.
- sgRNA length increased to 32 nucleotides.
- Optimized for identifying upstream gene targets for CRISPRtOE.
