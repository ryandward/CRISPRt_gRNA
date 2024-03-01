# CRISPRt_gRNA Design for CN PAMs

[![DOI](https://zenodo.org/badge/567437408.svg)](https://zenodo.org/badge/latestdoi/567437408)

## Overview
This repository enhances the `sgrna_design` toolkit for designing sgRNAs that target CN PAMs, extending to 32 nucleotides in length, specifically for CRISPRt Overexpression (CRISPRtOE) in bacterial systems. Building upon John S. Hawkins's foundational work and further modifications by Ryan D. Ward, this adaptation focuses on upstream gene target identification for precise CRISPRtOE applications.

## Getting Started

### Setup and Data Preparation

1. **Clone the Repository and Prepare the Environment**
   Begin by cloning this repository to your local machine and setting up the Conda environment to manage dependencies.
   ```bash
   git clone https://github.com/ryandward/CRISPRt_gRNA.git
   cd CRISPRt_gRNA
   conda env create -f environment.yml
   conda activate CRISPRt_gRNA
   ```

2. **Retrieve Genomic Data Using `seq_hunter`**
   Use the `seq_hunter` package to download the required genomic sequences from GenBank. This tool is part of the [ryandward/seq_hunter toolkit](https://github.com/ryandward/seq_hunter) repository, designed to streamline sequence retrieval.
   ```bash
   ❯ python seq_hunter.py --source GenBank GCF_000005845.2 GCA_003054575.1
   ```
   Follow the prompts to complete the download, which will store the sequences in the specified directory.

### sgRNA Design Execution

To generate sgRNAs, execute the design script with the input GenBank genome file and specify the output file for sgRNAs.
```bash
./CN32_design.py --input_genbank_genome_name ${ACC_NO}.gb --tsv_output_file ${ACC_NO}_sgrna.tsv
```

### Post-Design Processing

Filter and rank the sgRNAs using `awk`, focusing on specificity and limiting selection to the top ten sgRNAs per gene.
```bash
awk 'NR == 1 {print; next} $2 <= -95 && $9 == "sense" && $10 == 39 {print $0 | "sort -k1,1V -k2,2nr"}' ${ACC_NO}_sgrna.tsv | 
awk 'NR == 1 {print; next} {guides[$1]++} guides[$1] <= 10 {print}' > top_ten_sgRNA_outputs.tsv
```

## Key Modifications

- **CN PAM Compatibility:** Adapted for CN PAMs, broadening the toolkit's applicability.
- **Extended sgRNA Length:** Increases sgRNA length to 32 nucleotides to enhance target specificity.
- **Optimized for Upstream Gene Targets:** Focused on identifying upstream gene targets for CRISPRtOE, improving the precision of genetic modifications.

### Guide Intersection Analysis

Identify regions where CRISPRt might insert inside an upstream gene using the `GenBankParser` tool from the [ryandward/barcoder toolkit](https://github.com/ryandward/barcoder).
```bash
pipenv activate
python find_crisprt_sites.py
```

By integrating `seq_hunter` for genomic data retrieval, this documentation provides a comprehensive guide to utilizing the toolkit for sgRNA design with CN PAMs, ensuring users have the information needed to efficiently execute the design process and analyze potential guide RNA intersections with upstream genes.
