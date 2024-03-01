# CRISPRt_gRNA Design for CN PAMs

[![DOI](https://zenodo.org/badge/567437408.svg)](https://zenodo.org/badge/latestdoi/567437408)

## Overview
This repository enhances the `sgrna_design` toolkit for designing sgRNAs targeting CN PAMs, extending to 32 nucleotides in length, tailored for CRISPRt Overexpression (CRISPRtOE) in bacterial systems. Leveraging John S. Hawkins's foundational work with further modifications by Ryan D. Ward, this adaptation focuses on precise upstream gene target identification for CRISPRtOE applications.

## Getting Started

### Setup and Data Preparation

1. **Clone the Repository and Prepare the Pip Environment**

   Begin by cloning this repository to your local machine. Then, set up the project's environment using `pipenv` to manage dependencies.
   ```bash
   git clone https://github.com/ryandward/CRISPRt_gRNA.git
   cd CRISPRt_gRNA
   pipenv install
   ```

2. **Retrieve Genomic Data Using `seq_hunter`**

   To download the required genomic sequences from GenBank, use the `seq_hunter` package, housed within its separate repository at [ryandward/seq_hunter](https://github.com/ryandward/seq_hunter). Note that `seq_hunter` also utilizes `pipenv` for its environment, ensuring isolation and consistent dependency management across different projects.
   - First, ensure you have cloned the `seq_hunter` repository and prepared its environment:
     ```bash
     git clone https://github.com/ryandward/seq_hunter.git
     cd seq_hunter
     pipenv install
     ```
   - Next, run `seq_hunter` to download sequences:
     ```bash
     pipenv run python seq_hunter.py --source GenBank GCF_000005845.2 GCA_003054575.1
     ```
   This command initiates the sequence retrieval process, with the output detailing the downloaded sequences and their storage locations.

### sgRNA Design Execution

To generate sgRNAs, execute the design script within the `CRISPRt_gRNA` environment, specifying the GenBank genome file and the sgRNAs' output file.
```bash
pipenv run python ./CN32_design.py --input_genbank_genome_name ${ACC_NO}.gb --tsv_output_file ${ACC_NO}_sgrna.tsv
```

### Post-Design Processing

Filter and rank the sgRNAs using `awk`, limiting the selection to the top ten sgRNAs per gene based on specificity criteria.
```bash
awk 'NR == 1 {print; next} $2 <= -95 && $9 == "sense" && $10 == 39 {print $0 | "sort -k1,1V -k2,2nr"}' ${ACC_NO}_sgrna.tsv | 
awk 'NR == 1 {print; next} {guides[$1]++} guides[$1] <= 10 {print}' > top_ten_sgRNA_outputs.tsv
```

## Key Modifications

- **CN PAM Compatibility:** Specifically adapted for targeting CN PAMs.
- **Extended sgRNA Length:** The design now supports a 32-nucleotide sgRNA length, enhancing target specificity.
- **Upstream Gene Target Optimization:** Focused improvements on identifying upstream gene targets for CRISPRtOE.

### Guide Intersection Analysis

To determine potential insertion sites of CRISPRt within upstream genes, utilize the `GenBankParser` tool from the [barcoder toolkit](https://github.com/ryandward/barcoder).
```bash
pipenv run python find_crisprt_sites.py
```
