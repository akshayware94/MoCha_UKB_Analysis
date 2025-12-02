# *MoCha_UKB_Analysis*
## Overview
This repository contains workflows, scripts, and documentation for running the MoChA (Mosaic Chromosomal Alterations) WDL pipeline on UK Biobank genotype array data.
It includes:

1. Tools and instructions to reconstruct raw Affymetrix CEL-like files for each of the 106 UKB genotyping batches using ukb2txt.

2. A reproducible workflow for calling mosaic chromosomal alterations (including mLOY, mLOX, and autosomal mCAs) using the official MoChA WDL pipeline.

3. Supplementary scripts for preprocessing, QC, formatting input files, and downstream analysis.

This repository is designed as a complete end-to-end resource for large-scale mosaic event analysis in UK Biobank.

## Installation
Below are the core dependencies required to reconstruct UKB raw intensity files and run the MoChA WDL pipeline. <br />

## System Requirements

1. Unix-based OS (Linux recommended)
2. At least 6-8 TB free disk space for reconstruction and MoChA outputs
3. Access to UKB Project Space or local UKB genotype data

>
  ```
  1. JAVA
  2. Docker
  ```
