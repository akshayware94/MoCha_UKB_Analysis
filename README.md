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

### System Requirements

1. Unix-based OS (Linux recommended)
2. At least 6-8 TB free disk space for reconstruction and MoChA outputs
3. Access to UKB Project Space or local UKB genotype data

>
  ```
1. Java 8+        (needed for WDL and Cromwell)
2. Docker         (required for running WDL tasks)
3. Cromwell       (for executing WDL workflows)
4. ukb2txt        (to reconstruct Affymetrix raw intensity files)
5. bcftools       (VCF processing)
6. plink/plink2   (genotype file manipulation)
7. shapeit4       (phasing)
  ```
### Install Example (Ubuntu)

```
# Install Java
sudo apt-get update
sudo apt-get install -y default-jre

# Install Docker
sudo apt-get install -y docker.io
sudo usermod -aG docker $USER

# Install Cromwell
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar

# Install bcftools and plink
sudo apt-get install -y bcftools plink
```

## ukb2txt — Phase 1: Reconstructing Raw Affymetrix Files for UK Biobank

This module documents Phase 1 of the UK Biobank mosaic-chromosomal-alteration pipeline:
1. Reconstructing raw Affymetrix Axiom batch-specific files using ukb2txt for downstream use in the MoChA WDL pipeline.
2. UK Biobank provides genotype and intensity data in non-standard merged formats. The ukb2txt workflow reverses these transformations to recover the original Affymetrix structure (per-batch SNP posterior, calls, summary / intensity) required by MoChA.

Unlike other biobanks (FINNGEN, MVP) that distribute microarray batch files directly, UK Biobank:

- Merges intensities across all batches into .bin files
- Encodes genotypes using GRCh37 ref/alt alleles rather than A/B alleles
- Removes probesets such as JAK2 V617F, leaving intensities unavailable in standard output
- Provides LRR/BAF files with missing intensities for many SNPs

Because MoChA requires accurate A/B allele designations and batch-specific intensities, we reconstruct the raw Affymetrix files using the official manifests and UKBB binary sources.
This phase outputs 2.4 TiB of gzipped Affymetrix batch files that can be consumed directly by the MoChA WDL workflow.

### System Requirements
UK Biobank provides genotype and intensity resources in merged chromosome-wide binary formats:

```ukb_cal_chrN_v2.bed``` — genotype calls <br />
```ukb_int_chrN_v2.bin``` — probe intensities <br />
```ukb_sqc_v2.txt``` — sample tracker

These formats differ significantly from standard Affymetrix Axiom batch-level deliverables:

```
*.AxiomGT1.report.txt
*.AxiomGT1.snp-posteriors.txt
*.AxiomGT1.calls.txt
*.AxiomGT1.confidences.txt
*.AxiomGT1.summary.txt
```

To run MoChA (or any Axiom-based CNV/LOH pipeline), we must:
- Recover A/B allele orientation using array manifests
- Split merged binary files into 106 original UKBB batches
- Reconstruct all Affymetrix output files with correct formats

### Download Resources
You will need both private **UKBB resources** and public **Affymetrix/UKBB metadata**.

**Private resources (requires UKBB approval)**
```
ukb_sqc_v2.txt
ukb_cal_chr{1..22,X,Y,XY,MT}_v2.bed      # 92 GiB
ukb_int_chr{1..22,X,Y,XY,MT}_v2.bin      # 2.9 TiB
```
**Public resources**
Use the following commands:
```
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior_chrX_haploid.bim
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior.tar
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_bim.tar
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior.batch
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/ukb_genetic_data_description.txt
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/Array_BIL_34.zip
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/Array_UKB_34.zip
```
**Fix truncated Affymetrix manifest sequences**
For some reasons, the length of the flanking sequence field in the manifest files is capped at 250 characters. This causes the flanking sequences of three Affymetrix indels (Affx-89015252, Affx-92046163, Affx-92047500) to be truncated. The following commands can fix this issue
```
zcat Array_BIL_34.zip | \
  sed -e 's/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTAC/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTACTTTTTTGCATGCTATCAAAATCACTGCCAATAATATGAATGTCAGGTCAATTTCCATAGGTAAATCCGTTACCTTTTACCTCTTTAAAAGAAAAGTTTTGCAGAAGAGGGCTGAAAATTTCTTGAGCCATTTCAGCACAAAGAATGGAAGTTCATTTCTCACCATGATACAACTCTACCCTGCTGTCATCTTCATTGTGATGGTGGCAGAAGTTTAGCAGGGTGCAAGTGACCACTA\]AATGACATCTTTTCATGAACAATTGATAAATCTTT/' | \
  gzip > Axiom_UKBiLEVE.na34.annot.csv.gz
zcat Array_UKB_34.zip | \
  sed -e 's/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTAC/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTACTTTTTTGCATGCTATCAAAATCACTGCCAATAATATGAATGTCAGGTCAATTTCCATAGGTAAATCCGTTACCTTTTACCTCTTTAAAAGAAAAGTTTTGCAGAAGAGGGCTGAAAATTTCTTGAGCCATTTCAGCACAAAGAATGGAAGTTCATTTCTCACCATGATACAACTCTACCCTGCTGTCATCTTCATTGTGATGGTGGCAGAAGTTTAGCAGGGTGCAAGTGACCACTA\]AATGACATCTTTTCATGAACAATTGATAAATCTTT/' \
      -e 's/TATTGGCTCAATTTCTTGGGGAGGGGGTGCTGTCA\[-\/GAGATTGTTATCTGAGGATGTGACATAGATTTCTCAGGGCACAATTTCAACTACTTTTTCAGCTTTAGGGTTTTTAGATACGTTTGTACCACAATTGAGCATGGGAGGGAGAGGGGTGAGCCTAAGCAGTGATGGCTGATTTCTGTCATGTCTGTCATGTGTCCCCCAGTACCTCCAGAGGTAACTGTGCTCACAAACAGCCCTGTGGAACT/TATTGGCTCAATTTCTTGGGGAGGGGGTGCTGTCA\[-\/GAGATTGTTATCTGAGGATGTGACATAGATTTCTCAGGGCACAATTTCAACTACTTTTTCAGCTTTAGGGTTTTTAGATACGTTTGTACCACAATTGAGCATGGGAGGGAGAGGGGTGAGCCTAAGCAGTGATGGCTGATTTCTGTCATGTCTGTCATGTGTCCCCCAGTACCTCCAGAGGTAACTGTGCTCACAAACAGCCCTGTGGAACTGAGAGAGCCCAACGTCCTCATCTGTTTCATA\]GACAAGTTCACCCCACCAGTGGTCAATGTCACGTG/' \
      -e 's/GTAATAGGACTCCTTGAAGAGACTTTCGGCGGGGC\[-\/CTGGGGAAGACAGAGGGAGACACACTGAACAGTCTGGCTGTGACCTCCAGGGCCACTGTCACCCACACAGCTACCAGGTGGCCAGAGCCAGGATCTGAACCCAGGTCTGTGGGGGATCCACACCTGAATCCCATTCTTGGGGGAGTCTCATTGGCACCACAGCAGAGGAACCTCTAACCTAGGCCTTCGTTCAAGACTAGAACCTGCCCCCA/GTAATAGGACTCCTTGAAGAGACTTTCGGCGGGGC\[-\/CTGGGGAAGACAGAGGGAGACACACTGAACAGTCTGGCTGTGACCTCCAGGGCCACTGTCACCCACACAGCTACCAGGTGGCCAGAGCCAGGATCTGAACCCAGGTCTGTGGGGGATCCACACCTGAATCCCATTCTTGGGGGAGTCTCATTGGCACCACAGCAGAGGAACCTCTAACCTAGGCCTTCGTTCAAGACTAGAACCTGCCCCCACTCCAGTGCTGACCACTTCAGAGCAGAGGGGAGGCTGAAGAGGACACAGGGTCCTCAGTGTCCCAATGCCAGATCCCCACTCTCCTTGGTCACCAGCTTGTGAATCTGGGCAGTCGCCTGGCTCCTGCCTACTGTCCTGAGCCATGTTTCAGAGGGCAGGTAACAAATGAGAAGGGAAAAGTACAGCTCTAGTTCGGGGGGTGGGAGGCCGCTCTATCCTTTACTCTGAAGGCCTGGGGGAGGCTGACCTCCAGACCTGCAGCTGCCAGAAAACCCTGGGGCCCATCCACTGCTTAC\]CCATGGGGTCTGAGGAGTCAGTGATGATCACGTCG/' | \
  gzip > Axiom_UKB_WCSG.na34.annot.csv.gz
```
While Affymetrix orders marker's alleles according to an internal designation of alleles as A and B, the UK biobank has reordered the alleles for each marker in the genotype and SNP posterior (but not intensity) files as reference and alternate with respect to GRCh37. We will therefore recover information about which markers have A and B alleles swapped when compared to reference and alternate alleles and use this information to recover the original genotypes. The same list of markers were swapped for SNP posterior files, with the exception of ten Affymetrix indels (AX-82920347, AX-82999350, AX-83057578, AX-83106285, AX-83149193, AX-83166262, AX-83197070, AX-83253475, AX-83575267, AX-83587806) which, for unclear reasons, were inconsistently swapped between genotype files and SNP posterior files and incorrectly so in the SNP posterior files

**Build Auxiliary Tools**
Compile three lightweight C tools that enable high-speed processing of large UKBB files: <br />
```unpack``` — decode PLINK .bed into Affymetrix genotype structure <br />
```split``` — divide binary resources into batch files <br />
```dump``` — convert binary floats for SNP posterior reconstruction <br />
