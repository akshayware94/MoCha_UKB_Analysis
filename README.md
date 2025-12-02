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

**Private resources (requires UKBB approval)** <br />
```
ukb_sqc_v2.txt
ukb_cal_chr{1..22,X,Y,XY,MT}_v2.bed      # 92 GiB
ukb_int_chr{1..22,X,Y,XY,MT}_v2.bin      # 2.9 TiB
```
**Public resources** <br />
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
**Fix truncated Affymetrix manifest sequences** <br />
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

**Build Auxiliary Tools** <br />
Compile three lightweight C tools that enable high-speed processing of large UKBB files: <br />
```unpack``` — decode PLINK .bed into Affymetrix genotype structure <br />
```split``` — divide binary resources into batch files <br />
```dump``` — convert binary floats for SNP posterior reconstruction <br />

- Install dependencies:
  Make sure the GNU C compiler and C Library are available to compile three auxiliary tools (Debian/Ubuntu specific if you have admin privileges)
  
```
sudo apt install --no-install-recommends gcc libc6-dev
```
Generate a small ELF binary file that we will use to efficiently unpack PLINK binary files

```
echo '#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  int i, n, nread, swap;
  if (sscanf(argv[1], "%d", &n) != 1)
    return 1;
  char *in_buffer = (char *)malloc((n + 3) / 4 * sizeof(char));
  char *out_buffer = (char *)malloc((n + 3) / 4 * 4 * sizeof(char));
  FILE *in = fopen(argv[2], "r");
  while ((nread = fread(in_buffer, 1, (n + 3) / 4, stdin))) {
    if (nread != (n + 3) / 4)
      return 1;
    if (fscanf(in, "%d\n", &swap) != 1)
      return 1;
    for (i = 0; i < (n + 3) / 4; i++) {
      char x = in_buffer[i];
      char y = (swap ? (~x & 0xAA) : (x & 0x55) << 1) ^
               ((x & 0x55) ^ ((x & 0xAA) >> 1));
      int z = (int)(y & 0x03) ^ ((int)(y & 0x0C) << 6) ^
              ((int)(y & 0x30) << 12) ^ ((int)(y & 0xC0) << 18) ^ 0x30303030;
      memcpy(&out_buffer[4 * i], &z, 4);
    }
    fwrite(out_buffer, 1, n, stdout);
  }
  free(in_buffer);
  free(out_buffer);
  return 0;
}' > unpack.c
gcc -O3 -o unpack unpack.c
```
Generate a small ELF binary file that we will use to efficiently split binary resources by batches

```
echo '#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  int i, nread, n = 0, m = 4;
  int *offsets = (int *)malloc(m * sizeof(int));
  FILE **batches = (FILE **)malloc(m * sizeof(FILE *));
  FILE *in = fopen(argv[2], "r");
  char buffer[BUFSIZ];
  offsets[0] = 0;
  while (fscanf(in, "%d %s\n", &i, buffer) == 2) {
    sprintf(buffer + strlen(buffer), "%s", argv[1]);
    batches[n] = fopen(buffer, "w");
    n++;
    if (n == m) {
      m <<= 1;
      offsets = (int *)realloc(offsets, m * sizeof(int));
      batches = (FILE **)realloc(batches, m * sizeof(FILE *));
    }
    offsets[n] = offsets[n - 1] + i;
  }
  fclose(in);
  char *line = (char *)malloc(offsets[n] * sizeof(char *));
  while ((nread = fread(line, 1, offsets[n], stdin))) {
    if (nread != offsets[n])
      return 1;
    for (i = 0; i < n; i++)
      fwrite(line + offsets[i], 1, offsets[i + 1] - offsets[i], batches[i]);
  }
  for (i = 0; i < n; i++)
    fclose(batches[i]);
  free(offsets);
  free(batches);
  free(line);
  return 0;
}' > split.c
gcc -O3 -o split split.c
```
Generate a small ELF binary file that we will use to efficiently convert binary floats

```
echo '#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  int i, j, n, m, nread;
  if (sscanf(argv[1], "%d", &n) != 1)
    return 1;
  if (sscanf(argv[2], "%d", &m) != 1)
    return 1;
  char *buffer = (char *)malloc(n * m * sizeof(float));
  while ((nread = fread(buffer, sizeof(float), n * m, stdin))) {
    if (nread != n * m)
      return 1;
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        fprintf(stdout, "%g%c", *(float *)&buffer[4 * (n * j + i)],
                j + 1 == m ? 0x0A : 0x09);
  }
  free(buffer);
  return 0;
}' > dump.c
gcc -O3 -o dump dump.c
```
**Split by Batches**
Split sample report file into batch report files (~2 seconds)
```
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  awk 'NR==FNR {print "cel_files\tcomputed_gender\tcall_rate\tcn-probe-chrXY-ratio_gender_meanX\tcn-probe-chrXY-ratio_gender_meanY" > $1".AxiomGT1.report.txt"}
  NR>FNR {printf "%s\t%s\t%s\t%s\t%s\n",$1,$11,$7,$12,$13 > $4".AxiomGT1.report.txt"}' ukb_snp_posterior.batch -
sed 's/$/.AxiomGT1.report.txt/' ukb_snp_posterior.batch | tr '\n' '\0' | xargs -0 gzip --force --no-name
```
Split chromosome posterior files into batch posterior files (~1 minute)
```
tar xvf ukb_snp_posterior.tar 1>&2
cat ukb_snp_posterior_chr{{1..22},X,Y,XY,MT}.bin | ./split .snp-posteriors.bin <(sed 's/^/132 /' ukb_snp_posterior.batch)
```
Reconstruct Affymetrix SNP posterior files (~30 seconds per batch)
```
awk 'NR==FNR {x[$2]++} NR>FNR && FNR>1 {array=$9==0 || $9==2; print array"\t"$3;
  if ($1 in x) print array"\t"$3":1"}' ukb_snp_posterior_chrX_haploid.bim ukb_snp_qc.txt > snp-posteriors.UKBL.tsv
awk 'NR==FNR {x[$2]++} NR>FNR && FNR>1 {array=$9==1 || $9==2; print array"\t"$3;
  if ($1 in x) print array"\t"$3":1"}' ukb_snp_posterior_chrX_haploid.bim ukb_snp_qc.txt > snp-posteriors.UKBB.tsv
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq | while read array batch; do
  (echo -e "#%SnpPosteriorFormatVer=1\n#%data-order=meanX,varX,nObsMean,nObsVar,meanY,varY,covarXY\nid\tBB\tAB\tAA\tCV";
  cat $batch.snp-posteriors.bin | ./dump 1 33 | paste snp-posteriors.$array.tsv - | sed '/^0/d;s/^..//' | \
  awk '{fmt="%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
  if ($2<$16) printf fmt,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34;
         else printf fmt,$1,$16,$17,$18,$19,$20,$21,$22,$9,$10,$11,$12,$13,$14,$15,$2,$3,$4,$5,$6,$7,$8,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34}') | \
    gzip > $batch.AxiomGT1.snp-posteriors.txt.gz && /bin/rm $batch.snp-posteriors.bin
done
```
Split chromosome genotype calls files into batch genotype calls files (~30 minutes, could be batched by chromosomes)
```
(zcat Array_BIL_34.zip | grep -v ^# | tail -n+2;
 zcat Array_UKB_34.zip | grep -v ^# | tail -n+2) | tr -d '"' | \
  awk -F, '{if ($12==$14 && $13==$15) x=0; else if ($12==$15 && $13==$14) x=1; else x=0; print $1,x}' | \
  awk 'NR==FNR {x[$1]=$2} NR>FNR && FNR>1 {print x[$3]}' - ukb_snp_qc.txt > swap.lines
tail -qc+4 ukb_cal_chr{{1..22},X,Y,XY,MT}_v2.bed | ./unpack $(cat ukb_sqc_v2.txt | wc -l) swap.lines | \
  ./split .calls.bin <(sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | cut -d" " -f4 | uniq -c | sed 's/^ *//')
```
Reconstruct Affymetrix genotype calls files (~10 minutes per batch)
```
awk 'NR>1 {array=$9==0 || $9==2; print array"\t"$3}' ukb_snp_qc.txt > calls.UKBL.tsv
awk 'NR>1 {array=$9==1 || $9==2; print array"\t"$3}' ukb_snp_qc.txt > calls.UKBB.tsv
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq -c | while read n array batch; do
  (echo -en "#Calls: -1=NN, 0=AA, 1=AB, 2=BB\nprobeset_id\t";
  sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
    awk -v batch=$batch '$4==batch {print $1}' | tr '\n' '\t' | sed 's/\t$/\n/';
  fold -w$n $batch.calls.bin | tr '\n' '\r' | fold -w1 | tr '\n\r' '\t\n' | sed 's/3/-1/g' | \
  paste calls.$array.tsv - | sed '/^0/d;s/^..//') | \
    gzip > $batch.AxiomGT1.calls.txt.gz && /bin/rm $batch.calls.bin
done
```
Split chromosome intensity files into batch intensity files (~4 hours, could be batched by chromosomes)
```
cat ukb_int_chr{{1..22},X,Y,XY,MT}_v2.bin | ./split .summary.bin <(sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | cut -d" " -f4 | uniq -c | awk '{print 8*$1,$2}')
```
Reconstruct Affymetrix intensities summary files (~3.5 hours per batch)
```
awk 'NR>1 {array=$9==0 || $9==2; print array"\t"$3"-A"; print array"\t"$3"-B"}' ukb_snp_qc.txt > summary.UKBL.tsv
awk 'NR>1 {array=$9==1 || $9==2; print array"\t"$3"-A"; print array"\t"$3"-B"}' ukb_snp_qc.txt > summary.UKBB.tsv
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq -c | while read n array batch; do
  (echo -en "probeset_id\t";
  sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
    awk -v batch=$batch '$4==batch {print $1}' | tr '\n' '\t' | sed 's/\t$/\n/';
  cat $batch.summary.bin | ./dump 2 $n | paste summary.$array.tsv - | sed '/^0/d;s/^..//') | \
    gzip > $batch.AxiomGT1.summary.txt.gz && /bin/rm $batch.summary.bin
done
```
**Prepare input files required for MoChA WDL** <br />
Create sample table
```
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  awk 'BEGIN {print "sample_id\tbatch_id\tcel\tcomputed_gender\tcall_rate"}
  {printf "%s\t%s\t%s\t%s\t%s\n",$1,$4,$1,$11,$7/100}' > ukb.sample.tsv
```

Create batch table
```
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq -c | \
  awk 'BEGIN {csv["UKBB"]="Axiom_UKB_WCSG.na34.annot.csv.gz"; csv["UKBL"]="Axiom_UKBiLEVE.na34.annot.csv.gz"
  print "batch_id\tn_smpls\tcsv\tsnp\treport\tcalls\tsummary"}
  {printf "%s\t%s\t%s\t%s.AxiomGT1.snp-posteriors.txt.gz\t%s.AxiomGT1.report.txt.gz\t%s.AxiomGT1.calls.txt.gz\t%s.AxiomGT1.summary.txt.gz\n",
  $3,$1,csv[$2],$3,$3,$3,$3}' > ukb.batch.tsv
```
