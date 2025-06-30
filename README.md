# Tools

* [PLINK](https://www.cog-genomics.org/plink/)
* [PLINK2](https://www.cog-genomics.org/plink/2.0/)
* [bcftools](https://github.com/samtools/bcftools)
* [KING](https://www.kingrelatedness.com/)
* [ADMIXTURE](https://dalexander.github.io/admixture/)
* [REGENIE](https://github.com/rgcgithub/regenie)
* [SAIGE](https://github.com/saigegit/SAIGE)
* [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)
* [PRS-CSx](https://github.com/getian107/PRScsx)
* [fineSTRUCTURE v2 & GLOBETROTTER](https://people.maths.bris.ac.uk/~madjl/finestructure/index.html)
* [admixtools 2](https://uqrmaie1.github.io/admixtools/articles/recipes.html)

# Packages
## Python
* os
* pandas
* matplotlib
* numpy
* seaborn
* glob


## R
* future.apply
* data.table
* dplyr
* car
* sp
* pROC
* RColorBrewer

# Download data
1. [1000 Genomes phase 3](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg)
2. [PGS Catalog](https://www.pgscatalog.org/downloads/)


# Data Pre-processing
## Required files & format
- batchAAF.tsv

| CHROM | POS   | SNPID          | REF | ALT  | batch1   | batch2   | batch3 | batch4 | batch5 | batch6 | ... |
| ----- | ----- | -------------- | --- | ---- | -------- | -------- | ------ | ------ | ------ | ------ | --- |
| 1     | 13417 | rs777038595    | -   | GAGA | 0.000808 | 0        | 0      | 0      | 0      | 0      | ... |
| 1     | 30794 | Affx-474295932 | G   | A    | 0        | 0        | 0      | 0      | 0      | 0      | ... |
| 1     | 46434 | Affx-474295945 | A   | T    | 0.000403 | 0.000239 | 0      | 0      | 0      | 0      | ... |
| 1     | 74636 | rs62641290     | G   | A    | 1        | 1        | 1      | 1      | 1      | 1      | ... |
| 1     | 90081 | Affx-474296182 | A   | T    | 0        | 0        | 0      | 0      | 0      | 0      | ... |

- batchSampleSize.tsv

| CHROM | POS   | SNPID          | REF | ALT  | batch1 | batch2 | batch3 | batch4 | batch5 | batch6 | ... |
| ----- | ----- | -------------- | --- | ---- | ------ | ------ | ------ | ------ | ------ | ------ | --- |
| 1     | 13417 | rs777038595    | -   | GAGA | 2474   | 2094   | 2939   | 2389   | 3232   | 2457   | ... |
| 1     | 30794 | Affx-474295932 | G   | A    | 2480   | 2093   | 2939   | 2389   | 3233   | 2460   | ... |
| 1     | 46434 | Affx-474295945 | A   | T    | 2481   | 2094   | 2941   | 2387   | 3233   | 2459   | ... |
| 1     | 74636 | rs62641290     | G   | A    | 2478   | 2091   | 2939   | 2389   | 3232   | 2450   | ... |
| 1     | 90081 | Affx-474296182 | A   | T    | 2480   | 2091   | 2936   | 2391   | 3233   | 2457   | ... |

- g1k_EAS_unrl.afreq

| CHROM | POS   | ID           | REF | ALT | ALT_FREQS   | OBS_CT |
| ----- | ----- | ------------ | --- | --- | ----------- | ------ |
| 1     | 10397 | 1:10399:C:A  | C   | A   | 0.00396825  | 1008   |
| 1     | 10403 | 1:10405:A:AC | A   | AC  | 0           | 1008   |
| 1     | 10420 | 1:10420:A:C  | A   | C   | 0.000992063 | 1008   |
| 1     | 10438 | 1:10438:A:T  | A   | T   | 0.000992063 | 1008   |
| 1     | 10440 | 1:10440:C:A  | C   | A   | 0.0119048   | 1008   |
| 1     | 10444 | 1:10444:T:A  | T   | A   | 0.00198413  | 1008   |

- twb_af.tsv

| Chr   | Pos_start | Pos_end | Ref     | Obs | Freq   |
| ----- | --------- | ------- | ------- | --- | ------ |
| chr10 | 10319     | 10319   | C       | C   | 0.999  |
| chr10 | 10374     | 10380   | TTAACCC | T   | 0.0012 |
| chr10 | 10425     | 10425   | A       | C   | 0.0988 |
| chr10 | 10466     | 10466   | G       | C   | 0.0015 |
| chr10 | 10501     | 10501   | G       | T   | 0.014  |
| chr10 | 10505     | 10505   | C       | T   | 0.005  |


- batchID.clst

| FID     | IID     | Batch  |
| ------- | ------- | ------ |
| Sample1 | Sample1 | batch1 |
| Sample2 | Sample2 | batch1 |
| Sample3 | Sample3 | batch2 |
| Sample4 | Sample4 | batch2 |

- batchSNPs.zero

| SNP    | batch   |
| ------ | ------- |
| rs1111 | batch3  |
| rs2222 | batch4  |
| rs2222 | batch6  |
| rs2222 | batch10 |
| rs3333 | batch2  |
| rs4444 | batch6  |

- diagnoses_all.csv

| ID      | ICD10 CODE | DIAGNOSIS DATE | SOURCE     | DIVISION        |
| ------- | ---------- | -------------- | ---------- | --------------- |
| Sample1 | J00        | 2023/11/1      | Outpatient | Family Medicine |
| Sample1 | D23.9      | 2023/1/6       | Outpatient | Dermatology     |
| Sample1 | D23.9      | 2023/1/13      | Outpatient | Dermatology     |
| Sample1 | D23.9      | 2023/4/14      | Outpatient | Dermatology     |
| Sample1 | E66.9      | 2023/1/17      | Outpatient | Rheumatology    |
| Sample2 | E66.9      | 2023/2/21      | Outpatient | Rheumatology    |
| Sample2 | E66.9      | 2023/4/18      | Outpatient | Rheumatology    |

- lab_tests.csv

| ID      | Sampling Date | HB   | HB_UNIT | HB_REF    | HB_SPECIMEN_TYPE | HBA1C | HBA1C_UNIT | HBA1C_REF | HBA1C_SPECIMEN_TYPE | ... |
| ------- | ------------- | ---- | ------- | --------- | ---------------- | ----- | ---------- | --------- | ------------------- | --- |
| Sample1 | 2023-03-24    | 13.8 | g/dL    | 12.0-16.0 | Blood            | 6.2   | %          | 4.0-6.0   | Blood               | ... |
| Sample1 | 2023-04-06    |      | g/dL    | 12.0-16.0 | Blood            | 6.3   | %          | 4.0-6.0   | Blood               | ... |
| Sample1 | 2023-04-25    | 15.1 | g/dL    | 12.0-16.0 | Blood            | 6.8   | %          | 4.0-6.0   | Blood               | ... |
| Sample2 | 2021-01-21    | 14.2 | g/dL    | 12.0-16.0 | Blood            | 4.5   | %          | 4.0-6.0   | Blood               | ... |
| Sample2 | 2022-06-29    | 14.8 | g/dL    | 12.0-16.0 | Blood            |       | %          | 4.0-6.0   | Blood               | ... |


## Genetic Data
- Batch effect detection   
Detection of potential batch effect SNPs through analysis of between-batch frequency differences and comparison with other databases (Taiwan BioBank and 1000 Genomes Project).   
**[batchSNP.R](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Data%20Pre-processing/Genetic%20Data/batchSNP.R)**

- Batch SNP removal    
To maintain data quality by converting genotype calls to no-call status for SNPs that demonstrate potential batch effects.
```
./plink --bfile tpmi \
--within batchID.clst \
--zero-cluster batchSNPs.zero \
--make-bed \
--out tpmi_batchSNP
```

- Relatedness check
This step in the process focuses on recognizing the degree of kinship relationships, specifically identifying relationships within the third degree of kinship.
```
./king -b ForKing.bed \
--cpus 50 \
--related \
--degree 3 \
--prefix TPMI_Relatedness
```
## EMR Data
- Diagnosis reports   
For the diagnosis reports, we analyze the most frequent diagnostic codes in the database and calculate the age at first diagnosis.   
**[Diagnosis.py](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Data%20Pre-processing/EHR%20Data/Diagnosis.py)**

- Laboratory test results   
For laboratory measurements, we calculate the statistical distribution of data volume per patient and follow-up duration.   
**[Laboratory.py](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Data%20Pre-processing/EHR%20Data/Laboratory.py)**

# Population structure
## Principal Component Analysis(PCA)
- step1 : SNP pruning  
In this step, we aim to generate a pruned subset of variants for PCA. The command excludes the long-range linkage disequilibrium (LRLD) regions specified in **[LRLD.txt](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Data%20Pre-processing/Genetic%20Data/LRLD.txt)** and performs an independent pairwise pruning with a window size of 5000kb, step size of 1, and a linkage disequilibrium (LD) threshold of 0.2.

```
./plink2 --bfile tpmi \
--exclude range LRLD.txt \
--indep-pairwise 5000kb 1 0.2 \
--out tpmi_forPCA
```
- step2 : PCA of Birth Cohorts Before 1950  
In this step, we use participants born before 1950 (listed in G50.IDs) to perform PCA. The command extracts a pruned subset of variants from the file tpmi_forPCA.prune.in, keeps only the individuals listed in G50.IDs, and calculates allele frequencies with the --freq counts option. The PCA is then performed using the approximate method with biallelic variant weights, and the results are saved to the output file G50.
```
./plink2 --bfile tpmi \
--extract tpmi_forPCA.prune.in \
--keep G50.IDs \
--freq counts \
--pca approx biallelic-var-wts \
--out G50
```

- step3 : PCA projection  
In this step, we use the allele frequencies and principal component eigenvectors computed in the previous step to project all participants onto the PCA coordinates.
```
./plink2 --bfile tpmi \
--read-freq G50.acount \
--score G50.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize \
--score-col-nums 5-14 \
--out PCA_projection
```

## Admixture
- Reference samples
To select samples for Admixture analysis, the PC1 values on the X-axis and PC2 values on the Y-axis were divided into thirds from the 2nd to the 98th percentile. This division created a 3x3 grid, resulting in 9 blocks. **[Reference samples](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Population%20structure/Admixture/RefSample.png)** were then randomly selected from these 9 blocks.   
**[AdmixtureSample.py](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Population%20structure/Admixture/AdmixtureSample.py)**
- Run Admixture   
**[RunAdmixture.sh](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Population%20structure/Admixture/RunAdmixture.sh)**

- Projection
```
./admixture -P ForProjection.bed 10 -j10
```
## ChromoPainter & fineSTRUCTURE
-  step1 : Estimation of haplotypes (phasing)
```
shapeit --input-bed ${infile}.bed ${infile}.bim ${infile}.fam \
        -M ${GeneticMap}/chr${chr}.b38.gmap \
        -O ${outfile}.phased
```
- step 2 : Convert to ChromoPainter's PHASE and RECOMBFILES files
```
perl impute2chromopainter.pl -J ${infile} ${outfile}
```
- step 3 : Create recombination rate map file
```
perl makeuniformrecfile.pl ${infile} ${outfile}
```
- step 4 : Computational stage
```
./fs_linux_glibc2.3 ${cpfile}.cp -n -phasefiles ${phasefile} \
                                    -recombfiles ${recombfile} \
                                    -idfile ${idsfile} \
                                    -s1minsnps 5000 -s3iters 10000 -s4iters 10000 -go
```
## AdmixTools
**[ADMIXTOOLS_code.r](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Population%20structure/AdmixTools/ADMIXTOOLS_code.r)**

# Application
## GWAS   
- **[QC](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Application/GWAS/QC)**   
- **[Associaton](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Application/GWAS/ASSOC)**

## post-GWAS
- **[Polygenic risk score](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Application/post-GWAS/PRS)**

