# Tools

* [PLINK](https://www.cog-genomics.org/plink/)
* [PLINK2](https://www.cog-genomics.org/plink/2.0/)
* [KING](https://www.kingrelatedness.com/)
* [ADMIXTURE](https://dalexander.github.io/admixture/)
* [SAIGE](https://github.com/saigegit/SAIGE)
* [PRS-CSx](https://github.com/getian107/PRScsx)

# Packages
## Python
* os
* pandas
* matplotlib
* numpy
* seaborn


## R
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
1. Genetic Data
- Batch effect detection
- Batch SNP removal
## Required files & format
- batchID.clst

|FID	    |IID      |Batch   |
|---------|---------|----------|
|Sample1	|Sample1	|batch1|
|Sample2	|Sample2	|batch1|
|Sample3	|Sample3	|batch2|
|Sample4	|Sample4	|batch2|

- batchSNPs.zero 
|aa           |bb         |
|------------|----------|
| rs1111	| batch3  |
| rs2222	| batch4  |
| rs2222 	| batch6  |
| rs2222	| batch10|
| rs3333	| batch2  |
| rs4444	| batch6  |

```
./plink --bfile tpmi \
--within batchID.clst \
--zero-cluster batchSNPs.zero \
--make-bed \
--out tpmi_batchSNP
```

- Relatedness check
```
./king -b ForKing.bed \
--cpus 50 \
--related \
--degree 3 \
--prefix TPMI_Relatedness
```
2. EMR Data
- Diagnosis reports
- Laboratory test results
- Vital signs

# Population structure
## Principal Component Analysis(PCA)
- step1 : SNP pruning  
In this step, we aim to generate a pruned subset of variants for PCA. The command excludes the long-range linkage disequilibrium (LRLD) regions specified in LRLD.txt and performs an independent pairwise pruning with a window size of 5000kb, step size of 1, and a linkage disequilibrium (LD) threshold of 0.2.

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
To select samples for Admixture analysis, the PC1 values on the X-axis and PC2 values on the Y-axis were divided into thirds from the 2nd to the 98th percentile. This division created a 3x3 grid, resulting in 9 blocks. **[Reference samples](https://github.com/Jenn-Hwai/TPMI_Cohort/tree/main/Population%20structure/Admixture/AdmixtureSample.py)** were then randomly selected from these 9 blocks.
![Reference sample](Population%20structure/Admixture/RefSample.png)
- Run Admixture
```
./RunAdmixture.sh
```

- Projection
```
./admixture -P ForProjection.bed 10 -j10
```

# Application
## GWAS
- QC
- Associaton
## post-GWAS
- Polygenic risk score

