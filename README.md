# Tools

* [PLINK](https://www.cog-genomics.org/plink/)
* [PLINK2](https://www.cog-genomics.org/plink/2.0/)
* [KING](https://www.kingrelatedness.com/)
* [ADMIXTURE](https://dalexander.github.io/admixture/)

# Packages
## Python
* pandas
* matplotlib
* numpy

## R
* data.table
* dplyr
* car
* sp
* pROC
* RColorBrewer

# Download data
1. [1000 Genomes phase 3](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg)

# Required files & format
1. batchID.clst
```
FID	IID	Batch
Sample1	Sample1	batch1
Sample2	Sample2	batch1
Sample3	Sample3	batch2
Sample4	Sample4	batch2
```
2. batchSNPs.zero
```
rs1111	batch3
rs2222	batch4
rs2222 	batch6
rs2222	batch10
rs3333	batch2
rs4444	batch6
```
# Data Pre-processing
1. Genetic Data
- Batch effect detection
- Batch SNP removal
- Relatedness check
2. EHR Data
- Diagnosis data
- Laboratory test results
- Vital signs

# Population structure
## PCA
- SNP pruning
```
./plink2 --bfile infile \
--exclude range LRLD.txt \
--indep-pairwise 5000kb 1 0.2 \
--out outfile
```
- PCA of Birth Cohorts Before 1950
- PCA projection
2. Admixture
- Baseline Samples
- Run Admixture
- Projection

# Application
1. GWAS
2. Ploygenetic risk score

