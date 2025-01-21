
# Polygenic risk scores (PRS) calculation and modeling
This folder contains scripts for calculating polygenic risk scores (PRS) and modeling the calculated scores.

## Required files & format
This section lists the necessary files required before proceeding with the subsequent steps:
- phecov.txt

    | FID     | IID     | e11 | i10 | hba1c | sbp  | dbp  | hba1c_res_invnorm | age | sex | bmi  | PC1 | ... | PC10 | ...
    |---------|---------|-----|-----|-------|------|------|-------------------|-----|-----|------|-----|-----|------|-----
    | Sample1 | Sample1 |  1  |  0  | 7.6   | 117  | 71   |   -0.12           | 45  | 1   | 31.5 | 0.1 | ... | 0.7  |
    | Sample2 | Sample2 |  0  |  1  | 5.7   | 151  | 96   |   -0.77           | 50  | 2   | 24.0 | 0.2 | ... | 0.6  |
    | Sample3 | Sample3 |  0  |  0  | 5.3   | 111  | 69   |   -0.95           | 40  | 1   | 19.1 | 0.3 | ... | 0.8  |
    | Sample4 | Sample4 |  1  |  1  | 8.2   | 139  | 86   |    0.05           | 55  | 2   | 26.5 | 0.4 | ... | 0.5  |

## PRS Calculation
In this step, we obtain the variant weights in two ways:
- manually calculated variant weights by PRS-CSx
- externally sourced variant weights from PGS Catalog
```
./PRS.sh
```
## Modeling
In this step, we perform a logistic regression ananlysis to evaulate and compare the AUC of two models:
- $logit(P(Y=1)) = \beta_0 + \beta_1 \cdot \text{PRS}$
- $logit(P(Y=1)) = \beta_0 + \beta_1 \cdot \text{PRS} + \beta_2 \cdot \text{Age} + \beta_3 \cdot \text{Gender} + \beta_4 \cdot \text{BMI}$
```
R CMD BATCH --args --fname_prs=e11_DIAGRAM_prscsx.sscore --fname_phecov=phecov.txt --phe=e11 --covs=age,sex,bmi PRS.R
# R CMD BATCH --args --fname_prs=e11_PGS002308_prscsx.sscore --fname_phecov=phecov.txt --phe=e11 --covs=age,sex,bmi PRS.R

```

