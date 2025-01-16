
# PRS Calculation and Modeling
This folder contains scripts for calculating polygenic risk scores (PRS) and modeling the calculated scores.

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
Rscript PRS.R
```

