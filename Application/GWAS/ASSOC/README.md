
# Association test
This folder contains scripts for association tests on binary and quantitative traits for both independent and related sample datasets.

## Required files & format
This section lists the necessary files required before proceeding with the subsequent steps:
- phecov.txt

    | FID     | IID     | e11 | i10 | hba1c | sbp  | dbp  | hba1c_res_invnorm | age | sex | bmi  | PC1 | ... | PC10 | ...
    |---------|---------|-----|-----|-------|------|------|-------------------|-----|-----|------|-----|-----|------|-----
    | Sample1 | Sample1 |  1  |  0  | 7.6   | 117  | 71   |   -0.12           | 45  | 1   | 31.5 | 0.1 | ... | 0.7  |
    | Sample2 | Sample2 |  0  |  1  | 5.7   | 151  | 96   |   -0.77           | 50  | 2   | 24.0 | 0.2 | ... | 0.6  |
    | Sample3 | Sample3 |  0  |  0  | 5.3   | 111  | 69   |   -0.95           | 40  | 1   | 19.1 | 0.3 | ... | 0.8  |
    | Sample4 | Sample4 |  1  |  1  | 8.2   | 139  | 86   |    0.05           | 55  | 2   | 26.5 | 0.4 | ... | 0.5  |

## Association test for independent sample dataset
In this step, we conduct the association test for independent QC sample dataset by PLINK
- firth logistic regression for binary traits
- linear regression for quantitative traits
```
./ASSOC_sampleInd.sh ./QC/e11/tpmi_qcSampInd_qcVar phecov.txt e11 b
./ASSOC_sampleInd.sh ./QC/hba1c/tpmi_qcSampInd_qcVar phecov.txt hba1c_res_invnorm q
```

## Association test for related sample dataset
In this step, we conduct the association test for the related QC sample dataset by SAIGE
```
./ASSOC_sampleRel.sh ./QC/e11/tpmi_qcSampRel_qcVar phecov.txt e11 b
./ASSOC_sampleRel.sh ./QC/hba1c/tpmi_qcSampRel_qcVar phecov.txt hba1c q
```
