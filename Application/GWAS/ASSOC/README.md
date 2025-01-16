
## Association test for independent sample dataset
In this step, we conduct the association test for independent QC sample dataset by PLINK
- firth logistic regression for the binary trait
- linear regression for the quantitative trait
```
./ASSOC_sampleInd.sh ./QC/e11/tpmi_qcSampInd_qcVar example.phecov e11 b
./ASSOC_sampleInd.sh ./QC/hba1c/tpmi_qcSampInd_qcVar example.phecov hba1c q
```

## Association test for related sample dataset
In this step, we conduct the association test for the related QC sample dataset by SAIGE
```
./ASSOC_sampleRel.sh ./QC/e11/tpmi_qcSampRel_qcVar example.phecov e11 b
./ASSOC_sampleRel.sh ./QC/hba1c/tpmi_qcSampRel_qcVar example.phecov hba1c q
```
