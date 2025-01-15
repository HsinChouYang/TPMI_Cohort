
```
./QC_sample.sh tpmi example.tsv 1kgp? 1kgp?
```
fpath_geno=$1 # plink bfile name
fpath_sinf=$2 # tpmi037_486956_Sex_Dup_CR_HomoRate_relatedness.tsv
fpath_popu_dat=$3 # all_hg38_extract_tpm12_indep
fpath_popu_inf=$4 # relationships_w_pops.txt

```
./QC_marker.sh tpmi_qcSampInd example.tsv e11 b
./QC_marker.sh tpmi_qcSampInd example.tsv i10 b
./QC_marker.sh tpmi_qcSampInd example.tsv hba1c q
./QC_marker.sh tpmi_qcSampInd example.tsv sbp q
./QC_marker.sh tpmi_qcSampInd example.tsv dbp q
```

```
plink2 --bfile tpmi_qcSampInd --extract e11/tpmi_qcSampInd.snplist --make-bed --out e11/tpmi_qcSampInd_qcVar
plink2 --bfile tpmi_qcSampInd --extract i10/tpmi_qcSampInd.snplist --make-bed --out i10/tpmi_qcSampInd_qcVar
plink2 --bfile tpmi_qcSampInd --extract hba1c/tpmi_qcSampInd.snplist --make-bed --out hba1c/tpmi_qcSampInd_qcVar
plink2 --bfile tpmi_qcSampInd --extract sbp/tpmi_qcSampInd.snplist --make-bed --out sbp/tpmi_qcSampInd_qcVar
plink2 --bfile tpmi_qcSampInd --extract dbp/tpmi_qcSampInd.snplist --make-bed --out dbp/tpmi_qcSampInd_qcVar
```

