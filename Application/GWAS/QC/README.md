
## Sample QC
In this step, we perform the sample QC of the followings:
- inconsistent between EMR-gender and genetic gender
- inconsistent EMR records for the participants with highly similar genomes (might be the same ones in different hospitals)
- genotype call rate < 0.95
- excessive or reduced autosomal heterozygosity rate
- cryptic relatedness
- divergent ancestry

Finally, we will create two kinds of datasets, independent sample and related sample datasets, for the subsequent processing
```
./QC_sample.sh tpmi example.tsv 1kgp? 1kgp?
```
fpath_geno=$1 # plink bfile name
fpath_sinf=$2 # tpmi037_486956_Sex_Dup_CR_HomoRate_relatedness.tsv
fpath_popu_dat=$3 # all_hg38_extract_tpm12_indep
fpath_popu_inf=$4 # relationships_w_pops.txt

## Marker QC
In this step, we perform the marker QC for the indpendent QC samples according to different phenotypes of the followings:
- genotype call rate < 0.95
- failed the nonrandom missingness test (binary trait only)
- minor allele frequency < 0.01
- violation of Hardy-Weinberg equilibrium
```
./QC_marker.sh tpmi_qcSampInd example.phecov e11 b
./QC_marker.sh tpmi_qcSampInd example.phecov hba1c q
```

## QC data
In this step, we extract the QC data for each phenotype for the subsequent analyses
```
for phe in e11 hba1c; do
    plink2 --bfile tpmi_qcSampInd --extract ${phe}/tpmi_qcSampInd.snplist --make-bed --out ${phe}/tpmi_qcSampInd_qcVar
    plink2 --bfile tpmi_qcSampRel --extract ${phe}/tpmi_qcSampInd.snplist --make-bed --out ${phe}/tpmi_qcSampRel_qcVar
done
```

