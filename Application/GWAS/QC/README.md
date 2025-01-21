
# Qualtiy control (QC)
This folder contains scripts for quality control (QC) of genotype data, including steps for sample QC and marker QC.

## Required files & format
This section lists the necessary files required before proceeding with the subsequent QC steps:
- info.txt (the order of the first 5 columns must not be altered)

    | FID     | IID     | SexCheck | RecordCheck | RelativeCheck | ...
    |---------|---------|----------|-------------|---------------|-----
    | Sample1 | Sample1 |          |             |  remove       |
    | Sample2 | Sample2 |          |             |               |
    | Sample3 | Sample3 | fail     | fail        |               |
    | Sample4 | Sample4 |          | fail        |               |

- all_hg38.bed, .bim, .fam (downloaded from [here](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg))

    for the downloaded files, we extract the variants typed in the TPM arrays

    ```
    awk '{print $1,$4,$4,$2}' tpmi.bim > lst_var.txt

    plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
    plink2 --pfile all_hg38 vzs --allow-extra-chr --extract bed1 lst_var.txt --snps-only --max-alleles 2 --set-all-var-ids @:# --rm-dup exclude-all --make-bed --out all_hg38_extract 
    ```

- phecov.txt

    | FID     | IID     | e11 | i10 | hba1c | sbp  | dbp  | hba1c_res_invnorm | age | sex | bmi  | PC1 | ... | PC10 | ...
    |---------|---------|-----|-----|-------|------|------|-------------------|-----|-----|------|-----|-----|------|-----
    | Sample1 | Sample1 |  1  |  0  | 7.6   | 117  | 71   |   -0.12           | 45  | 1   | 31.5 | 0.1 | ... | 0.7  |
    | Sample2 | Sample2 |  0  |  1  | 5.7   | 151  | 96   |   -0.77           | 50  | 2   | 24.0 | 0.2 | ... | 0.6  |
    | Sample3 | Sample3 |  0  |  0  | 5.3   | 111  | 69   |   -0.95           | 40  | 1   | 19.1 | 0.3 | ... | 0.8  |
    | Sample4 | Sample4 |  1  |  1  | 8.2   | 139  | 86   |    0.05           | 55  | 2   | 26.5 | 0.4 | ... | 0.5  |

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
./QC_sample.sh tpmi info.txt all_hg38_extract all_hg38.psamp
```

## Marker QC
In this step, we perform the marker QC for the **indpendent QC samples** according to different phenotypes of the followings:
- genotype call rate < 0.95
- failed the nonrandom missingness test (binary trait only)
- minor allele frequency < 0.01
- violation of Hardy-Weinberg equilibrium
```
./QC_marker.sh tpmi_qcSampInd phecov.txt e11 b
./QC_marker.sh tpmi_qcSampInd phecov.txt hba1c q
```

## QC data
In this step, we extract the QC data for each phenotype to use in the following analyses
```
for phe in e11 hba1c; do
    plink2 --bfile tpmi_qcSampInd --extract ${phe}/tpmi_qcSampInd.snplist --make-bed --out ${phe}/tpmi_qcSampInd_qcVar
    plink2 --bfile tpmi_qcSampRel --extract ${phe}/tpmi_qcSampInd.snplist --make-bed --out ${phe}/tpmi_qcSampRel_qcVar
done
```

