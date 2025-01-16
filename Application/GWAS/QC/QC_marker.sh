
#!/bin/bash

fpath_geno=$1 # *_qcSampInd
fpath_phecov=$2 # tpmi037_pcJH.fpath_phecov
phe=$3 # e11, colname in $2
ttype=$4 # b (binary) or q (quantitative)

fname=$(basename $fpath_geno)

## (7,8,final) GCR, MAF & HWE
mkdir ${phe}

plink2 --bfile ${fpath_geno} --pheno ${fpath_phecov} --pheno-name $phe --1 --require-pheno --write-samples --out ${phe}/keepInd_phe.txt # --1 has no effect on quantitative trait
if [[ ${ttype} == "b" ]]; then # binary
	plink2 --bfile ${fpath_geno} --keep ${phe}/keepInd_phe.txt --geno 0.05 --test-missing --out ${phe}/${fname}_6
	n=$(wc -l < ${fpath_geno}.bim)
	thres=$(bc -l <<< 0.05/$n)
	awk -v thres="$thres" '$5<thres {print $2}' ${phe}/${fname}_6.missing > ${phe}/rmSNP_missingTest.txt
	plink2 --bfile ${fpath_geno} --keep ${phe}/keepInd_phe.txt --geno 0.05 --exclude ${phe}/rmSNP_missingTest.txt --write-snplist --out ${phe}/${fname}_7

	plink2 --bfile ${fpath_geno} --keep ${phe}/keepInd_phe.txt --extract ${phe}/${fname}_7.snplist --maf 0.01 --write-snplist --out ${phe}/${fname}_8

	n=$(wc -l < ${phe}/${fname}_8.snplist)
	thres=$(bc -l <<< 0.05/$n)
	plink2 --bfile ${fpath_geno} --pheno ${fpath_phecov} --keep-if "${phe}==0" --extract ${phe}/${fname}_8.snplist --hwe $thres --write-snplist --out ${phe}/${fname}_9
else # quantitative
	plink2 --bfile ${fpath_geno} --keep ${phe}/keepInd_phe.txt --geno 0.05 --write-snplist --out ${phe}/${fname}_7

	plink2 --bfile ${fpath_geno} --keep ${phe}/keepInd_phe.txt --extract ${phe}/${fname}_7.snplist --maf 0.01 --write-snplist --out ${phe}/${fname}_8

	n=$(wc -l < ${phe}/${fname}_8.snplist)
	thres=$(bc -l <<< 0.05/$n)
	plink2 --bfile ${fpath_geno} --pheno ${fpath_phecov} --pheno-name $phe --require-pheno --extract ${phe}/${fname}_8.snplist --hwe $thres --write-snplist --out ${phe}/${fname}_9
fi

# plink2 --bfile ${fpath_geno} --keep ${fpath_samp} --extract ${phe}/${fname}_9.snplist --make-bed --out ${phe}/${fname}_qcSampInd_qcVar
# plink2 --bfile ${fpath}_qcSampRel --keep ${fpath_samp} --extract ${phe}/${fname}_9.snplist --make-bed --out ${phe}/${fname}_qcSampRel_qcVar


