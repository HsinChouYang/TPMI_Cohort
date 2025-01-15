
fpath_geno=$1 # *_qcSampInd
fpath_phecov=$2 # tpmi037_pcJH.fpath_phecov

fname=$(basename $fpath_geno)

## (7,8,final) GCR, MAF & HWE
for phe in e11 i10 hba1c sbp dbp; do
	mkdir ${phe}

	plink2 --bfile ${fpath}_qcSampInd --pheno ${fpath_phecov} --pheno-name $phe --1 --require-pheno --write-samples --out ${phe}/keepInd_phe.txt # --1 has no effect on quantitative trait
	if [[ ${phe} =~ ^(e11|i10)$ ]]; then # binary
		plink2 --bfile ${fpath}_qcSampInd --keep ${phe}/keepInd_phe.txt --geno 0.05 --test-missing --out ${phe}/${fname}_6
		n=$(wc -l < ${fpath}_qcSampInd.bim)
		thres=$(bc -l <<< 0.05/$n)
		awk -v thres="$thres" '$5<thres {print $2}' ${phe}/${fname}_6.missing > ${phe}/rmSNP_missingTest.txt
		plink2 --bfile ${fpath}_qcSampInd --keep ${phe}/keepInd_phe.txt --geno 0.05 --exclude ${phe}/rmSNP_missingTest.txt --write-snplist --out ${phe}/${fname}_7

		plink2 --bfile ${fpath}_qcSampInd --keep ${phe}/keepInd_phe.txt --extract ${phe}/${fname}_7.snplist --maf 0.01 --write-snplist --out ${phe}/${fname}_8

		n=$(wc -l < ${phe}/${fname}_8.snplist)
		thres=$(bc -l <<< 0.05/$n)
		plink2 --bfile ${fpath}_qcSampInd --pheno ${fpath_phecov} --keep-if "${phe}==0" --extract ${phe}/${fname}_8.snplist --hwe $thres --write-snplist --out ${phe}/${fname}_9
	else # quantitative
		plink2 --bfile ${fpath}_qcSampInd --keep ${phe}/keepInd_phe.txt --geno 0.05 --write-snplist --out ${phe}/${fname}_7

		plink2 --bfile ${fpath}_qcSampInd --keep ${phe}/keepInd_phe.txt --extract ${phe}/${fname}_7.snplist --maf 0.01 --write-snplist --out ${phe}/${fname}_8

		n=$(wc -l < ${phe}/${fname}_8.snplist)
		thres=$(bc -l <<< 0.05/$n)
		plink2 --bfile ${fpath}_qcSampInd --pheno ${fpath_phecov} --pheno-name $phe --require-pheno --extract ${phe}/${fname}_8.snplist --hwe $thres --write-snplist --out ${phe}/${fname}_9
	fi

	plink --bfile ${fpath}_qcSampInd --keep ${fpath_samp} --extract ${phe}/${fname}_9.snplist --make-bed --out ${phe}/${fname}_qcSampInd_qcVar
	plink --bfile ${fpath}_qcSampRel --keep ${fpath_samp} --extract ${phe}/${fname}_9.snplist --make-bed --out ${phe}/${fname}_qcSampRel_qcVar
done


