
#!/bin/bash

fpath_geno=$1 # *_qcSampInd
fpath_phecov=$2 # tpmi037_pcJH.fpath_phecov
phe=$3 # e11, colname in $2
ttype=$4 # b (binary) or q (quantitative)

if [[ ${ttype} == "b" ]]; then # binary
	plink2 --bfile ${fpath_geno} \
	       --pheno ${fpath_phecov} --1 --pheno-name $phe \
	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out assoc_sampInd_${phe}
else # quantitative
	# original vlaues
#	plink2 --bfile ${fpath_geno} \
#	       --pheno ${fpath_phecov} --pheno-name ${phe} \
#	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 16 \
#	       --out assoc_sampInd_${phe}_ori

	# inverse normal transformation on original vlaues
#	plink2 --bfile ${fpath_geno} \
#	       --pheno ${fpath_phecov} --pheno-name ${phe} \
#	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 16 \
#	       --out assoc_sampInd_${phe}_inv

	# inverse normal transformation on residuals
	plink2 --bfile ${fpath_geno} \
	       --pheno ${fpath_phecov} --pheno-name ${phe} \
	       --glm no-x-sex allow-no-covars cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out assoc_sampInd_${phe}_resInv
fi







