
fpath_geno=$1 # *_qcSampInd
fpath_phecov=$2 # tpmi037_pcJH.fpath_phecov
phe=$3 # e11, colname in $2

## unrelated samples
# binary trait
for phe in e11 i10; do
	plink2 --bfile ${fpath_geno} \
	       --pheno ${fpath_phecov} --1 --pheno-name $phe \
	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out assoc_sampInd_${phe}
done

# quantitative trait
for phe in hba1c sbp dbp; do
	# original vlaues
#	plink2 --bfile ${fpath_geno} \
#	       --pheno ${fpath_phecov} --pheno-name $phe \
#	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 16 \
#	       --out assoc_sampInd_${phe}_ori

	# inverse normal transformation on original vlaues
#	plink2 --bfile ${fpath_geno} \
#	       --pheno ${fpath_phecov} --pheno-name ${phe}_invnorm \
#	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 16 \
#	       --out assoc_sampInd_${phe}_inv

	# inverse normal transformation on residuals
	plink2 --bfile ${fpath_geno} \
	       --pheno ${fpath_phecov} --pheno-name ${phe}_res_invnorm \
	       --glm no-x-sex allow-no-covars cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out assoc_sampInd_${phe}_resInv
done


## related samples
for phe in e11 i10 hba1c sbp dbp; do
	plink2 --bfile ${fpath_geno} --export bgen-1.2 bits=8 --out ${fpath_geno}
	awk 'NR>2' ${fpath_geno}.sample > tmp && mv tmp ${fpath_geno}.sample
done

for phe in e11 i10 hba1c sbp dbp; do
	if [[ ${phe} =~ ^(e11|i10)$ ]]; then # binary
		ptype=binary
	else # quantitative
		ptype=quantitative
	fi

	# step 1
	Rscript step1_fitNULLGLMM.R \
	 --plinkFile=${fpath_geno} \
	 --phenoFile=${fpath_phecov} \
	 --phenoCol=${phe} \
	 --traitType=${ptype} \
	 --covarColList=age,sex,bmi,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	 --sampleIDColinphenoFile=IID \
	 --LOCO=FALSE \
	 --outputPrefix=assoc_sampRel_${phe} \
	 --nThreads=32 \
	 --IsOverwriteVarianceRatioFile=TRUE

	# step 2
	Rscript step2_SPAtests.R \
		--bgenFile=${fpath_geno}.bgen \
		--sampleFile=${fpath_geno}.sample \
		--minMAF=0.0001 \
		--minMAC=1 \
		--LOCO=FALSE \
		--GMMATmodelFile=assoc_sampRel_${phe}.rda \
		--varianceRatioFile=assoc_sampRel_${phe}.varianceRatio.txt \
		--SAIGEOutputFile=assoc_sampRel_${phe}.SAIGE.results_step2.txt \
		--IsOutputAFinCaseCtrl=TRUE
done













