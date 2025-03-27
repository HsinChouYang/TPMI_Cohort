
fpath_geno_i=$1 # independent QC sample dataset, e.g., *_qcSampInd
fpath_geno_r=$2 # related QC sample dataset, e.g., *_qcSampRel
fpath_phecov=$3 # phenotype + covariate file, phecov.txt

## unrelated samples
# binary trait
for phe in e11 i10; do
	plink2 --bfile ${fpath_geno_i} \
	       --pheno ${fpath_phecov} --1 --pheno-name $phe \
	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out assoc_sampInd_${phe}
done

# quantitative trait
for phe in hba1c sbp dbp; do
	# original vlaues
#	plink2 --bfile ${fpath_geno_i} \
#	       --pheno ${fpath_phecov} --pheno-name $phe \
#	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 16 \
#	       --out assoc_sampInd_${phe}_ori

	# inverse normal transformation on original vlaues
#	plink2 --bfile ${fpath_geno_i} \
#	       --pheno ${fpath_phecov} --pheno-name ${phe}_invnorm \
#	       --covar ${fpath_phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 16 \
#	       --out assoc_sampInd_${phe}_inv

	# inverse normal transformation on residuals
	plink2 --bfile ${fpath_geno_i} \
	       --pheno ${fpath_phecov} --pheno-name ${phe}_res_invnorm \
	       --glm no-x-sex allow-no-covars cols=chrom,pos,ref,alt,a1freq,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out assoc_sampInd_${phe}_resInv
done

## related samples
for phe in e11 i10 hba1c sbp dbp; do
	plink2 --bfile ${fpath_geno_r} --export bgen-1.2 bits=8 --out ${fpath_geno_r}
	awk 'NR>2' ${fpath_geno_r}.sample > tmp && mv tmp ${fpath_geno_r}.sample
done

# REGENIE
for phe in e11 i10 hba1c sbp dbp; do
	if [[ ${phe} =~ ^(e11|i10)$ ]]; then # binary
		# step 1
		./regenie \
		--step 1 \
		--bed ${fpath_geno_r} \
		--phenoFile ${fpath_phecov} \
		--phenoCol ${phe} \
		--covarFile ${fpath_phecov} \
		--covarColList age,sex,bmi,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--bsize 100 \
		--bt \
		--out assoc_sampRel_${phe}.regenie.step1

		# step 2
		./regenie \
		--step 2 \
		--bgen ${fpath_geno_r}.bgen \
		--phenoFile ${fpath_phecov} \
		--phenoCol ${phe} \
		--covarFile ${fpath_phecov} \
		--covarColList age,sex,bmi,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--bsize 200 \
		--bt \
		--firth --approx \
		--pThresh 0.01 \
		--pred assoc_sampRel_${phe}.regenie.step1_pred.list \
		--out assoc_sampRel_${phe}.regenie.step2
	else # quantitative
		# step 1
		./regenie \
		--step 1 \
		--bed ${fpath_geno_r} \
		--phenoFile ${fpath_phecov} \
		--phenoCol ${phe} \
		--covarFile ${fpath_phecov} \
		--covarColList age,sex,bmi,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--bsize 100 \
		--qt --apply-rint \
		--out assoc_sampRel_${phe}.regenie.step1

		# step 2
		./regenie \
		--step 2 \
		--bgen ${fpath_geno_r}.bgen \
		--phenoFile ${fpath_phecov} \
		--phenoCol ${phe} \
		--covarFile ${fpath_phecov} \
		--covarColList age,sex,bmi,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--bsize 200 \
		--qt \
		--apply-rint \
		--pThresh 0.01 \
		--pred assoc_sampRel_${phe}.regenie.step1_pred.list \
		--out assoc_sampRel_${phe}.regenie.step2
	fi
done

# SAIGE
for phe in e11 i10 hba1c sbp dbp; do
	if [[ ${phe} =~ ^(e11|i10)$ ]]; then # binary
		ptype=binary
	else # quantitative
		ptype=quantitative
	fi

	# step 1
	Rscript step1_fitNULLGLMM.R \
	 --plinkFile=${fpath_geno_r} \
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
		--bgenFile=${fpath_geno_r}.bgen \
		--sampleFile=${fpath_geno_r}.sample \
		--minMAF=0.0001 \
		--minMAC=1 \
		--LOCO=FALSE \
		--GMMATmodelFile=assoc_sampRel_${phe}.rda \
		--varianceRatioFile=assoc_sampRel_${phe}.varianceRatio.txt \
		--SAIGEOutputFile=assoc_sampRel_${phe}.SAIGE.results_step2.txt \
		--IsOutputAFinCaseCtrl=TRUE
done













