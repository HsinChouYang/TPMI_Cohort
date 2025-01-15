
phecov= # phenotype and covariate file

## unrelated samples
# binary trait
for phe in e11 i10_sbp120 i10_sbp130 i10_sbp140; do
	plink2 --bfile ${phe}/${fname}_qcSampInd_qcVar \
	       --pheno ${phecov} --1 --pheno-name $phe \
	       --covar ${phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
	       --threads 16 \
	       --out ${out}/${fname}_qcSampInd_qcVar_${phe}
done

# quantitative trait
for phe in hba1c sbp dbp; do
	# original vlaues
#	plink2 --bfile ${phe}/${fname}_qcSampInd_qcVar \
#	       --pheno ${phecov} --pheno-name $phe \
#	       --covar ${phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 8 \
#	       --out ${out}/${fname}_qcSampInd_qcVar_${phe}

	# inverse normal transformation on original vlaues
#	plink2 --bfile ${phe}/${fname}_qcSampInd_qcVar \
#	       --pheno ${phecov} --pheno-name ${phe}_invnorm \
#	       --covar ${phecov} --covar-name age,sex,bmi,PC1-PC10 --covar-variance-standardize \
#	       --glm no-x-sex hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
#	       --threads 8 \
#	       --out ${out}/${fname}_qcSampInd_qcVar_${phe}_invnorm

	# inverse normal transformation on residuals
	plink2 --bfile ${phe}/${fname}_qcSampInd_qcVar \
	       --pheno ${phecov} --pheno-name ${phe}_res_invnorm \
	       --glm no-x-sex allow-no-covars cols=chrom,pos,ref,alt,a1freq,a1freqcc,gcountcc,nobs,orbeta,se,tz,p,firth,err \
	       --threads 8 \
	       --out ${out}/${fname}_qcSampInd_qcVar_${phe}_res_invnorm
done


## related samples
for phe in e11 i10 hba1c sbp dbp; do
	plink2 --bfile ${phe}/tpm1_tpm2_qcSampRel_qcVar --export bgen-1.2 bits=8 --out ${phe}/tpm1_tpm2_qcSampRel_qcVar
	awk 'NR>2' ${phe}/tpm1_tpm2_qcSampRel_qcVar.sample > tmp && mv tmp ${phe}/tpm1_tpm2_qcSampRel_qcVar.sample
done

for phe in e11 i10 hba1c sbp dbp; do
	if [[ ${phe} =~ ^(e11|i10)$ ]]; then # binary
		ptype=binary
	else # quantitative
		ptype=quantitative
	fi

	# step 1
	Rscript step1_fitNULLGLMM.R \
	 --plinkFile=${phe}/tpm1_tpm2_qcSampRel_qcVar \
	 --phenoFile=${pcname} \
	 --phenoCol=${phe} \
	 --traitType=${ptype} \
	 --covarColList=age,sex,bmi,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	 --sampleIDColinphenoFile=IID \
	 --LOCO=FALSE \
	 --outputPrefix=${out}/tpm1_tpm2_qcSampRel_qcVar_${phe} \
	 --nThreads=32 \
	 --IsOverwriteVarianceRatioFile=TRUE

	# step 2
	Rscript step2_SPAtests.R \
	 --bgenFile=/home/hsinchou@pandora.sinica/jiawei/TPMI/tpmi037/marker_paper/QC/${phe}/tpm1_tpm2_qcSampRel_qcVar.bgen \
	 --sampleFile=/home/hsinchou@pandora.sinica/jiawei/TPMI/tpmi037/marker_paper/QC/${phe}/tpm1_tpm2_qcSampRel_qcVar.sample \
	 --minMAF=0.0001 \
	 --minMAC=1 \
	 --LOCO=FALSE \
	 --GMMATmodelFile=${out}/tpm1_tpm2_qcSampRel_qcVar_${phe}.rda \
	 --varianceRatioFile=${out}/tpm1_tpm2_qcSampRel_qcVar_${phe}.varianceRatio.txt \
	 --SAIGEOutputFile=${out}/tpm1_tpm2_qcSampRel_qcVar_${phe}.SAIGE.results_step2.txt \
	 --IsOutputAFinCaseCtrl=TRUE
done














