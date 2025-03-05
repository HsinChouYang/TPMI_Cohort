
#!/bin/bash

fpath_geno=$1 # related QC sample dataset, e.g., *_qcSampRel
fpath_phecov=$2 # phenotype + covariate file, phecov.txt
phe=$3 # colname in $2, e.g., e11
ttype=$4 # b (binary) or q (quantitative)

# create bgen files for the step 2 in REGENIE/SAIGE
plink2 --bfile ${fpath_geno} --export bgen-1.2 bits=8 --out ${fpath_geno}
awk 'NR>2' ${fpath_geno}.sample > tmp && mv tmp ${fpath_geno}.sample # for SAIGE

# REGENIE
if [[ ${ttype} == "b" ]]; then # binary
	# step 1
	./regenie \
	--step 1 \
	--bed ${fpath_geno} \
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
	--bgen ${fpath_geno}.bgen \
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
	--bed ${fpath_geno} \
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
	--bgen ${fpath_geno}.bgen \
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

# SAIGE
if [[ ${ttype} == "b" ]]; then # binary
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














