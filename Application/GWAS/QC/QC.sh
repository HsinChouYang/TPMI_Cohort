

fpath_geno=$1 # plink bfile name
fpath_sinfo=$2 # tpmi037_486956_Sex_Dup_CR_HomoRate_relatedness.tsv

fname=$(basename $fpath_geno)

fpath_popu=/home/jiawei@pandora.sinica/Desktop/1KGP/all_hg38_extract_tpm12_indep
fpath_clut=/home/jiawei@pandora.sinica/Desktop/1KGP/all_hg38.clust
fpath_info=/home/jiawei@pandora.sinica/Desktop/1KGP/relationships_w_pops.txt
data_name=TPMI

# (1) sex check -- 705
awk '$10=="fail" {print $2,$2}' ${fpath_sinfo} > rmInd_sex.txt
plink2 --bfile ${fpath_geno} --remove rmInd_sex.txt --write-samples --out ${fname}_1

# (2) inconsistent record --372
awk '$11=="fail" {print $2,$2}' ${fpath_sinfo} > rmInd_wrong.txt
plink2 --bfile ${fpath_geno} --keep ${fname}_1.id --remove rmInd_wrong.txt --write-samples --out ${fname}_2

# (3) GCR (0.95) check --19
# (4) homozygosity check (Mean rule) --8123
# plink --bfile ${fpath_geno} --keep ${fname}_2.id --het --missing --out ${fname}_2
plink2 --bfile ${fpath_geno} --keep ${fname}_2.id --het --missing --out ${fname}_2

R CMD BATCH --args --process=2 --imiss_thre=0.05 --het_thre=1 QC.R
plink2 --bfile ${fpath_geno} --keep ${fname}_2.id --remove rmInd_missing_het.txt --write-samples --out ${fname}_3_4

# (5) cryptic relatedness check --111489
awk '$15=="remove" {print $2,$2}' ${fpath_sinfo} > rmInd_rel.txt
plink2 --bfile ${fpath_geno} --keep ${fname}_3_4.id --remove rmInd_rel.txt --write-samples --out ${fname}_5
# R CMD BATCH --args --process=3 /home/hsinchou@pandora.sinica/jiawei/TPMI/tpmi037/marker_paper/QC.R

# (6) divergent ancestry check (# using 1KGP)
plink2 --bfile ${fpath_geno} --set-all-var-ids @:# --rm-dup exclude-all --make-bed --out ${fname}_
plink --bfile ${fname}_ --bmerge $fpath_popu --make-bed --out tmp
plink --bfile $fpath_popu --exclude tmp-merge.missnp --make-bed --out popu_rm3alleles
plink --bfile ${fname}_ --bmerge popu_rm3alleles --exclude tmp-merge.missnp --maf 0.01 --make-bed --out ${fname}__
plink2 --bfile ${fname}__ \
       --read-freq ${fpath_popu}.acount \
       --score ${fpath_popu}.eigenvec.var 2 4 header-read variance-standardize no-mean-imputation \
       --score-col-nums 5-14 \
       --out ${fname}__projAll
R CMD BATCH --args --process=4 --fname_id_ind=${fname}_5.id --fname_id_rel=${fname}_3_4.id --fname_popu=$fpath_popu --ancestry_info=$fpath_info --data_name=$data_name --confid=0.9999 QC.R

plink2 --bfile ${fpath_geno} --keep ${fname}_5.id --remove rmInd_divAncestry.txt --write-samples --out ${fname}_6ind
plink2 --bfile ${fpath_geno} --keep ${fname}_3_4.id --remove rmInd_divAncestry.txt --write-samples --out ${fname}_6rel

plink --bfile ${fpath_geno} --keep ${fname}_6ind.id --make-bed --out ${fpath_geno}_qcSampInd
plink --bfile ${fpath_geno} --keep ${fname}_6rel.id --make-bed --out ${fpath_geno}_qcSampRel
# causion!! when using plink2 to generate bfiles, ${fpath_geno}_qcSampInd had problems on some markers that would have zero MAFs (but ${fpath_geno}_qcSampRel was fine), thus, using plink1.9 to generate bfiles
# chr15 Affx-81719549 etc., over 85000 markers


## (7,8,final) GCR, MAF & HWE
phecov= # phenotype and covariate file

# for phe in e11 i10_sbp120 i10_sbp130 i10_sbp140 hba1c sbp dbp; do
for phe in e11 i10 hba1c sbp dbp; do
	mkdir ${phe}

	fpath_samp=/home/hsinchou@pandora.sinica/jiawei/TPMI/tpmi037/marker_paper/review/lst_sample_${phe}.txt
	# if [[ ${phe} =~ ^(e11|i10_sbp120|i10_sbp130|i10_sbp140)$ ]]; then # binary
	if [[ ${phe} =~ ^(e11|i10)$ ]]; then # binary
		fpath_ctrl=/home/hsinchou@pandora.sinica/jiawei/TPMI/tpmi037/marker_paper/review/lst_health_${phe}.txt

		plink --bfile ${fpath}_qcSampInd --keep ${fpath_samp} --geno 0.05 --pheno ${phecov} --pheno-name $phe --1 --test-missing --out ${phe}/${fname}_6
		n=$(wc -l < ${fpath}_qcSampInd.bim)
		thres=$(bc -l <<< 0.05/$n)
		awk -v thres="$thres" '$5<thres {print $2}' ${phe}/${fname}_6.missing > ${phe}/rmSNP_missingTest.txt
		plink --bfile ${fpath}_qcSampInd --keep ${fpath_samp} --geno 0.05 --exclude ${phe}/rmSNP_missingTest.txt --write-snplist --out ${phe}/${fname}_7
	else # quantitative
		fpath_ctrl=/home/hsinchou@pandora.sinica/jiawei/TPMI/tpmi037/marker_paper/review/lst_sample_${phe}.txt

		plink --bfile ${fpath}_qcSampInd --keep ${fpath_samp} --geno 0.05 --write-snplist --out ${phe}/${fname}_7
	fi

	plink --bfile ${fpath}_qcSampInd --keep ${fpath_samp} --extract ${phe}/${fname}_7.snplist --maf 0.01 --write-snplist --out ${phe}/${fname}_8

	n=$(wc -l < ${phe}/${fname}_8.snplist)
	thres=$(bc -l <<< 0.05/$n)
	plink --bfile ${fpath}_qcSampInd --keep ${fpath_ctrl} --extract ${phe}/${fname}_8.snplist --hwe $thres --write-snplist --out ${phe}/${fname}_9

	plink --bfile ${fpath}_qcSampInd --keep ${fpath_samp} --extract ${phe}/${fname}_9.snplist --make-bed --out ${phe}/${fname}_qcSampInd_qcVar
	plink --bfile ${fpath}_qcSampRel --keep ${fpath_samp} --extract ${phe}/${fname}_9.snplist --make-bed --out ${phe}/${fname}_qcSampRel_qcVar
done

















