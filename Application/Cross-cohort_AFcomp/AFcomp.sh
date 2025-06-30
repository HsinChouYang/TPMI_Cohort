
#!/bin/bash

fpath_geno=$1 # plink bfile name
fpath_sinf=$2 # sample information file, e.g., info.txt
fpath_popu_dat=$3 # plink bfile name of 1000 genomes project, e.g., all_hg38_extract_tpm12_indep
fpath_popu_inf=$4 # sample information of 1000 genomes project, e.g., all_hg38.psamp

mkdir -p {ChinaMAP,WBBC}

## download
cd ChinaMAP
bcftools index mbiobank_ChinaMAP.phase1.vcf.gz (needs to login for download)
cd ..

cd WBBC # https://wbbc.westlake.edu.cn/downloads_38.html
for chr in `seq 1 22`; do wget https://wbbc.westlake.edu.cn/data/WBBC.chr${chr}.GRCh38.vcf.gz; done
for chr in `seq 1 22`; do wget https://wbbc.westlake.edu.cn/data/WBBC.chr${chr}.GRCh38.vcf.gz.tbi; done
cd ..

## extract AF information
# TPMI imputation data (PLINK format, by chromosomes)
for chr in `seq 1 22`; do 
       plink2 --bfile TPM1_TPM2_Han_imputed_INFO07_MAF001_rmEffects_chr${chr} --freq --out tpmi.chr${chr}
done
awk 'FNR==1 && NR!=1 {next} {print}' tpmi.chr*.afreq > tpmi_af.txt
rm tpmi.chr*.afreq

# ChinaMAP (vcf format)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AN\n' ChinaMAP/mbiobank_ChinaMAP.phase1.vcf.gz > chinamap_af.txt

# WBBC (vcf format, by chromosomes)
for chr in `seq 1 22`; do 
 bcftools view -h WBBC/WBBC.chr${chr}.GRCh38.vcf.gz > header.txt
 awk '/^##contig/ {
  print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
  print "##INFO=<ID=South_AF,Number=A,Type=Float,Description=\"South population allele frequency\">"
  print "##INFO=<ID=AN,Number=A,Type=Integer,Description=\"Allele frequency\">"
  print "##INFO=<ID=South_AN,Number=A,Type=Integer,Description=\"South population allele frequency\">"
 }
 { print }' header.txt > tmp && mv tmp header.txt # add tags before the ##contig line
 bcftools reheader -h header.txt -o tmp.vcf.gz WBBC/WBBC.chr${chr}.GRCh38.vcf.gz

 bcftools query --allow-undef-tags -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AN\n' tmp.vcf.gz > wbbc.chr${chr}.txt
done
cat $(ls wbbc.chr*.txt | sort -V) > wbbc_af.txt
rm tmp.vcf.gz wbbc.chr*.txt

## run R script to compare AFs