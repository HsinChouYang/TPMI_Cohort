
#!/bin/bash

## Manually calculated variant weights
# DIAGRAM
python3 PRScsx.py --ref_dir=LD_ref \
--bim_prefix=/data/home/jiawei/prs/dm/tpm1_imput \
--sst_file=DIAMANTE-EAS.sumstat_ok.txt,DIAMANTE-EUR.sumstat_ok.txt,DIAMANTE-SAS.sumstat_ok.txt \
--n_gwas=283423,933970,49492 \
--pop=EAS,EUR,SAS \
--seed=1 --meta=TRUE \
--out_dir=e11 --out_name=e11_prscsx

cat e11/t2d_prscsx_META_*.txt > e11/t2d_prscsx_META.txt 

plink2 --bfile tpmi --score e11/t2d_prscsx_META.txt 2 4 6 cols=fid,nallele,scoresums,scoreavgs --out e11_DIAGRAM_prscsx

## Externally sourced variant weights
# PGS Catalog: PGS002308 (https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002308/ScoringFiles/Harmonized/PGS002308_hmPOS_GRCh38.txt.gz)
plink2 --bfile tpmi --score PGS002308_hmPOS_GRCh38.txt 1 4 6 cols=fid,nallele,scoresums,scoreavgs --out e11_PGS002308_prscsx

## note: if possible, try to use chr_pos to represent for the variant ids when using PLINK to calculate scores