
#!/bin/bash

fpath_list=$1 # annovar input (refer to https://annovar.openbioinformatics.org/en/latest/user-guide/input/)

# download reference database
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avdblist humandb/ 
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug humandb/ # take time
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp151 humandb/ # take time
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47c humandb/ # take time
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240611 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar icgc28 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 humandb/

./annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
./annotate_variation.pl -buildver hg38 -downdb gwasCatalog humandb/
./annotate_variation.pl -buildver hg38 -downdb gtexGeneModelV8 humandb/

# annotate
./table_annovar.pl $fpath_list humandb/ -buildver hg38 -out ${fpath_list%.*}
 -remove -protocol refGene,cytoBand,EAS.sites.2015_08,exac03,avsnp151,dbnsfp47c,gwasCatalog,clinvar_20240611,icgc28,cosmic70 -operation gx,r,f,f,f,f,r,f,f,f -nastring . -polish -xref example/gene_fullxref.txt
