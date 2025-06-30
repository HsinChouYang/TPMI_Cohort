
library(data.table)

af <- fread("tpmi_af.txt")
 af <- af[,c("#CHROM","POS"):=list(paste0("chr",`#CHROM`),as.integer(gsub("chr[0-9]+_([0-9]+)_[ACGT-]+_[ACGT1]+","\\1",ID)))] # ID is in the format of 'chr_pos_ref_alt'
 
af_chinamap <- fread("chinamap_af.txt")
 names(af_chinamap) <- c("#CHROM","POS","REF","ALT","AF.CMAP","AN.CMAP") 

af_wbbc <- fread("wbbc_af.txt")
 names(af_wbbc) <- c("#CHROM","POS","REF","ALT","AF.WBBC","AN.WBBC") 

af_type <- lapply(paste0("chr",1:22),function(chr){
 dt <- merge(af[`#CHROM`==chr,c("POS","ID","REF","ALT","ALT_FREQS")],
             merge(af_chinamap[`#CHROM`==chr],af_wbbc[`#CHROM`==chr],by=c("POS","REF","ALT"),all=T),by=c("POS","REF","ALT"),all=T)
 dt[,c(.SD[,.(POS,REF,ALT)],lapply(.SD,function(x) fifelse(x>.05,"C",fifelse(x>.01,"L",fifelse(x>1e-3,"R","U"))))),.SDcols=c("ALT_FREQS","AF.CMAP","AF.WBBC")]
})
# af_type <- rbindlist(af_type,idcol="chr")

af_type_comp <- rbindlist(lapply(af_type,function(dt){
 dcast(melt(dt,id.var=c("POS","ALT_FREQS"),measure.vars=c("AF.CMAP","AF.WBBC")),ALT_FREQS+value~variable)
}),idcol="chr")

