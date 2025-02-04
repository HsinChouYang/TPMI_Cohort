############################################################################################
## This program was written to detect SNPs in specific batches of participants with       ##
## significantly different allele frequencies compared to other batches. Each array       ##
## (TPMv1/TPMv2) has the corresponding batchAAF.tsv and batchSampleSize.tsv for the work. ##
############################################################################################

## step1: identify SNP candidates with different AFs across batches (using Fisher's exact test)
rm(list=ls())
library(future.apply)

dtin <- read.table("batchAAF.tsv",header=T,sep="\t",stringsAsFactors=FALSE)
dtin_msng <- read.table("batchSampleSize.tsv",header=T,sep="\t",stringsAsFactors=FALSE)
colAND <- intersect(colnames(dtin),colnames(dtin_msng))
dtin <- dtin[,colAND]; dtin_msng <- dtin_msng[,colAND]

dtin <- dtin[dtin[,1] %in% c(1:23),]
dtin_bfreq <- dtin[,6:ncol(dtin)]; dtin_bfreq <- as.matrix(dtin_bfreq)
ind_nz <- which(apply(dtin_bfreq,1,function(x){any(x>=0)}))
length(ind_nz)==nrow(dtin_bfreq)
if(length(ind_nz)!=nrow(dtin_bfreq)){dtin <- dtin[ind_nz,] ;dtin_bfreq <- dtin_bfreq[ind_nz,]} 

dtin_sz <- dtin_msng[,6:ncol(dtin_msng)]; dtin_sz <- as.matrix(dtin_sz); dtin_sz <- 2*dtin_sz
identical(colnames(dtin_sz),colnames(dtin_bfreq))

dtin_idtmp <- apply(dtin,1,function(x){paste(x[1:5],collapse="_")})
dtinmsng_idtmp <- apply(dtin_msng,1,function(x){paste(x[1:5],collapse="_")})
dtin_sz_m <- dtin_sz[match(dtin_idtmp,dtinmsng_idtmp),]
any(is.na(dtin_sz_m[,1]))

fun_mktb <- function(x){ 
ind_nn <- which((dtin_bfreq[x,]>=0) & (dtin_sz_m[x,]>0))
tmpout <- matrix(round(c(as.numeric(dtin_sz_m[x,ind_nn]*dtin_bfreq[x,ind_nn]),
                      as.numeric(dtin_sz_m[x,ind_nn]*(1-dtin_bfreq[x,ind_nn])))),nrow=2,byrow=TRUE)
tmpout}

plan(multicore)
options(future.globals.maxSize=Inf)
st1 <- proc.time()
tmplst <- future_lapply(as.list(1:nrow(dtin_sz_m)), fun_mktb)
ed1 <- proc.time()

plan(multicore,workers=16)
options(future.globals.maxSize=Inf)					  
st2 <- proc.time()
outmp <- future_lapply(tmplst,function(x){tryCatch(fisher.test(x,simulate.p.value=TRUE,B=10^5)$p.value,error=function(e) NA)},future.seed=TRUE)
ed2 <- proc.time()
outmp_m <- unlist(outmp)

out_fisher_p <- cbind(dtin[,1:5],fisher_p=outmp_m)

write.table(out_fisher_p,"batch_freq_Fisherp_e5.txt",sep="\t",quote=F,col.name=T,row.name=F)
write.table(dtin,"batch_freq_Freq.txt",sep="\t",quote=F,col.name=T,row.name=F)

## check if the pvalue_NA came from a row of zeros
## TRUE  
tmplst_tmp <- tmplst[which(is.na(outmp_m))]
tmpcheck <- lapply(tmplst_tmp,function(x){any(rowSums(x)==0)})
tmpcheck_m <- unlist(tmpcheck); all(tmpcheck_m)

if(!all(tmpcheck_m)){
 if(all(sapply(tmplst_tmp[which(!tmpcheck_m)],ncol)==1)){ ## check if the pvalue_NA came from a 2x1 table
  oneBatch <- apply(dtin[which(is.na(outmp_m))[which(!tmpcheck_m)],][,-(1:5)],1,function(x) which(x!=-1))
  write.table(cbind(dtin[names(oneBatch),"SNPID"],names(dtin)[-(1:5)][oneBatch]),"Call_by_oneBatch.txt",sep="\t",quote=F,col.name=c("SNPID","Batch"),row.name=F)
 }
}


## step2: for each SNP candidate, identifying batches with AFs away from AFs of (1) TPM array (median AF of batches) or (2) TWB(NGS) and EAS(1000 genomes project) datasets
## (pos twb2: hg38)
rm(list=ls())

dtin_bfreqtmp <- read.table("batch_freq_Freq.txt",header=T,sep="\t")
dtin_bptmp <- read.table("batch_freq_Fisherp_e5.txt",header=T,sep="\t")
identical(paste(dtin_bfreqtmp[,1:5],collapse="_"), paste(dtin_bptmp[,1:5],collapse="_"))

dtin_bfreq <- dtin_bfreqtmp[,6:ncol(dtin_bfreqtmp)]; dtin_bfreq <- as.matrix(dtin_bfreq)
dtin_p <- as.vector(dtin_bptmp[,6])

diff_a <- 0.1

tmp1 <- apply(dtin_bfreq,1,function(x){diff(range(x[x>=0],na.rm=TRUE)) > diff_a})
tmp2 <- dtin_p < 10^-5
tmp3 <- which(tmp1 & tmp2)
tmp3_id <- dtin_bfreqtmp[tmp3,3]

af_wgs <- read.table("twb_af.tsv",header=T,sep="\t",stringsAsFactors=FALSE)
af_wgs_idtmp <- paste(af_wgs[,1],af_wgs[,2],af_wgs[,4],af_wgs[,5],sep="_") 
dtin_idtmp <- paste(paste("chr",dtin_bfreqtmp[,1],sep=""),dtin_bfreqtmp[,2],dtin_bfreqtmp[,4],dtin_bfreqtmp[,5],sep="_") # ann_refalt_m[,13]
af_wgs_m <- af_wgs[match(dtin_idtmp,af_wgs_idtmp),]

af_eas <- read.table("g1k_EAS_unrl.afreq",header=T,sep="\t",stringsAsFactors=FALSE,comment.char="")
af_eas_idtmp <- paste(af_eas[,1],af_eas[,2],af_eas[,4],af_eas[,5],sep="_")
dtin_idtmp2 <- paste(dtin_bfreqtmp[,1],dtin_bfreqtmp[,2],dtin_bfreqtmp[,4],dtin_bfreqtmp[,5],sep="_")
af_eas_m <- af_eas[match(dtin_idtmp2,af_eas_idtmp),] 

af_eas_m[,6] <- as.numeric(af_eas_m[,6]) 
# af_eas_m[rev_tmp,6] <- 1-af_eas_m[rev_tmp,6]

af_pub_tmp <- cbind(as.numeric(af_wgs_m[,6]),as.numeric(af_eas_m[,6]))
af_pub <- apply(af_pub_tmp,1,function(x){ if(!is.na(x[1])){x[1]}else if(!is.na(x[2])){x[2]}else{NA}})

st3 <- proc.time()
zcls <- c()
for(i in tmp3){
   dtin_bfreq_i <- dtin_bfreq[i,]
   dtin_bfreq_i[dtin_bfreq_i < 0]=NA

   srt_x_i <- factor(1:ncol(dtin_bfreq),levels=order(dtin_bfreq_i))

   isambg <- paste(dtin_bfreqtmp[i,4],dtin_bfreqtmp[i,5],sep="_") %in% c("A_T","T_A","C_G","G_C")

   fmed_thr_l <- pmax(0,median(dtin_bfreq_i,na.rm=T)-diff_a)
   fmed_thr_u <- pmin(1,median(dtin_bfreq_i,na.rm=T)+diff_a)
   iout_med <- (dtin_bfreq_i < fmed_thr_l) | (dtin_bfreq_i > fmed_thr_u)
   iin_med <- (dtin_bfreq_i >= fmed_thr_l) & (dtin_bfreq_i <= fmed_thr_u)

   if((!is.na(af_pub[i])) & (!isambg)){
      fpub_thr_l <- pmax(0,af_pub[i]-diff_a)
      fpub_thr_u <- pmin(1,af_pub[i]+diff_a)
      iout_pub <- (dtin_bfreq_i < fpub_thr_l) | (dtin_bfreq_i > fpub_thr_u)
      iin_pub <- (dtin_bfreq_i >= fpub_thr_l) & (dtin_bfreq_i <= fpub_thr_u)

      if(any(iout_med & iout_pub,na.rm=TRUE)){zclstmp <- cbind(dtin_bfreqtmp[i,3],colnames(dtin_bfreq)[which(iout_med & iout_pub)],"3_1");  # paste("b",which(iout_med & iout_pub),sep="")
                                    zcls <- rbind(zcls,zclstmp)}
      if((sum(iin_pub,na.rm=TRUE)>1) & (!any(iin_med & iin_pub,na.rm=TRUE)) & ( (sum(iin_med | iin_pub,na.rm=TRUE)/sum(!is.na(iin_med | iin_pub),na.rm=TRUE)) >= 0.9)){
                                    zclstmp <- cbind(dtin_bfreqtmp[i,3],colnames(dtin_bfreq)[which(iout_pub)],"3_2");  # paste("b",which(iout_pub),sep="")
                                    zcls <- rbind(zcls,zclstmp)}
   }else{
      if(any(iout_med,na.rm=TRUE)){
                  zclstmp <- cbind(dtin_bfreqtmp[i,3],colnames(dtin_bfreq)[which(iout_med)],"3_3");  # paste("b",which(iout_med),sep="")
                        zcls <- rbind(zcls,zclstmp)}
   }


   tiff(paste("snp_ind_",i,".tiff",sep=""),width=800,height=600)
   #par(mfrow=c(4,1))
   plot(dtin_bfreq_i~srt_x_i,ylim=c(0,1),
      main=paste("chr", dtin_bfreqtmp[i,1],"; snp",dtin_bfreqtmp[i,3], "; range:", diff(range(dtin_bfreq_i,na.rm=TRUE)), "; af_pub:", af_pub[i], sep=" "))
   abline(h=c(pmax(0,median(dtin_bfreq_i,na.rm=T)-diff_a),pmin(1,median(dtin_bfreq_i,na.rm=T)+diff_a)),col="red")
   if(!is.na(af_pub[i])){abline(h=c(pmax(0,af_pub[i]-diff_a),pmin(1,af_pub[i]+diff_a)),col="orange")}

   dev.off()
}
ed3 <- proc.time()

write.table(zcls,"batchSNPs.zero",col.name=F,row.name=F,quote=F,sep="\t")
write.table(tmp3_id,"ref_batchSNPs.txt",col.name=F,row.name=F,quote=F,sep="\t")


