
args <- commandArgs(T)
progPath <- commandArgs(F)

# parse arguments & evaluate statements
args <- sapply(args,function(x){
 x <- gsub("-","",x)
 x <- gsub("(.*)=(.*)","\\1='\\2'",x)
})
for(i in 1:length(args))
eval(parse(text=args[[i]]))

# default parameter settings
th.imiss <- .05
th.het <- 1 # 1 = MEAN +/- 3SD, 2 = MEDIAN +/- 3IQR

for(a in gsub("--(.*?)=.*","\\1",names(args))){
 switch(a,process = check <- as.integer(process),
          icr_thre = th.imiss <- as.numeric(imiss_thre),
          het_thre = th.het <- as.integer(het_thre),
          ancestry_info = f.ancestry <- ancestry_info,
          data_name = data_name <- data_name,
          confid = level <- as.numeric(confid),

          fname_a = fn_a <- fname_a,
          fname_p = fn_p <- fname_p)
}

pkgs <- .packages(all.available=TRUE)
if(all(pkgs!="data.table")) install.packages("data.table",repos="https://cloud.r-project.org/")
if(all(pkgs!="car")) install.packages("car",repos="https://cloud.r-project.org/")
if(all(pkgs!="sp")) install.packages("sp",repos="https://cloud.r-project.org/")
if(all(pkgs!="RColorBrewer")) install.packages("RColorBrewer",repos="https://cloud.r-project.org/")

library(data.table)
library(car)
library(sp)
library(RColorBrewer)

wd <- getwd()

f.sex <- dir(wd,pattern="\\.sexcheck",full.name=T)
#f.imiss <- dir(wd,pattern="\\.imiss",full.name=T)
f.smiss <- dir(wd,pattern="\\.smiss",full.name=T)
f.het <- dir(wd,pattern="\\.het",full.name=T)
f.genome <- dir(wd,pattern="\\.genome",full.name=T)
f.eigen1 <- dir(wd,pattern="\\.eigenval",full.name=T)
f.eigen2 <- dir(wd,pattern="\\.eigenvec",full.name=T)

#if(check==1){
# if(length(f.sex)>0) {
#  sex <- fread(f.sex)
#  fwrite(sex[STATUS%in%"PROBLEM",c("FID","IID")],sprintf("%s/rmInd_sex.txt",wd),sep=" ",col.name=F)
# } else {
#  cat("\nERROR! Please provide PLINK output file: *.sexcheck\n\n")
# }
#}

if(check==2){
 if(length(f.smiss)>0 & length(f.het)>0) {
 # if(length(f.imiss)>0 & length(f.het)>0) {
  # imiss <- fread(f.imiss) # plink1.9 
  imiss <- fread(f.smiss) # plink2
  het <- fread(f.het)

  mr <- imiss$F_MISS
  # hr <- het[,(`N(NM)`-`O(HOM)`)/`N(NM)`] # plink1.9 
  hr <- het[,(OBS_CT-`O(HOM)`)/OBS_CT] # plink2

  pdf("cr_hr_check.pdf")
   xlim <- range(pretty(hr))
   ylim <- range(pretty(mr))
   plot(hr,mr,main="Heterozygosity Rate vs. Missing Rate",xlab="Heterozygosity Rate",ylab="Missing Rate",xlim=xlim,ylim=ylim)
    # imiss
     abline(h=th.imiss,lty=2,col="red")
    # het
     a <- ifelse(th.het==1,mean(hr),median(hr))
     w <- ifelse(th.het==1,3*sd(hr),3*diff(quantile(hr,c(.25,.75))))
     abline(v=c(a-w,a+w),lty=2,col="indianred")
     arrows(a-w,rev(axTicks(2))[2],a+w,rev(axTicks(2))[2],length=.1,code=3,col="indianred")
     text(a,rev(axTicks(2))[2],ifelse(th.het==1,"MEAN +/- 3SD","MEDIAN +/- 3IQR"),pos=3,col="indianred")

    idx <- which(mr>th.imiss | (hr<a-w|hr>a+w))
    points(hr[idx],mr[idx],pch=16,col="red")
    text(hr[idx],mr[idx],imiss$IID[idx],pos=3,font=2,col="red",xpd=NA)
  dev.off()
  fwrite(cbind(imiss[idx,c("FID","IID")],mr[idx],hr[idx]),sprintf("%s/rmInd_missing_het.txt",wd),sep=" ",col.name=F)
 } else {
  cat("\nERROR! Please provide PLINK output files: *.imiss and *.het\n\n")
 }
}

#if(check==3){
# if(length(f.genome)>0 & length(f.imiss)>0){
# # data import
#  imiss <- read.table(f.imiss,header=T,stringsAsFactors=F)
#  ibd <- read.table(f.genome,header=T,stringsAsFactors=F)
#
#  ibd.bad <- ibd[ibd$PI_HAT>.1875,]
#  if(nrow(ibd.bad)>0){
#   ibd.bad <- cbind(ibd.bad,t(apply(ibd.bad[,c("IID1","IID2")],1,function(x) imiss$F_MISS[match(x,imiss$IID)])))
#    names(ibd.bad)[(ncol(ibd.bad)-1):ncol(ibd.bad)] <- c("Miss_rate1","Miss_rate2")
#
#   idx <- ifelse(ibd.bad[,"Miss_rate1"]>ibd.bad[,"Miss_rate2"],1,2)
#    idx <- cbind(idx,3-idx)
#
#   write.table(cbind(t(sapply(1:nrow(ibd.bad),function(x) ibd.bad[x,paste0(rep(c("FID","IID","Miss_rate"),2),rep(idx[x,],c(3,3)))])),ibd.bad$PI_HAT),"rmInd_relate.txt",sep=" ",quote=F,row.name=F,col.name=F)
#  } else {
#   cat("\nERROR! Please provide PLINK output files: *.genome and *.imiss\n\n")
#  }
# }
#}

if(check==4){
 if(length(f.eigen1)>0 & length(f.eigen2)>0 & length(f.ancestry)>0){
  # pc scores
   var <- scan(f.eigen1) 
   pc <- fread(f.eigen2,header=F)
   popu <- fread(f.ancestry) # IID, population

## 1000 genomes
   col.popu <- setNames(c(brewer.pal(9,"Oranges")[-(1:2)],brewer.pal(9,"Purples")[c(3,5,7,9)],brewer.pal(9,"Blues")[c(3,4,5,7,9)],brewer.pal(9,"Greens")[c(3,4,5,7,9)],brewer.pal(9,"Greys")[c(3,4,5,7,9)]),
                        c("MSL","ESN","GWD","ACB","ASW","LWK","YRI", "PUR","PEL","CLM","MXL", "GBR","FIN","IBS","TSI","CEU", "KHV","CDX","CHS","CHB","JPT", "PJL","BEB","STU","ITU","GIH"))
   col.popu.mat <- cbind(col.popu[1:7],c(col.popu[8:11],rep(NA,3)),c(col.popu[12:16],rep(NA,2)),c(col.popu[17:21],"red",rep(NA,1)),c(col.popu[22:26],rep(NA,2)))
   col.name.mat <- names(col.popu)[c(1:7,8:11,rep(27,3),12:16,rep(27,2),17:21,rep(27,2),22:26,rep(27,2))]; col.name.mat[27] <- data_name; col.name.mat <- matrix(col.name.mat,7)

   col <- col.popu[toupper(popu[match(pc[,V2],popu[,IID]),population])]
    col[is.na(col)] <- "red"

  pdf("divAncestry_check.pdf")
   plot(pc[,V3],pc[,V4],main="Ancestry check",xlab=sprintf("PC1 (%.2f%%)",var[1]/sum(var)*100),ylab=sprintf("PC2 (%.2f%%)",var[2]/sum(var)*100),pch=16,cex=.5,col=col)
   legend("bottom",col.name.mat,text.col=col.popu.mat,fill=col.popu.mat,border=col.popu.mat,text.font=2,cex=.75,ncol=5,bty="n")

  plot(pc[,V3],pc[,V4],main="Ancestry check",xlab="PC1",ylab="PC2",pch=16,cex=.8,col=col,xlim=range(pretty(pc[col%in%"red",V3]),na.rm=T),ylim=range(pretty(pc[col%in%"red",V4]),na.rm=T))
   xy <- dataEllipse(pc[col%in%"red",V3],pc[col%in%"red",V4],
                     #labels=pc[col%in%"red",V2],id.col="red",id.n=0,
                     levels=level,center.pch=FALSE,add=TRUE,col="red",grid=FALSE,plot.points=FALSE,lty=2)

   idx <- which(point.in.polygon(pc[col%in%"red",V3],pc[col%in%"red",V4],xy[,1],xy[,2])==0)

   dataEllipse(pc[col%in%"red",V3],pc[col%in%"red",V4],
               #labels=pc[col%in%"red",2],id.n=length(idx),id.col="red",
               levels=level,center.pch=FALSE,add=TRUE,col="red",grid=FALSE,plot.points=FALSE,lty=2)
   text(pc[col%in%"red"][idx,3:4],pc[col%in%"red",][idx,V2],col="red",cex=.8,font=2,pos=3,xpd=NA)
  dev.off()

  fwrite(pc[col=="red",1:2][idx,],"rmInd_divAncestry.txt",sep=" ",col.name=F)
 } else {
  cat("\nERROR! Please provide PLINK output file: *.eigenval and *.eigenvec\n")
 }
}








