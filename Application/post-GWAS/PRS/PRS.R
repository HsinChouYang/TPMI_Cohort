
pkgs <- .packages(all.available=TRUE)
for(pkg in c("data.table","dplyr","ggpubr","pROC","Metrics","cutpointr"))
 if(all(pkgs!=pkg)) install.packages(pkg)
 
library(data.table)
library(dplyr) # ntile
library(ggpubr)
library(pROC)
library(Metrics)
library(cutpointr)


## import calculated scores & phenotype/covariate file
path <- # path of score files (calculated by PLINK)
prs <- sapply(dir(path,pattern="score$",full.name=T),fread,simplify=F)
 names(prs) <- gsub("\\.score","",basename(names(prs)))
 prs <- dcast(rbindlist(prs,idcol="type"),sample_id~type)

path <- # path of phenotype & covariate file
phecov <- fread(path)

phecov <- merge(phecov[,c("GID","e11","i10","age","sex","bmi")],prs,by.x="GID",by.y="sample_id")
setnames(phecov,"GID","IID")


## modeling
mod <- sapply(names(prs)[-1],function(x){
 sapply(c("e11","i10"),function(d){
  dat <- copy(phecov)[,prs:=phecov[[x]]][!is.na(phecov[[d]])]

  # ~ age + gender + BMI + PRS
  mod <- glm(formula(sprintf("%s~age+sex+bmi+prs",d)),data=dat[IID%in%unrel],family="binomial")
  roc_mod <- pROC::roc(mod$y,mod$fitted)
  roc_pred <- pROC::roc(dat[!IID%in%unrel][[d]],predict(mod,dat[!IID%in%unrel,c("age","sex","bmi","prs")]))
  mod_asbp <- list(mod=mod,roc_mod=roc_mod,roc_pred=roc_pred,n_mod=dat[,sum(IID%in%unrel)],n_pred=dat[,sum(!IID%in%unrel)])

  # ~ PRS
  mod <- glm(formula(sprintf("%s~prs",d)),data=dat[IID%in%unrel],family="binomial")
  roc_mod <- pROC::roc(mod$y,mod$fitted)
  roc_pred <- pROC::roc(dat[!IID%in%unrel][[d]],predict(mod,dat[!IID%in%unrel,c("age","sex","bmi","prs")]))
  mod_p <- list(mod=mod,roc_mod=roc_mod,roc_pred=roc_pred,n_mod=dat[,sum(IID%in%unrel)],n_pred=dat[,sum(!IID%in%unrel)])

  list(mod_asbp=mod_asbp,mod_p=mod_p)
 },simplify=F)
},simplify=F)


## plot
prs_decile <- function(mod,level){
 level_ <- (1-level)/2

 decile <- setDT(data.frame(coef(summary(mod))),keep.rownames=T)
 decile[,c("OR","OR_se_delta"):=list(exp(Estimate),abs(exp(Estimate))*`Std..Error`)]
 decile[,c("OR_lb_mle","OR_ub_mle","OR_lb_delta","OR_ub_delta"):=list(exp(Estimate-qnorm(level_,lower=F)*`Std..Error`),
                                                                      exp(Estimate+qnorm(level_,lower=F)*`Std..Error`),
                                                                      OR-qnorm(level_,lower=F)*OR_se_delta,
                                                                      OR+qnorm(level_,lower=F)*OR_se_delta)]

 decile
}

prs_decile_plot <- function(decile,col="black",adj=0){
 plot((1:nrow(decile))+adj,c(1,decile[-1,OR]),axes=F,pch=16,xlab="PRS decile",ylab="OR",xlim=c(0,nrow(decile)+1),ylim=c(0,10),col=col)
# plot(1:nrow(decile)+adj,decile[,OR],axes=F,pch=c(1,rep(16,nrow(decile)-1)),xlab="PRS decile",ylab="OR",xlim=c(0,nrow(decile)+1),ylim=c(0,10),col=col)
 arrows(x0=(2:nrow(decile))+adj,y0=decile[-1,OR_lb_mle],y1=decile[-1,OR_ub_mle],code=3,angle=90,length=.1,col=col)
 text((1:nrow(decile))+adj,c(1,decile[-1,OR]),round(c(1,decile[-1,OR]),2),col=col,cex=.8,font=2,pos=c(3,rep(4,nrow(decile)-1)))
 box()
 axis(1,1:10)
 axis(2,las=1)
 abline(h=1,lty=2)
# text(1,decile[1,OR],"Baseline\nOdds",pos=3)
}

for(s in c("dm_PGS002308","dm_tpm1_prscsx_meta")){
 dat <- copy(phecov)[,prs:=phecov[[s]]][!is.na(e11)][IID%in%unrel]
  dat[,ntile_dm:=ntile(prs,10)]
  dat[,ntile_dm_nm:=sprintf("%s%%~%s%%",c(0:4,4,6:9)*10,c(1:4,6,6:10)*10)[ntile_dm]]
  dat[,ntile_dm_nm:=factor(ntile_dm_nm)]
  dat[,ntile_dm_nm:=relevel(ntile_dm_nm,ref="40%~60%")]
 
  glm_dm_p <- dat[,glm(e11~ntile_dm_nm,family="binomial")]
  glm_dm_asbp <- dat[,glm(e11~age+sex+bmi+ntile_dm_nm,family="binomial")]
 
  prs_decile_dm_p <- prs_decile(glm_dm_p,.95)
  prs_decile_dm_asbp <- prs_decile(glm_dm_asbp,.95)
 
 pdf(sprintf("dm_PRS_%s.pdf",s),16,8)
  layout(matrix(1:2,1,))

  # ROC curve
  plot(mod[[s]]$e11$mod_p$roc_pred,main="",col="#E41A1C",las=1,mar=c(6.1,4.1,4.1,2.1))
  lines(mod[[s]]$e11$mod_asbp$roc_pred,col="#377EB8")
  labels <- c(as.expression(bquote(bold(AUC[PRS]:~.(sprintf("%.4f",pROC::auc(mod[[s]]$e11$mod_p$roc_pred)))))),
              as.expression(bquote(bold(AUC[PRS+age+sex+BMI]:~.(sprintf("%.4f",pROC::auc(mod[[s]]$e11$mod_asbp$roc_pred)))))))
  legend("bottomright",labels,text.col=c("#E41A1C","#377EB8"),bty="n",cex=1.25)
  text(par('usr')[1],par('usr')[4],"(A)",font=2,cex=1.5,xpd=NA,pos=3,offset=1.25)
 
  # decile plot 
  ylim <- range(pretty(c(prs_decile_dm_p[,OR],prs_decile_dm_asbp[,OR])))
  par(mar=c(6.1,4.1,4.1,2.1))
  adj <- 0
  decile <- prs_decile_dm_p
  col <- "#E31A1C"
  plot(c(1:4,5,6:9)+adj,c(decile[2:5,OR],1,decile[6:9,OR]),type="b",axes=F,pch=16,xlab="",ylab="OR",xlim=c(0,9+1),ylim=ylim,col=col)
  arrows(x0=c(1:4,6:9)+adj,y0=decile[-1,OR_lb_mle],y1=decile[-1,OR_ub_mle],code=3,angle=90,length=.1,col=col)
  text(c(1:4,6:9)+adj,decile[-1,OR],round(decile[-1,OR],2),col=col,cex=.8,font=2,pos=c(rep(2,4),rep(4,4)))
  box(); axis(1,1:9,sprintf("%s%%-%s%%",c(0:4,6:9)*10,c(1:4,6,7:10)*10),las=2); axis(2,las=1); abline(h=1,lty=2)
 
  par(new=T)
  adj <- 0
  col <- "#377EB8"
 # decile <- prs_decile_dm_asp
 # plot(c(1:4,5,6:9)+adj,c(decile[4:7,OR],1,decile[8:11,OR]),type="b",axes=F,pch=16,xlab="",ylab="OR",xlim=c(0,9+1),ylim=ylim,col=col)
 # arrows(x0=c(1:4,6:9)+adj,y0=decile[-(1:3),OR_lb_mle],y1=decile[-(1:3),OR_ub_mle],code=3,angle=90,length=.1,col=col)
 # text(c(1:4,6:9)+adj,decile[-(1:3),OR],round(decile[-(1:3),OR],2),col=col,cex=.8,font=2,pos=c(rep(4,4),rep(2,4)))
  decile <- prs_decile_dm_asbp
  plot(c(1:4,5,6:9)+adj,c(decile[5:8,OR],1,decile[9:12,OR]),type="b",axes=F,pch=16,xlab="",ylab="OR",xlim=c(0,9+1),ylim=ylim,col=col)
  arrows(x0=c(1:4,6:9)+adj,y0=decile[-(1:4),OR_lb_mle],y1=decile[-(1:4),OR_ub_mle],code=3,angle=90,length=.1,col=col)
  text(c(1:4,6:9)+adj,decile[-(1:4),OR],round(decile[-(1:4),OR],2),col=col,cex=.8,font=2,pos=c(rep(4,4),rep(2,4)))
  box(); axis(1,1:9,sprintf("%s%%-%s%%",c(0:4,6:9)*10,c(1:4,6,7:10)*10),las=2); axis(2,las=1); abline(h=1,lty=2)
 
  legend("topleft",c("PRS","PRS + age + sex + BMI"),col=c("#E31A1C","#377EB8"),text.col=c("#E31A1C","#377EB8"),text.font=2,pch=16,cex=1.25,bty="n")
  text(par('usr')[1],par('usr')[4],"(B)",font=2,cex=1.5,xpd=NA,pos=3,offset=1.25)
 dev.off()
}




