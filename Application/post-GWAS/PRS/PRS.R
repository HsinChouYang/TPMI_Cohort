
args <- commandArgs(T)
progPath <- commandArgs(F)

# parse arguments & evaluate statements
args <- sapply(args,function(x){
 x <- gsub("-","",x)
 x <- gsub("(.*)=(.*)","\\1='\\2'",x)
})
for(i in 1:length(args))
eval(parse(text=args[[i]]))

# default parameter settings (none)

for(a in gsub("--(.*?)=.*","\\1",names(args))){
 switch(a,fname_prs = fn_prs <- fname_prs,
          fname_phecov = fn_pc <- fnfname_phecovame_p,
          phe = phe <- phe
          covs = covs <- covs)
}
covs <- gsub(",","+",covs)

pkgs <- .packages(all.available=TRUE)
for(pkg in c("data.table","dplyr","pROC"))
 if(all(pkgs!=pkg)) install.packages(pkg,repos="https://cloud.r-project.org/")
 
library(data.table)
library(dplyr) # ntile
library(pROC)

## import calculated scores & phenotype/covariate file
prs <- fread(fn_prs)
 setnames(prs,"SCORE1_SUM","prs") # check the score name
phecov <- fread(fn_pc)
 phecov <- merge(phecov,prs[,c("IID","prs")],by="IID")

## modeling & prediction
dat <- phecov[!is.na(phecov[[phe]])]

# ~ PRS
mod <- glm(formula(sprintf("%s~prs",phe)),data=dat,family="binomial")
roc_mod <- pROC::roc(mod$y,mod$fitted)
# roc_pred <- pROC::roc(dat_new,predict(mod,dat_new)) # make sure 'dat_new' involved the same column names in mod
mod_p <- list(mod=mod,
              roc_mod=roc_mod,auc_mod=pROC::auc(roc_mod),
#              roc_pred=roc_pred,auc_pred=pROC::auc(roc_pred)
             )

# ~ age + gender + BMI + PRS
mod <- glm(formula(sprintf("%s~%s+prs",phe,covs)),data=dat,family="binomial")
roc_mod <- pROC::roc(mod$y,mod$fitted)
roc_pred <- pROC::roc(dat_new,predict(mod,dat_new)) # make sure 'dat_new' involved the same column names in mod
mod_asbp <- list(mod=mod,
                 roc_mod=roc_mod,auc_mod=pROC::auc(roc_mod),
#                 roc_pred=roc_pred,auc_pred=pROC::auc(roc_pred)
                )


## decile model (for training data)
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

dat[,ntile:=ntile(prs,10)]
dat[,ntile_nm:=sprintf("%s%%~%s%%",c(0:4,4,6:9)*10,c(1:4,6,6:10)*10)[ntile_dm]]
dat[,ntile_nm:=factor(ntile_nm)]
dat[,ntile_nm:=relevel(ntile_nm,ref="40%~60%")]

glm_p <- glm(formula(sprintf("%s~ntile_nm",phe)),data=dat,family="binomial")
glm_asbp <- glm(formula(sprintf("%s~age+sex+bmi+ntile_nm",phe)),data=dat,family="binomial")

prs_decile_p <- prs_decile(glm_p,.95)
prs_decile_asbp <- prs_decile(glm_asbp,.95)

save(mod_p,mod_asbp,prs_decile_p,prs_decile_asbp,file=sprintf("%s_%s_mod.RData",phe,gsub("\\..*","",basename(fname_prs))))


