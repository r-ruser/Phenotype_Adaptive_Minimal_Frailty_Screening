if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(metafor)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

out<-score_source
p<-read.csv(file.path(out,"FI_item_number_performance.csv"));p<-p[p$strategy=="Universal"&p$k>=1,]
coh<-c("HRS","ELSA","SHARE","CHARLS","MHAS");meta<-list()
for(outcome in c("concurrent","followup_frail","incident_adl"))for(k in 1:16){
 z<-p[p$outcome==outcome&p$k==k,];z<-z[match(coh,z$cohort),]
 se<-(z$AUROC_high-z$AUROC_low)/(2*1.96);se[!is.finite(se)|se<=0]<-.005
 yi<-qlogis(pmin(pmax(z$AUROC,.501),.999));sei<-se/(z$AUROC*(1-z$AUROC))
 fit<-tryCatch(metafor::rma.uni(yi=yi,sei=sei,method="REML"),error=function(e)NULL)
 if(!is.null(fit))meta[[paste(outcome,k,"AUROC")]]<-data.frame(outcome=outcome,k=k,metric="AUROC",pooled=plogis(as.numeric(fit$b)),low=plogis(fit$ci.lb),high=plogis(fit$ci.ub),tau2=fit$tau2,I2=fit$I2,n_cohorts=nrow(z))
 if(outcome=="concurrent"){
  zrho<-atanh(pmin(pmax(z$Spearman,-.999),.999));seso<-1/sqrt(pmax(z$n-3,1))
  rf<-tryCatch(metafor::rma.uni(yi=zrho,sei=seso,method="REML"),error=function(e)NULL)
  if(!is.null(rf))meta[[paste(outcome,k,"Spearman")]]<-data.frame(outcome=outcome,k=k,metric="Spearman",pooled=tanh(as.numeric(rf$b)),low=tanh(rf$ci.lb),high=tanh(rf$ci.ub),tau2=rf$tau2,I2=rf$I2,n_cohorts=nrow(z))
 }
}
m<-do.call(rbind,meta);write.csv(m,file.path(out,"FI_reduction_meta_summary.csv"),row.names=FALSE)

file.copy(file.path(out,"FI_reduction_meta_summary.csv"),file.path(publication_dir,"figures/Figure_2_pooled_source.csv"),overwrite=TRUE)
