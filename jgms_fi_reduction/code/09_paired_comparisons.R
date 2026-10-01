options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
source(file.path(package_dir,"code/02_functions.R"),encoding="UTF-8")
folds<-readRDS(file.path(score_source,"FI_reduction_locked_scores.rds"));prior<-read.csv(file.path(score_source,"FI5_SHAP_locked_predictions.csv"));opt<-read.csv(file.path(out,"FI_optimal_item_number_summary.csv"))
chosen<-unique(c(5L,8L,opt$min_k_4_of_5[is.finite(opt$min_k_4_of_5)]));chosen<-chosen[chosen>=1&chosen<=16]
set.seed(20261219);rows<-list()
one<-function(coh,k,outcome,comparison,y,a,u,full,risk,B=500){ix<-which(risk&is.finite(a)&is.finite(u));if(length(ix)<100||length(unique(y[ix]))<2)return(NULL)
 yy<-y[ix];aa<-a[ix];uu<-u[ix];delta<-auc_fast(yy,aa)-auc_fast(yy,uu)
 boot<-replicate(B,{z<-sample.int(length(ix),length(ix),replace=TRUE);auc_fast(yy[z],aa[z])-auc_fast(yy[z],uu[z])})
 data.frame(cohort=coh,k=k,outcome=outcome,comparison=comparison,n=length(ix),events=sum(yy),AUROC_A=auc_fast(yy,aa),AUROC_universal=auc_fast(yy,uu),delta_AUROC=delta,delta_low=unname(quantile(boot,.025,na.rm=TRUE)),delta_high=unname(quantile(boot,.975,na.rm=TRUE)),Spearman_A=if(outcome=="concurrent")suppressWarnings(cor(aa,full[ix],method="spearman")) else NA_real_,Spearman_universal=if(outcome=="concurrent")suppressWarnings(cor(uu,full[ix],method="spearman")) else NA_real_,MAE_A=if(outcome=="concurrent")mean(abs(aa-full[ix])) else NA_real_,MAE_universal=if(outcome=="concurrent")mean(abs(uu-full[ix])) else NA_real_)
}
for(coh in names(folds)){
 f<-folds[[coh]];old<-prior[match(f$uid,prior$uid),]
 common<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
 outcomes<-list(concurrent=as.numeric(f$complete_FI>=.25),followup_frail=f$followup_frail,incident_adl=f$incident_adl)
 for(outcome in names(outcomes)){
  y<-outcomes[[outcome]];risk<-common&!is.na(y)
  for(k in chosen){rows[[paste(coh,outcome,k,"Phenotype")]]<-one(coh,k,outcome,"Phenotype_minus_Universal",y,f$phenotype_score[,k],f$universal[,k],f$complete_FI,risk)
   if(k==5)rows[[paste(coh,outcome,k,"Forced")]]<-one(coh,k,outcome,"ForcedDistinct_minus_Universal",y,old$forced_unique_FI5,f$universal[,k],f$complete_FI,risk)
  }
 }
 cat("PAIRED",coh,"\n");flush.console()
}
write.csv(do.call(rbind,Filter(Negate(is.null),rows)),file.path(out,"FI_universal_vs_phenotype.csv"),row.names=FALSE)
cat("PASS: paired person-bootstrap universal versus phenotype comparisons at",paste(chosen,collapse=","),"items\n")
