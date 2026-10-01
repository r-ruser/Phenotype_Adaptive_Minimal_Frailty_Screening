options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
source(file.path(package_dir,"code/02_functions.R"),encoding="UTF-8")
folds<-readRDS(file.path(score_source,"FI_reduction_locked_scores.rds"))
ks<-1:16;set.seed(20261201)
auc_ci<-function(y,s,B=150){a<-auc_fast(y,s);if(length(unique(y))<2||B<=0)return(c(a,NA_real_,NA_real_));n<-length(y);v<-replicate(B,{ix<-sample.int(n,n,replace=TRUE);auc_fast(y[ix],s[ix])});c(a,unname(quantile(v,.025,na.rm=TRUE)),unname(quantile(v,.975,na.rm=TRUE)))}
one<-function(coh,k,strategy,outcome,y,s,risk,available,full_score,threshold=NA_real_,B=150,rule="80pct"){
 z<-which(risk&is.finite(s));yy<-y[z];ss<-s[z]
 if(length(z)<30||length(unique(yy))<2)return(NULL)
 a<-auc_ci(yy,ss,B);rho<-if(outcome=="concurrent")suppressWarnings(cor(ss,full_score[z],method="spearman")) else NA_real_
 mae<-if(outcome=="concurrent")mean(abs(ss-full_score[z])) else NA_real_;rmse<-if(outcome=="concurrent")sqrt(mean((ss-full_score[z])^2)) else NA_real_
 sens<-if(outcome=="concurrent"&&is.finite(threshold))mean(ss[yy==1]>=threshold) else NA_real_
 spec<-if(outcome=="concurrent"&&is.finite(threshold))mean(ss[yy==0]<threshold) else NA_real_
 data.frame(cohort=coh,k=k,strategy=strategy,outcome=outcome,rule=rule,n=length(z),events=sum(yy),eligible_n=length(y),calculable_n=available,calculable_fraction=available/length(y),AUROC=a[1],AUROC_low=a[2],AUROC_high=a[3],Spearman=rho,MAE=mae,RMSE=rmse,threshold85=threshold,sensitivity=sens,specificity=spec)
}
rows<-list();ccrows<-list();riskrows<-list()
for(coh in names(folds)){
 f<-folds[[coh]];n<-length(f$uid);full<-f$complete_FI
 common<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
 outcomes<-list(concurrent=as.numeric(full>=.25),followup_frail=f$followup_frail,incident_adl=f$incident_adl)
 for(outcome in names(outcomes)){
  y<-outcomes[[outcome]];risk<-common&!is.na(y)
  riskrows[[paste(coh,outcome)]]<-data.frame(cohort=coh,outcome=outcome,baseline_n=n,common_all_k_n=sum(common),analysis_n=sum(risk),events=sum(y[risk]))
  for(strategy in c("Universal","Phenotype")){
   m<-if(strategy=="Universal")f$universal else f$phenotype_score
   mc<-if(strategy=="Universal")f$universal_complete_case else f$phenotype_complete_case
   th<-if(strategy=="Universal")f$threshold_universal else f$threshold_phenotype
   for(k in ks){
    key<-paste(coh,outcome,strategy,k)
    rows[[key]]<-one(coh,k,strategy,outcome,y,m[,k],risk,sum(is.finite(m[,k])),full,th[k])
    ccrisk<-risk&is.finite(mc[,k]);ccrows[[key]]<-one(coh,k,strategy,outcome,y,mc[,k],ccrisk,sum(is.finite(mc[,k])),full,th[k],B=0,rule="complete_items")
   }
  }
  if(outcome!="concurrent"){
   rows[[paste(coh,outcome,"FullFI")]]<-one(coh,0,"Full_FI",outcome,y,full,risk,sum(is.finite(full)),full)
  }
 }
 cat("PERFORMANCE",coh,"\n");flush.console()
}
perf<-do.call(rbind,Filter(Negate(is.null),rows));cc<-do.call(rbind,Filter(Negate(is.null),ccrows))
write.csv(perf,file.path(out,"FI_item_number_performance.csv"),row.names=FALSE)
write.csv(cc,file.path(out,"FI_item_number_complete_case.csv"),row.names=FALSE)
write.csv(do.call(rbind,riskrows),file.path(out,"FI_reduction_risk_sets.csv"),row.names=FALSE)
cat("PASS:",nrow(perf),"common-risk-set item-number performance rows and complete-item sensitivity\n")
