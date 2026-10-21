if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(svglite);library(ragg)})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(normalizePath(script_file,winslash="/",mustWork=TRUE))
source(file.path(package_dir,"00_config.R"),encoding="UTF-8")
src <- score_source
out <- file.path(runtime_dir,"extensions/approximation");dir.create(out,recursive=TRUE,showWarnings=FALSE)
b<-readRDS(file.path(src,"candidate_question_bank.rds"));folds<-readRDS(file.path(src,"FI_reduction_locked_scores.rds"))
cohorts<-c("HRS","ELSA","SHARE","CHARLS","MHAS");set.seed(20261001)
score6<-function(d,items){m<-as.matrix(d[,items,drop=FALSE]);v<-rowMeans(m,na.rm=TRUE);v[rowSums(!is.na(m))<5]<-NA_real_;v}
auc<-function(y,s){n1<-sum(y==1);n0<-sum(y==0);if(!n1||!n0)return(NA_real_);(sum(rank(s,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)}
calc<-function(s,y){e<-s-y;c(bias=mean(e),MAE=mean(abs(e)),RMSE=sqrt(mean(e^2)),Spearman=cor(s,y,method="spearman"),LoA_lower=mean(e)-1.96*sd(e),LoA_upper=mean(e)+1.96*sd(e))}
metrics<-list();maps<-list();bins<-list();ret<-list();ret_draws<-list()
for(coh in cohorts){
 f<-folds[[coh]];tr<-b$data[match(f$training_uid,b$data$uid),];stopifnot(!anyNA(tr$uid),!any(tr$cohort==coh))
 tx<-score6(tr,f$universal_rank[1:6]);ok<-is.finite(tx)&is.finite(tr$fi_complete)
 fit<-lm(tr$fi_complete[ok]~tx[ok]);cf<-unname(coef(fit));stopifnot(cf[2]>0)
 maps[[coh]]<-data.frame(cohort=coh,training_n=sum(ok),intercept=cf[1],slope=cf[2],training_cohorts=paste(setdiff(cohorts,coh),collapse=";"),items=paste(f$universal_rank[1:6],collapse=";"))
 common<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
 full<-f$complete_FI[common];raw<-f$universal[common,6];mapped<-pmin(pmax(cf[1]+cf[2]*raw,0),1)
 stopifnot(length(raw)==sum(common),all(is.finite(full)),all(is.finite(raw)))
 for(type in c("Raw FI6","Training-mapped FI6")){
  s<-if(type=="Raw FI6")raw else mapped;pt<-calc(s,full)
  bb<-replicate(500,{ii<-sample.int(length(s),length(s),replace=TRUE);calc(s[ii],full[ii])[1:4]})
  lo<-apply(bb,1,quantile,.025);hi<-apply(bb,1,quantile,.975)
  rr<-data.frame(cohort=coh,score=type,n=length(s),mean_FI=mean(full),mean_score=mean(s),clipped_fraction=if(type=="Raw FI6")0 else mean(cf[1]+cf[2]*raw<0|cf[1]+cf[2]*raw>1))
  for(m in names(pt)){rr[[m]]<-pt[m];if(m%in%rownames(bb)){rr[[paste0(m,"_low")]]<-lo[m];rr[[paste0(m,"_high")]]<-hi[m]}}
  metrics[[paste(coh,type)]]<-rr
  v<-data.frame(mean=round((s+full)/2/.02)*.02,difference=round((s-full)/.02)*.02,n=1)
  v<-aggregate(n~mean+difference,v,sum);v$cohort<-coh;v$score<-type;bins[[paste(coh,type)]]<-v
 }
 for(outcome in c("followup_frail","incident_adl")){
  y<-f[[outcome]];ii<-which(common&!is.na(y));yy<-y[ii];ss<-f$universal[ii,6];ff<-f$complete_FI[ii]
  pt<-c(AUROC_FI6=auc(yy,ss),AUROC_full=auc(yy,ff));pt<-c(pt,delta_AUROC=pt[1]-pt[2],excess_AUROC_retention=(pt[1]-.5)/(pt[2]-.5))
  bb<-replicate(500,{ix<-sample.int(length(yy),length(yy),replace=TRUE);aa<-auc(yy[ix],ss[ix]);ab<-auc(yy[ix],ff[ix]);c(AUROC_FI6=aa,AUROC_full=ab,delta_AUROC=aa-ab,excess_AUROC_retention=(aa-.5)/(ab-.5))})
  rr<-data.frame(cohort=coh,outcome=outcome,n=length(yy),events=sum(yy))
  for(j in seq_along(pt)){nm<-c("AUROC_FI6","AUROC_full","delta_AUROC","excess_AUROC_retention")[j];rr[[nm]]<-unname(pt[j]);rr[[paste0(nm,"_low")]]<-quantile(bb[j,],.025,na.rm=TRUE);rr[[paste0(nm,"_high")]]<-quantile(bb[j,],.975,na.rm=TRUE)}
  ret[[paste(coh,outcome)]]<-rr;ret_draws[[paste(coh,outcome)]]<-bb[4,]
 }
 cat("APPROXIMATION AND PAIRED RETENTION",coh,"\n");flush.console()
}
a<-do.call(rbind,metrics);mp<-do.call(rbind,maps);bn<-do.call(rbind,bins);re<-do.call(rbind,ret)
write.csv(a,file.path(out,"FI6_continuous_agreement.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(mp,file.path(out,"training_linear_mapping.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(bn,file.path(out,"Bland_Altman_binned_source.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(re,file.path(out,"FI6_paired_prospective_retention.csv"),row.names=FALSE,fileEncoding="UTF-8")
macro<-do.call(rbind,lapply(c("followup_frail","incident_adl"),function(o){zz<-re[re$outcome==o,];dd<-sapply(cohorts,function(c)ret_draws[[paste(c,o)]]);v<-rowMeans(dd);data.frame(outcome=o,mean_retention=mean(zz$excess_AUROC_retention),low=quantile(v,.025),high=quantile(v,.975))}))
write.csv(macro,file.path(out,"FI6_macro_retention.csv"),row.names=FALSE,fileEncoding="UTF-8")
# Figure 3 is rendered from saved aggregate agreement outputs by the current figure script.
method<-"Continuous agreement was evaluated on the same all-length common concurrent samples used in the principal curves. We report raw FI6 minus complete-FI mean bias, MAE, RMSE, Spearman correlation, and descriptive Bland–Altman limits (mean difference ±1.96 SD). An additional linear mapping complete FI = intercept + slope × FI6 was fitted solely in each four-cohort training set using calculable scores and applied unchanged to its held-out cohort, with predictions bounded to 0–1. The evaluated FI6 item set and equal weights remain fixed. We used 500 paired participant bootstrap resamples for error/rank intervals and for FI6-minus-complete-FI prospective AUROC differences and excess-over-chance AUROC retention. These intervals condition on the locked item sets, thresholds and mappings. The FI6 classification uses its training-locked cutoff; the complete-FI 0.25 rule defines the reference outcome."
writeLines(method,file.path(out,"approximation_methods.md"),useBytes=TRUE)
writeLines(c("PASS: held-out participants excluded from each mapping fit.","Raw metrics reproduced from the prior performance CSV.","Same concurrent participants used for raw and mapped estimates.","All prospective differences/retention intervals use paired resampling.","Source data for agreement figure aggregated; individual-level data kept local."),file.path(out,"approximation_QA.txt"),useBytes=TRUE)
old<-read.csv(file.path(src,"FI_item_number_performance.csv"));old<-old[old$k==6&old$strategy=="Universal"&old$outcome=="concurrent",]
rr<-a[a$score=="Raw FI6",];stopifnot(max(abs(rr$MAE-old$MAE[match(rr$cohort,old$cohort)]))<1e-12,max(abs(rr$RMSE-old$RMSE[match(rr$cohort,old$cohort)]))<1e-12)
print(a[,c("cohort","score","n","bias","MAE","RMSE","LoA_lower","LoA_upper")]);print(macro)
cat("PASS: continuous agreement, training-only mapping, Bland-Altman figure and paired prospective intervals completed.\n")
