# Six-item item-selection benchmarks; source files remain read-only.
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file,winslash="/",mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"),encoding="UTF-8")
basedir <- score_source
outdir <- file.path(runtime_dir,"extensions/benchmarks")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
logcon <- file(file.path(outdir,'benchmark_run.log'),open='wt',encoding='UTF-8')
sink(logcon,split=TRUE);sink(logcon,type='message')
cat('Started:',format(Sys.time()),'\n')
b <- readRDS(file.path(basedir,'candidate_question_bank.rds'))
folds <- readRDS(file.path(basedir,'FI_reduction_locked_scores.rds'))
d <- b$data;bank <- b$bank;cohorts <- c('HRS','ELSA','SHARE','CHARLS','MHAS')
stopifnot(nrow(d)==34954L,length(bank)==16L,!anyDuplicated(d$uid),all(is.finite(d$fi_complete)))
num <- function(x) as.numeric(x)
score6 <- function(x,items){m<-as.matrix(x[,items,drop=FALSE]);n<-rowSums(is.finite(m));v<-rowMeans(m,na.rm=TRUE);v[n<5L]<-NA_real_;v}
auc <- function(y,s){n1<-sum(y==1);n0<-sum(y==0);if(!n1||!n0)return(NA_real_);r<-rank(s,ties.method='average');(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}
metric <- function(y,s,full,outcome){a<-auc(y,s);if(outcome=='concurrent')c(AUROC=a,Spearman=suppressWarnings(cor(s,full,method='spearman')),MAE=mean(abs(s-full)),RMSE=sqrt(mean((s-full)^2))) else c(AUROC=a)}
all_combinations <- t(combn(seq_along(bank),6L))
set.seed(20261001)
selected_ids <- sample.int(nrow(all_combinations),1000L,replace=FALSE)
combos <- all_combinations[selected_ids,,drop=FALSE]
stopifnot(nrow(combos)==1000L,!anyDuplicated(apply(combos,1,paste,collapse='+')))
combo_tbl<-data.frame(combination_id=seq_len(1000L),lexicographic_combination_id=selected_ids,seed=20261001,k=6,items=apply(combos,1,function(i)paste(bank[i],collapse='+')))
for(j in seq_len(6L))combo_tbl[[paste0('item_',j)]]<-bank[combos[,j]]
write.csv(combo_tbl,file.path(outdir,'random_six_item_combinations.csv'),row.names=FALSE)
performance <- list();paired <- list();ranks <- list();flows <- list();choices <- list();training_audit <- list();uid_flow <- list();audit_rows <- list()
perf_i <- pair_i <- flow_i <- 0L
for(coh in cohorts){
 cat('Start cohort',coh,format(Sys.time()),'\n');flush.console()
 f<-folds[[coh]];tr<-d[d$cohort!=coh,,drop=FALSE];te<-d[match(f$uid,d$uid),,drop=FALSE]
 stopifnot(identical(f$uid,te$uid),setequal(f$training_uid,tr$uid),!any(te$uid%in%tr$uid),all(te$cohort==coh))
 rr<-vapply(bank,function(item)suppressWarnings(cor(tr[[item]],tr$fi_complete,method='spearman',use='complete.obs')),numeric(1))
 nn<-vapply(bank,function(item)sum(is.finite(tr[[item]])&is.finite(tr$fi_complete)),integer(1))
 ord<-order(-abs(rr),match(bank,bank));ranked<-bank[ord];chosen<-ranked[1:6]
 ranks[[coh]]<-data.frame(held_out=coh,item=ranked,correlation=rr[ord],absolute_correlation=abs(rr[ord]),rank=seq_along(bank),pairwise_n=nn[ord],training_n=nrow(tr),selected_top6=seq_along(bank)<=6)
 choices[[coh]]<-data.frame(held_out=coh,k=6,SHAP_items=paste(f$universal_rank[1:6],collapse='+'),simple_correlation_items=paste(chosen,collapse='+'),training_n=nrow(tr),heldout_n=nrow(te))
 training_audit[[coh]]<-data.frame(held_out=coh,training_cohort=cohorts[cohorts!=coh],training_n=as.integer(table(factor(tr$cohort,levels=cohorts[cohorts!=coh]))),heldout_UID_overlap=0)
 common<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
 observed<-rowSums(is.finite(as.matrix(te[,bank,drop=FALSE])))
 fair<-common&observed>=15L
 s_shap<-f$universal[,6];s_simple<-score6(te,chosen)
 stopifnot(max(abs(score6(te,f$universal_rank[1:6])[is.finite(s_shap)]-s_shap[is.finite(s_shap)]))<1e-12)
 uid_flow[[coh]]<-data.frame(cohort=coh,uid=te$uid,common_all_k=common,observed_candidate_count=observed,benchmark_base_risk=fair,followup_outcome_observed=is.finite(f$followup_frail),incident_ADL_outcome_observed=is.finite(f$incident_adl))
 outcomes<-list(concurrent=as.numeric(f$complete_FI>=.25),followup_frail=f$followup_frail,incident_adl=f$incident_adl)
 sm<-matrix(NA_real_,sum(fair),1000L)
 te_fair<-te[fair,,drop=FALSE]
 for(j in seq_len(1000L))sm[,j]<-score6(te_fair,bank[combos[j,]])
 stopifnot(all(is.finite(sm)),all(is.finite(s_shap[fair])),all(is.finite(s_simple[fair])))
 for(outcome in names(outcomes)){
  y<-outcomes[[outcome]];risk<-fair&is.finite(y);yf<-y[risk];full<-f$complete_FI[risk];a<-s_shap[risk];c<-s_simple[risk];rand<-sm[is.finite(y[fair]),,drop=FALSE]
  stopifnot(nrow(rand)==length(yf),all(is.finite(rand)),length(unique(yf))==2L)
  flow_i<-flow_i+1L;flows[[flow_i]]<-data.frame(cohort=coh,outcome=outcome,baseline_n=nrow(te),common_all_k_n=sum(common),candidate_15_of_16_n=sum(observed>=15L),fair_base_n=sum(fair),fair_analysis_n=sum(risk),events=sum(yf),common_all_k_outcome_n=sum(common&is.finite(y)),removed_by_candidate_rule=sum(common&is.finite(y)&observed<15L))
  ma<-metric(yf,a,full,outcome);mc<-metric(yf,c,full,outcome)
  perf_i<-perf_i+1L;performance[[perf_i]]<-data.frame(cohort=coh,outcome=outcome,method='Locked_SHAP',combination_id=NA_integer_,n=length(yf),events=sum(yf),AUROC=ma['AUROC'],Spearman=if(outcome=='concurrent')ma['Spearman'] else NA_real_,MAE=if(outcome=='concurrent')ma['MAE'] else NA_real_,RMSE=if(outcome=='concurrent')ma['RMSE'] else NA_real_)
  perf_i<-perf_i+1L;performance[[perf_i]]<-data.frame(cohort=coh,outcome=outcome,method='Simple_correlation',combination_id=NA_integer_,n=length(yf),events=sum(yf),AUROC=mc['AUROC'],Spearman=if(outcome=='concurrent')mc['Spearman'] else NA_real_,MAE=if(outcome=='concurrent')mc['MAE'] else NA_real_,RMSE=if(outcome=='concurrent')mc['RMSE'] else NA_real_)
  rm<-vapply(seq_len(1000L),function(j)metric(yf,rand[,j],full,outcome),numeric(length(ma)))
  if(is.null(dim(rm)))rm<-matrix(rm,nrow=length(ma),dimnames=list(names(ma),NULL))
  perf_i<-perf_i+1L;performance[[perf_i]]<-data.frame(cohort=coh,outcome=outcome,method='Random',combination_id=seq_len(1000L),n=length(yf),events=sum(yf),AUROC=rm['AUROC',],Spearman=if(outcome=='concurrent')rm['Spearman',] else NA_real_,MAE=if(outcome=='concurrent')rm['MAE',] else NA_real_,RMSE=if(outcome=='concurrent')rm['RMSE',] else NA_real_)
  set.seed(20261001+match(coh,cohorts)*100L+match(outcome,names(outcomes)))
  if(max(abs(a-c))<1e-12)boot<-replicate(500L,{ix<-sample.int(length(yf),length(yf),replace=TRUE);stopifnot(all(abs(a[ix]-c[ix])<1e-12));setNames(rep(0,length(ma)),names(ma))}) else boot<-replicate(500L,{ix<-sample.int(length(yf),length(yf),replace=TRUE);metric(yf[ix],a[ix],full[ix],outcome)-metric(yf[ix],c[ix],full[ix],outcome)})
  if(is.null(dim(boot)))boot<-matrix(boot,nrow=length(ma),dimnames=list(names(ma),NULL))
  for(m in names(ma)){
   pair_i<-pair_i+1L;paired[[pair_i]]<-data.frame(cohort=coh,outcome=outcome,metric=m,n=length(yf),SHAP=unname(ma[m]),simple_correlation=unname(mc[m]),difference_SHAP_minus_simple=unname(ma[m]-mc[m]),CI_low=unname(quantile(boot[m,],.025,na.rm=TRUE)),CI_high=unname(quantile(boot[m,],.975,na.rm=TRUE)),bootstrap_replicates=500L,comparison_direction=if(m%in%c('MAE','RMSE'))'Lower is better' else 'Higher is better')
  }
  cat('Done',coh,outcome,'n',length(yf),format(Sys.time()),'\n');flush.console()
 }
 audit_rows[[coh]]<-data.frame(cohort=coh,UID_alignment='PASS',training_heldout_separation='PASS',recomputed_SHAP_score='PASS',all_random_scores_available='PASS',random_count=1000L)
 saveRDS(list(performance=performance,paired=paired,ranks=ranks,flows=flows,choices=choices,training_audit=training_audit),file.path(outdir,'benchmark_progress.rds'))
 rm(sm,te_fair);gc(verbose=FALSE)
}
perf<-do.call(rbind,performance);pair<-do.call(rbind,paired);flow<-do.call(rbind,flows)
write.csv(perf,file.path(outdir,'six_item_benchmark_performance.csv'),row.names=FALSE)
write.csv(pair,file.path(outdir,'SHAP_vs_simple_correlation_paired_bootstrap.csv'),row.names=FALSE)
write.csv(flow,file.path(outdir,'benchmark_sample_flow.csv'),row.names=FALSE)
write.csv(do.call(rbind,ranks),file.path(outdir,'simple_correlation_training_rankings.csv'),row.names=FALSE)
write.csv(do.call(rbind,choices),file.path(outdir,'selected_six_item_sets.csv'),row.names=FALSE)
write.csv(do.call(rbind,training_audit),file.path(outdir,'training_heldout_separation_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,uid_flow),file.path(outdir,'benchmark_participant_risk_sets.csv'),row.names=FALSE)
write.csv(do.call(rbind,audit_rows),file.path(outdir,'benchmark_computational_QA.csv'),row.names=FALSE)
summary_rows<-list();summ_i<-0L;macro_rows<-list();macro_i<-0L
for(outcome in c('concurrent','followup_frail','incident_adl')){
 metrics<-if(outcome=='concurrent')c('AUROC','Spearman','MAE','RMSE') else 'AUROC'
 for(m in metrics){
  for(coh in c(cohorts,'Cohort_mean')){
   z<-perf[perf$outcome==outcome,,drop=FALSE]
   if(coh=='Cohort_mean'){
    ra<-vapply(1:1000,function(j)mean(z[z$method=='Random'&z$combination_id==j,m]),numeric(1));sha<-mean(z[z$method=='Locked_SHAP',m]);sim<-mean(z[z$method=='Simple_correlation',m]);n<-sum(flow$fair_analysis_n[flow$outcome==outcome])
    macro_i<-macro_i+1L;macro_rows[[macro_i]]<-data.frame(outcome=outcome,metric=m,combination_id=1:1000,cohort_mean=ra)
   }else{
    z<-z[z$cohort==coh,,drop=FALSE];ra<-z[z$method=='Random',m];sha<-z[z$method=='Locked_SHAP',m];sim<-z[z$method=='Simple_correlation',m];n<-unique(z$n)
   }
   higher<-m%in%c('AUROC','Spearman');tol<-1e-12
   pct_sha<-100*if(higher)mean(ra<sha-tol)+.5*mean(abs(ra-sha)<=tol) else mean(ra>sha+tol)+.5*mean(abs(ra-sha)<=tol)
   pct_sim<-100*if(higher)mean(ra<sim-tol)+.5*mean(abs(ra-sim)<=tol) else mean(ra>sim+tol)+.5*mean(abs(ra-sim)<=tol)
   summ_i<-summ_i+1L;summary_rows[[summ_i]]<-data.frame(cohort=coh,outcome=outcome,metric=m,n=n,SHAP=sha,simple_correlation=sim,random_median=median(ra),random_2p5=unname(quantile(ra,.025)),random_97p5=unname(quantile(ra,.975)),SHAP_better_than_random_percentile=pct_sha,simple_better_than_random_percentile=pct_sim,random_sets=1000L,interval_type='2.5th-97.5th percentiles across random item sets',better_direction=if(higher)'Higher' else 'Lower')
  }
 }
}
summ<-do.call(rbind,summary_rows)
write.csv(summ,file.path(outdir,'random_combination_distribution_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,macro_rows),file.path(outdir,'random_combination_cohort_mean_distribution.csv'),row.names=FALSE)
stopifnot(nrow(perf)==15030L,nrow(pair)==30L,nrow(summ)==36L,nrow(flow)==15L,all(flow$fair_analysis_n<=flow$common_all_k_outcome_n))
saveRDS(list(performance=perf,paired_bootstrap=pair,random_summary=summ,sample_flow=flow,seed=20261001,random_combinations=combo_tbl,rankings=do.call(rbind,ranks)),file.path(outdir,'six_item_benchmark_results.rds'))
cat('\nKEY RESULTS\n');print(summ[summ$cohort=='Cohort_mean',]);print(pair[pair$metric=='AUROC',]);print(flow)
writeLines(capture.output(sessionInfo()),file.path(outdir,'R_sessionInfo.txt'),useBytes=TRUE)
cat('PASS: 1,000 distinct random six-item combinations; 5 held-out cohorts; paired bootstrap 500; all methods share outcome-specific participants. Completed:',format(Sys.time()),'\n')
sink(type='message');sink();close(logcon)
