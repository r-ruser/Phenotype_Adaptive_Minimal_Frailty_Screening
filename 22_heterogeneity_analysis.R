if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE)
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
suppressPackageStartupMessages({library(metafor);library(grid);library(ragg);library(svglite)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

base <- score_source
out <- file.path(runtime_dir,'extensions/heterogeneity')
src <- score_source
cohorts <- c('HRS','ELSA','SHARE','CHARLS','MHAS')
cols <- c(HRS='#4E79A7',ELSA='#59A14F',SHARE='#E15759',CHARLS='#B07AA1',MHAS='#F28E2B')
folds <- readRDS(file.path(base,'FI_reduction_locked_scores.rds'))
map <- read.csv(file.path(runtime_dir,'extensions/approximation/training_linear_mapping.csv'))
orig <- read.csv(file.path(runtime_dir,'extensions/approximation/FI6_continuous_agreement.csv'))
pros <- read.csv(file.path(runtime_dir,'extensions/approximation/FI6_paired_prospective_retention.csv'))
auc <- function(y,s){n1<-sum(y==1);n0<-sum(y==0);if(!n1||!n0)return(NA_real_);(sum(rank(s)[y==1])-n1*(n1+1)/2)/(n1*n0)}
B <- 1000L;set.seed(20261019);rows<-list();draws<-list()
add <- function(coh,id,n,point,bb){
 rows[[paste(coh,id)]] <<- data.frame(cohort=coh,metric=id,n=n,estimate=point,se=sd(bb),lower=unname(quantile(bb,.025)),upper=unname(quantile(bb,.975)),B=B)
 draws[[paste(coh,id)]] <<- bb
}
if (!file.exists(file.path(out,'cohort_effects.csv'))) {
for(coh in cohorts){
 f<-folds[[coh]]; common<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
 x<-f$universal[common,6];full<-f$complete_FI[common];mp<-map[map$cohort==coh,];mapped<-pmin(pmax(mp$intercept+mp$slope*x,0),1)
 y<-as.integer(full>=.25);n<-length(x)
 pt<-c(mean(x-full),mean(mapped-full),auc(y,x))
 bb<-replicate(B,{ii<-sample.int(n,n,TRUE);c(mean(x[ii]-full[ii]),mean(mapped[ii]-full[ii]),auc(y[ii],x[ii]))})
 for(j in 1:3)add(coh,c('raw_bias','mapped_bias','concurrent_auc')[j],n,pt[j],bb[j,])
 stopifnot(abs(pt[1]-orig$bias[orig$cohort==coh&orig$score=='Raw FI6'])<1e-10,abs(pt[2]-orig$bias[orig$cohort==coh&orig$score=='Training-mapped FI6'])<1e-10)
 for(o in c('followup_frail','incident_adl')){
  ii<-which(common&!is.na(f[[o]]));yy<-f[[o]][ii];ss<-f$universal[ii,6];ff<-f$complete_FI[ii];nn<-length(ii)
  p<-c(auc(yy,ss),auc(yy,ss)-auc(yy,ff))
  z<-replicate(B,{ix<-sample.int(nn,nn,TRUE);aa<-auc(yy[ix],ss[ix]);c(aa,aa-auc(yy[ix],ff[ix]))})
  prefix<-if(o=='followup_frail')'future_fi' else 'adl'
  add(coh,paste0(prefix,'_auc'),nn,p[1],z[1,]);add(coh,paste0(prefix,'_delta'),nn,p[2],z[2,])
  ref<-pros[pros$cohort==coh&pros$outcome==o,];stopifnot(abs(p[1]-ref$AUROC_FI6)<1e-10,abs(p[2]-ref$delta_AUROC)<1e-10)
 }
 cat('BOOTSTRAP COMPLETE',coh,'\n');flush.console()
}
d<-do.call(rbind,rows);rownames(d)<-NULL
write.csv(d,file.path(out,'cohort_effects.csv'),row.names=FALSE)
if(length(draws))saveRDS(draws,file.path(out,'bootstrap_draws_aggregate.rds'))
} else { d<-read.csv(file.path(out,'cohort_effects.csv')); stopifnot(all(d$B==B), nrow(d)==35) }
fits<-list();sums<-list();loo<-list();sens<-list()
for(id in unique(d$metric)){
 z<-d[d$metric==id,];m<-rma.uni(yi=z$estimate,sei=z$se,method='REML',test='knha');ci<-confint(m);pr<-predict(m)
 fits[[id]]<-m
 sums[[id]]<-data.frame(metric=id,k=m$k,estimate=as.numeric(m$b),se=m$se,lower=m$ci.lb,upper=m$ci.ub,p_mean=if(grepl('_auc$',id))2*pt(-abs((as.numeric(m$b)-.5)/m$se),df=m$k-1) else m$pval,mean_null=if(grepl('_auc$',id)).5 else 0,I2=m$I2,I2_lower=ci$random['I^2(%)','ci.lb'],I2_upper=ci$random['I^2(%)','ci.ub'],tau2=m$tau2,Q=m$QE,Q_df=m$k-1,p_Q=m$QEp,PI_lower=pr$pi.lb,PI_upper=pr$pi.ub)
 d$weight[d$metric==id]<-as.numeric(weights(m))
 for(j in seq_len(nrow(z))){mm<-rma.uni(yi=z$estimate[-j],sei=z$se[-j],method='REML',test='knha');loo[[paste(id,j)]]<-data.frame(metric=id,omitted=z$cohort[j],estimate=as.numeric(mm$b),lower=mm$ci.lb,upper=mm$ci.ub,I2=mm$I2,tau2=mm$tau2)}
 # Sensitivity: modified Knapp-Hartung CI prevents a scale factor below one.
 scale<-max(1,m$s2w);ss<-sqrt(scale/sum(1/(z$se^2+m$tau2)));qq<-qt(.975,m$k-1)
 sens[[id]]<-data.frame(metric=id,modified_HK_se=ss,lower=as.numeric(m$b)-qq*ss,upper=as.numeric(m$b)+qq*ss)
}
s<-do.call(rbind,sums);rownames(s)<-NULL
write.csv(d,file.path(out,'cohort_effects.csv'),row.names=FALSE)
write.csv(s,file.path(out,'heterogeneity_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,loo),file.path(out,'leave_one_cohort_out.csv'),row.names=FALSE)
write.csv(do.call(rbind,sens),file.path(out,'modified_HK_sensitivity.csv'),row.names=FALSE)
saveRDS(fits,file.path(out,'meta_models.rds'));capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
# Forest-density quantitative grid; curves are normal uncertainty approximations.
labels<-c(raw_bias='Raw FI6 minus complete FI: mean bias',mapped_bias='Mapped FI6 minus complete FI: mean bias',future_fi_delta='Follow-up FI: AUROC(FI6) - AUROC(complete FI)',adl_delta='Incident ADL: AUROC(FI6) - AUROC(complete FI)',concurrent_auc='Concurrent FI: universal six-item AUROC',future_fi_auc='Follow-up FI: universal six-item AUROC',adl_auc='Incident ADL: universal six-item AUROC')

print(s)
