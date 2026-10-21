if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE,warn=1)
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(svglite)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

BASE <- score_source
OUT<-file.path(runtime_dir,'extensions/training_bootstrap_500')
files<-file.path(OUT,'checkpoints',sprintf('rep_%03d.rds',1:500));stopifnot(all(file.exists(files)))
a<-lapply(files,readRDS);stopifnot(identical(as.integer(vapply(a,function(z)z$rep,numeric(1))),1:500),all(vapply(a,function(z)all(unlist(z$QA)),logical(1))))
m<-do.call(rbind,lapply(a,'[[','metrics'));rk<-do.call(rbind,lapply(a,'[[','rankings'));mp<-do.call(rbind,lapply(a,'[[','mapping'));path<-do.call(rbind,lapply(a,'[[','length_curves'));avail<-do.call(rbind,lapply(a,'[[','availability'))
stopifnot(nrow(m)==75000,nrow(rk)==40000,nrow(mp)==2500,nrow(path)==480000,all(is.finite(m$estimate)))
write.csv(m,file.path(OUT,'all_bootstrap_metrics.csv'),row.names=FALSE);write.csv(mp,file.path(OUT,'all_mapping_threshold_draws.csv'),row.names=FALSE)
write.csv(avail,file.path(OUT,'score_availability.csv'),row.names=FALSE)
raw<-read.csv(file.path(runtime_dir,'extensions/approximation/FI6_continuous_agreement.csv'))
pros<-read.csv(file.path(runtime_dir,'extensions/approximation/FI6_paired_prospective_retention.csv'))
class<-read.csv(file.path(score_source,'FI6_locked_threshold_classification.csv'))
perf<-read.csv(file.path(score_source,'FI_item_number_performance.csv'))
cohorts<-c('HRS','ELSA','SHARE','CHARLS','MHAS')
point<-list()
for(c in cohorts){
 ra<-raw[raw$cohort==c&raw$score=='Raw FI6',];ma<-raw[raw$cohort==c&raw$score=='Training-mapped FI6',];cl<-class[class$cohort==c,];pe<-perf[perf$cohort==c&perf$strategy=='Universal'&perf$k==6&perf$outcome=='concurrent',]
 pt<-c(raw_bias=ra$bias,mapped_bias=ma$bias,concurrent_auc=pe$AUROC,raw_MAE=ra$MAE,mapped_MAE=ma$MAE,raw_RMSE=ra$RMSE,mapped_RMSE=ma$RMSE,sensitivity=cl$sensitivity,specificity=cl$specificity)
 for(o in c('followup_frail','incident_adl')){z<-pros[pros$cohort==c&pros$outcome==o,];pp<-if(o=='followup_frail')'future_fi' else 'adl';v<-c(z$AUROC_FI6,z$delta_AUROC,z$excess_AUROC_retention);names(v)<-paste0(pp,c('_auc','_delta','_retention'));pt<-c(pt,v)}
 point[[c]]<-data.frame(cohort=c,metric=names(pt),point=unname(pt))
}
point<-do.call(rbind,point)
sumrows<-lapply(split(m,interaction(m$cohort,m$mode,m$metric,drop=TRUE)),function(z){pt<-point$point[point$cohort==z$cohort[1]&point$metric==z$metric[1]];data.frame(cohort=z$cohort[1],mode=z$mode[1],metric=z$metric[1],point=pt,B=nrow(z),bootstrap_mean=mean(z$estimate),bootstrap_sd=sd(z$estimate),lower=quantile(z$estimate,.025,names=FALSE),upper=quantile(z$estimate,.975,names=FALSE),min_n=min(z$n),max_n=max(z$n))})
summary<-do.call(rbind,sumrows);rownames(summary)<-NULL;stopifnot(nrow(summary)==150)
write.csv(summary,file.path(OUT,'bootstrap_cohort_summary.csv'),row.names=FALSE)
relative<-list()
for(coh in cohorts)for(md in c('training_only','joint'))for(er in c('MAE','RMSE')){
 z<-m[m$cohort==coh&m$mode==md,];r<-z[z$metric==paste0('raw_',er),];mm<-z[z$metric==paste0('mapped_',er),];r<-r[order(r$rep),];mm<-mm[order(mm$rep),];vv<-1-mm$estimate/r$estimate
 pt<-1-point$point[point$cohort==coh&point$metric==paste0('mapped_',er)]/point$point[point$cohort==coh&point$metric==paste0('raw_',er)]
 relative[[paste(coh,md,er)]]<-data.frame(cohort=coh,mode=md,error_metric=er,point=pt,bootstrap_mean=mean(vv),sd=sd(vv),lower=quantile(vv,.025,names=FALSE),upper=quantile(vv,.975,names=FALSE))
}
write.csv(do.call(rbind,relative),file.path(OUT,'relative_error_improvement_full_bootstrap.csv'),row.names=FALSE)

selected<-rk[rk$rank<=6,];freq<-aggregate(I(rank<=6)~cohort+item,rk,mean);names(freq)[3]<-'selection_fraction'
# Full rank summaries retain all sixteen questions, including zero selection rates.
freq<-merge(freq,aggregate(rank~cohort+item,rk,mean),by=c('cohort','item'));names(freq)[4]<-'mean_rank'
write.csv(freq,file.path(OUT,'item_selection_stability.csv'),row.names=FALSE)
oldrank<-read.csv(file.path(score_source,'FI_universal_SHAP_rankings.csv'))
sets<-do.call(rbind,lapply(cohorts,function(c){original<-oldrank$item[oldrank$held_out==c&oldrank$SHAP_rank<=6];zz<-selected[selected$cohort==c,];same<-vapply(split(zz$item,zz$rep),function(x)setequal(x,original),logical(1));data.frame(cohort=c,rep=as.integer(names(same)),same_six_set=same)}))
write.csv(sets,file.path(OUT,'same_six_set_by_replicate.csv'),row.names=FALSE)
mp_summary<-do.call(rbind,lapply(cohorts,function(c)do.call(rbind,lapply(c('intercept','slope','threshold'),function(v){x<-mp[mp$cohort==c,v];data.frame(cohort=c,parameter=v,mean=mean(x),sd=sd(x),lower=quantile(x,.025,names=FALSE),upper=quantile(x,.975,names=FALSE),distinct_values=length(unique(x)))}))))
write.csv(mp_summary,file.path(OUT,'mapping_threshold_summary.csv'),row.names=FALSE)
threshold_freq<-aggregate(rep~cohort+threshold,mp,length);names(threshold_freq)[3]<-'count';threshold_freq$fraction<-threshold_freq$count/500;write.csv(threshold_freq,file.path(OUT,'threshold_distribution.csv'),row.names=FALSE)
curve<-do.call(rbind,lapply(split(path,interaction(path$cohort,path$mode,path$k,path$metric,drop=TRUE)),function(z)data.frame(cohort=z$cohort[1],mode=z$mode[1],k=z$k[1],metric=z$metric[1],bootstrap_mean=mean(z$estimate),sd=sd(z$estimate),lower=quantile(z$estimate,.025,names=FALSE),upper=quantile(z$estimate,.975,names=FALSE))))
write.csv(curve,file.path(OUT,'all_length_bootstrap_summary.csv'),row.names=FALSE)
# Coupled draws retain cross-cohort covariance from shared training samples.
old<-read.csv(file.path(runtime_dir,'extensions/heterogeneity/cohort_effects.csv'))
pooled<-list();covrows<-list()
for(mode in c('training_only','joint'))for(metric in unique(old$metric)){
 z<-m[m$mode==mode&m$metric==metric,];mat<-sapply(cohorts,function(c){v<-z[z$cohort==c,];v$estimate[order(v$rep)]})
 w<-old$weight[match(paste(cohorts,metric),paste(old$cohort,old$metric))]/100
 stopifnot(abs(sum(w)-1)<1e-9)
 draws<-as.vector(mat%*%w);cv<-cov(mat)
 pooled[[paste(mode,metric)]]<-data.frame(mode=mode,metric=metric,point=sum(w*point$point[match(paste(cohorts,metric),paste(point$cohort,point$metric))]),bootstrap_mean=mean(draws),bootstrap_sd=sd(draws),lower=quantile(draws,.025,names=FALSE),upper=quantile(draws,.975,names=FALSE),variance_with_covariance=as.numeric(t(w)%*%cv%*%w),variance_ignoring_covariance=sum(w*w*diag(cv)))
 covrows[[paste(mode,metric)]]<-transform(expand.grid(cohort1=cohorts,cohort2=cohorts),mode=mode,metric=metric,covariance=as.vector(cv))
}
po<-do.call(rbind,pooled);write.csv(po,file.path(OUT,'coupled_fixed_weight_pooled_summary.csv'),row.names=FALSE);write.csv(do.call(rbind,covrows),file.path(OUT,'cross_cohort_covariance.csv'),row.names=FALSE)
# Four-panel uncertainty comparison, keeping the established cohort palette.
cols<-c(HRS='#4E79A7',ELSA='#59A14F',SHARE='#E15759',CHARLS='#B07AA1',MHAS='#F28E2B')
ids<-c('concurrent_auc','future_fi_auc','adl_auc','mapped_bias')
ax<-c(concurrent_auc='Concurrent FI: AUROC',future_fi_auc='Follow-up FI: AUROC',adl_auc='Incident ADL: AUROC',mapped_bias='Mapped FI6 minus complete FI: mean bias')
old$mode<-'Fixed rules';old$point<-old$estimate
new<-summary[summary$mode=='joint',];new$mode<-'Training + evaluation'
pd<-rbind(old[,c('cohort','metric','mode','point','lower','upper')],new[,c('cohort','metric','mode','point','lower','upper')]);pd$cohort<-factor(pd$cohort,levels=rev(cohorts));pd$mode<-factor(pd$mode,levels=c('Fixed rules','Training + evaluation'))
plots<-lapply(ids,function(id){z<-pd[pd$metric==id,];ggplot(z,aes(x=point,y=cohort,colour=cohort,shape=mode,linetype=mode))+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',position=position_dodge(width=.42),width=.08,linewidth=.4)+geom_point(position=position_dodge(width=.42),size=2.4,stroke=.6)+scale_colour_manual(values=cols,guide='none')+scale_shape_manual(values=c(1,16),name=NULL)+scale_linetype_manual(values=c('dashed','solid'),name=NULL)+labs(x=ax[id],y=NULL)+theme_classic(base_family='Arial',base_size=9)+theme(axis.text=element_text(colour='#333333'),legend.position='bottom',plot.margin=margin(4,8,4,4))})
p<-(plots[[1]]|plots[[2]])/(plots[[3]]|plots[[4]])+plot_layout(guides='collect')+plot_annotation(tag_levels='a')&theme(legend.position='bottom')
for(fmt in c('png','tiff','svg')){
 file<-file.path(OUT,paste0('Figure_S7_training_bootstrap_500.',fmt))
 if(fmt=='png')agg_png(file,width=183,height=154,units='mm',res=600,background='white')
 if(fmt=='tiff')agg_tiff(file,width=183,height=154,units='mm',res=600,compression='lzw',background='white')
 if(fmt=='svg')svglite(file,width=183/25.4,height=154/25.4)
 print(p);dev.off()
}
agg_png(file.path(OUT,'Figure_S7_training_bootstrap_500_preview.png'),width=183,height=154,units='mm',res=170,background='white');print(p);dev.off()
qa<-data.frame(check=c('500_complete_replicates','7500_models_retrained','same_uid_kept_in_inner_fold','heldout_excluded_from_training','150_summary_rows','all_original_risk_sets_recorded'),pass=TRUE)
write.csv(qa,file.path(OUT,'analysis_QA.csv'),row.names=FALSE);capture.output(sessionInfo(),file=file.path(OUT,'R_sessionInfo.txt'))
print(aggregate(same_six_set~cohort,sets,mean));print(mp_summary);print(po[po$mode=='joint',])
