if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(svglite)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

BASE <- score_source
OUT<-file.path(runtime_dir,'extensions/training_bootstrap_500/sources');dir.create(OUT,showWarnings=FALSE)
b<-readRDS(file.path(BASE,'candidate_question_bank.rds'));d<-b$data;bank<-b$bank;folds<-readRDS(file.path(BASE,'FI_reduction_locked_scores.rds'))
cohorts<-c('HRS','ELSA','SHARE','CHARLS','MHAS');cols<-c(HRS='#4E79A7',ELSA='#59A14F',SHARE='#E15759',CHARLS='#B07AA1',MHAS='#F28E2B')
mp<-read.csv(file.path(runtime_dir,'extensions/approximation/training_linear_mapping.csv'))
rawtab<-read.csv(file.path(runtime_dir,'extensions/approximation/FI6_continuous_agreement.csv'))
score<-function(data,items){m<-as.matrix(data[,items,drop=FALSE]);v<-rowMeans(m,na.rm=TRUE);v[rowSums(is.finite(m))<5]<-NA_real_;v}
panels<-list();pop<-list();qdist<-list();resbin<-list();gains<-list();set.seed(202610191)
for(c in cohorts){
 f<-folds[[c]];dd<-d[match(f$uid,d$uid),];common<-complete.cases(f$universal)&complete.cases(f$phenotype_score);raw<-f$universal[,6];cf<-mp[mp$cohort==c,];mapped<-pmin(pmax(cf$intercept+cf$slope*raw,0),1)
 for(which in c('baseline_FI','common_scores')){
  zz<-if(which=='baseline_FI')dd else dd[common,]
  fi<-zz$fi_complete
  row<-data.frame(cohort=c,population=which,n=nrow(zz),FI_mean=mean(fi),FI_sd=sd(fi),FI_p10=quantile(fi,.1),FI_p25=quantile(fi,.25),FI_median=median(fi),FI_p75=quantile(fi,.75),FI_p90=quantile(fi,.9),frailty_prevalence=mean(fi>=.25),age_mean=mean(zz$age,na.rm=TRUE))
  for(v in c('hypertension','diabetes','heart_disease','cancer','stroke','arthritis','adl','iadl')){row[[paste0(v,'_n')]]<-sum(!is.na(zz[[v]]));row[[paste0(v,'_positive')]]<-mean(zz[[v]]>0,na.rm=TRUE)}
  pop[[paste(c,which)]]<-row
  for(v in bank){x<-zz[[v]];qdist[[paste(c,which,v)]]<-data.frame(cohort=c,population=which,item=v,n=nrow(zz),observed_n=sum(is.finite(x)),missing_fraction=mean(!is.finite(x)),mean_deficit=mean(x,na.rm=TRUE),any_deficit=mean(x>0,na.rm=TRUE))}
 }
 z<-data.frame(cohort=c,FI=dd$fi_complete[common],raw=raw[common],mapped=mapped[common]);z$residual<-z$mapped-z$FI;z$bin<-cut(z$FI,breaks=c(seq(0,.8,.1),Inf),include.lowest=TRUE,right=FALSE)
 panels[[c]]<-z
 for(g in levels(z$bin)){
  zz<-z[z$bin==g,];if(nrow(zz)==0)next
  resbin[[paste(c,g)]]<-data.frame(cohort=c,FI_bin=g,n=nrow(zz),mean_FI=mean(zz$FI),mean_residual=mean(zz$residual),lower=mean(zz$residual)-1.96*sd(zz$residual)/sqrt(nrow(zz)),upper=mean(zz$residual)+1.96*sd(zz$residual)/sqrt(nrow(zz)),MAE=mean(abs(zz$residual)),RMSE=sqrt(mean(zz$residual^2)))
 }
 r<-z$raw-z$FI;m<-z$mapped-z$FI
 pt<-c(MAE_raw=mean(abs(r)),MAE_mapped=mean(abs(m)),MAE_relative_reduction=1-mean(abs(m))/mean(abs(r)),RMSE_raw=sqrt(mean(r*r)),RMSE_mapped=sqrt(mean(m*m)),RMSE_relative_reduction=1-sqrt(mean(m*m)/mean(r*r)))
 bt<-replicate(500,{ii<-sample.int(nrow(z),nrow(z),TRUE);c(1-mean(abs(m[ii]))/mean(abs(r[ii])),1-sqrt(mean(m[ii]^2)/mean(r[ii]^2)))})
 gains[[c]]<-data.frame(cohort=c,n=nrow(z),t(pt),MAE_reduction_lower=quantile(bt[1,],.025),MAE_reduction_upper=quantile(bt[1,],.975),RMSE_reduction_lower=quantile(bt[2,],.025),RMSE_reduction_upper=quantile(bt[2,],.975))
 rr<-rawtab[rawtab$cohort==c&rawtab$score=='Training-mapped FI6',];stopifnot(abs(mean(abs(m))-rr$MAE)<1e-10,abs(sqrt(mean(m*m))-rr$RMSE)<1e-10)
}
pop<-do.call(rbind,pop);qdist<-do.call(rbind,qdist);bins<-do.call(rbind,resbin);gains<-do.call(rbind,gains);plots<-do.call(rbind,panels)
write.csv(pop,file.path(OUT,'baseline_and_common_distribution.csv'),row.names=FALSE)
write.csv(qdist,file.path(OUT,'candidate_item_distribution.csv'),row.names=FALSE)
write.csv(bins,file.path(OUT,'mapped_residual_by_FI.csv'),row.names=FALSE)
write.csv(gains,file.path(OUT,'relative_error_improvement.csv'),row.names=FALSE)
# MHAS versus the exact four-cohort training population of its scale mapping.
cf<-mp[mp$cohort=='MHAS',];items<-strsplit(cf$items,';',fixed=TRUE)[[1]];train<-d[d$cohort!='MHAS',];ts<-score(train,items);tr<-data.frame(raw=ts,FI=train$fi_complete);tr<-tr[is.finite(tr$raw),];stopifnot(nrow(tr)==cf$training_n)
test<-panels[['MHAS']];tr$key<-sprintf('%.10f',tr$raw);test$key<-sprintf('%.10f',test$raw)
grp<-function(z){meanfi<-aggregate(FI~key,z,mean);n<-aggregate(FI~key,z,length);names(n)[2]<-'n';ans<-merge(meanfi,n,by='key');ans$raw<-as.numeric(ans$key);ans}
gtr<-grp(tr);gte<-grp(test);names(gtr)[names(gtr)=='FI']<-'training_mean_FI';names(gtr)[names(gtr)=='n']<-'training_n';names(gte)[names(gte)=='FI']<-'MHAS_mean_FI';names(gte)[names(gte)=='n']<-'MHAS_n'
compare<-merge(gte,gtr,by=c('key','raw'),all.x=TRUE);compare$mapped<-pmin(pmax(cf$intercept+cf$slope*compare$raw,0),1);support<-!is.na(compare$training_mean_FI)
cc<-compare[support,];w<-cc$MHAS_n/sum(cc$MHAS_n)
decomp<-data.frame(target='MHAS',training_n=nrow(tr),target_n=nrow(test),common_support_n=sum(cc$MHAS_n),support_fraction=sum(cc$MHAS_n)/nrow(test),training_mean_raw=mean(tr$raw),target_mean_raw=mean(test$raw),training_mean_FI=mean(tr$FI),target_mean_FI=mean(test$FI),target_mean_mapped=mean(test$mapped),target_bias=mean(test$residual),score_distribution_and_linear_fit_component=sum(w*(cc$mapped-cc$training_mean_FI)),conditional_reference_FI_component=sum(w*(cc$training_mean_FI-cc$MHAS_mean_FI)))
stopifnot(abs(sum(decomp[1,c('score_distribution_and_linear_fit_component','conditional_reference_FI_component')])-sum(w*(cc$mapped-cc$MHAS_mean_FI)))<1e-10)
write.csv(compare,file.path(OUT,'MHAS_exact_score_comparison.csv'),row.names=FALSE);write.csv(decomp,file.path(OUT,'MHAS_bias_decomposition.csv'),row.names=FALSE)
# Shared-score standardization is descriptive and separates score mix from conditional FI differences.
base_theme<-theme_classic(base_family='Arial',base_size=9)+theme(axis.text=element_text(colour='#303030'),legend.position='bottom',plot.margin=margin(6,8,6,5))
a<-ggplot(plots,aes(FI,colour=cohort))+stat_ecdf(linewidth=.7)+scale_colour_manual(values=cols)+coord_cartesian(xlim=c(0,.8))+labs(x='Complete baseline FI',y='Cumulative proportion',colour=NULL)+base_theme
q<-qdist[qdist$population=='common_scores',];q$item<-factor(q$item,levels=rev(bank));q$cohort<-factor(q$cohort,levels=cohorts)
bb<-ggplot(q,aes(cohort,item,fill=mean_deficit))+geom_tile(colour='white',linewidth=.3)+scale_fill_gradient(low='white',high='#176B8A',limits=c(0,1),name='Mean deficit')+scale_y_discrete(labels=function(x)sub('q_','',x))+labs(x=NULL,y=NULL)+base_theme+theme(axis.ticks=element_blank())
zz<-bins[bins$n>=30,]
c<-ggplot(zz,aes(mean_FI,mean_residual,colour=cohort,group=cohort))+geom_hline(yintercept=0,linetype='dashed',colour='#999999',linewidth=.35)+geom_errorbar(aes(ymin=lower,ymax=upper),width=.008,linewidth=.3)+geom_line(linewidth=.6)+geom_point(size=1.2)+scale_colour_manual(values=cols)+labs(x='Mean complete FI within bin',y='Mapped FI6 minus complete FI',colour=NULL)+base_theme
cp<-compare[compare$MHAS_n>=30&compare$training_n>=30,];long<-rbind(data.frame(raw=cp$raw,FI=cp$MHAS_mean_FI,group='MHAS'),data.frame(raw=cp$raw,FI=cp$training_mean_FI,group='Four-cohort training'),data.frame(raw=cp$raw,FI=cp$mapped,group='Training mapping'))
e<-ggplot(long,aes(raw,FI,colour=group,linetype=group))+geom_line(linewidth=.65)+geom_point(size=1.3)+scale_colour_manual(values=c('MHAS'=unname(cols['MHAS']),'Four-cohort training'='#566873','Training mapping'='#176B8A'))+scale_linetype_manual(values=c('MHAS'='solid','Four-cohort training'='solid','Training mapping'='dashed'))+labs(x='Raw six-item score',y='Mean complete FI / mapped value',colour=NULL,linetype=NULL)+base_theme+theme(legend.text=element_text(size=8))
p<-(a|bb)/(c|e)+plot_annotation(tag_levels='a')
for(fmt in c('png','tiff','svg')){file<-file.path(OUT,paste0('Figure_S8_heterogeneity_sources.',fmt));if(fmt=='png')agg_png(file,width=183,height=175,units='mm',res=600,background='white');if(fmt=='tiff')agg_tiff(file,width=183,height=175,units='mm',res=600,compression='lzw',background='white');if(fmt=='svg')svglite(file,width=183/25.4,height=175/25.4);print(p);dev.off()}
agg_png(file.path(OUT,'Figure_S8_heterogeneity_sources_preview.png'),width=183,height=175,units='mm',res=170,background='white');print(p);dev.off()
print(decomp);print(gains)
