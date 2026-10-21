if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# Identical linear training-only calibration for every universal score length.
# Paired 500 evaluation bootstraps condition on original ranks and fitted maps;
# these intervals are distinct from the 500 fully retrained six-item analyses.
options(stringsAsFactors=FALSE)
fpath<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);code<-dirname(normalizePath(fpath,winslash="/"))
source(file.path(code,"00_config.R"),encoding="UTF-8")
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
out<-file.path(runtime_dir,"extensions/all_length_calibration");dir.create(out,recursive=TRUE,showWarnings=FALSE)
b<-readRDS(file.path(score_source,"candidate_question_bank.rds"));folds<-readRDS(file.path(score_source,"FI_reduction_locked_scores.rds"))
cohorts<-c("HRS","ELSA","SHARE","CHARLS","MHAS");colors<-c(HRS="#4E79A7",ELSA="#59A14F",SHARE="#E15759",CHARLS="#B07AA1",MHAS="#F28E2B")
score<-function(d,items){x<-as.matrix(d[,items,drop=FALSE]);s<-rowMeans(x,na.rm=TRUE);s[rowSums(!is.na(x))<ceiling(.8*length(items))]<-NA_real_;s}
metrics<-list();maps<-list();contrasts<-list();set.seed(20261021)
for(coh in cohorts){
 f<-folds[[coh]];tr<-b$data[match(f$training_uid,b$data$uid),];stopifnot(!anyNA(tr$uid),!any(tr$cohort==coh))
 common<-complete.cases(f$universal)&complete.cases(f$phenotype_score);y<-f$complete_FI[common];raw<-f$universal[common,,drop=FALSE];pred<-raw
 for(k in 1:16){
  x<-score(tr,f$universal_rank[1:k]);ok<-is.finite(x)&is.finite(tr$fi_complete);fit<-lm(tr$fi_complete[ok]~x[ok]);cf<-unname(coef(fit));stopifnot(all(is.finite(cf)),cf[2]>0)
  pred[,k]<-pmin(1,pmax(0,cf[1]+cf[2]*raw[,k]))
  maps[[paste(coh,k)]]<-data.frame(cohort=coh,k=k,training_n=sum(ok),intercept=cf[1],slope=cf[2],items=paste(f$universal_rank[1:k],collapse=";"))
 }
 # One paired evaluation draw is used for all lengths and both scoring scales.
 errors<-list(Raw=sweep(raw,1,y,"-"),Mapped=sweep(pred,1,y,"-"));draws<-list()
 for(type in names(errors)){
  e<-errors[[type]];pt<-rbind(bias=colMeans(e),MAE=colMeans(abs(e)),RMSE=sqrt(colMeans(e^2)));dd<-array(NA_real_,c(500,3,16))
  # Restore seed per scale so paired raw/mapped draws use identical participants.
  set.seed(20261021+match(coh,cohorts))
  for(j in 1:500){ix<-sample.int(length(y),length(y),TRUE);z<-e[ix,,drop=FALSE];dd[j,,]<-rbind(colMeans(z),colMeans(abs(z)),sqrt(colMeans(z^2)))}
  draws[[type]]<-dd
  for(k in 1:16)for(m in 1:3)metrics[[paste(coh,type,k,m)]]<-data.frame(cohort=coh,k=k,scale=type,metric=rownames(pt)[m],n=length(y),estimate=pt[m,k],low=quantile(dd[,m,k],.025),high=quantile(dd[,m,k],.975))
 }
 for(m in 2:3){v<-draws$Mapped[,m,6]-draws$Mapped[,m,16];e<-errors$Mapped;val<-function(k)if(m==2)mean(abs(e[,k]))else sqrt(mean(e[,k]^2));contrasts[[paste(coh,m)]]<-data.frame(cohort=coh,metric=c("bias","MAE","RMSE")[m],n=length(y),six=val(6),sixteen=val(16),decrease_6_to_16=val(6)-val(16),relative_decrease_percent=100*(1-val(16)/val(6)),low=quantile(v,.025),high=quantile(v,.975))}
 cat("ALL-LENGTH CALIBRATION",coh,length(y),"\n");flush.console()
}
a<-do.call(rbind,metrics);mp<-do.call(rbind,maps);ct<-do.call(rbind,contrasts)
write.csv(a,file.path(out,"all_length_calibrated_errors.csv"),row.names=FALSE);write.csv(mp,file.path(out,"all_length_training_mappings.csv"),row.names=FALSE);write.csv(ct,file.path(out,"mapped_six_vs_sixteen_paired.csv"),row.names=FALSE)
macro<-aggregate(estimate~k+scale+metric,a,mean);write.csv(macro,file.path(out,"all_length_cohort_mean.csv"),row.names=FALSE)
# Six-item calibration must exactly reproduce the established six-item estimates.
old<-read.csv(file.path(runtime_dir,"extensions/approximation/FI6_continuous_agreement.csv"));qa<-c()
for(coh in cohorts)for(m in c("bias","MAE","RMSE")){v<-a$estimate[a$cohort==coh&a$k==6&a$scale=="Mapped"&a$metric==m];ref<-old[old$cohort==coh&old$score=="Training-mapped FI6",m];qa<-c(qa,abs(v-ref))}
stopifnot(max(qa)<1e-12,nrow(mp)==80,nrow(a)==480)
plots<-list();i<-0
for(m in c("MAE","RMSE"))for(type in c("Raw","Mapped")){
 i<-i+1;z<-a[a$metric==m&a$scale==type,];z$cohort<-factor(z$cohort,cohorts)
 plots[[i]]<-ggplot(z,aes(k,estimate,color=cohort,fill=cohort))+geom_ribbon(aes(ymin=low,ymax=high),alpha=.08,colour=NA)+geom_line(linewidth=.5)+geom_point(size=.85)+geom_vline(xintercept=6,linetype=3,colour="grey50",linewidth=.3)+scale_color_manual(values=colors)+scale_fill_manual(values=colors)+scale_x_continuous(breaks=c(1,4,6,8,12,16))+labs(x="Number of questions",y=paste(if(type=="Raw")"Raw score"else"Training-mapped score",m),tag=letters[i])+theme_classic(base_size=9,base_family="Arial")+theme(legend.position="bottom",legend.title=element_blank(),plot.tag=element_text(face="bold"))
}
fig<-wrap_plots(plots,ncol=2,guides="collect")&theme(legend.position="bottom")
for(ext in c("png","tiff","svg")){file<-file.path(out,paste0("Figure_S9_all_length_calibration.",ext));if(ext=="svg")ggsave(file,fig,width=7.1,height=5.4,device=svglite::svglite)else ggsave(file,fig,width=7.1,height=5.4,dpi=600,device=if(ext=="png")ragg::agg_png else ragg::agg_tiff)}
writeLines(c("PASS: all five held-out cohorts excluded from calibration training.","PASS: 80 length-specific training mappings and 480 error rows.","PASS: identical evaluation participants across 16 lengths and raw/mapped scales.",paste("PASS: maximum difference from original six-item mapping metrics",max(qa)),"Intervals: 500 paired evaluation bootstraps, fixed original rankings and maps; separate from retraining bootstrap."),file.path(out,"all_length_calibration_QA.txt"))
print(ct);print(macro[macro$k%in%c(6,16)&macro$scale=="Mapped",])
