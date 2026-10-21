if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE)
script_file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
package_dir<-dirname(normalizePath(script_file,winslash="/",mustWork=TRUE))
source(file.path(package_dir,"00_config.R"),encoding="UTF-8")

f<-readRDS(file.path(score_source,"FI_reduction_locked_scores.rds"))
wilson<-function(k,n){z<-1.96;p<-k/n;d<-1+z*z/n;mid<-(p+z*z/(2*n))/d;half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d;c(low=mid-half,high=mid+half)}
rows<-lapply(names(f),function(c){x<-f[[c]];common<-complete.cases(x$universal)&complete.cases(x$phenotype_score);sc<-x$universal[common,6];y<-x$complete_FI[common]>=.25;cut<-x$threshold_universal[6];pred<-sc>=cut
 TP<-sum(pred&y);FN<-sum(!pred&y);FP<-sum(pred&!y);TN<-sum(!pred&!y)
 z<-data.frame(cohort=c,threshold=cut,n=length(y),TP=TP,FN=FN,FP=FP,TN=TN,sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),PPV=TP/(TP+FP),NPV=TN/(TN+FN),positive_fraction=mean(pred),frail_fraction=mean(y),FI6_calculable_n=sum(is.finite(x$universal[,6])),baseline_FI_n=length(x$uid))
 for(nm in c("sensitivity","specificity","PPV","NPV")){kn<-switch(nm,sensitivity=c(TP,TP+FN),specificity=c(TN,TN+FP),PPV=c(TP,TP+FP),NPV=c(TN,TN+FN));ci<-wilson(kn[1],kn[2]);z[[paste0(nm,"_low")]]<-ci[1];z[[paste0(nm,"_high")]]<-ci[2]};z})
write.csv(do.call(rbind,rows),file.path(score_source,"FI6_locked_threshold_classification.csv"),row.names=FALSE)
