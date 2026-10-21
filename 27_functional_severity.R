if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

BASE <- score_source
OUT<-file.path(runtime_dir,'extensions/training_bootstrap_500/sources')
d<-readRDS(file.path(BASE,'candidate_question_bank.rds'))$data;ff<-readRDS(file.path(BASE,'FI_reduction_locked_scores.rds'))
r<-list();hist<-list()
for(c in names(ff)){
 f<-ff[[c]];dd<-d[match(f$uid,d$uid),];ok<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
 for(pop in c('baseline_FI','common_scores')){
  z<-if(pop=='baseline_FI')dd else dd[ok,]
  a<-data.frame(cohort=c,population=pop,n=nrow(z))
  for(v in c('adl','iadl')){
   x<-z[[v]];x<-x[is.finite(x)];stopifnot(all(x>=0),all(x==as.integer(x)))
   a[[paste0(v,'_n')]]<-length(x);a[[paste0(v,'_mean')]]<-mean(x);a[[paste0(v,'_sd')]]<-sd(x);a[[paste0(v,'_median')]]<-median(x);a[[paste0(v,'_two_or_more')]]<-mean(x>=2)
   tt<-table(x);hist[[paste(c,pop,v)]]<-data.frame(cohort=c,population=pop,measure=v,count=as.integer(names(tt)),n=as.integer(tt),fraction=as.numeric(tt)/sum(tt))
  }
  r[[paste(c,pop)]]<-a
 }
}
write.csv(do.call(rbind,r),file.path(OUT,'functional_limitation_severity.csv'),row.names=FALSE)
write.csv(do.call(rbind,hist),file.path(OUT,'functional_limitation_count_distribution.csv'),row.names=FALSE)
print(do.call(rbind,r))
