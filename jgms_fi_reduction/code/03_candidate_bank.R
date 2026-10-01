options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
z<-readRDS(file.path(base,"prepared_data.rds"));d0<-z$main
d<-d0[!is.na(d0$fi_complete),]
cohorts<-c("HRS","ELSA","SHARE","CHARLS","MHAS")
core<-c("shlt","chaira","climsa","stoopa","lifta","dimea","armsa","dressa","batha","eata","beda","moneya","medsa","toilta","mealsa","shopa")
stopifnot(setequal(unique(d$cohort),cohorts),all(vapply(z$panels[cohorts],function(x)all(core%in%colnames(x$deficits)),logical(1))))
bank<-paste0("q_",core);q<-matrix(NA_real_,nrow(d),length(bank),dimnames=list(NULL,bank))
for(coh in cohorts){
 target<-which(d$cohort==coh);allcoh<-which(d0$cohort==coh)
 panel<-z$panels[[coh]]$deficits
 stopifnot(length(allcoh)==nrow(panel),!anyDuplicated(d0$uid[allcoh]))
 ix<-match(d$uid[target],d0$uid[allcoh]);stopifnot(!anyNA(ix))
 q[target,paste0("q_",core)]<-panel[ix,core,drop=FALSE]
}
stopifnot(all(q[!is.na(q)]>=0&q[!is.na(q)]<=1),all(colSums(!is.na(q))>0))
dat<-cbind(d,as.data.frame(q))
coverage<-do.call(rbind,lapply(cohorts,function(coh){x<-dat[dat$cohort==coh,bank,drop=FALSE];data.frame(cohort=coh,item=bank,n=nrow(x),observed=colSums(!is.na(x)),availability=colMeans(!is.na(x)),unique_observed=vapply(x,function(v)length(unique(v[!is.na(v)])),integer(1)))}))
mapping<-read.csv(file.path(base,"reference_item_mapping.csv"));mapping<-mapping[mapping$cohort%in%cohorts&mapping$item%in%core,c("cohort","item","source","baseline","label","rule","domain")]
mapping$item<-paste0("q_",mapping$item)
mapping<-mapping[order(mapping$cohort,mapping$item),]
write.csv(coverage,file.path(out,"candidate_question_coverage.csv"),row.names=FALSE)
write.csv(mapping,file.path(out,"candidate_question_mapping.csv"),row.names=FALSE)
saveRDS(list(data=dat,bank=bank,cohorts=cohorts),file.path(out,"candidate_question_bank.rds"))
cat("PASS:",length(bank),"single-question candidates across",length(cohorts),"cohorts; main n",nrow(dat),"\n")
