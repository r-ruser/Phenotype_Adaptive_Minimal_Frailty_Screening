options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
source(file.path(package_dir,"code/02_functions.R"),encoding="UTF-8")
d<-readRDS(file.path(base,"prepared_data.rds"));scores<-readRDS(file.path(out,"FI_reduction_locked_scores.rds"))
bank<-c("shlt","chaira","climsa","stoopa","lifta","dimea","armsa","dressa","batha","eata","beda","moneya","medsa","toilta","mealsa","shopa")
ans<-list()
for(coh in names(scores)){
 f<-scores[[coh]];p<-d$panels[[coh]]$deficits
 base_coh<-d$main[d$main$cohort==coh,];ix<-match(f$uid,base_coh$uid);stopifnot(!anyNA(ix))
 cols<-setdiff(colnames(p),bank);z<-p[ix,cols,drop=FALSE]
 n<-rowSums(!is.na(z));alt<-rowSums(z,na.rm=TRUE)/n;alt[n<ceiling(.8*ncol(z))]<-NA_real_
 s<-f$universal[,6];full<-f$complete_FI;ok<-is.finite(s)&is.finite(full)&is.finite(alt)
 for(target in c("complete_FI","no_candidate_FI")){
  y<-if(target=="complete_FI")full else alt
  ans[[paste(coh,target)]]<-data.frame(cohort=coh,reference=target,reference_items=if(target=="complete_FI")ncol(p) else ncol(z),n=sum(ok),frail_fraction=mean(y[ok]>=.25),AUROC=auc_fast(as.numeric(y[ok]>=.25),s[ok]),Spearman=suppressWarnings(cor(s[ok],y[ok],method="spearman")))
 }
}
res<-do.call(rbind,ans)
write.csv(res,file.path(out,"JGMS_reference_overlap_sensitivity.csv"),row.names=FALSE)
print(res,row.names=FALSE)
cat("PASS: same-participant no-candidate-deficit reference sensitivity\n")
