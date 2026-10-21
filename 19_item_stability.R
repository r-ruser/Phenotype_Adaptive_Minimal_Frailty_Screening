if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

out<-score_source
p<-read.csv(file.path(out,"FI5_SHAP_importance.csv"));u<-read.csv(file.path(out,"FI_universal_SHAP_rankings.csv"))
stopifnot(nrow(p)==29*16,nrow(u)==5*16)
summarize<-function(x,type,group){
 do.call(rbind,lapply(split(x,x$item),function(z){
  v<-z$SHAP_rank
  data.frame(level=type,group=group,item=z$item[1],n_rankings=length(v),top3=sum(v<=3),top5=sum(v<=5),top8=sum(v<=8),top3_fraction=mean(v<=3),top5_fraction=mean(v<=5),top8_fraction=mean(v<=8),mean_rank=mean(v),median_rank=median(v),rank_q25=unname(quantile(v,.25)),rank_q75=unname(quantile(v,.75)))
 }))
}
rows<-list(summarize(p,"phenotype_all","all_29"),summarize(u,"universal_all","all_5"))
for(hold in unique(p$held_out))rows[[length(rows)+1]]<-summarize(p[p$held_out==hold,],"heldout_fold",hold)
for(hold in unique(u$held_out))rows[[length(rows)+1]]<-summarize(u[u$held_out==hold,],"universal_fold",hold)
for(hold in unique(p$held_out))for(k in sort(unique(p$phenotype[p$held_out==hold])))rows[[length(rows)+1]]<-summarize(p[p$held_out==hold&p$phenotype==k,],"phenotype",paste0(hold,"_P",k))
s<-do.call(rbind,rows);write.csv(s,file.path(out,"FI_item_stability.csv"),row.names=FALSE)
write.csv(p[,c("held_out","phenotype","item","SHAP_rank","mean_abs_SHAP")],file.path(out,"FI_item_stability_matrix.csv"),row.names=FALSE)
cat("PASS: top-3/5/8 frequency and SHAP-rank stability for 29 phenotypes and five universal folds\n")
