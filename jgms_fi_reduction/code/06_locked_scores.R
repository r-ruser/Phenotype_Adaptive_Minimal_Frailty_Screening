options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
source(file.path(package_dir,"code/02_functions.R"),encoding="UTF-8")
b<-readRDS(file.path(out,"candidate_question_bank.rds"));d<-b$data;bank<-b$bank
rank<-read.csv(file.path(out,"FI_universal_SHAP_rankings.csv"));pi<-read.csv(file.path(out,"FI5_SHAP_importance.csv"));prior<-read.csv(file.path(out,"FI5_SHAP_locked_predictions.csv"))
cohorts<-c("HRS","ELSA","SHARE","CHARLS","MHAS");ks<-1:16
score<-function(x,items,minfrac=.8,complete=FALSE){m<-as.matrix(x[,items,drop=FALSE]);storage.mode(m)<-"double";n<-rowSums(!is.na(m));v<-rowMeans(m,na.rm=TRUE);v[n<if(complete)length(items) else ceiling(minfrac*length(items))]<-NA_real_;v}
make_scores<-function(x,ranked,group=NULL,complete=FALSE){ans<-matrix(NA_real_,nrow(x),length(ks),dimnames=list(NULL,paste0("FI",ks)))
 if(is.null(group)){for(k in ks)ans[,k]<-score(x,ranked[1:k],complete=complete)}else{
  for(g in sort(unique(group))){ix<-which(group==g);rr<-ranked[[as.character(g)]];stopifnot(length(rr)==16);for(k in ks)ans[ix,k]<-score(x[ix,,drop=FALSE],rr[1:k],complete=complete)}
 }
 ans}
folds<-list();choices<-list()
for(hold in cohorts){
 tr<-d[d$cohort!=hold,,drop=FALSE];te<-d[d$cohort==hold,,drop=FALSE]
 lca<-readRDS(file.path(out,paste0("FI5_SHAP_holdout_",hold,".rds")))$lca
 stopifnot(!any(tr$uid%in%te$uid))
 tc<-latent_predict(as.matrix(tr[,context,drop=FALSE]),lca)$class
 ec<-prior$phenotype[match(te$uid,prior$uid)];stopifnot(!anyNA(ec))
 ur<-rank$item[rank$held_out==hold][order(rank$SHAP_rank[rank$held_out==hold])];stopifnot(length(ur)==16)
 pr<-setNames(lapply(seq_len(lca$k),function(g){v<-pi[pi$held_out==hold&pi$phenotype==g,];v$item[order(v$SHAP_rank)]}),as.character(seq_len(lca$k)))
 u_train<-make_scores(tr,ur);p_train<-make_scores(tr,pr,tc)
 u_test<-make_scores(te,ur);p_test<-make_scores(te,pr,ec)
 u_cc<-make_scores(te,ur,complete=TRUE);p_cc<-make_scores(te,pr,ec,complete=TRUE)
 ytr<-as.numeric(tr$fi_complete>=.25)
 tu<-vapply(ks,function(k)threshold(ytr,u_train[,k],.85),numeric(1));tp<-vapply(ks,function(k)threshold(ytr,p_train[,k],.85),numeric(1))
 pprev<-prior[match(te$uid,prior$uid),];stopifnot(identical(pprev$uid,te$uid))
 folds[[hold]]<-list(cohort=hold,uid=te$uid,phenotype=ec,complete_FI=te$fi_complete,full_FI_items=unique(te$fi_n_items),followup_FI=te$fu_fi_complete,followup_frail=pprev$followup_frail,incident_adl=pprev$incident_adl,universal=u_test,phenotype_score=p_test,universal_complete_case=u_cc,phenotype_complete_case=p_cc,threshold_universal=tu,threshold_phenotype=tp,training_uid=tr$uid,universal_rank=ur,phenotype_rank=pr)
 choices[[hold]]<-data.frame(held_out=hold,k=ks,universal_items=vapply(ks,function(k)paste(ur[1:k],collapse="+"),character(1)),threshold85_universal=tu,threshold85_phenotype=tp,training_n=nrow(tr),heldout_n=nrow(te))
 cat("SCORED",hold,nrow(te),"\n");flush.console()
}
saveRDS(folds,file.path(out,"FI_reduction_locked_scores.rds"))
write.csv(do.call(rbind,choices),file.path(out,"FI_reduction_training_choices.csv"),row.names=FALSE)
cat("PASS: universal and phenotype FI1–FI16 scores, training-locked thresholds, 80% and complete-item rules\n")
