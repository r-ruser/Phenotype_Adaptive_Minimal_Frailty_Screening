options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
suppressPackageStartupMessages(library(xgboost))
source(file.path(package_dir,"code/02_functions.R"),encoding="UTF-8")
b<-readRDS(file.path(out,"candidate_question_bank.rds"));d<-b$data;bank<-b$bank
cohorts<-c("HRS","ELSA","SHARE","CHARLS","MHAS")
importance<-list();beeswarm<-list()
for(hold in cohorts){
 tr<-d[d$cohort!=hold,,drop=FALSE];seed<-20261120+match(hold,cohorts)*100+9
 fold<-folds_stratified(as.numeric(tr$fi_complete>=.25),tr$cohort,3,seed)
 sums<-setNames(numeric(length(bank)),bank);nshap<-0L
 for(j in 1:3){
  itr<-which(fold!=j);iva<-which(fold==j)
  xt<-as.matrix(tr[itr,bank,drop=FALSE]);xv<-as.matrix(tr[iva,bank,drop=FALSE]);storage.mode(xt)<-"double";storage.mode(xv)<-"double"
  set.seed(seed+j)
  fit<-xgb.train(params=list(objective="reg:squarederror",max_depth=3,eta=.06,subsample=.8,colsample_bytree=.8,min_child_weight=15,nthread=2),data=xgb.DMatrix(xt,label=tr$fi_complete[itr],missing=NA_real_),nrounds=100,verbose=0)
  sh<-as.matrix(predict(fit,xgb.DMatrix(xv,missing=NA_real_),predcontrib=TRUE))
  stopifnot(all(bank%in%colnames(sh)))
  sums<-sums+colSums(abs(sh[,bank,drop=FALSE]));nshap<-nshap+nrow(sh)
  set.seed(seed+1000+j);take<-sample.int(nrow(sh),min(350,nrow(sh)))
  beeswarm[[paste(hold,j)]]<-do.call(rbind,lapply(bank,function(item)data.frame(held_out=hold,inner_fold=j,item=item,SHAP=sh[take,item],feature_value=xv[take,item])))
 }
 imp<-sums/nshap;ranking<-names(sort(imp,decreasing=TRUE))
 expected<-readRDS(file.path(out,paste0("FI5_SHAP_holdout_",hold,".rds")))$global_items
 stopifnot(identical(ranking[1:5],expected))
 importance[[hold]]<-data.frame(held_out=hold,item=ranking,SHAP_rank=seq_along(bank),mean_abs_SHAP=as.numeric(imp[ranking]),training_n=nrow(tr))
 cat("UNIVERSAL SHAP",hold,paste(ranking[1:5],collapse="+"),"\n");flush.console()
}
write.csv(do.call(rbind,importance),file.path(out,"FI_universal_SHAP_rankings.csv"),row.names=FALSE)
write.csv(do.call(rbind,beeswarm),file.path(out,"FI_universal_SHAP_beeswarm_sample.csv"),row.names=FALSE,na="")
cat("PASS: all 16 universal ranks match locked models; beeswarm contributions are cross-fitted on outer training cohorts\n")
