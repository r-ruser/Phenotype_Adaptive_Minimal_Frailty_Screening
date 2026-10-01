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
stopifnot(setequal(unique(d$cohort),cohorts),length(bank)==16,!anyDuplicated(d$uid))
fi5<-function(data,items){x<-as.matrix(data[,items,drop=FALSE]);storage.mode(x)<-"double";n<-rowSums(!is.na(x));z<-rowMeans(x,na.rm=TRUE);z[n<4]<-NA_real_;z}
spearman<-function(x,y){ok<-is.finite(x)&is.finite(y);if(sum(ok)<20)return(NA_real_);suppressWarnings(cor(x[ok],y[ok],method="spearman"))}
select_shap<-function(data,seed){
 y<-data$fi_complete;stopifnot(nrow(data)>=300,all(is.finite(y)))
 fold<-folds_stratified(as.numeric(y>=.25),data$cohort,3,seed)
 sums<-setNames(numeric(length(bank)),bank);counts<-setNames(numeric(length(bank)),bank)
 for(j in 1:3){
  tr<-which(fold!=j);va<-which(fold==j)
  xtr<-as.matrix(data[tr,bank,drop=FALSE]);xva<-as.matrix(data[va,bank,drop=FALSE])
  storage.mode(xtr)<-"double";storage.mode(xva)<-"double"
  set.seed(seed+j)
  fit<-xgb.train(params=list(objective="reg:squarederror",max_depth=3,eta=.06,subsample=.8,colsample_bytree=.8,min_child_weight=15,nthread=2),data=xgb.DMatrix(xtr,label=y[tr],missing=NA_real_),nrounds=100,verbose=0)
  sh<-predict(fit,xgb.DMatrix(xva,missing=NA_real_),predcontrib=TRUE)
  sh<-as.matrix(sh);stopifnot(all(bank%in%colnames(sh)))
  sums<-sums+colSums(abs(sh[,bank,drop=FALSE]));counts<-counts+length(va)
 }
 importance<-sums/counts;ranked<-names(sort(importance,decreasing=TRUE))
 list(importance=importance,ranked=ranked,raw5=ranked[1:5])
}
choose_unique<-function(selected,ranked){key<-paste(sort(selected),collapse="+");if(!key%in%ranked$used)return(list(items=selected,adjusted=FALSE))
 for(pos in 5:1)for(candidate in ranked$ranking){z<-selected;z[pos]<-candidate;if(length(unique(z))!=5)next;k<-paste(sort(z),collapse="+");if(!k%in%ranked$used)return(list(items=z,adjusted=TRUE))}
 stop("Unable to find distinct five-question set")
}
pred<-list();choices<-list();importance_rows<-list();architecture<-list();models<-list()
for(hold in cohorts){
 tr<-d[d$cohort!=hold,,drop=FALSE];te<-d[d$cohort==hold,,drop=FALSE]
 seed<-20261120+match(hold,cohorts)*100
 lca<-latent_fit(as.matrix(tr[,context,drop=FALSE]),seed,starts=5)
 cltr<-latent_predict(as.matrix(tr[,context,drop=FALSE]),lca)$class
 clte<-latent_predict(as.matrix(te[,context,drop=FALSE]),lca)
 global<-select_shap(tr,seed+9)
 selected<-vector("list",lca$k);distinct<-vector("list",lca$k);raw<-vector("list",lca$k);adjusted<-logical(lca$k);used<-character()
 for(k in seq_len(lca$k)){
  sub<-tr[cltr==k,,drop=FALSE]
  imp<-select_shap(sub,seed+100+k)
  raw[[k]]<-imp$raw5
  choice<-choose_unique(imp$raw5,list(used=used,ranking=imp$ranked))
  selected[[k]]<-imp$raw5;distinct[[k]]<-choice$items;adjusted[k]<-choice$adjusted
  used<-c(used,paste(sort(choice$items),collapse="+"))
  importance_rows[[paste(hold,k)]]<-data.frame(held_out=hold,phenotype=k,training_n=nrow(sub),item=bank,mean_abs_SHAP=as.numeric(imp$importance[bank]),SHAP_rank=match(bank,imp$ranked),raw_top5=bank%in%imp$raw5,forced_unique_selected=bank%in%choice$items,unique_adjusted=choice$adjusted)
  choices[[paste(hold,k)]]<-data.frame(held_out=hold,phenotype=k,training_n=nrow(sub),position=1:5,item=imp$raw5,forced_unique_item=choice$items,unique_adjusted=choice$adjusted)
 }
 rows<-vector("list",lca$k)
 for(k in seq_len(lca$k)){
  ix<-which(clte$class==k);if(!length(ix))next
  z<-te[ix,,drop=FALSE];sc<-fi5(z,selected[[k]]);alt<-fi5(z,distinct[[k]]);un<-fi5(z,global$raw5)
  rows[[k]]<-data.frame(uid=z$uid,cohort=hold,phenotype=k,posterior_max=clte$confidence[ix],fi_complete=z$fi_complete,fi_frail=as.numeric(z$fi_complete>=.25),FI5=sc,forced_unique_FI5=alt,universal_FI5=un,FI5_items=paste(selected[[k]],collapse="+"),forced_unique_items=paste(distinct[[k]],collapse="+"),universal_items=paste(global$raw5,collapse="+"),followup_fi=z$fu_fi_complete,followup_frail=ifelse(is.na(z$fu_fi_complete),NA_real_,as.numeric(z$fu_fi_complete>=.25)),incident_adl=ifelse(is.na(z$adl)|z$adl>=1|is.na(z$fu_adl),NA_real_,as.numeric(z$fu_adl>=1)))
 }
 pred[[hold]]<-do.call(rbind,Filter(Negate(is.null),rows))
 architecture[[hold]]<-data.frame(held_out=hold,phenotypes=lca$k,training_n=nrow(tr),heldout_n=nrow(te),raw_distinct_sets=length(unique(vapply(raw,function(x)paste(sort(x),collapse="+"),character(1)))),forced_distinct_sets=length(unique(vapply(distinct,function(x)paste(sort(x),collapse="+"),character(1)))),unique_adjusted=sum(adjusted),universal_items=paste(global$raw5,collapse="+"))
 models[[hold]]<-list(lca=lca,selected=selected,forced_unique=distinct,global_items=global$raw5,training_uid=tr$uid,held_out=hold)
 saveRDS(models[[hold]],file.path(out,paste0("FI5_SHAP_holdout_",hold,".rds")))
 write.csv(do.call(rbind,pred),file.path(out,"FI5_SHAP_locked_predictions.partial.csv"),row.names=FALSE,na="")
 cat("FI5 SHAP",hold,"classes",lca$k,"raw distinct",architecture[[hold]]$raw_distinct_sets,"forced",sum(adjusted),"\n");flush.console()
}
p<-do.call(rbind,pred);stopifnot(nrow(p)==nrow(d),!anyDuplicated(p$uid),setequal(p$uid,d$uid))
write.csv(p,file.path(out,"FI5_SHAP_locked_predictions.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,choices),file.path(out,"FI5_SHAP_selected_questions.csv"),row.names=FALSE)
write.csv(do.call(rbind,importance_rows),file.path(out,"FI5_SHAP_importance.csv"),row.names=FALSE)
write.csv(do.call(rbind,architecture),file.path(out,"FI5_SHAP_architecture.csv"),row.names=FALSE)
cat("PASS: training-only latent phenotypes, cross-fitted XGBoost SHAP raw top-five FI5 and separate forced-unique sensitivity scores\n")
