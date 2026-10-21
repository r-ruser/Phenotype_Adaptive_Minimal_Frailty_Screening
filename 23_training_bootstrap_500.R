if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE,warn=1)
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
Sys.setenv(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages(library(xgboost))
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

BASE <- score_source
OUT <- file.path(runtime_dir,'extensions/training_bootstrap_500')
dir.create(file.path(OUT,'checkpoints'),recursive=TRUE,showWarnings=FALSE)
args<-commandArgs(TRUE);mode<-if(length(args))args[1] else 'run'
b<-readRDS(file.path(BASE,'candidate_question_bank.rds'));locked<-readRDS(file.path(BASE,'FI_reduction_locked_scores.rds'))
cohorts<-c('HRS','ELSA','SHARE','CHARLS','MHAS');bank<-b$bank;d<-b$data
X<-as.matrix(d[,bank]);storage.mode(X)<-'double';Y<-d$fi_complete;C<-d$cohort
cohort_ix<-setNames(lapply(cohorts,function(c)which(C==c)),cohorts)
F<-setNames(lapply(cohorts,function(c){f<-locked[[c]];ix<-match(f$uid,d$uid);stopifnot(!anyNA(ix),all(C[ix]==c));list(ix=ix,common=complete.cases(f$universal)&complete.cases(f$phenotype_score),followup_frail=f$followup_frail,incident_adl=f$incident_adl)}),cohorts)
rm(b,d,locked);gc()
auc<-function(y,s){a<-sum(y==1);z<-sum(y==0);if(!a||!z)return(NA_real_);(sum(rank(s)[y==1])-a*(a+1)/2)/(a*z)}
scoreall<-function(m,r){z<-m[,r,drop=FALSE];ans<-matrix(NA_real_,nrow(m),16);nn<-numeric(nrow(m));ss<-numeric(nrow(m));for(k in 1:16){ok<-is.finite(z[,k]);nn<-nn+ok;ss<-ss+replace(z[,k],!ok,0);v<-ss/nn;v[nn<ceiling(.8*k)]<-NA_real_;ans[,k]<-v};ans}

getfold<-function(ix,seed){u<-unique(ix);set.seed(seed);ff<-integer(length(u));for(gg in split(seq_along(u),interaction(C[u],Y[u]>=.25,drop=TRUE)))ff[gg]<-sample(rep(1:3,length.out=length(gg)));ff[match(ix,u)]}
# All copies of a bootstrap participant stay in the same inner fold.
runrep<-function(repno,save=TRUE){
 file<-file.path(OUT,'checkpoints',sprintf('rep_%03d.rds',repno));if(save&&file.exists(file))return(list(rep=repno,status='cached'))
 start<-proc.time()[3];set.seed(202610190+repno)
 samples<-lapply(cohort_ix,function(ix)sample(ix,length(ix),TRUE))
 # Shared cohort draws couple the overlapping training sets across outer folds.
 res<-list();ranks<-list();maps<-list();paths<-list();quality<-list()
 for(h in cohorts){
  train<-unlist(samples[setdiff(cohorts,h)],use.names=FALSE)
  stopifnot(all(C[train]!=h));seed<-202610190+repno*1000+match(h,cohorts)*10
  fi<-getfold(train,seed);imp<-numeric(16)
  for(j in 1:3){
   it<-train[fi!=j];iv<-train[fi==j];stopifnot(length(intersect(unique(it),unique(iv)))==0)
   set.seed(seed+j)
   fit<-xgb.train(params=list(objective='reg:squarederror',max_depth=3,eta=.06,subsample=.8,colsample_bytree=.8,min_child_weight=15,nthread=2),data=xgb.DMatrix(X[it,,drop=FALSE],label=Y[it],missing=NA_real_,nthread=2),nrounds=100,verbose=0)
   shap<-as.matrix(predict(fit,xgb.DMatrix(X[iv,,drop=FALSE],missing=NA_real_,nthread=2),predcontrib=TRUE))
   imp<-imp+colSums(abs(shap[,bank,drop=FALSE]))
  }
  ord<-order(-imp,seq_along(imp));rank<-bank[ord]
  trscore<-scoreall(X[train,,drop=FALSE],rank);tr6<-trscore[,6];good<-is.finite(tr6)
  lmfit<-lm(Y[train][good]~tr6[good]);cf<-unname(coef(lmfit));stopifnot(length(cf)==2,all(is.finite(cf)))
  cutoff<-unname(quantile(tr6[Y[train]>=.25],.15,type=1,na.rm=TRUE))
  f<-F[[h]];test<-f$ix;testscore<-scoreall(X[test,,drop=FALSE],rank);raw<-testscore[,6];mapped<-pmin(pmax(cf[1]+cf[2]*raw,0),1);full<-Y[test]
  ranks[[h]]<-data.frame(rep=repno,cohort=h,item=rank,rank=1:16,mean_abs_SHAP=imp[ord]/length(train))
  maps[[h]]<-data.frame(rep=repno,cohort=h,intercept=cf[1],slope=cf[2],threshold=cutoff,training_n=length(train),training_unique=length(unique(train)))
  for(evalmode in c('training_only','joint')){
   evalix<-if(evalmode=='training_only')seq_along(test) else match(samples[[h]],test)
   stopifnot(!anyNA(evalix))
   ok<-evalix[f$common[evalix]&is.finite(raw[evalix])]
   rr<-raw[ok];mm<-mapped[ok];ff<-full[ok];yy<-as.integer(ff>=.25)
   metric<-c(raw_bias=mean(rr-ff),mapped_bias=mean(mm-ff),concurrent_auc=auc(yy,rr),raw_MAE=mean(abs(rr-ff)),mapped_MAE=mean(abs(mm-ff)),raw_RMSE=sqrt(mean((rr-ff)^2)),mapped_RMSE=sqrt(mean((mm-ff)^2)),sensitivity=mean(rr[yy==1]>=cutoff),specificity=mean(rr[yy==0]<cutoff))
   res[[paste(h,evalmode,'concurrent')]]<-data.frame(rep=repno,cohort=h,mode=evalmode,metric=names(metric),estimate=unname(metric),n=length(ok))
   for(o in c('followup_frail','incident_adl')){
    ii<-ok[!is.na(f[[o]][ok])];y<-f[[o]][ii];a<-auc(y,raw[ii]);af<-auc(y,full[ii]);prefix<-if(o=='followup_frail')'future_fi' else 'adl'
    v<-c(a,a-af,(a-.5)/(af-.5));names(v)<-paste0(prefix,c('_auc','_delta','_retention'))
    res[[paste(h,evalmode,o)]]<-data.frame(rep=repno,cohort=h,mode=evalmode,metric=names(v),estimate=unname(v),n=length(ii))
   }
   # Re-evaluate the universal 1–16-item sequence on its original risk sets.
   for(k in 1:16){
    ss<-testscore[,k];ii<-evalix[f$common[evalix]&is.finite(ss[evalix])]
    vv<-c(concurrent_auc=auc(as.integer(full[ii]>=.25),ss[ii]),raw_MAE=mean(abs(ss[ii]-full[ii])),raw_RMSE=sqrt(mean((ss[ii]-full[ii])^2)),Spearman=suppressWarnings(cor(ss[ii],full[ii],method='spearman')))
    for(o in c('followup_frail','incident_adl')){jj<-ii[!is.na(f[[o]][ii])];vv<-c(vv,setNames(auc(f[[o]][jj],ss[jj]),if(o=='followup_frail')'future_fi_auc' else 'adl_auc'))}
    paths[[paste(h,evalmode,k)]]<-data.frame(rep=repno,cohort=h,mode=evalmode,k=k,metric=names(vv),estimate=unname(vv),n=length(ii))
   }
   quality[[paste(h,evalmode)]]<-data.frame(rep=repno,cohort=h,mode=evalmode,reference_common_n=sum(f$common[evalix]),score6_available_n=length(ok))
  }
 }
 ans<-list(rep=repno,metrics=do.call(rbind,res),rankings=do.call(rbind,ranks),mapping=do.call(rbind,maps),length_curves=do.call(rbind,paths),availability=do.call(rbind,quality),elapsed_seconds=unname(proc.time()[3]-start),QA=list(training_test_separated=TRUE,duplicate_participants_grouped=TRUE))
 if(save){tmp<-paste0(file,'.tmp');saveRDS(ans,tmp);stopifnot(file.rename(tmp,file))}
 list(rep=repno,status='success',seconds=ans$elapsed_seconds)
}
if(mode=='pilot'){
 print(system.time(ans<-runrep(1)));print(ans)
}else{
 todo<-which(!file.exists(file.path(OUT,'checkpoints',sprintf('rep_%03d.rds',1:500))))
 cat('Remaining',length(todo),'\n');flush.console()
 if(!length(todo)){cat('All 500 checkpoints already complete\n');quit(status=0)}
 cl<-parallel::makePSOCKcluster(8,outfile=file.path(Sys.getenv('TEMP'),'ms5_training_bootstrap_workers_20261019.log'));invisible(NULL)
 parallel::clusterEvalQ(cl,{Sys.setenv(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8');library(xgboost);NULL})
 parallel::clusterExport(cl,c('OUT','X','Y','C','cohorts','cohort_ix','F','bank','auc','scoreall','getfold','runrep'),envir=.GlobalEnv)
 for(i in seq(1,length(todo),by=8)){
  batch<-todo[i:min(i+7,length(todo))];zz<-parallel::parLapply(cl,batch,function(j)tryCatch(runrep(j),error=function(e)list(rep=j,status='ERROR',message=conditionMessage(e))))
  for(v in zz){cat(paste(unlist(v),collapse=' | '),'\n');if(identical(v$status,'ERROR'))saveRDS(v,file.path(OUT,sprintf('error_%03d.rds',v$rep)))}
  cat('COMPLETED',sum(file.exists(file.path(OUT,'checkpoints',sprintf('rep_%03d.rds',1:500)))),'OF 500\n');flush.console()
 }
 parallel::stopCluster(cl)
}
