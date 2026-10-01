suppressPackageStartupMessages(library(pROC))
context <- c("hypertension","diabetes","heart_disease","cancer","stroke","arthritis","not_married","low_education")
lambda_grid <- c(1,.1,.01,.001)
auc_fast <- function(y,p) {
 ok<-!is.na(y)&!is.na(p);y<-y[ok];p<-p[ok];n1<-sum(y==1);n0<-sum(y==0)
 if(!n1||!n0)return(NA_real_)
 (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
softmax <- function(z) {e<-exp(z-apply(z,1,max));e/rowSums(e)}
latent_predict <- function(z,fit) {
 n<-nrow(z);k<-nrow(fit$profile);lp<-matrix(log(fit$mix),n,k,byrow=TRUE)
 for(j in seq_len(ncol(z))) {obs<-which(!is.na(z[,j]));if(length(obs))for(cl in seq_len(k))lp[obs,cl]<-lp[obs,cl]+z[obs,j]*log(fit$profile[cl,j])+(1-z[obs,j])*log(1-fit$profile[cl,j])}
 p<-softmax(lp);list(class=max.col(p,ties.method="first"),confidence=apply(p,1,max),posterior=p)
}
latent_fit <- function(z,seed,starts=5) {
 z<-as.matrix(z[complete.cases(z),,drop=FALSE]);if(nrow(z)<300)stop("Insufficient context sample")
 key<-apply(z,1,paste0,collapse="");tab<-table(key);u<-do.call(rbind,strsplit(names(tab),""));storage.mode(u)<-"numeric";w<-as.numeric(tab);n<-sum(w);best_models<-list();stats<-list()
 for(k in 2:6) {
  best<-NULL;bestll<--Inf
  for(st in seq_len(starts)) {
   set.seed(seed+100*k+st);pr<-matrix(runif(k*ncol(u),.1,.9),k,ncol(u));mix<-rep(1/k,k);prev<--Inf
   for(iter in 1:500) {
    lp<-u %*% t(log(pr))+(1-u) %*% t(log(1-pr))+matrix(log(mix),nrow(u),k,byrow=TRUE)
    post<-softmax(lp);ll<-sum(w*(apply(lp,1,max)+log(rowSums(exp(lp-apply(lp,1,max))))))
    if(is.finite(prev)&&abs(ll-prev)<1e-7*(1+abs(prev)))break
    wr<-post*w;mass<-colSums(wr);mix<-pmax(mass/n,1e-8);mix<-mix/sum(mix)
    pr<-pmin(pmax(t(wr) %*% u/pmax(mass,1e-8),1e-6),1-1e-6);prev<-ll
   }
   if(ll>bestll){bestll<-ll;best<-list(profile=pr,mix=mix,ll=ll,k=k,iterations=iter,converged=iter<500)}
  }
  pp<-latent_predict(u,best)$posterior;cls<-max.col(pp);share<-sapply(1:k,function(v)sum(w[cls==v])/n)
  ent<-1+sum(w*rowSums(pp*log(pp+1e-12)))/(n*log(k))
  bic<--2*bestll+(k*ncol(u)+k-1)*log(n)
  stats[[k-1]]<-data.frame(k=k,BIC=bic,entropy=ent,min_class_fraction=min(share),converged=best$converged)
  best_models[[as.character(k)]]<-best
 }
 st<-do.call(rbind,stats);valid<-st$min_class_fraction>=.03 & st$entropy>=.5 & st$converged
 if(!any(valid))valid<-st$min_class_fraction>=.02
 if(!any(valid))valid<-rep(TRUE,nrow(st))
 idx<-which(valid)[which.min(st$BIC[valid])];ans<-best_models[[as.character(st$k[idx])]];ans$comparison<-st;ans
}
folds_stratified <- function(y,cohort,folds=3,seed=1) {
 set.seed(seed);f<-integer(length(y));for(ix in split(seq_along(y),interaction(y,cohort,drop=TRUE)))f[ix]<-sample(rep(seq_len(folds),length.out=length(ix)));f
}

threshold <- function(y,p,target=.85) {unname(quantile(p[y==1],1-target,type=1,na.rm=TRUE))}

