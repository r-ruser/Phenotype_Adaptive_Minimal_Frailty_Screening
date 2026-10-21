if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE)
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

BASE <- score_source
OUT<-file.path(runtime_dir,'extensions/training_bootstrap_500/sources')
p<-readRDS(file.path(dirname(BASE),'prepared_data.rds'));b<-readRDS(file.path(BASE,'candidate_question_bank.rds'));d<-b$data;folds<-readRDS(file.path(BASE,'FI_reduction_locked_scores.rds'))
cohorts<-c('HRS','ELSA','SHARE','CHARLS','MHAS');domains<-c('disease','daily_function','mobility','affect','other_health')
contribution<-list();summary<-list();items<-list()
for(c in cohorts){
 all<-p$main[p$main$cohort==c,];panel<-p$panels[[c]];z<-panel$deficits;ix<-match(d$uid[d$cohort==c],all$uid);stopifnot(!anyNA(ix));z<-z[ix,,drop=FALSE];zz<-d[d$cohort==c,];stopifnot(max(abs(rowMeans(z,na.rm=TRUE)-zz$fi_complete))<1e-10)
 a<-panel$items[match(colnames(z),panel$items$item),];stopifnot(!anyNA(a$domain));den<-rowSums(is.finite(z));sums<-sapply(domains,function(dom)rowSums(z[,a$domain==dom,drop=FALSE],na.rm=TRUE)/den)
 stopifnot(max(abs(rowSums(sums)-zz$fi_complete))<1e-10)
 f<-folds[[c]];common<-complete.cases(f$universal)&complete.cases(f$phenotype_score);dd<-data.frame(uid=zz$uid,cohort=c,sums)
 contribution[[c]]<-dd
 for(dom in domains){x<-z[common,a$domain==dom,drop=FALSE];n<-rowSums(is.finite(x));v<-rowMeans(x,na.rm=TRUE);v[n<ceiling(.8*ncol(x))]<-NA_real_;summary[[paste(c,dom)]]<-data.frame(cohort=c,domain=dom,items=ncol(x),common_n=sum(common),domain_calculable_n=sum(is.finite(v)),mean_domain_deficit=mean(v,na.rm=TRUE),mean_contribution_to_FI=mean(sums[common,dom]))}
 for(j in 1:ncol(z))items[[paste(c,j)]]<-data.frame(cohort=c,item=colnames(z)[j],domain=a$domain[j],observed_n=sum(is.finite(z[common,j])),mean_deficit=mean(z[common,j],na.rm=TRUE))
}
write.csv(do.call(rbind,summary),file.path(OUT,'FI_domain_distribution.csv'),row.names=FALSE);write.csv(do.call(rbind,items),file.path(OUT,'reference_item_distribution.csv'),row.names=FALSE)
cc<-do.call(rbind,contribution);cc<-cc[match(d$uid,cc$uid),];stopifnot(identical(cc$uid,d$uid))
f<-folds[['MHAS']];it<-match(f$training_uid,d$uid);te<-match(f$uid,d$uid);common<-complete.cases(f$universal)&complete.cases(f$phenotype_score)
score<-function(ix){x<-as.matrix(d[ix,f$universal_rank[1:6]]);v<-rowMeans(x,na.rm=TRUE);v[rowSums(is.finite(x))<5]<-NA_real_;v}
ts<-score(it);ok<-is.finite(ts);it<-it[ok];ts<-ts[ok];te<-te[common];es<-f$universal[common,6]
tr<-data.frame(key=sprintf('%.10f',ts),cc[it,domains]);target<-data.frame(key=sprintf('%.10f',es),cc[te,domains]);n<-table(target$key)
trm<-aggregate(.~key,tr,mean);tem<-aggregate(.~key,target,mean);mt<-match(tem$key,trm$key);stopifnot(!anyNA(mt));w<-as.numeric(n[tem$key])/sum(n)
comp<-data.frame(domain=domains,training_conditional_contribution=colSums(trm[mt,domains]*w),MHAS_contribution=colSums(tem[,domains]*w))
comp$training_minus_MHAS<-comp$training_conditional_contribution-comp$MHAS_contribution
write.csv(comp,file.path(OUT,'MHAS_score_standardized_domain_gap.csv'),row.names=FALSE)
ref<-read.csv(file.path(OUT,'MHAS_bias_decomposition.csv'))
stopifnot(abs(sum(comp$training_minus_MHAS)-ref$conditional_reference_FI_component)<1e-10)
print(comp)
