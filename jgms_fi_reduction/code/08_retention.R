options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
p<-read.csv(file.path(out,"FI_item_number_performance.csv"))
ref<-read.csv(file.path(base,"reference_summary.csv"))
ref<-ref[match(c("HRS","ELSA","SHARE","CHARLS","MHAS"),ref$cohort),]
stopifnot(all(!is.na(ref$reference_items)))
rows<-list()
for(coh in ref$cohort)for(strategy in c("Universal","Phenotype"))for(k in 1:16){
 concurrent<-p[p$cohort==coh&p$strategy==strategy&p$k==k&p$outcome=="concurrent",]
 maxrow<-p[p$cohort==coh&p$strategy==strategy&p$k==16&p$outcome=="concurrent",]
 for(outcome in c("followup_frail","incident_adl")){
  z<-p[p$cohort==coh&p$strategy==strategy&p$k==k&p$outcome==outcome,]
  full<-p[p$cohort==coh&p$strategy=="Full_FI"&p$outcome==outcome,]
  if(!nrow(z)||!nrow(full))next
  rows[[length(rows)+1]]<-data.frame(cohort=coh,strategy=strategy,k=k,full_FI_items=ref$reference_items[ref$cohort==coh],item_reduction=1-k/ref$reference_items[ref$cohort==coh],outcome=outcome,n=z$n,short_AUROC=z$AUROC,full_FI_AUROC=full$AUROC,discrimination_retention=if(full$AUROC>.5)(z$AUROC-.5)/(full$AUROC-.5) else NA_real_,concurrent_AUROC=concurrent$AUROC,Spearman=concurrent$Spearman,correlation_retention=if(maxrow$Spearman>0)concurrent$Spearman/maxrow$Spearman else NA_real_)
 }
}
r<-do.call(rbind,rows);write.csv(r,file.path(out,"FI_information_retention.csv"),row.names=FALSE)
rules<-list();selected<-list()
for(cutoff in c(.90,.95,.975))for(definition in c("future_FI","future_FI_and_ADL")){
 for(k in 1:16){
  ok<-logical(length(ref$cohort));names(ok)<-ref$cohort
  for(coh in ref$cohort){
   a<-r[r$cohort==coh&r$strategy=="Universal"&r$k==k,]
   curve<-p[p$cohort==coh&p$strategy=="Universal"&p$outcome=="concurrent",]
   cur<-a[a$outcome=="followup_frail",];adl<-a[a$outcome=="incident_adl",]
   ok[coh]<-nrow(cur)==1&&is.finite(cur$discrimination_retention)&&cur$concurrent_AUROC>=cutoff*max(curve$AUROC,na.rm=TRUE)&&cur$Spearman>=.80&&cur$discrimination_retention>=.85&&(definition=="future_FI"||(nrow(adl)==1&&is.finite(adl$discrimination_retention)&&adl$discrimination_retention>=.85))
  }
  rules[[length(rules)+1]]<-data.frame(AUROC_fraction_of_max=cutoff,prospective_rule=definition,k=k,cohorts_meeting=sum(ok),cohorts=paste(names(ok)[ok],collapse="+"),meets_4_of_5=sum(ok)>=4,meets_5_of_5=all(ok))
 }
 z<-do.call(rbind,rules);z<-z[z$AUROC_fraction_of_max==cutoff&z$prospective_rule==definition,]
 kk<-z$k[z$meets_4_of_5];selected[[length(selected)+1]]<-data.frame(AUROC_fraction_of_max=cutoff,prospective_rule=definition,min_k_4_of_5=if(length(kk))min(kk) else NA_integer_,min_k_5_of_5=if(any(z$meets_5_of_5))min(z$k[z$meets_5_of_5]) else NA_integer_)
}
write.csv(do.call(rbind,rules),file.path(out,"FI_optimal_item_number.csv"),row.names=FALSE)
write.csv(do.call(rbind,selected),file.path(out,"FI_optimal_item_number_summary.csv"),row.names=FALSE)
cat("PASS: prospective discrimination retention and descriptive multi-cohort parsimony criteria\n")
