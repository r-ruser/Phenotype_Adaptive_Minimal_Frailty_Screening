options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
suppressPackageStartupMessages(library(haven))
old <- readRDS(eligible_input)
meta <- read.csv(file.path(package_dir,"mappings/baseline_variable_metadata.csv"),fileEncoding="UTF-8")
configs <- list(HRS=list(id="hhidpn",b="r14",f="r15"), ELSA=list(id="idauniq",b="r9",f="r10"),
 SHARE=list(id="mergeid",b="r2",f="r4"),CHARLS=list(id="ID",b="r3",f="r4"),MHAS=list(id="unhhidnp",b="r5",f="r6"))
binary <- c("hibpe","diabe","hearte","cancre","stroke","arthre","lunge")
function_items <- c("walkra","dressa","batha","eata","beda","toilta","phonea","moneya","medsa","shopa","mealsa","housewka","mapa",
 "joga","walksa","walk1a","walk1kma","walk100a","sita","chaira","climsa","clim1a","stoopa","lifta","dimea","armsa","pusha")
cesd <- c("depres","effort","sleepr","whappy","flone","fsad","going","enlife")
euro <- c("depress","pessim","suicid","guilt","sleep","intrst","irritb","appett","fatig","concnt","enjoym","tearfl")
charls_cesd <- c("depresl","effortl","sleeprl","whappyl","flonel","botherl","goingl","mindtsl","fhopel","fearll")
panels <- list(); audit <- list(); clean <- list(); missing_codes <- list(); ancillary <- list()
num <- function(x) {v <- as.numeric(x); v[!is.finite(v)|v<0] <- NA_real_; v}
idstr <- function(x) if(is.numeric(x)) format(as.numeric(x),scientific=FALSE,trim=TRUE) else trimws(as.character(x))
score <- function(x,type) {
 v <- num(x)
 if(type %in% c("binary","reverse")) {v[!v %in% 0:1] <- NA; if(type=="reverse")v<-1-v}
 if(type=="difficulty") {v[!v %in% 0:2] <- NA; v <- as.numeric(v>0)}
 if(type %in% c("ordinal4","reverse4")) {v[!v %in% 1:4] <- NA; v<-(v-1)/3; if(type=="reverse4")v<-1-v}
 if(type=="ordinal5") {v[!v %in% 1:5] <- NA; v<-(v-1)/4}
 if(type=="sensory") {v[!v %in% 1:6] <- NA; v<-pmin((v-1)/4,1)}
 if(type=="bmi") {v[v<10|v>80] <- NA; v<-ifelse(v<18.5,.75,ifelse(v<25,0,ifelse(v<30,.25,ifelse(v<35,.5,ifelse(v<40,.75,1)))))}
 v
}
fi <- function(z,minobs=.8) {n<-rowSums(!is.na(z));v<-rowSums(z,na.rm=TRUE)/n;v[n<ceiling(ncol(z)*minobs)]<-NA;v}
for(coh in names(configs)) {
 cat("Building full reference:",coh,"\n")
 cf<-configs[[coh]];m<-meta[meta$cohort==coh,]; m<-m[!duplicated(m$variable),]
 suff<-sub(paste0("^",cf$b),"",m$variable)
 affect<-if(coh=="SHARE")euro else if(coh=="CHARLS")charls_cesd else if(coh=="MHAS")c(cesd,"ftired","energ") else cesd
 basic<-c(binary,function_items,affect,"shlt","sight","dsight","nsight","hearing","fall","fall_s","painfr","pain_s","urinai","urinai6m","bmi","mbmi")
 if("sight" %in% suff)basic<-setdiff(basic,c("dsight","nsight"))
 if("bmi" %in% suff)basic<-setdiff(basic,"mbmi")
 selected<-m[suff %in% basic,]; selected$suffix<-sub(paste0("^",cf$b),"",selected$variable)
 dd<-old[old$cohort==coh,]; zz<-matrix(NA_real_,nrow(dd),nrow(selected));ff<-zz
 colnames(zz)<-colnames(ff)<-selected$suffix
 dates<-data.frame(uid=dd$uid,year=NA_real_,month=NA_real_,fu_year=NA_real_,fu_month=NA_real_,fu_iwstat=NA_real_)
 for(relative in unique(selected$source)) {
  path <- source_file(relative)
  s<-selected[selected$source==relative,];vars<-c(s$variable,sub(cf$b,cf$f,s$variable,fixed=TRUE),paste0(cf$b,c("iwy","iwm","iwendy","iwendm","mstat")),paste0(cf$f,c("iwy","iwm","iwendy","iwendm","iwstat")))
  raw<-read_dta(path,col_select=any_of(c(cf$id,vars)))
  stopifnot(!anyDuplicated(idstr(raw[[cf$id]])))
  ix<-match(dd$id,idstr(raw[[cf$id]]))
  for(j in seq_len(nrow(s))) {
   v<-s$variable[j];su<-s$suffix[j];dom<-if(su %in% binary)"disease" else if(su %in% function_items)if(su %in% function_items[1:13])"daily_function" else "mobility" else if(su %in% affect)"affect" else "other_health"
   typ<-if(su %in% function_items)"difficulty" else if(su %in% c("shlt"))"ordinal5" else if(su %in% c("sight","dsight","nsight","hearing"))"sensory" else if(su %in% c("bmi","mbmi"))"bmi" else if(su %in% charls_cesd)if(su %in% c("whappyl","fhopel"))"reverse4" else "ordinal4" else if(su %in% c("whappy","enlife","energ"))"reverse" else "binary"
   # MHAS ftired and energy are distinct CESD symptoms; no aggregate score counted.
   zz[,su]<-score(raw[[v]][ix],typ)
   vf<-sub(cf$b,cf$f,v,fixed=TRUE);if(vf %in% names(raw))ff[,su]<-score(raw[[vf]][ix],typ)
   ov<-su %in% c("shlt","walkra",if(coh=="CHARLS")"walk100a", "effort","effortl","fatig","ftired","dressa","batha","eata","beda","toilta","phonea","moneya","medsa","shopa","mealsa","housewka","mapa")
   audit[[length(audit)+1]]<-data.frame(cohort=coh,item=su,source=relative,baseline=v,followup=if(vf %in% names(raw))vf else "NONE",label=attr(raw[[v]],"label"),rule=typ,domain=dom,screen_overlap=ov,background_overlap=su %in% setdiff(binary,"lunge"),observed_n=sum(!is.na(zz[,su])),eligible_n=nrow(dd),availability=mean(!is.na(zz[,su])),value_labels=paste(names(attr(raw[[v]],"labels")),attr(raw[[v]],"labels"),collapse="; "))
   vals<-as.numeric(raw[[v]][ix]);bad<-vals[!is.na(vals)&(vals<0|is.na(zz[,su]))];if(length(bad))missing_codes[[length(missing_codes)+1]]<-data.frame(cohort=coh,variable=v,code=names(table(bad)),n=as.integer(table(bad)))
  }
  for(nm in c("year","month","fu_year","fu_month","fu_iwstat")) {
   vv<-if(coh=="HRS")switch(nm,year=paste0(cf$b,"iwendy"),month=paste0(cf$b,"iwendm"),fu_year=paste0(cf$f,"iwendy"),fu_month=paste0(cf$f,"iwendm"),fu_iwstat=paste0(cf$f,"iwstat")) else switch(nm,year=paste0(cf$b,"iwy"),month=paste0(cf$b,"iwm"),fu_year=paste0(cf$f,"iwy"),fu_month=paste0(cf$f,"iwm"),fu_iwstat=paste0(cf$f,"iwstat"))
   if(vv %in% names(raw))dates[[nm]]<-num(raw[[vv]][ix])
  }
 }
 # Availability selection is an outcome-blind fixed rule; no selection by discrimination.
 keep<-colMeans(!is.na(zz))>=.5 & apply(zz,2,function(v)length(unique(v[!is.na(v)]))>1)
 zz<-zz[,keep,drop=FALSE];ff<-ff[,keep,drop=FALSE]
 a<-do.call(rbind,audit);a<-a[a$cohort==coh,];a<-a[match(colnames(zz),a$item),]
 dd$scr_srh<-ifelse(is.na(dd$srh),NA_real_,as.numeric(dd$srh>=4))
 dd$fi_complete<-fi(zz);dd$fu_fi_complete<-fi(ff);dd$fi_n_observed<-rowSums(!is.na(zz));dd$fi_n_items<-ncol(zz)
 dd$fi_no_screen<-fi(zz[,!a$screen_overlap,drop=FALSE]);dd$fi_no_context<-fi(zz[,!a$background_overlap,drop=FALSE]);dd$fi_no_both<-fi(zz[,!a$screen_overlap & !a$background_overlap,drop=FALSE])
 doms<-sapply(unique(a$domain),function(dom)fi(zz[,a$domain==dom,drop=FALSE]))
 dd$fi_domain_equal<-fi(doms)
 dd$interval_years<-(dates$fu_year+ifelse(is.na(dates$fu_month),6,dates$fu_month)/12)-(dates$year+ifelse(is.na(dates$month),6,dates$month)/12)
 dd$fu_status<-dates$fu_iwstat
 dd$fu_dead<-dates$fu_iwstat %in% if(coh=="HRS")c(2,3,5,6) else if(coh=="ELSA")c(3,5,6) else c(5,6)
 dd$fu_dead[is.na(dates$fu_iwstat)]<-NA
 panels[[coh]]<-list(deficits=zz,followup_deficits=ff,items=a,dates=dates)
 clean[[coh]]<-dd
 cat(coh,"items",ncol(zz),"FI available",sum(!is.na(dd$fi_complete)),"\n")
}


all_audit<-do.call(rbind,audit);all_audit$included<-FALSE
for(coh in names(clean))all_audit$included[all_audit$cohort==coh & all_audit$item %in% panels[[coh]]$items$item]<-TRUE
d<-do.call(rbind,clean);rownames(d)<-NULL
stopifnot(!anyDuplicated(d$uid),all(d$fi_complete>=0&d$fi_complete<=1,na.rm=TRUE))
saveRDS(list(main=d,panels=panels),file.path(base,"prepared_data.rds"))
write.csv(all_audit,file.path(base,"reference_item_mapping.csv"),row.names=FALSE,na="")
if(length(missing_codes))write.csv(do.call(rbind,missing_codes),file.path(base,"special_missing_codes.csv"),row.names=FALSE)
summary<-do.call(rbind,lapply(names(clean),function(coh){z<-clean[[coh]];p<-panels[[coh]];data.frame(cohort=coh,eligible_n=nrow(z),reference_items=ncol(p$deficits),reference_available=sum(!is.na(z$fi_complete)),followup_reference_available=sum(!is.na(z$fu_fi_complete)))}))
write.csv(summary,file.path(base,"reference_summary.csv"),row.names=FALSE)
print(table(d$cohort[!is.na(d$fi_complete)]))
