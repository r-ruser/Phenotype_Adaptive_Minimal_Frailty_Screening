if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(svglite)})
out<-file.path(publication_dir,"figures");paired<-read.csv(file.path(score_source,"FI_universal_vs_phenotype.csv"))
th<-theme_classic(base_size=10,base_family="Arial")
save3<-function(fig,stem,w,h){ggsave(file.path(out,paste0(stem,".png")),fig,width=w,height=h,dpi=600,device=ragg::agg_png,bg="white");ggsave(file.path(out,paste0(stem,".svg")),fig,width=w,height=h,device=svglite::svglite,bg="white")}
z<-paired[paired$comparison=="Phenotype_minus_Universal"&paired$k==6,]
z$cohort<-factor(z$cohort,levels=rev(c("HRS","ELSA","SHARE","CHARLS","MHAS")))
z$outcome<-factor(z$outcome,levels=c("concurrent","followup_frail","incident_adl"),labels=c("Concurrent full-FI frailty","Follow-up FI frailty","Incident ADL limitation"))
fig4<-ggplot(z,aes(x=delta_AUROC,y=cohort))+geom_vline(xintercept=0,lty=2,colour="#555555")+geom_segment(aes(x=delta_low,xend=delta_high,yend=cohort),linewidth=.8,colour="#176B8A")+geom_point(size=2.1,colour="#176B8A")+facet_wrap(~outcome,nrow=1,scales="free_x")+labs(title=NULL,subtitle=NULL,x="Phenotype-specific minus universal AUROC",y=NULL)+th+theme(strip.text=element_text(face="bold",size=9))

fig4j<-fig4+facet_wrap(~outcome,ncol=1,scales="free_x")+theme(strip.text=element_text(size=9,face="bold"))
save3(fig4j,"Figure_S3",6.5,8)
