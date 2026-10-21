if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors=FALSE)
Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
suppressPackageStartupMessages({library(grid);library(ragg);library(svglite)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

base <- file.path(runtime_dir,'extensions/heterogeneity')
out <- file.path(base,'forest_points');dir.create(out,showWarnings=FALSE)
d <- read.csv(file.path(base,'cohort_effects.csv'))
s <- read.csv(file.path(base,'heterogeneity_summary.csv'))
fits <- readRDS(file.path(base,'meta_models.rds'))
cohorts <- c('HRS','ELSA','SHARE','CHARLS','MHAS')
cols <- c(HRS='#4E79A7',ELSA='#59A14F',SHARE='#E15759',CHARLS='#B07AA1',MHAS='#F28E2B')
stopifnot(nrow(d)==35,nrow(s)==7)
# Weighted-circle forest plots from the unchanged analysis results.
labels<-c(raw_bias='Raw FI6 minus complete FI: mean bias',mapped_bias='Mapped FI6 minus complete FI: mean bias',future_fi_delta='Follow-up FI: AUROC(FI6) - AUROC(complete FI)',adl_delta='Incident ADL: AUROC(FI6) - AUROC(complete FI)',concurrent_auc='Concurrent FI: universal six-item AUROC',future_fi_auc='Follow-up FI: universal six-item AUROC',adl_auc='Incident ADL: universal six-item AUROC')
txt<-function(x,y,label,size=8.4,bold=FALSE,col='#303030',just='left')grid.text(label,x,y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
line<-function(x,y,col='#888888',lwd=.6,lty=1)grid.lines(x,y,gp=gpar(col=col,lwd=lwd,lty=lty))
pfmt<-function(x)if(x<.001)'< 0.001' else paste0('= ',formatC(x,3,format='f'))
plotpanel<-function(id,tag){
 z<-d[d$metric==id,];z<-z[match(cohorts,z$cohort),];sm<-s[s$metric==id,];m<-fits[[id]]
 lim<-range(z$estimate-3.5*z$se,z$estimate+3.5*z$se,sm$lower,sm$upper,sm$PI_lower,sm$PI_upper)
 if(!grepl('_auc$',id))lim<-range(lim,0)
 lim<-lim+c(-1,1)*diff(lim)*.06
 xx<-function(v).30+(v-lim[1])/diff(lim)*.265
 line(c(.015,.985),c(.985,.985),lwd=.55)
 txt(.02,.93,tag,9,TRUE)
 null<-if(grepl('_auc$',id)).5 else 0
 if(null>=lim[1] && null<=lim[2])line(rep(xx(null),2),c(.17,.87),col='#A7BBCB',lty=2)
 ys<-seq(.84,.44,length.out=5)
 for(j in 1:5){
  yy<-ys[j];co<-z$cohort[j];mu<-z$estimate[j];se<-z$se[j]
  txt(.02,yy+.015,co,8.5,TRUE);txt(.16,yy+.015,paste0('N=',format(z$n[j],big.mark=',',scientific=FALSE)),7.7)
  line(xx(c(z$lower[j],z$upper[j])),rep(yy+.012,2),col=cols[co],lwd=.65)
  # Circle area is directly proportional to the random-effects weight.
  grid.circle(x=xx(mu),y=yy+.012,r=unit(1.8*sqrt(z$weight[j]/20),'mm'),gp=gpar(fill=cols[co],col=NA))
  txt(.59,yy+.025,sprintf('%.4f [%.4f, %.4f]',mu,z$lower[j],z$upper[j]),7.7)
  txt(.59,yy-.026,sprintf('w=%.1f%%',z$weight[j]),7.2,col='#777777')
 }
 yy<-.30;txt(.02,yy,'Random-effects mean',8,TRUE)
 grid.polygon(xx(c(sm$lower,sm$estimate,sm$upper,sm$estimate)),c(yy,yy+.020,yy,yy-.020),gp=gpar(fill='#303030',col=NA))
 txt(.59,yy,sprintf('%.4f [%.4f, %.4f]',sm$estimate,sm$lower,sm$upper),7.7,TRUE)
 txt(.02,.235,'95% prediction interval',7.6,col='#777777')
 line(xx(c(sm$PI_lower,sm$PI_upper)),c(.235,.235),col='#89969E',lwd=1)
 for(v in c(sm$PI_lower,sm$PI_upper))line(rep(xx(v),2),c(.224,.246),col='#89969E')
 txt(.84,.81,sprintf('I\u00b2 = %.1f%%',sm$I2),8.3,TRUE)
 txt(.84,.75,sprintf('95%% CI %.1f\u2013%.1f',sm$I2_lower,sm$I2_upper),6.9,col='#777777')
 txt(.84,.66,sprintf('tau\u00b2 = %.5f',sm$tau2),7.6)
 txt(.84,.57,sprintf('Q(%d) = %.2f',sm$Q_df,sm$Q),7.6)
 txt(.84,.48,paste('P(Q)',pfmt(sm$p_Q)),7.6)
 txt(.84,.39,paste('P(mean)',pfmt(sm$p_mean)),7.6)
 ticks<-pretty(lim,3);ticks<-ticks[ticks>=lim[1]&ticks<=lim[2]]
 line(c(.30,.565),c(.15,.15),col='#555555')
 for(v in ticks){line(rep(xx(v),2),c(.15,.139),col='#555555');txt(xx(v),.109,formatC(v,digits=if(diff(lim)<.06)3 else 2,format='f'),7.3,just='center')}
 txt(.43,.038,labels[id],7.4,just='center')
}
draw<-function(ids){grid.newpage();pushViewport(viewport(layout=grid.layout(length(ids),1)));for(j in seq_along(ids)){pushViewport(viewport(layout.pos.row=j));plotpanel(ids[j],letters[j]);popViewport()};popViewport()}
exports<-function(ids,stem,height){
 ragg::agg_png(file.path(out,paste0(stem,'.png')),width=183,height=height,units='mm',res=600,background='white');draw(ids);dev.off()
 ragg::agg_png(file.path(out,paste0(stem,'_preview.png')),width=183,height=height,units='mm',res=170,background='white');draw(ids);dev.off()
 ragg::agg_tiff(file.path(out,paste0(stem,'.tiff')),width=183,height=height,units='mm',res=600,compression='lzw',background='white');draw(ids);dev.off()
 svglite::svglite(file.path(out,paste0(stem,'.svg')),width=183/25.4,height=height/25.4);draw(ids);dev.off()
}
exports(c('raw_bias','mapped_bias','future_fi_delta','adl_delta'),'Figure_S5_heterogeneity_calibration_and_loss',230)
exports(c('concurrent_auc','future_fi_auc','adl_auc'),'Figure_S6_heterogeneity_discrimination',177)
print(s)

qa<-transform(d,radius_mm=1.8*sqrt(weight/20),area_mm2=pi*(1.8*sqrt(weight/20))^2)
stopifnot(max(qa$area_mm2/qa$weight)-min(qa$area_mm2/qa$weight)<1e-10)
write.csv(qa,file.path(out,'point_size_QA.csv'),row.names=FALSE)
