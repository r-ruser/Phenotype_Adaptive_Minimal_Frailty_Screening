# Figure S3: quantitative grid of random six-item selection benchmarks.
# Core claim: the trained six-item set provides strong ranking/discrimination,
# while raw numerical agreement and SHAP's gain over simple ranking are separate evidence.
# Panel map: a concurrent AUROC; b follow-up FI AUROC; c incident ADL AUROC;
# d raw concurrent MAE. Backend: R exclusively; 183 x 155 mm; editable vectors.
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid);library(svglite);library(ragg);library(xml2);library(pdftools)})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file,winslash="/",mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"),encoding="UTF-8")
out <- file.path(runtime_dir,"extensions/benchmarks")
d <- read.csv(file.path(out,'random_combination_cohort_mean_distribution.csv'))
s <- read.csv(file.path(out,'random_combination_distribution_summary.csv'))
s <- s[s$cohort=='Cohort_mean',]
ink <- '#203642';signal <- '#176B8A';random_fill <- '#83AEBB';neutral <- '#757D82';family <- 'Arial'
width_mm <- 183;height_mm <- 155
contract <- data.frame(panel=letters[1:4],outcome=c('concurrent','followup_frail','incident_adl','concurrent'),metric=c('AUROC','AUROC','AUROC','MAE'),panel_title=c('Concurrent FI classification','Follow-up FI classification','New ADL limitation','Raw FI approximation error'),xlabel=c('AUROC (cohort mean)','AUROC (cohort mean)','AUROC (cohort mean)','Raw-score MAE (cohort mean)'),xmin=c(.70,.62,.52,.08),xmax=c(.96,.85,.73,.21),binwidth=c(.01,.01,.01,.005))
plot_list <- list();bins <- list();source <- list();positions <- list()
theme_nature <- theme_classic(base_family=family,base_size=7)+theme(axis.line=element_line(colour=ink,linewidth=.32),axis.ticks=element_line(colour=ink,linewidth=.28),axis.ticks.length=unit(1.2,'mm'),axis.text=element_text(size=6.5,colour=ink),axis.title=element_text(size=7,colour=ink),axis.title.x=element_text(margin=margin(t=4)),axis.title.y=element_text(margin=margin(r=4)),plot.title=element_text(size=8,face='bold',colour=ink,hjust=0),plot.subtitle=element_text(size=6.5,colour=neutral,margin=margin(b=4)),plot.tag=element_text(size=9,face='bold',colour=ink),plot.tag.position=c(0,1),plot.margin=margin(7,8,5,7,'pt'),panel.grid=element_blank(),legend.position='none',plot.background=element_rect(fill='white',colour=NA))
for(i in seq_len(nrow(contract))){
 z<-contract[i,];vals<-d[d$outcome==z$outcome&d$metric==z$metric,];sumrow<-s[s$outcome==z$outcome&s$metric==z$metric,]
 stopifnot(nrow(vals)==1000L,nrow(sumrow)==1L,abs(sumrow$SHAP-sumrow$simple_correlation)<1e-12)
 stopifnot(all(vals$cohort_mean>=z$xmin&vals$cohort_mean<=z$xmax))
 br<-seq(z$xmin,z$xmax+1e-9,by=z$binwidth);hh<-hist(vals$cohort_mean,breaks=br,plot=FALSE,include.lowest=TRUE,right=TRUE)
 hd<-data.frame(panel=z$panel,bin_left=hh$breaks[-length(hh$breaks)],bin_right=hh$breaks[-1],count=hh$counts);stopifnot(sum(hd$count)==1000L)
 ymax<-max(hd$count)*1.28
 src<-data.frame(panel=z$panel,outcome=z$outcome,metric=z$metric,combination_id=vals$combination_id,cohort_mean=vals$cohort_mean);source[[i]]<-src;bins[[i]]<-hd
 bounds<-data.frame(value=c(sumrow$random_2p5,sumrow$random_97p5))
 p<-ggplot(hd)+geom_rect(aes(xmin=bin_left,xmax=bin_right,ymin=0,ymax=count),fill=random_fill,colour='white',linewidth=.18)+geom_vline(data=bounds,aes(xintercept=value),colour=neutral,linewidth=.35,linetype='22')+geom_vline(xintercept=sumrow$SHAP,colour=signal,linewidth=.65)+annotate('label',x=sumrow$SHAP-.004,y=ymax*.94,label=paste0('Six-item = ',formatC(sumrow$SHAP,format='f',digits=4)),colour=signal,fill='white',linewidth=0,label.padding=unit(.65,'mm'),size=6.5/2.845,fontface='bold',family=family,hjust=1)+scale_x_continuous(limits=c(z$xmin,z$xmax),breaks=if(i==4)seq(.08,.20,.03) else pretty(c(z$xmin,z$xmax),n=5),labels=function(v)formatC(v,format='f',digits=2),expand=expansion(mult=c(0,0)))+scale_y_continuous(limits=c(0,ymax),breaks=pretty(c(0,max(hd$count)),n=4),expand=expansion(mult=c(0,0)))+labs(title=z$panel_title,subtitle=paste0('n_sets = 1,000; n = ',format(sumrow$n,big.mark=',',trim=TRUE)),x=z$xlabel,y='Random item sets',tag=z$panel)+theme_nature
 if(i==4){p<-p+annotate('segment',x=.125,xend=.088,y=ymax*.855,yend=ymax*.855,colour=neutral,linewidth=.35,arrow=arrow(length=unit(1.6,'mm'),type='closed'))+annotate('text',x=.1065,y=ymax*.94,label='Lower error',colour=neutral,size=6.5/2.845,family=family)}
 plot_list[[i]]<-p
 positions[[i]]<-data.frame(panel=z$panel,outcome=z$outcome,metric=z$metric,n=sumrow$n,n_sets=1000L,selected_six_item=sumrow$SHAP,random_median=sumrow$random_median,random_2p5=sumrow$random_2p5,random_97p5=sumrow$random_97p5,xmin=z$xmin,xmax=z$xmax,bins=sum(hd$count))
}
fig <- (plot_list[[1]]|plot_list[[2]])/(plot_list[[3]]|plot_list[[4]])
draw <- function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA));pushViewport(viewport(x=.5,y=.55,width=1,height=.90));grid.draw(patchwork::patchworkGrob(fig));popViewport()
 grid.rect(x=unit(12,'mm'),y=unit(11,'mm'),width=unit(4,'mm'),height=unit(2.5,'mm'),gp=gpar(fill=random_fill,col=NA))
 grid.text('Random six-item combinations',x=unit(16,'mm'),y=unit(11,'mm'),just='left',gp=gpar(fontfamily=family,fontsize=6.8,col=ink))
 grid.lines(x=unit(c(84,91),'mm'),y=unit(c(11,11),'mm'),gp=gpar(col=signal,lwd=1.8))
 grid.text('Selected six-item set (SHAP / simple correlation)',x=unit(94,'mm'),y=unit(11,'mm'),just='left',gp=gpar(fontfamily=family,fontsize=6.8,col=ink))
 grid.lines(x=unit(c(12,19),'mm'),y=unit(c(5,5),'mm'),gp=gpar(col=neutral,lwd=1,lty='22'))
 grid.text('Random-set 2.5th–97.5th percentiles',x=unit(22,'mm'),y=unit(5,'mm'),just='left',gp=gpar(fontfamily=family,fontsize=6.8,col=ink))
}
stem <- file.path(out,'Supplementary_Figure_S3_item_selection_benchmarks')
svglite::svglite(paste0(stem,'.svg'),width=width_mm/25.4,height=height_mm/25.4,bg='white');draw();dev.off()
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family=family,bg='white');draw();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=width_mm,height=height_mm,units='mm',res=600,background='white');draw();dev.off()
ragg::agg_tiff(paste0(stem,'.tiff'),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw',background='white');draw();dev.off()
ragg::agg_png(paste0(stem,'_preview.png'),width=width_mm,height=height_mm,units='mm',res=220,background='white');draw();dev.off()
write.csv(do.call(rbind,source),file.path(out,'Figure_S3_source_data.csv'),row.names=FALSE)
write.csv(do.call(rbind,bins),file.path(out,'Figure_S3_histogram_bins.csv'),row.names=FALSE)
write.csv(do.call(rbind,positions),file.path(out,'Figure_S3_reference_positions.csv'),row.names=FALSE)
legend <- 'Supplementary Figure S3. Item-selection benchmarks for the six-item candidate score. Histograms show the distributions of equally weighted mean performance across the five held-out cohorts for 1,000 distinct six-item combinations sampled uniformly without replacement from 8,008 possible combinations (seed 20261001). Panels a–c show AUROC for concurrent complete FI ≥0.25, follow-up FI ≥0.25, and new ADL limitation; panel d shows the mean absolute error (MAE) of the raw six-item mean relative to concurrent continuous complete FI, with lower values indicating closer numerical agreement. The solid teal line marks the trained six-item score. SHAP and training-only absolute-Spearman ranking selected the same six items in all five outer folds, yielding identical equal-weight scores. Dashed grey lines mark the 2.5th and 97.5th percentiles across random item combinations and characterize variability across combinations. Within each cohort and outcome, every method uses the same participants. The shared benchmark risk set requires all original FI1–FI16 scores to be available and at least 15 of 16 candidate items to be observed; at least five items are therefore observed for every six-item combination. The combined participant counts are 34,126 for panel a, 21,818 for panel b, 18,472 for panel c, and 34,126 for panel d. The benchmarks form part of the descriptive, exploratory length analysis.'
writeLines(legend,file.path(out,'Figure_S3_legend.txt'),useBytes=TRUE)
# R-only PDF/font and editable-SVG checks; ASCII temporary copy for pdftools.
svg<-xml2::read_xml(paste0(stem,'.svg'));texts<-xml2::xml_find_all(svg,'.//*[local-name()="text"]');stopifnot(length(texts)>40)
qa_pdf<-file.path(tempdir(),'MS5_Figure_S3_QA.pdf');file.copy(paste0(stem,'.pdf'),qa_pdf,overwrite=TRUE);fonts<-pdftools::pdf_fonts(qa_pdf);info<-pdftools::pdf_info(qa_pdf);stopifnot(info$pages==1L,all(fonts$embedded))
write.csv(fonts,file.path(out,'Figure_S3_PDF_font_QA.csv'),row.names=FALSE)
writeLines(c('PASS: four panels, 1,000 random combinations per panel, all histogram counts sum to 1,000.', 'PASS: selected six-item reference equals SHAP and simple-correlation macro performance.', 'PASS: all random performance values are displayed within plotting limits.', 'PASS: editable SVG text and embedded PDF fonts; white background; R-only export.', 'Visual QA: preview and PDF rendering remain for inspection.', paste0('Figure size: ',width_mm,' x ',height_mm,' mm; PNG/TIFF 600 dpi; preview 220 dpi.')),file.path(out,'Figure_S3_QA.txt'),useBytes=TRUE)
cat('PASS: Figure S3 exported in PDF/SVG/600-dpi PNG/TIFF and preview.\n')

