if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# Figure 4: quantitative grid of random six-item selection benchmarks.
# Core claim: the trained six-item set provides strong ranking/discrimination,
# while raw numerical agreement and SHAP's gain over simple ranking are separate evidence.
# Panel map: a concurrent AUROC; b follow-up FI AUROC; c incident ADL AUROC;
# d raw concurrent MAE. Backend: R exclusively; 183 x 155 mm; editable vectors.
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid);library(svglite);library(ragg);library(xml2)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

root <- publication_dir
data_dir <- file.path(runtime_dir,'extensions/benchmarks')
out <- file.path(root,'figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_dir <- file.path(root,'figure_source');dir.create(source_dir,recursive=TRUE,showWarnings=FALSE)
combinations <- read.csv(file.path(data_dir,'random_six_item_combinations.csv'))
analysis_seed <- unique(combinations$seed)
stopifnot(nrow(combinations)==1000L,length(unique(combinations$combination_id))==1000L,
          length(analysis_seed)==1L,analysis_seed==20261001L,all(combinations$k==6L))
d <- read.csv(file.path(data_dir,'random_combination_cohort_mean_distribution.csv'))
s <- read.csv(file.path(data_dir,'random_combination_distribution_summary.csv'))
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
stem <- file.path(out,'Figure_4')
svglite::svglite(paste0(stem,'.svg'),width=width_mm/25.4,height=height_mm/25.4,bg='white');draw();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=width_mm,height=height_mm,units='mm',res=600,background='white');draw();dev.off()
ragg::agg_tiff(paste0(stem,'.tiff'),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw',background='white');draw();dev.off()
ragg::agg_png(paste0(stem,'_preview.png'),width=width_mm,height=height_mm,units='mm',res=220,background='white');draw();dev.off()
write.csv(do.call(rbind,source),file.path(source_dir,'Figure_4_source_data.csv'),row.names=FALSE)
write.csv(do.call(rbind,bins),file.path(source_dir,'Figure_4_histogram_bins.csv'),row.names=FALSE)
write.csv(do.call(rbind,positions),file.path(source_dir,'Figure_4_reference_positions.csv'),row.names=FALSE)
legend <- 'Figure 4. Item-selection benchmarks for the six-item candidate score. Histograms show the distributions of equally weighted mean performance across the five held-out cohorts for 1,000 distinct six-item combinations sampled uniformly without replacement from 8,008 possible combinations (ANALYSIS_SEED_TOKEN). Panels a–c show AUROC for concurrent complete FI ≥0.25, follow-up FI ≥0.25, and new ADL limitation; panel d shows the mean absolute error (MAE) of the raw six-item mean relative to concurrent continuous complete FI, with lower values indicating closer numerical agreement. The solid teal line marks the trained six-item score. SHAP and training-only absolute-Spearman ranking selected the same six items in all five outer folds, yielding identical equal-weight scores. Dashed grey lines mark the 2.5th and 97.5th percentiles across random item combinations and characterize variability across combinations. Within each cohort and outcome, every method uses the same participants. The shared benchmark risk set requires all original FI1–FI16 scores to be available and at least 15 of 16 candidate items to be observed; at least five items are therefore observed for every six-item combination. The combined participant counts are 34,126 for panel a, 21,818 for panel b, 18,472 for panel c, and 34,126 for panel d. The benchmarks form part of the descriptive, exploratory length analysis.'
legend <- gsub('ANALYSIS_SEED_TOKEN',paste0('seed ',analysis_seed),legend,fixed=TRUE)
writeLines(legend,file.path(out,'Figure_4_legend.txt'),useBytes=TRUE)
# Export verification uses editable SVG structure; every figure device is R-only.
svg <- xml2::read_xml(paste0(stem,'.svg'))
texts <- xml2::xml_find_all(svg,'.//*[local-name()="text"]')
stopifnot(length(texts)>40L)
writeLines(c('PASS: four panels; each histogram sums to 1,000 fixed random combinations.',
 'PASS: seed read from saved combinations is 20261001.',
 'PASS: solid line marks identical SHAP and simple-correlation six-item scores.',
 'PASS: all random performance values are displayed within plotting limits.',
 'PASS: editable SVG text and R-only PNG/TIFF 600-dpi exports.',
 'PASS: dashed reference bounds are combination percentiles; participant CIs are separate.',
 paste0('Figure size: ',width_mm,' x ',height_mm,' mm; PNG/TIFF 600 dpi; preview 220 dpi.')),
 file.path(out,'Figure_4_QA.txt'),useBytes=TRUE)
cat('PASS: Figure 4 exported in SVG/600-dpi PNG/TIFF and PNG preview.\n')
# Official caption is shared with the submitted DOCX.
canonical_path <- file.path(root,"figure_source/figure_legends_master.csv")
if(file.exists(canonical_path)) {
  canonical <- read.csv(canonical_path,fileEncoding="UTF-8",stringsAsFactors=FALSE)
  writeLines(canonical$legend[canonical$figure_id=="4"],file.path(root,"figures/Figure_4_legend.txt"),useBytes=TRUE)
}
