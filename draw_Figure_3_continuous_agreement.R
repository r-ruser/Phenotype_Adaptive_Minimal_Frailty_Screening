if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# Figure 3 contract: training-only scale mapping improves overall agreement,
# while remaining value-dependent differences define the measurement boundary.
# Archetype: quantitative grid, all five held-out cohorts in five rows x two columns.
# Run: Rscript scripts/draw_Figure_3_continuous_agreement.R [revision_root]
# Existing aggregated source data only; no model refitting or participant reanalysis.
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(grid);library(svglite);library(ragg);library(xml2)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

root <- publication_dir
out <- file.path(root,'figures'); dir.create(out,showWarnings=FALSE,recursive=TRUE)
src <- file.path(root,'figure_source'); dir.create(src,showWarnings=FALSE,recursive=TRUE)
bn <- read.csv(file.path(runtime_dir,'extensions/approximation/Bland_Altman_binned_source.csv'))
a <- read.csv(file.path(runtime_dir,'extensions/approximation/FI6_continuous_agreement.csv'))
cohorts <- c('HRS','ELSA','SHARE','CHARLS','MHAS')
scores <- c('Raw FI6','Training-mapped FI6')
stopifnot(nrow(a)==10L,all(a$n>0),all(is.finite(a$bias)),all(is.finite(a$MAE)),all(is.finite(a$RMSE)))
count <- aggregate(n~cohort+score,bn,sum)
idx <- match(paste(a$cohort,a$score),paste(count$cohort,count$score))
stopifnot(!anyNA(idx),all(a$n==count$n[idx]))
stopifnot(all(bn$mean>=0 & bn$mean<=1),all(bn$difference>=-.45 & bn$difference<=.75))
stopifnot(all(a$MAE[a$score==scores[2]]<a$MAE[a$score==scores[1]]))
bn$cohort <- factor(bn$cohort,levels=cohorts);bn$score <- factor(bn$score,levels=scores)
a$cohort <- factor(a$cohort,levels=cohorts);a$score <- factor(a$score,levels=scores)
hr <- rbind(data.frame(a[,c('cohort','score')],value=a$bias,kind='Mean bias'),
            data.frame(a[,c('cohort','score')],value=a$LoA_lower,kind='Descriptive 95% limits'),
            data.frame(a[,c('cohort','score')],value=a$LoA_upper,kind='Descriptive 95% limits'))
lab <- a
lab$label <- paste0('n = ',format(lab$n,big.mark=',',trim=TRUE),'   bias = ',sprintf('%+.3f',lab$bias),
                    '\nMAE = ',sprintf('%.3f',lab$MAE),'   RMSE = ',sprintf('%.3f',lab$RMSE))
ink <- '#203642';signal <- '#176B8A';family <- 'Arial'
width_mm <- 183;height_mm <- 230
p <- ggplot(bn,aes(mean,difference))+
  geom_hline(yintercept=0,colour='#BAC7CD',linewidth=.28)+
  geom_point(aes(size=n),shape=16,colour=signal,alpha=.30)+
  geom_hline(data=hr,aes(yintercept=value,linetype=kind),linewidth=.40,colour='#3F5964')+
  geom_text(data=lab,aes(x=.025,y=.72,label=label),inherit.aes=FALSE,hjust=0,vjust=1,
            size=6.6/2.845,family=family,colour=ink,lineheight=1.10)+
  facet_grid(cohort~score)+scale_size_area(max_size=2.4,guide='none')+
  scale_linetype_manual(values=c('Mean bias'='solid','Descriptive 95% limits'='dashed'),
                        breaks=c('Mean bias','Descriptive 95% limits'))+
  scale_x_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),
                     labels=function(x)sprintf('%.2f',x),expand=expansion(mult=c(0,.01)))+
  scale_y_continuous(limits=c(-.45,.75),breaks=c(-.4,0,.4),expand=expansion(mult=c(0,0)))+
  labs(x='Mean of short score and complete FI',y='Short score minus complete FI',linetype=NULL)+
  theme_classic(base_size=7,base_family=family)+
  theme(strip.background=element_rect(fill='#EDF4F7',colour=NA),
        strip.text=element_text(size=7.5,face='bold',colour=ink),
        legend.position='bottom',legend.text=element_text(size=6.9,colour=ink),
        legend.key.width=unit(8,'mm'),legend.margin=margin(1,1,1,1),
        panel.spacing=unit(2,'mm'),panel.spacing.x=unit(5,'mm'),axis.line=element_line(linewidth=.32,colour=ink),
        axis.ticks=element_line(linewidth=.28,colour=ink),axis.ticks.length=unit(1,'mm'),
        axis.text=element_text(size=6.5,colour=ink),axis.title=element_text(size=7,colour=ink),
        axis.title.x=element_text(margin=margin(t=4)),axis.title.y=element_text(margin=margin(r=5)),
        plot.margin=margin(3,3,3,3),plot.background=element_rect(fill='white',colour=NA))
stem <- file.path(out,'Figure_3')
ggsave(paste0(stem,'.svg'),p,width=width_mm,height=height_mm,units='mm',device=svglite::svglite,bg='white')
ggsave(paste0(stem,'.png'),p,width=width_mm,height=height_mm,units='mm',dpi=600,device=ragg::agg_png,bg='white')
ggsave(paste0(stem,'.tiff'),p,width=width_mm,height=height_mm,units='mm',dpi=600,device=ragg::agg_tiff,compression='lzw',bg='white')
ggsave(paste0(stem,'_preview.png'),p,width=width_mm,height=height_mm,units='mm',dpi=220,device=ragg::agg_png,bg='white')
write.csv(bn,file.path(src,'Figure_3_binned_source.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(a,file.path(src,'Figure_3_agreement_metrics.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(hr,file.path(src,'Figure_3_horizontal_references.csv'),row.names=FALSE,fileEncoding='UTF-8')
legend <- paste0('Figure 3. Continuous agreement before and after training-only FI-scale mapping. ',
  'Rows show the five held-out cohorts, and columns show the raw equal-weight six-item score and its linear FI-scale mapping estimated exclusively in the other four cohorts. ',
  'Within each cohort, both columns use the same participants. The horizontal coordinate is the mean of the short score and complete FI, and the vertical coordinate is short score minus complete FI. ',
  'Points aggregate observations into 0.02-by-0.02 bins, with area proportional to participant count. ',
  'Solid lines show mean bias; dashed lines show descriptive limits of agreement (mean difference ±1.96 SD), which describe the spread of participant differences. ',
  'These limits are distinct from confidence intervals for a mean. Each panel reports participant n, mean bias, MAE and RMSE. Mapping predictions are bounded to 0–1. ',
  'The value-dependent residual differences and MHAS mapped mean bias of −0.049 characterize the limits of numerical agreement. ',
  'Individual-level interchangeability requires prespecified acceptable error and independent validation. Supplementary Table 6 supplies participant-bootstrap intervals for the error estimates.')
writeLines(legend,file.path(out,'Figure_3_legend.txt'),useBytes=TRUE)
svg <- xml2::read_xml(paste0(stem,'.svg'));texts <- xml2::xml_find_all(svg,'.//*[local-name()="text"]')
stopifnot(length(texts)>50L)
writeLines(c('PASS: all five cohorts, both score types, ten panels.',
 'PASS: aggregated participant counts exactly match the unchanged agreement metrics.',
 'PASS: common plotting axes include every saved binned observation.',
 'PASS: n/bias/MAE/RMSE annotations originate in the saved result CSV.',
 'PASS: all mapped MAEs are lower than raw MAEs; residual bias is preserved.',
 'PASS: editable SVG text and R-only PNG/TIFF 600-dpi exports.',
 'Export size: 183 x 230 mm; preview 220 dpi.'),file.path(out,'Figure_3_QA.txt'),useBytes=TRUE)
cat('PASS: Figure 3 exports and source checks completed.\n')

# Official caption is shared with the submitted DOCX.
canonical_path <- file.path(root,"figure_source/figure_legends_master.csv")
if(file.exists(canonical_path)) {
  canonical <- read.csv(canonical_path,fileEncoding="UTF-8",stringsAsFactors=FALSE)
  writeLines(canonical$legend[canonical$figure_id=="3"],file.path(root,"figures/Figure_3_legend.txt"),useBytes=TRUE)
}
