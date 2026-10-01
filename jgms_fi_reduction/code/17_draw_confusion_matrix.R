# R-only four-quadrant confusion matrices with manuscript FI6 colours.
options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
suppressPackageStartupMessages({library(grid);library(svglite);library(ragg);library(xml2)})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file,winslash="/",mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"),encoding="UTF-8")
out<-file.path(runtime_dir,"extensions/confusion_matrix");dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(file.path(package_dir,"results/FI6_locked_threshold_classification.csv"))
cohorts<-c("HRS","ELSA","SHARE","CHARLS","MHAS")
d<-d[match(cohorts,d$cohort),]
stopifnot(nrow(d)==5,all(rowSums(d[,c("TP","FN","FP","TN")])==d$n),sum(d$n)==34202)
stopifnot(max(abs(d$sensitivity-d$TP/(d$TP+d$FN)))<1e-12,max(abs(d$specificity-d$TN/(d$TN+d$FP)))<1e-12,max(abs(d$PPV-d$TP/(d$TP+d$FP)))<1e-12,max(abs(d$NPV-d$TN/(d$TN+d$FN)))<1e-12)
d$threshold_fraction<-ifelse(abs(d$threshold-5/12)<1e-12,"5/12","11/24")
stopifnot(max(abs(d$threshold-ifelse(d$threshold_fraction=="5/12",5/12,11/24)))<1e-12)
# Same deep/light teal family used by the manuscript's FI6 comparison figures.
width_mm<-183;height_mm<-230;family<-"Arial"
ink<-"#203642";muted<-"#566873";line_col<-"#607E8D"
positive_header<-"#176B8A";agreement_fill<-"#83AEBB"
negative_header<-"#DBE8ED";discordance_fill<-"#DBE8ED"
pct<-function(v)paste0(formatC(100*v,format="f",digits=1),"%")
num<-function(v)format(v,big.mark=",",trim=TRUE,scientific=FALSE)
qa<-list();cells<-list()
text_mm<-function(label,x,y,size=7.2,face="plain",color=ink,just="center",rot=0,max_width=NULL){
 gp<-gpar(fontfamily=family,fontsize=size,fontface=face,col=color)
 grid.text(label,x=unit(x,"mm"),y=unit(y,"mm"),just=just,rot=rot,gp=gp)
 if(!is.null(max_width)&&rot==0){ww<-convertWidth(grobWidth(textGrob(label,gp=gp)),"mm",valueOnly=TRUE);qa[[length(qa)+1]]<<-data.frame(text=label,size_pt=size,width_mm=ww,available_mm=max_width,fits=ww<=max_width)}
}
rect_mm<-function(x,y,w,h,fill,col=NA,lwd=.5){grid.rect(x=unit(x,"mm"),y=unit(y,"mm"),width=unit(w,"mm"),height=unit(h,"mm"),gp=gpar(fill=fill,col=col,lwd=lwd))}
seg_mm<-function(x1,y1,x2,y2,lwd=.7){grid.lines(x=unit(c(x1,x2),"mm"),y=unit(c(y1,y2),"mm"),gp=gpar(col=line_col,lwd=lwd))}
sector<-function(code,count,total,cx,cy,angles,cohort,side=22){
 prop<-count/total;r<-side*sqrt(prop)
 a<-seq(angles[1],angles[2],length.out=181)*pi/180
 xx<-c(cx,cx+r*cos(a),cx);yy<-c(cy,cy+r*sin(a),cy)
 grid.polygon(x=unit(xx,"mm"),y=unit(yy,"mm"),gp=gpar(fill=if(code%in%c("TP","TN"))agreement_fill else discordance_fill,col=NA))
 cells[[length(cells)+1]]<<-data.frame(cohort=cohort,cell=code,n=count,reference_category_n=total,reference_percent=100*prop,rate_label=c(TP="TPR",FP="FPR",FN="FNR",TN="TNR")[[code]],radius_mm=r,sector_area_mm2=pi*r*r/4)
}
cohort_panel<-function(i,x0,top){
 z<-d[i,];text_mm(letters[i],x0,top,9,"bold",just="left",max_width=4)
 text_mm(z$cohort,x0+5,top,9.5,"bold",just="left",max_width=70)
 text_mm(paste0("n = ",num(z$n),"   |   FI6 cutoff = ",z$threshold_fraction),x0,top-5,7,color=muted,just="left",max_width=79)
 l<-x0+22;r<-l+44;cx<-l+22;gt<-top-20;gb<-gt-44;cy<-gt-22
 text_mm("Actual (complete FI)",cx,top-10,7.2,"bold",max_width=58)
 rect_mm(l+11,top-15,22,6,positive_header)
 rect_mm(l+33,top-15,22,6,negative_header)
 text_mm("Positive",l+11,top-15,7,"bold",color="white",max_width=20)
 text_mm("Negative",l+33,top-15,7,"bold",max_width=20)
 rect_mm(l-5,gt-11,6,22,positive_header)
 rect_mm(l-5,gt-33,6,22,negative_header)
 text_mm("Positive",l-5,gt-11,7,"bold",color="white",rot=90)
 text_mm("Negative",l-5,gt-33,7,"bold",rot=90)
 text_mm("Predicted (FI6)",l-11,cy,7.2,"bold",rot=90)
 pos<-z$TP+z$FN;neg<-z$FP+z$TN
 sector("TP",z$TP,pos,cx,cy,c(90,180),z$cohort)
 sector("FP",z$FP,neg,cx,cy,c(0,90),z$cohort)
 sector("FN",z$FN,pos,cx,cy,c(180,270),z$cohort)
 sector("TN",z$TN,neg,cx,cy,c(270,360),z$cohort)
 rect_mm(cx,cy,44,44,NA,line_col,.9)
 seg_mm(cx,gb,cx,gt,.75);seg_mm(l,cy,r,cy,.75)
 # Rates remain visible at each outer corner, outside the quarter-disc.
 text_mm(pct(z$sensitivity),l+1.8,gt-3.2,8,"bold",just="left",max_width=18)
 text_mm(pct(1-z$specificity),r-1.8,gt-3.2,8,"bold",just="right",max_width=18)
 text_mm(pct(1-z$sensitivity),l+1.8,gb+3.2,8,"bold",just="left",max_width=18)
 text_mm(pct(z$specificity),r-1.8,gb+3.2,8,"bold",just="right",max_width=18)
 # Cell labels identify rate and count; area encodes the reference-column percentage.
 for(j in 1:4){
  xx<-c(l+11,l+33,l+11,l+33)[j];yy<-c(gt-12,gt-12,gt-32,gt-32)[j]
  codes<-c("TPR","FPR","FNR","TNR");counts<-c(z$TP,z$FP,z$FN,z$TN)
  text_mm(codes[j],xx,yy+1.5,8.2,"bold",max_width=20)
  text_mm(paste0("n = ",num(counts[j])),xx,yy-2.7,6.8,max_width=20)
 }
 text_mm(paste0("PPV ",pct(z$PPV),"   ·   NPV ",pct(z$NPV)),cx,gb-5,7,max_width=70)
}
draw<-function(){
 grid.newpage();grid.rect(gp=gpar(fill="white",col=NA))
 cohort_panel(1,7,225);cohort_panel(2,95,225)
 cohort_panel(3,7,149);cohort_panel(4,95,149)
 cohort_panel(5,51,73)
}
stem<-file.path(out,"Supplementary_Figure_S1_FI6_confusion_matrix")
svglite::svglite(paste0(stem,".svg"),width=width_mm/25.4,height=height_mm/25.4,bg="white");draw();dev.off()
grDevices::cairo_pdf(paste0(stem,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family=family,bg="white");draw();dev.off()
ragg::agg_png(paste0(stem,".png"),width=width_mm,height=height_mm,units="mm",res=600,background="white");draw();dev.off()
ragg::agg_tiff(paste0(stem,".tiff"),width=width_mm,height=height_mm,units="mm",res=600,compression="lzw",background="white");draw();dev.off()
qa<-list();cells<-list()
ragg::agg_png(paste0(stem,"_preview.png"),width=width_mm,height=height_mm,units="mm",res=220,background="white");draw();dev.off()
fit<-do.call(rbind,qa);source_cells<-do.call(rbind,cells)
stopifnot(all(fit$fits),nrow(source_cells)==20)
stopifnot(all(abs(aggregate(reference_percent~cohort+reference_category_n,source_cells,sum)$reference_percent-100)<1e-10))
stopifnot(max(abs(source_cells$sector_area_mm2/(pi*22^2/4)-source_cells$reference_percent/100))<1e-12)
write.csv(source_cells,file.path(out,"Figure_S1_cell_source_data.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(d,file.path(out,"Figure_S1_classification_metrics.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(fit,file.path(out,"Figure_S1_text_fit_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
legend<-paste0("Supplementary Figure S1. Cohort-specific confusion matrices for the universal six-item score. Panels a–e show HRS, ELSA, SHARE, CHARLS, and MHAS. Columns indicate actual complete-FI reference status (positive, FI ≥0.25; negative, FI <0.25), and rows indicate the FI6 classification. TPR, FPR, FNR, and TNR denote true-positive, false-positive, false-negative, and true-negative rates, respectively. Percentages and quarter-disc areas are normalized within each reference-status column; TPR + FNR = 100% and FPR + TNR = 100%. The displayed counts identify participants in each classification cell. Darker teal sectors indicate concordant classifications; lighter teal sectors indicate discordant classifications. Cutoffs were selected in the four training cohorts to target 85% sensitivity and applied unchanged to the held-out cohort: 5/12 (approximately 0.4167) for ELSA and SHARE and 11/24 (approximately 0.4583) for HRS, CHARLS, and MHAS. Analyses use the common concurrent samples (total n = 34,202). PPV and NPV denote positive and negative predictive values. Classification metrics with Wilson 95% intervals are supplied in the accompanying source table.")
writeLines(legend,file.path(out,"Figure_S1_legend.txt"),useBytes=TRUE)
svg<-read_xml(paste0(stem,".svg"));labels<-xml_text(xml_find_all(svg,"//*[local-name()='text']"))
stopifnot(any(grepl("82.0%",labels,fixed=TRUE)),any(grepl("92.3%",labels,fixed=TRUE)),!any(grepl("Reading the matrices",labels,fixed=TRUE)),!any(grepl("Sensitivity =",labels,fixed=TRUE)),any(grepl("11/24",labels,fixed=TRUE)),any(grepl("5/12",labels,fixed=TRUE)))
writeLines(c("Figure S1 numerical and export QA: PASS.","R-only quantitative grid: five four-quadrant matrices, 20 classification cells, total n = 34,202.","Columns: actual complete-FI reference; rows: predicted FI6 classification.","Reference-column percentages: TPR+FNR and FPR+TNR each sum to 100%.","Quarter-disc area is proportional to the displayed rate; radius = 22 mm * sqrt(rate).","Counts and all metrics match the locked classification aggregate table.","Exact cutoff fractions shown: 5/12 and 11/24.","230 x 183 mm; PNG/TIFF 600 dpi; editable SVG; vector PDF.","Manuscript FI6 colours: #176B8A and #83AEBB with a pale tint.","Explanatory reading block removed.",paste("All measured text labels fit:",nrow(fit)),"Sensitivity range at one decimal: 82.0%-92.3%."),file.path(out,"Figure_S1_QA.txt"),useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"),useBytes=TRUE)
cat("PASS: R four-quadrant matrices exported; counts, rates, exact cutoffs, sector areas and text-fit checks passed.\n")

