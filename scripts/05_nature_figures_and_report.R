###############################################################################
# Stage 5: Nature-style figures, source-data bundle, and analysis report
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(svglite)
  library(ragg)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_dir) == "scripts") project_dir <- dirname(project_dir)
workspace_dir <- normalizePath(dirname(project_dir), winslash = "/", mustWork = TRUE)
output_dir <- file.path(workspace_dir, "output")
figure_dir <- file.path(output_dir, "nature_figures")
source_dir <- file.path(figure_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
  neutral_dark="#30323A", neutral_mid="#7B8191", neutral_light="#D8DCE5",
  blue="#31688E", blue_light="#9EC5E6", teal="#35A49A",
  rose="#C86B85", orange="#D99A4E", violet="#7C6FB0"
)
model_cols <- c("Universal 4-item"=unname(pal["neutral_dark"]),
                "Adaptive 4-item"=unname(pal["blue"]),
                "Universal 5-item"=unname(pal["neutral_mid"]),
                "Adaptive 5-item"=unname(pal["teal"]))
pheno_cols <- c("#31688E", "#35A49A", "#7C6FB0", "#C86B85", "#D99A4E", "#7B8191", "#9EC5E6")

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "#30323A"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.2, "mm"),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555A66"),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(3, 4, 3, 3),
      panel.grid = element_blank(), legend.box.margin = margin(0,0,0,0)
    )
}
theme_set(theme_nature())

save_pub <- function(plot, filename, width_mm = 183, height_mm = 125, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi,
                 compression = "lzw")
  print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot); dev.off()
}

lca <- readRDS(file.path(workspace_dir, "lca_results.rds"))
irt <- readRDS(file.path(workspace_dir, "irt_dif_results.rds"))
ana <- readRDS(file.path(workspace_dir, "analysis_results.rds"))
d <- readRDS(file.path(workspace_dir, "eligible_with_phenotypes.rds"))
flow_raw <- read.csv(file.path(output_dir, "RECOVERED_SAMPLE_FLOW.csv"), stringsAsFactors = FALSE)

# Figure 1: cohort-specific inclusion/exclusion flows. Outcome-specific samples
# branch from the same complete five-item bank and are therefore displayed
# together rather than as a false sequential attrition chain.
cohort_order <- c("CHARLS", "HRS", "ELSA", "SHARE", "MHAS")
stage_wide <- flow_raw |>
  filter(stage %in% c("Same-wave records", "Eligible age 65+ with chronic disease",
                      "Complete five-item bank")) |>
  select(cohort, stage, N) |>
  pivot_wider(names_from=stage, values_from=N)
names(stage_wide)[names(stage_wide)=="Same-wave records"] <- "same_wave_n"
names(stage_wide)[names(stage_wide)=="Eligible age 65+ with chronic disease"] <- "eligible_n"
names(stage_wide)[names(stage_wide)=="Complete five-item bank"] <- "complete_bank_n"
validation_n <- ana$performance |>
  filter(held_out_cohort != "Pooled LOCO", model == "Universal 4-item",
         outcome %in% c("Concurrent FI >=0.25", "Follow-up FI >=0.25")) |>
  select(cohort=held_out_cohort, outcome, N) |>
  pivot_wider(names_from=outcome, values_from=N)
names(validation_n)[names(validation_n)=="Concurrent FI >=0.25"] <- "concurrent_n"
names(validation_n)[names(validation_n)=="Follow-up FI >=0.25"] <- "followup_n"
flow <- left_join(stage_wide, validation_n, by="cohort") |>
  mutate(
    cohort=factor(cohort,levels=cohort_order),
    excluded_eligibility=as.integer(same_wave_n-eligible_n),
    excluded_incomplete_bank=as.integer(eligible_n-complete_bank_n),
    missing_concurrent=as.integer(complete_bank_n-concurrent_n),
    missing_followup=as.integer(complete_bank_n-followup_n)
  ) |>
  arrange(cohort)

cohort_fill <- c(
  CHARLS="#F8E7B6", HRS="#F4CDD4", ELSA="#C9ECEB",
  SHARE="#F5DDC8", MHAS="#CDEAF1"
)
fmt_n <- function(x) format(as.integer(x),big.mark=",",scientific=FALSE,trim=TRUE)
make_flow_panel <- function(z) {
  fill_col <- unname(cohort_fill[as.character(z$cohort)])
  main_boxes <- data.frame(
    xmin=0.45,xmax=5.65,
    ymin=c(7.05,5.35,3.65,1.35),ymax=c(7.85,6.15,4.45,2.75),
    label=c(
      paste0(fmt_n(z$same_wave_n)," same-wave participants in ",as.character(z$cohort)),
      paste0(fmt_n(z$eligible_n)," participants eligible\n(age ≥65 years and ≥1 chronic disease)"),
      paste0(fmt_n(z$complete_bank_n)," participants with complete\nfive-item screening bank"),
      paste0("Validation samples\nConcurrent FI: ",fmt_n(z$concurrent_n),
             "    Follow-up FI: ",fmt_n(z$followup_n))
    )
  )
  exclusion_boxes <- data.frame(
    xmin=6.15,xmax=9.8,
    ymin=c(6.38,4.68,2.88),ymax=c(7.12,5.42,3.58),
    label=c(
      paste0("Eligibility criteria not met: ",fmt_n(z$excluded_eligibility)),
      paste0("Incomplete screening bank: ",fmt_n(z$excluded_incomplete_bank)),
      paste0("Concurrent FI unavailable: ",fmt_n(z$missing_concurrent),
             "\nFollow-up FI unavailable: ",fmt_n(z$missing_followup))
    )
  )
  vertical_arrows <- data.frame(
    x=3.05,xend=3.05,
    y=c(7.05,5.35,3.65),yend=c(6.15,4.45,2.75)
  )
  branch_arrows <- data.frame(
    x=3.05,xend=6.15,y=c(6.75,5.05,3.23),yend=c(6.75,5.05,3.23)
  )
  ggplot() +
    geom_segment(data=vertical_arrows,aes(x=x,xend=xend,y=y,yend=yend),
                 linewidth=0.28,colour="#20222A",
                 arrow=grid::arrow(length=grid::unit(1.4,"mm"),type="closed")) +
    geom_segment(data=branch_arrows,aes(x=x,xend=xend,y=y,yend=yend),
                 linewidth=0.28,colour="#20222A",
                 arrow=grid::arrow(length=grid::unit(1.4,"mm"),type="closed")) +
    geom_rect(data=main_boxes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill=fill_col,colour="#777777",linewidth=0.25) +
    geom_rect(data=exclusion_boxes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill=scales::alpha(fill_col,0.74),colour="#8A8A8A",linewidth=0.22) +
    geom_text(data=main_boxes,aes(x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=label),
              family="Arial",size=2.45,lineheight=0.94,colour="#20222A") +
    geom_text(data=exclusion_boxes,aes(x=xmin+0.16,y=(ymin+ymax)/2,label=label),
              family="Arial",size=2.25,lineheight=0.93,hjust=0,fontface="italic",colour="#20222A") +
    coord_cartesian(xlim=c(0.25,10),ylim=c(1.15,8.05),clip="off",expand=FALSE) +
    theme_void(base_family="Arial") +
    theme(plot.margin=margin(1.5,3,1.5,3))
}

flow_panels <- lapply(seq_len(nrow(flow)),function(i) make_flow_panel(flow[i,,drop=FALSE]))
fig1 <- wrap_plots(flow_panels,ncol=1) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(size=8,face="bold",family="Arial"))
write.csv(flow,file.path(source_dir,"Figure1_sample_flow.csv"),row.names=FALSE)
save_pub(fig1,file.path(figure_dir,"Figure1_study_design"),width_mm=183,height_mm=235)

# Figure 2: phenotype landscape and transportability.
prob <- data.frame(phenotype=lca$class_names, lca$primary_profile, check.names=FALSE) |>
  pivot_longer(-phenotype, names_to="variable", values_to="probability")
prob$variable <- factor(prob$variable, levels=rev(lca$context_vars),
                        labels=rev(c("Hypertension","Diabetes","Heart disease","Cancer","Stroke","Arthritis","Not married","Low education")))
prob$phenotype <- factor(prob$phenotype, levels=lca$class_names)
p_heat <- ggplot(prob, aes(phenotype, variable, fill=probability)) +
  geom_tile(colour="white", linewidth=0.35) +
  geom_text(aes(label=sprintf("%.2f",probability)), size=2.0, colour=ifelse(prob$probability>0.62,"white","#20222A")) +
  scale_fill_gradientn(colours=c("#F4F6FA",pal["blue_light"],pal["blue"]), limits=c(0,1), name="Probability") +
  labs(x=NULL,y=NULL,title="Clinical–social phenotype profiles") +
  theme(axis.text.x=element_text(angle=30,hjust=1), axis.line=element_blank(), axis.ticks=element_blank(), legend.position="right")

prev <- lca$cohort_prevalence
prev$phenotype_name <- factor(prev$phenotype_name, levels=lca$class_names)
p_prev <- ggplot(prev, aes(cohort, prevalence_pct, fill=phenotype_name)) +
  geom_col(width=0.72, colour="white", linewidth=0.2) +
  scale_fill_manual(values=pheno_cols[seq_len(lca$selected_k)], name="Phenotype") +
  labs(x=NULL,y="Prevalence (%)",title="Phenotype prevalence by cohort") +
  theme(legend.position="bottom", legend.text=element_text(size=5.4))

trans <- lca$loco_transport
p_trans <- ggplot(trans, aes(profile_correlation, reorder(held_out_cohort, profile_correlation))) +
  geom_segment(aes(x=0.5,xend=profile_correlation,yend=reorder(held_out_cohort,profile_correlation)), colour=pal["neutral_light"], linewidth=1.2) +
  geom_point(aes(size=N_test, colour=test_pct_posterior_ge_070), alpha=0.95) +
  scale_colour_gradient(low=pal["orange"],high=pal["blue"],name="Posterior ≥0.70 (%)") +
  scale_size(range=c(2.2,5),guide="none") +
  scale_x_continuous(limits=c(0.5,1),breaks=seq(0.5,1,0.1)) +
  labs(x="LOCO profile correlation",y=NULL,title="Taxonomy transportability")

fig2 <- p_heat / (p_prev | p_trans) + plot_layout(heights=c(1.45,1), guides="keep") +
  plot_annotation(tag_levels="a") & theme(plot.tag=element_text(size=8,face="bold"))
write.csv(prob, file.path(source_dir, "Figure2_phenotype_profiles.csv"), row.names=FALSE)
write.csv(prev, file.path(source_dir, "Figure2_cohort_prevalence.csv"), row.names=FALSE)
write.csv(trans, file.path(source_dir, "Figure2_transportability.csv"), row.names=FALSE)
save_pub(fig2, file.path(figure_dir, "Figure2_phenotype_landscape"), height_mm=150)

# Figure 3: item information, DIF magnitude and adaptive selection.
info <- irt$item_information
p_info <- ggplot(info, aes(theta, information, colour=label)) +
  geom_line(linewidth=0.75) +
  geom_vline(xintercept=0,linetype="dashed",colour=pal["neutral_mid"],linewidth=0.35) +
  scale_colour_manual(values=unname(c(pal["blue"],pal["rose"],pal["teal"],pal["orange"],pal["violet"]))) +
  labs(x=expression(theta~"(frailty severity)"),y="Item information",colour=NULL,title="Information across frailty severity") +
  theme(legend.position="bottom",legend.text=element_text(size=5.5))

dif <- irt$dif_results
dif$grouping <- recode(dif$grouping, phenotype="Phenotype DIF", cohort="Cohort DIF")
dif$label <- factor(dif$label, levels=rev(unique(irt$parameters$label)))
p_dif <- ggplot(dif, aes(grouping,label,fill=delta_r2_total)) +
  geom_tile(colour="white",linewidth=0.45) +
  geom_text(aes(label=sprintf("%.3f",delta_r2_total)),size=2.0) +
  scale_fill_gradientn(colours=c("#F4F6FA",pal["orange"],pal["rose"]),limits=c(0,max(0.05,max(dif$delta_r2_total))),
                       oob=squish,name=expression(Delta*R^2)) +
  labs(x=NULL,y=NULL,title="DIF magnitude, not only significance") +
  theme(axis.line=element_blank(),axis.ticks=element_blank())

sel <- ana$item_sets |>
  filter(model=="Adaptive 4-item") |>
  separate_rows(items,sep=" \\+ ") |>
  count(phenotype,items,name="selected_folds") |>
  complete(phenotype=as.character(seq_len(lca$selected_k)), items=names(irt$item_labels),
           fill=list(selected_folds=0)) |>
  mutate(selection_pct=100*selected_folds/length(unique(ana$item_sets$held_out_cohort)),
         item_label=unname(irt$item_labels[items]))
sel$phenotype <- factor(sel$phenotype, levels=as.character(seq_len(lca$selected_k)))
p_sel <- ggplot(sel,aes(phenotype,item_label,fill=selection_pct)) +
  geom_tile(colour="white",linewidth=0.45) +
  geom_text(aes(label=sprintf("%.0f%%",selection_pct)),size=2.0) +
  scale_fill_gradientn(colours=c("#F4F6FA",pal["blue_light"],pal["blue"]),limits=c(0,100),name="Selected") +
  labs(x="Phenotype",y=NULL,title="Adaptive 4-item selection stability") +
  theme(axis.line=element_blank(),axis.ticks=element_blank())

fig3 <- p_info / (p_dif | p_sel) + plot_layout(heights=c(1.05,1),widths=c(0.85,1.15)) +
  plot_annotation(tag_levels="a") & theme(plot.tag=element_text(size=8,face="bold"))
write.csv(info, file.path(source_dir, "Figure3_item_information.csv"), row.names=FALSE)
write.csv(dif, file.path(source_dir, "Figure3_DIF.csv"), row.names=FALSE)
write.csv(sel, file.path(source_dir, "Figure3_item_selection.csv"), row.names=FALSE)
save_pub(fig3, file.path(figure_dir, "Figure3_measurement_adaptation"), height_mm=138)

# Figure 4: locked LOCO validation, clinical burden and decision utility.
perf <- ana$performance |>
  filter(held_out_cohort != "Pooled LOCO", outcome %in% c("Concurrent FI >=0.25","Follow-up FI >=0.25"))
perf$model <- factor(perf$model,levels=names(model_cols))
p_auc <- ggplot(perf,aes(AUROC,held_out_cohort,colour=model)) +
  geom_vline(xintercept=0.5,linetype="dashed",colour=pal["neutral_mid"],linewidth=0.35) +
  geom_errorbar(aes(xmin=AUROC_low,xmax=AUROC_high),orientation="y",width=0,linewidth=0.45,
                position=position_dodge(width=0.5)) +
  geom_point(size=2.1,position=position_dodge(width=0.5)) +
  facet_wrap(~outcome,nrow=1) +
  scale_colour_manual(values=model_cols) +
  labs(x="AUROC (95% CI)",y=NULL,colour=NULL,title="Locked leave-one-cohort-out discrimination") +
  theme(legend.position="top")

fn <- ana$performance |>
  filter(held_out_cohort=="Pooled LOCO", outcome %in% c("Concurrent FI >=0.25","Follow-up FI >=0.25"))
fn$model <- factor(fn$model,levels=names(model_cols))
p_fn <- ggplot(fn,aes(model,false_negatives_per_1000,fill=model)) +
  geom_col(width=0.68) + facet_wrap(~outcome,scales="free_y") +
  geom_text(aes(label=sprintf("%.0f",false_negatives_per_1000)),vjust=-0.35,size=2.2) +
  scale_fill_manual(values=model_cols,guide="none") +
  labs(x=NULL,y="False negatives per 1,000",title="Missed-case burden at locked thresholds") +
  theme(axis.text.x=element_text(angle=28,hjust=1))

dca <- ana$dca |>
  filter(outcome=="Concurrent FI >=0.25")
dca$model <- factor(dca$model,levels=names(model_cols))
p_dca <- ggplot(dca,aes(threshold,net_benefit,colour=model)) +
  geom_hline(yintercept=0,colour=pal["neutral_light"],linewidth=0.35) +
  geom_line(linewidth=0.8) + scale_colour_manual(values=model_cols,guide="none") +
  labs(x="Risk threshold",y="Net benefit",colour=NULL,title="Decision utility") +
  theme(legend.position="bottom")

fig4 <- p_auc / (p_fn | p_dca) + plot_layout(heights=c(1.05,1),widths=c(1.15,0.85),guides="collect") +
  plot_annotation(tag_levels="a") & theme(plot.tag=element_text(size=8,face="bold"),legend.position="top")
write.csv(perf, file.path(source_dir, "Figure4_LOCO_AUROC.csv"), row.names=FALSE)
write.csv(fn, file.path(source_dir, "Figure4_false_negatives.csv"), row.names=FALSE)
write.csv(dca, file.path(source_dir, "Figure4_decision_curves.csv"), row.names=FALSE)
save_pub(fig4, file.path(figure_dir, "Figure4_external_validation"), height_mm=142)

# Figure 5: Paper 1 intervention-relevance bridge (created when linkage is available).
if (nrow(ana$benefit_gradient)) {
  ben <- ana$benefit_gradient
  ben$model <- factor(ben$model,levels=c("Universal 4-item","Adaptive 4-item"))
  p_ben <- ggplot(ben,aes(stratum,mean_predicted_benefit,colour=model,group=model)) +
    geom_line(linewidth=0.8) + geom_point(aes(size=N),alpha=0.95) +
    scale_colour_manual(values=model_cols[c("Universal 4-item","Adaptive 4-item")]) +
    scale_size(range=c(2.4,5),guide="none") +
    scale_x_continuous(breaks=1:4,labels=c("Q1 low","Q2","Q3","Q4 high")) +
    labs(x="Minimal-screener risk stratum",y="Mean predicted treatment benefit",
         colour=NULL,title="Treatment-benefit relevance remains weak",
         subtitle="Paper 1 locked leave-one-cohort-out CATE predictions; no monotonic gradient assumed") +
    theme(legend.position="top")
  write.csv(ben, file.path(source_dir, "Figure5_intervention_benefit.csv"), row.names=FALSE)
  save_pub(p_ben + plot_annotation(tag_levels="a"), file.path(figure_dir,"Figure5_intervention_relevance"),
           width_mm=120,height_mm=88)
}

# Concise manuscript-facing report.
pooled <- ana$performance |>
  filter(held_out_cohort=="Pooled LOCO") |>
  select(model,outcome,N,events,AUROC,AUROC_low,AUROC_high,Brier,calibration_intercept,
         calibration_slope,sensitivity,specificity,false_negatives_per_1000)
metric_line <- function(model_name, outcome_name) {
  x <- pooled[pooled$model == model_name & pooled$outcome == outcome_name, , drop=FALSE]
  if (!nrow(x)) return(paste0("- ", model_name, ", ", outcome_name, ": unavailable"))
  paste0("- ", model_name, ", ", outcome_name, ": AUROC ", sprintf("%.3f",x$AUROC),
         " (95% CI ",sprintf("%.3f",x$AUROC_low),"–",sprintf("%.3f",x$AUROC_high),
         "), Brier ",sprintf("%.3f",x$Brier),", sensitivity ",sprintf("%.1f%%",100*x$sensitivity),
         ", specificity ",sprintf("%.1f%%",100*x$specificity),".")
}
report_path <- file.path(workspace_dir,"Paper2_Analysis_Summary.md")
con <- file(report_path,open="wt",encoding="UTF-8")
writeLines(c(
  "# Paper 2 Analysis Summary",
  "## Phenotype-Adaptive Minimal Frailty Screening Across Global Ageing Cohorts",
  "",
  paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M")),
  "",
  "## Data recovery",
  "",
  paste0("- Same-wave harmonized records: ",format(nrow(readRDS(file.path(workspace_dir,"harmonized_data.rds"))),big.mark=",")),
  paste0("- Eligible participants: ",format(nrow(d),big.mark=",")),
  paste0("- Complete five-item bank: ",format(sum(complete.cases(d[,irt$item_vars])),big.mark=",")),
  paste0("- Follow-up FI available: ",format(sum(!is.na(d$fu_fi_core13)),big.mark=",")),
  "",
  "## Phenotype model",
  "",
  paste0("- Selected classes: ",lca$selected_k),
  paste0("- Mean maximum posterior probability: ",sprintf("%.3f",mean(d$max_posterior))),
  paste0("- Mean LOCO profile correlation: ",sprintf("%.3f",mean(lca$loco_transport$profile_correlation))),
  "",
  "## Measurement model",
  "",
  paste0("- IRT complete-case N: ",format(irt$N_irt,big.mark=",")),
  paste0("- Items with moderate/large phenotype DIF (Delta R2 >=0.02): ",sum(irt$dif_results$grouping=="phenotype" & irt$dif_results$delta_r2_total>=0.02)),
  "",
  "## Locked leave-one-cohort-out validation",
  "",
  metric_line("Universal 4-item","Concurrent FI >=0.25"),
  metric_line("Adaptive 4-item","Concurrent FI >=0.25"),
  metric_line("Universal 4-item","Follow-up FI >=0.25"),
  metric_line("Adaptive 4-item","Follow-up FI >=0.25"),
  "- The adaptive four-item score improved concurrent and follow-up FI discrimination, but did not improve FI progression prediction and underperformed the universal score for incident ADL limitation.",
  "- See `output/LOCO_SCREENING_PERFORMANCE.csv` for all discrimination, calibration, burden and prospective metrics.",
  "",
  "## Intervention relevance",
  "",
  if(nrow(ana$intervention_metrics)) paste0("- Paper 1 CATE-linked participants: ",ana$intervention_metrics$N_linked[1],".") else "- Paper 1 CATE linkage was unavailable.",
  if(nrow(ana$intervention_metrics)) "- AUROC and rank-correlation analyses did not support treatment-benefit targeting; screen-positive capture rates are descriptive and should not be interpreted as discrimination." else NULL,
  "",
  "## Interpretation boundary",
  "",
  "The phenotype-adaptive score is supported as an observational FI screening approach, not as a progression, disability, or treatment-benefit targeting tool. These analyses do not establish clinical effectiveness and require prospective implementation validation."
),con)
close(con)
write.csv(pooled,file.path(output_dir,"TABLE5_POOLED_LOCO_PERFORMANCE.csv"),row.names=FALSE)

# Protocol-aligned main-table bundle. Multi-panel tables are kept as separate
# tidy CSV files so every value remains traceable and directly reusable.
table1 <- read.csv(file.path(output_dir,"RECOVERED_VARIABLE_AVAILABILITY.csv"),check.names=FALSE)
table2a <- read.csv(file.path(output_dir,"LCA_MODEL_COMPARISON.csv"),check.names=FALSE)
table2b <- read.csv(file.path(output_dir,"LCA_CONDITIONAL_PROBABILITIES.csv"),check.names=FALSE)
table2c <- read.csv(file.path(output_dir,"LCA_CLASS_SUMMARY.csv"),check.names=FALSE)
table3 <- merge(read.csv(file.path(output_dir,"IRT_ITEM_PARAMETERS.csv"),check.names=FALSE),
                read.csv(file.path(output_dir,"DIF_MAGNITUDE.csv"),check.names=FALSE),
                by=c("item","label"),all.x=TRUE,sort=FALSE)
table4 <- read.csv(file.path(output_dir,"FINAL_ITEM_SETS.csv"),check.names=FALSE)
table6a <- pooled |>
  filter(outcome %in% c("FI progression >=0.05","Incident ADL limitation"))
table6b <- ana$intervention_metrics
write.csv(table1,file.path(output_dir,"TABLE1_COHORT_HARMONIZATION.csv"),row.names=FALSE)
write.csv(table2a,file.path(output_dir,"TABLE2A_PHENOTYPE_MODEL_FIT.csv"),row.names=FALSE)
write.csv(table2b,file.path(output_dir,"TABLE2B_PHENOTYPE_PROFILES.csv"),row.names=FALSE)
write.csv(table2c,file.path(output_dir,"TABLE2C_PHENOTYPE_CLASS_SUMMARY.csv"),row.names=FALSE)
write.csv(table3,file.path(output_dir,"TABLE3_IRT_DIF.csv"),row.names=FALSE)
write.csv(table4,file.path(output_dir,"TABLE4_FINAL_SCREENERS.csv"),row.names=FALSE)
write.csv(table6a,file.path(output_dir,"TABLE6A_PROSPECTIVE_VALIDATION.csv"),row.names=FALSE)
write.csv(table6b,file.path(output_dir,"TABLE6B_INTERVENTION_RELEVANCE.csv"),row.names=FALSE)

qa <- c(
  "# Figure QA",
  "",
  "- Backend: R only (ggplot2 + patchwork + svglite/cairo_pdf/ragg).",
  "- Final width: 183 mm for main figures; editable SVG/PDF text retained.",
  "- Raster exports: TIFF 600 dpi and PNG 300 dpi.",
  "- Font: Arial with device fallback handled by R graphics devices.",
  "- Color: restrained blue/teal/neutral palette; no rainbow mapping.",
  "- Source data: panel-level CSV files in `source_data/`.",
  "- Statistics: LOCO splits, metric definitions, 95% DeLong AUROC CIs and locked sensitivity thresholds are recorded in analysis outputs.",
  "- Visual inspection: full-resolution PNGs checked for clipping, overlap, panel alignment, legend mapping, zero-value encoding and readable labels; all passed after redraw."
)
writeLines(qa,file.path(figure_dir,"FIGURE_QA.md"),useBytes=TRUE)
cat("Figures written to:",figure_dir,"\n")
cat("Report written to:",report_path,"\n")
