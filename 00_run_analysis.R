if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# FULL ANALYSIS AND FIGURE ENTRY POINT (R files only).
# Recorded training runtime: R 4.4.3, xgboost 1.7.11.1; current versions are
# written to private run_logs/software_versions.csv and sessionInfo.txt.
# Required packages: haven, dplyr, xgboost, ggplot2, patchwork, metafor, ragg,
# svglite, systemfonts, ggbeeswarm, xml2, png, pROC. Graphics use Arial.
# Binary deficits and graded items use 0-1 coding; see 02_prepare_complete_FI.R
# and 03_candidate_bank.R for valid response codes, missingness and completeness.
# Complete FI mappings are generated from baseline_metadata() in 00_config.R.
# No external CSV mapping or bundled analysis result is needed to start from
# licensed inputs. All intermediate objects and result files are generated locally.
# Required environment variables:
# MS5_DATA_ROOT: directory holding HRS_USA, ELSA_UK, SHARE_Europe,
# CHARLS_China, MHAS_Mexico (official licensed harmonized source releases).
# MS5_RUNTIME_DIR: private generated-data/results directory, outside this repository.
# Usage: Rscript 00_run_analysis.R plan | all | figures | <exact stage basename>
# plan writes the run order, outputs and current installed software versions locally.
# all refits the primary and phenotype models plus 500 coupled training bootstraps.
# Expected reference population: age >=65; >=1 of six diagnoses; N=34,954.
# Fixed-wave prefixes: HRS14/15, ELSA9/10, SHARE2/4, CHARLS3/4, MHAS5/6.
# Seeds/iteration counts remain in their analysis scripts. The six-item length stays
# exploratory and fixed in the 500-repeat bootstrap; phenotype uncertainty is separate.
# Outputs below are generated locally and must not be committed to the R-only release.
options(stringsAsFactors=FALSE)
script_file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
package_dir<-dirname(normalizePath(script_file,winslash="/",mustWork=TRUE))
source(file.path(package_dir,"00_config.R"),encoding="UTF-8")

stages<-data.frame(script=c(
"01_recover_baseline.R","02_prepare_complete_FI.R","03_candidate_bank.R","04_training_phenotypes_SHAP.R","05_universal_SHAP.R","06_locked_scores.R","07_performance.R","08_retention.R","09_paired_comparisons.R","10_overlap_sensitivity.R","11_continuous_agreement.R","12_six_item_benchmarks.R","14_SHARE_country_source.R","15_attrition_IPW.R","16_shared_five_ADL_sensitivity.R","18_classification_results.R","19_item_stability.R","20_pooled_curves.R","21_length_errors_and_sensitivity.R","22_heterogeneity_analysis.R","23_training_bootstrap_500.R","24_summarize_training_bootstrap.R","25_heterogeneity_sources.R","26_FI_domain_sources.R","27_functional_severity.R","draw_Figure_1_framework.R","draw_Figure_2_length_tradeoffs.R","draw_Figure_3_continuous_agreement.R","draw_Figure_4_benchmarks.R","draw_Figure_5_prospective.R","17_draw_confusion_matrix.R","draw_Figure_S2_SHAP.R","draw_Figure_S3_phenotype_comparison.R","draw_Figures_S5_S6_weighted_forest.R"),
output=c("eligible_baseline.rds","analysis/prepared_data.rds; reference_item_mapping.csv; reference_summary.csv","analysis/results/candidate_question_bank.rds","FI5_SHAP_holdout_*.rds; FI5_SHAP_*.csv","FI_universal_SHAP_rankings.csv; FI_universal_SHAP_beeswarm_sample.csv","FI_reduction_locked_scores.rds; FI_reduction_training_choices.csv","FI_item_number_performance.csv; FI_item_number_complete_case.csv; FI_reduction_risk_sets.csv","FI_optimal_item_number*.csv; retention results","FI_universal_vs_phenotype.csv","JGMS_reference_overlap_sensitivity.csv","extensions/approximation/*.csv (Table6/7)","extensions/benchmarks/*.csv (Table8)","extensions/attrition/SHARE_country_source.rds","extensions/attrition/*.csv (Tables9-11)","extensions/attrition/shared_ADL*.csv (Table12)","FI6_locked_threshold_classification.csv (S1)","FI_item_stability*.csv","FI_reduction_meta_summary.csv; publication/figures/Figure_2_pooled_source.csv","analysis/results/length_curve/*.csv; Figure_S4* (Tables13-17)","extensions/heterogeneity/*.csv (Table18)","extensions/training_bootstrap_500/checkpoints/rep_001.rds through rep_500.rds (private; 7500 fits)","extensions/training_bootstrap_500/*.csv; Figure_S7* (Table19/20)","extensions/training_bootstrap_500/sources/*.csv; Figure_S8* (Table21)","extensions/training_bootstrap_500/sources/FI_domain_distribution.csv; MHAS_score_standardized_domain_gap.csv","extensions/training_bootstrap_500/sources/functional_limitation*.csv","publication/figures/Figure_1*","publication/figures/Figure_2*","publication/figures/Figure_3*","publication/figures/Figure_4*","publication/figures/Figure_5*","extensions/confusion_matrix/* (S1)","analysis/results/Figure_S2_SHAP*","publication/figures/Figure_S3*","extensions/heterogeneity/forest_points/Figure_S5* and Figure_S6*"))
stages$order<-seq_len(nrow(stages))
logs<-file.path(runtime_dir,"run_logs");dir.create(logs,showWarnings=FALSE,recursive=TRUE)
write.csv(stages,file.path(logs,"analysis_stage_output_map.csv"),row.names=FALSE)
packages<-c("haven","dplyr","xgboost","ggplot2","patchwork","metafor","ragg","svglite","systemfonts","ggbeeswarm","xml2","png","pROC")
versions<-data.frame(package=c("R",packages),version=c(as.character(getRversion()),vapply(packages,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else "MISSING",character(1))))
write.csv(versions,file.path(logs,"software_versions.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(logs,"sessionInfo.txt"))
args<-commandArgs(TRUE);mode<-if(length(args))args[1] else "plan"
if(mode=="plan"){print(stages,row.names=FALSE);quit(status=0)}
todo<-if(mode=="all")stages$script else if(mode=="figures")stages$script[grepl("^draw_|^17_",stages$script)] else mode
stopifnot(all(todo %in% stages$script))
Sys.setenv(MS5_CODE_DIR=package_dir,MS5_PACKAGE_DIR=package_dir)
for(f in todo){
 cat("RUN",f,"\n");flush.console()
 log<-file.path(logs,paste0(f,".log"))
 status<-system2(file.path(R.home("bin"),"Rscript"),shQuote(file.path(package_dir,f)),stdout=log,stderr=log)
 if(status!=0L)stop("Stage failed: ",f,"; see ",log)
}
