# Source files and participant-level outputs remain in user-controlled locations.
data_root <- Sys.getenv("MS5_DATA_ROOT", unset="")
runtime_dir <- Sys.getenv("MS5_RUNTIME_DIR", unset="")
if (!nzchar(runtime_dir)) stop("Set MS5_RUNTIME_DIR to a private output directory before running analysis")
runtime_dir <- normalizePath(runtime_dir, winslash="/", mustWork=FALSE)
base <- file.path(runtime_dir,"analysis")
out <- file.path(base,"results")
score_source <- out
dir.create(out,recursive=TRUE,showWarnings=FALSE)
eligible_input <- file.path(runtime_dir,"eligible_baseline.rds")
source_file <- function(relative) { if(!nzchar(data_root))stop("Set MS5_DATA_ROOT to the parent folder containing the five cohort folders"); file.path(data_root,relative) }
