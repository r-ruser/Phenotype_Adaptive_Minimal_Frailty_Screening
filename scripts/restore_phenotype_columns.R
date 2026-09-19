# Utility used only when Stage 1 is re-run after a completed Stage 2 model.
options(stringsAsFactors = FALSE)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_dir) == "scripts") project_dir <- dirname(project_dir)
workspace_dir <- normalizePath(dirname(project_dir), winslash = "/", mustWork = TRUE)
new <- readRDS(file.path(workspace_dir, "eligible_with_bmi.rds"))
old <- readRDS(file.path(workspace_dir, "eligible_with_phenotypes.tmp.rds"))
keep <- grep("^phenotype|^max_posterior$", names(old), value = TRUE)
idx <- match(new$uid, old$uid)
stopifnot(!anyNA(idx))
for (v in keep) new[[v]] <- old[[v]][idx]
saveRDS(new, file.path(workspace_dir, "eligible_with_phenotypes.rds"))
cat("Restored", length(keep), "phenotype columns into", nrow(new), "rows.\n")
