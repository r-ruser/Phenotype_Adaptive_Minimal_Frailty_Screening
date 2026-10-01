# Reader-facing scoring function. Inputs are already coded 0-1 deficits.
# Required columns: q_lifta, q_climsa, q_shlt, q_chaira, q_stoopa, q_shopa.
score_FI6 <- function(data) {
 items <- c("q_lifta","q_climsa","q_shlt","q_chaira","q_stoopa","q_shopa")
 stopifnot(all(items %in% names(data)))
 x <- as.matrix(data[,items,drop=FALSE]);storage.mode(x)<-"double"
 stopifnot(all(x >= 0 & x <= 1,na.rm=TRUE))
 n <- rowSums(!is.na(x));score <- rowMeans(x,na.rm=TRUE)
 score[n < 5] <- NA_real_
 data.frame(FI6=score,observed_items=n)
}
classify_FI6 <- function(FI6,cohort) {
 stopifnot(all(cohort%in%c("HRS","ELSA","SHARE","CHARLS","MHAS")))
 cutoff <- ifelse(cohort%in%c("ELSA","SHARE"),5/12,11/24)
 data.frame(locked_cutoff=cutoff,positive=FI6>=cutoff)
}
