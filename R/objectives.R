# Clinical objectives: fixed penalty profiles, evaluation metrics, and the
# scalarized loss used to learn the combiner. All internal.

# Fixed (lambda1, lambda2, lambda3) profile for each objective. Because the
# threshold is recalibrated to the target specificity for every candidate weight
# vector, the specificity-floor term (lambda1) is essentially inactive and the
# referral term (lambda2 * Ref, with Ref = 1 - Spec) is nearly constant across
# candidates; lambda3 is therefore the operative lead-time lever.
.lambda_profile <- function(objective) {
  switch(objective,
    sens     = c(lam1 = 50, lam2 = 0, lam3 = 0.0),
    leadtime = c(lam1 = 50, lam2 = 0, lam3 = 0.7),
    combined = c(lam1 = 50, lam2 = 1, lam3 = 0.5),
    stop("Unknown objective: ", objective))
}

# Metrics after the PEB decision at operating specificity `alpha`.
# Sensitivity and lead time are per patient (a case is detected if ANY of its
# pre-diagnostic visits crosses the threshold); specificity is per visit. When
# `window` is finite (in months) detection is restricted to that pre-diagnostic
# horizon; the default `window = Inf` scores over the whole trajectory.
# `ctrl_score` must be the control-VISIT scores, so the threshold is a per-visit
# quantile.
.metrics <- function(score, D, ids, ttd, ctrl_score, alpha, window = Inf) {
  thr <- stats::quantile(ctrl_score, alpha, na.rm = TRUE, type = 8)
  wdays <- if (is.finite(window)) window * 30.4375 else Inf
  cids <- unique(ids[D == 1])
  lead <- vapply(cids, function(id) {
    v <- which(ids == id & ttd > 0 & ttd <= wdays & score > thr)
    if (length(v) == 0L) NA_real_ else max(ttd[v])
  }, numeric(1))
  SensW <- mean(!is.na(lead))
  LT <- if (all(is.na(lead))) 0 else stats::median(lead, na.rm = TRUE) / 365.25
  Spec <- mean(score[D == 0] <= thr, na.rm = TRUE)
  c(SensW = SensW, LT = LT, Spec = Spec, Ref = 1 - Spec)
}

# Scalarized loss to MINIMIZE (manuscript form). lam3 > 0 rewards lead time.
.loss <- function(score, D, ids, ttd, ctrl_score, alpha, window = Inf,
                  lam1 = 50, lam2 = 0, lam3 = 0, LT_ref = 2) {
  mt <- .metrics(score, D, ids, ttd, ctrl_score, alpha, window)
  as.numeric(-mt["SensW"] + lam1 * max(alpha - mt["Spec"], 0) +
               lam2 * mt["Ref"] - lam3 * (mt["LT"] / LT_ref))
}

# Rank-based AUC (Mann-Whitney); returns NA if a class is absent.
.auc <- function(s, y) {
  ok <- is.finite(s) & !is.na(y)
  s <- s[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Subject-level AUC using each subject's maximum score across visits.
.auc_subject <- function(score, D, ids) {
  ss <- tapply(score, ids, function(z) max(z, na.rm = TRUE))
  yy <- tapply(D, ids, function(z) max(z, na.rm = TRUE))
  .auc(as.numeric(ss), as.numeric(yy))
}
