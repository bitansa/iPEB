# Shared data generator and evaluation helpers for the iPEB simulation
# reproducibility scripts. The data-generating process is identical to the one
# used in the paper, so the scenario scripts reproduce the reported numbers. All
# iPEB and PEB scoring is done through the installed package (library(iPEB)); the
# helpers below are only the simulation scaffolding (data generation, ROC/AUC,
# and the per-patient sensitivity / lead-time metrics).

suppressMessages({ library(iPEB); library(parallel) })
RNGkind("L'Ecuyer-CMRG")                          # reproducible parallel RNG
.NCORES <- max(1L, parallel::detectCores() - 1L)

# Parallel lapply where available (mclapply on Unix/macOS; serial on Windows).
# Each cohort self-seeds, so results are identical serial or parallel.
par_lapply <- function(X, FUN) {
  if (.Platform$OS.type == "windows") lapply(X, FUN)
  else parallel::mclapply(X, FUN, mc.cores = .NCORES)
}

## ---- data generator --------------------------------------------------------
## Long data frame with columns: IDvar, D, t, draw_to_dx_days, M1, M2, composite.
## Controls have random baselines (a fraction constitutively high on M2); cases
## are heterogeneous (a fraction driven by an M1 rise, the rest by M2), each over
## its own lead window. `composite` is the fixed comparator panel scored through
## PEB: a fitted-then-frozen logistic combination of the markers (or a supplied
## external panel), the simulation analogue of a deployed fixed panel. (The
## "JAMA 4MP" name in the paper refers only to the published four-marker panel in
## the real lung analysis, not to this simulated logistic composite.)
simulate_cohort <- function(n_ctrl = 600, n_case = 200, K = 8, dmax = 1080, dmin = 90,
                            mu = 10, tau = 2, sd_noise = 1,
                            trend_sd = 0, ar_phi = 0,
                            hi_base_frac = 0.20, hi_base_shift = 8,
                            pA = 0.60, delta1 = 6, delta2 = 6,
                            leadA = 365, leadB = 365,
                            both_markers = FALSE, panel_weights = NULL,
                            dev_panel = FALSE, seed = 1) {
  set.seed(seed)
  ramp <- function(d2dx, lead) pmax(0, 1 - d2dx / lead)
  d2 <- round(seq(dmax, dmin, length.out = K))
  yr <- (max(d2) - d2) / 365.25
  ar_noise <- function() {
    if (ar_phi <= 0) return(rnorm(K, 0, sd_noise))
    e <- rnorm(K, 0, sd_noise * sqrt(1 - ar_phi^2)); x <- numeric(K); x[1] <- rnorm(1, 0, sd_noise)
    for (k in 2:K) x[k] <- ar_phi * x[k - 1] + e[k]; x
  }
  mk <- function(id, D) {
    b1 <- rnorm(1, mu, tau); b2 <- rnorm(1, mu, tau)
    if (runif(1) < hi_base_frac) b2 <- b2 + hi_base_shift
    s1 <- rnorm(1, 0, trend_sd); s2 <- rnorm(1, 0, trend_sd)
    M1 <- b1 + s1 * yr + ar_noise(); M2 <- b2 + s2 * yr + ar_noise()
    if (D == 1) {
      if (both_markers) {
        M1 <- M1 + delta1 * ramp(d2, leadA); M2 <- M2 + delta2 * ramp(d2, leadB)
      } else if (runif(1) < pA) M1 <- M1 + delta1 * ramp(d2, leadA)
      else                      M2 <- M2 + delta2 * ramp(d2, leadB)
    }
    data.frame(IDvar = id, D = D, draw_to_dx_days = d2, t = -d2, M1 = M1, M2 = M2)
  }
  df <- do.call(rbind, c(lapply(seq_len(n_ctrl),        function(i) mk(sprintf("c%04d", i), 0)),
                         lapply(seq_len(n_case) + n_ctrl, function(i) mk(sprintf("k%04d", i), 1))))
  df$obs_number <- ave(seq_len(nrow(df)), df$IDvar, FUN = seq_along)
  if (dev_panel) {
    dev <- simulate_cohort(n_ctrl = n_ctrl, n_case = n_case, K = K, dmax = dmax, dmin = dmin,
                           mu = mu, tau = tau, sd_noise = sd_noise, hi_base_frac = hi_base_frac,
                           hi_base_shift = hi_base_shift, both_markers = both_markers,
                           delta1 = 0, delta2 = delta2, leadA = leadA, leadB = leadB,
                           panel_weights = NULL, dev_panel = FALSE, seed = seed + 100000L)
    g <- glm(D ~ M1 + M2, family = binomial, data = dev[dev$obs_number == K, ])
    df$composite <- as.numeric(predict(g, newdata = df))
  } else if (is.null(panel_weights)) {
    g <- glm(D ~ M1 + M2, family = binomial, data = df[df$obs_number == K, ])
    df$composite <- as.numeric(predict(g, newdata = df))
  } else {
    df$composite <- panel_weights[1] * df$M1 + panel_weights[2] * df$M2
  }
  df
}

## ---- evaluation helpers ----------------------------------------------------
auc_trap <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  if (x[1] > 0) { x <- c(0, x); y <- c(0, y) }
  if (utils::tail(x, 1) < 1) { x <- c(x, 1); y <- c(y, 1) }
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}
roc_curve <- function(score, D, ids, g = seq(0, 1, 0.02)) {
  pmx <- tapply(score, ids, max, na.rm = TRUE); pl <- tapply(D, ids, max)
  case <- as.numeric(pmx[names(pl)[pl == 1]]); ctrl <- score[D == 0]
  thr <- quantile(ctrl, 1 - g, na.rm = TRUE, type = 8)
  sapply(thr, function(c) mean(case > c, na.rm = TRUE))
}
roc_auc <- function(score, D, ids) { g <- seq(0, 1, 0.02); auc_trap(g, roc_curve(score, D, ids, g)) }

## Per-patient sensitivity and median lead (whole trajectory), with a per-visit
## threshold calibrated on the supplied control visits.
metrics_subj <- function(score, D, ids, ttd, ctrl_score, alpha) {
  thr <- quantile(ctrl_score, alpha, na.rm = TRUE, type = 8)
  cids <- unique(ids[D == 1])
  lead <- vapply(cids, function(id) {
    v <- which(ids == id & ttd > 0 & score > thr)
    if (length(v) == 0) NA_real_ else max(ttd[v])
  }, numeric(1))
  c(SensW = mean(!is.na(lead)),
    LT = if (all(is.na(lead))) 0 else median(lead, na.rm = TRUE) / 365.25,
    Spec = mean(score[D == 0] <= thr, na.rm = TRUE))
}

## Evaluate any score vector on the test rows: AUC + Sens_W + lead time.
## ctrl_score = the score on TRAIN control visits (per-visit calibration).
eval_score <- function(score, D, ids, ttd, ctrl_score, alpha) {
  e <- metrics_subj(score, D, ids, ttd, ctrl_score, alpha)
  c(AUC = roc_auc(score, D, ids), SensW = unname(e["SensW"]), LT = unname(e["LT"]))
}

## Split subjects into train/test (stratified by case status).
split_ids <- function(df, train_frac = 0.7, seed = 1) {
  set.seed(seed)
  isC <- tapply(df$D, df$IDvar, max)
  cs <- names(isC)[isC == 1]; ct <- names(isC)[isC == 0]
  c(sample(cs, round(train_frac * length(cs))), sample(ct, round(train_frac * length(ct))))
}
