# Scenario 4 (irregular visit spacing): the value of time-gap awareness.
#
# Companion to Scenario 3, but on an IRREGULAR visit grid. Healthy subjects carry
# a random slope and continuous-time (Ornstein-Uhlenbeck) serial correlation, and
# visits are spaced unevenly (gaps from a few months up to about two years). A
# gap-ignorant rule mistakes the drift and the long gaps for signal and detects
# few cases, whereas the time-gap-aware iPEB layer scales the innovation variance
# by the elapsed gap and detects nearly all of them. Cases and controls share the
# same visit-spacing distribution, so spacing alone is uninformative about status.
#
# Outputs: the PEB-vs-iPEB row (AUC / Sensitivity / Lead time) for the main-text
# mechanism table, and one 300-dpi sensitivity figure for the Web supplement.
#
# Run:  source("reproduce_S4.R")

source("sim_helpers.R")

## ---- scenario configuration ------------------------------------------------
gen <- list(n_ctrl = 600, n_case = 200, K = 8, dmin = 90,
            gap_pool = c(90, 180, 365, 730),          # 3 mo, 6 mo, 1 yr, 2 yr
            mu = 10, tau = 2, sd_noise = 1, trend_sd = 1.5, ar_phi = 0.6,
            hi_base_frac = 0.20, hi_base_shift = 8,
            pA = 0.60, delta1 = 6, delta2 = 6, leadA = 365, leadB = 365)
R    <- 200                                           # Monte-Carlo cohorts
SPEC <- 0.95                                          # target per-visit specificity

## ---- irregular-spacing data generator --------------------------------------
## Two-marker structure matching Scenarios 1-3, but each subject is measured on
## its own irregular schedule and residuals follow a continuous-time OU process,
## so the elapsed gap between visits varies within and across subjects.
simulate_cohort_irregular <- function(n_ctrl, n_case, K, dmin, gap_pool,
                                      mu, tau, sd_noise, trend_sd, ar_phi,
                                      hi_base_frac, hi_base_shift,
                                      pA, delta1, delta2, leadA, leadB, seed = 1) {
  set.seed(seed)
  ramp <- function(d2dx, lead) pmax(0, 1 - d2dx / lead)

  # Continuous-time OU residual: correlation ar_phi^Delta with the gap in years.
  ou_noise <- function(d2) {
    n <- length(d2); x <- numeric(n); x[1] <- rnorm(1, 0, sd_noise)
    for (k in 2:n) {
      rho  <- ar_phi ^ (abs(d2[k - 1] - d2[k]) / 365.25)
      x[k] <- rho * x[k - 1] + rnorm(1, 0, sd_noise * sqrt(1 - rho^2))
    }
    x
  }

  make_subject <- function(id, D) {
    gaps <- sample(gap_pool, K - 1, replace = TRUE)
    d2   <- dmin + rev(cumsum(c(0, gaps)))            # days to diagnosis, decreasing to dmin
    yr   <- (max(d2) - d2) / 365.25
    b1 <- rnorm(1, mu, tau); b2 <- rnorm(1, mu, tau)
    if (runif(1) < hi_base_frac) b2 <- b2 + hi_base_shift
    s1 <- rnorm(1, 0, trend_sd); s2 <- rnorm(1, 0, trend_sd)
    M1 <- b1 + s1 * yr + ou_noise(d2)
    M2 <- b2 + s2 * yr + ou_noise(d2)
    if (D == 1) {
      if (runif(1) < pA) M1 <- M1 + delta1 * ramp(d2, leadA)
      else               M2 <- M2 + delta2 * ramp(d2, leadB)
    }
    data.frame(IDvar = id, D = D, draw_to_dx_days = d2, t = -d2, M1 = M1, M2 = M2)
  }

  df <- do.call(rbind, c(
    lapply(seq_len(n_ctrl),          function(i) make_subject(sprintf("c%04d", i), 0)),
    lapply(seq_len(n_case) + n_ctrl, function(i) make_subject(sprintf("k%04d", i), 1))))
  df$obs_number <- ave(seq_len(nrow(df)), df$IDvar, FUN = seq_along)

  # Fixed comparator composite: a frozen logistic combination, fit on the visit
  # closest to diagnosis (the simulation analogue of a deployed fixed panel).
  fit_panel    <- glm(D ~ M1 + M2, family = binomial, data = df[df$obs_number == K, ])
  df$composite <- as.numeric(predict(fit_panel, newdata = df))
  df
}

## ---- one Monte-Carlo cohort ------------------------------------------------
one_cohort <- function(r) {
  df  <- do.call(simulate_cohort_irregular, c(gen, seed = r))
  tr  <- split_ids(df, train_frac = 0.7, seed = 100 + r)
  te  <- !(df$IDvar %in% tr)
  ctr <- (df$IDvar %in% tr) & df$D == 0

  # Gap-ignorant PEB: the frozen composite through an intercept-only, i.i.d. layer.
  peb     <- ipeb_innovations(df, "composite", reference_ids = tr, id = "IDvar",
                              case = "D", time = "t", slope = "off", innovation = "iid")[, 1]
  peb_row <- eval_score(peb[te], df$D[te], df$IDvar[te], df$draw_to_dx_days[te], peb[ctr], SPEC)

  # Gap-aware iPEB: random slope with gap-scaled AR(1)/OU innovations.
  fit    <- ipeb(df[df$IDvar %in% tr, ], markers = c("M1", "M2"),
                 id = "IDvar", case = "D", time = "t", time_to_dx = "draw_to_dx_days",
                 objective = "sensitivity", alpha = SPEC, slope = "on", innovation = "ar1")
  ip     <- predict(fit, df[te, ])
  ip_row <- eval_score(ip, df$D[te], df$IDvar[te], df$draw_to_dx_days[te], fit$train_ctrl_scores, SPEC)

  rbind(`PEB (gap-ignorant)` = peb_row, `iPEB (gap-aware)` = ip_row)
}

## ---- run and summarize -----------------------------------------------------
acc   <- par_lapply(seq_len(R), function(r) tryCatch(one_cohort(r), error = function(e) NULL))
acc   <- acc[!vapply(acc, is.null, logical(1))]
arr   <- simplify2array(acc)
means <- apply(arr, c(1, 2), mean, na.rm = TRUE)
sds   <- apply(arr, c(1, 2), sd,   na.rm = TRUE)

cat(sprintf("\nScenario 4 (irregular visit spacing): mean (sd) over %d cohorts\n", length(acc)))
print(matrix(sprintf("%.3f (%.3f)", means, sds), nrow(means), dimnames = dimnames(means)), quote = FALSE)
s4_summary <- data.frame(
  method    = rownames(means),
  AUC_mean  = means[, "AUC"],   AUC_sd  = sds[, "AUC"],
  Sens_mean = means[, "SensW"], Sens_sd = sds[, "SensW"],
  Lead_mean = means[, "LT"],    Lead_sd = sds[, "LT"],
  row.names = NULL)
utils::write.csv(s4_summary, "reproduce_S4_means.csv", row.names = FALSE)

## ---- Figure: sensitivity (PEB vs iPEB) -------------------------------------
png("reproduce_S4_sensitivity.png", width = 2600, height = 2000, res = 300)
graphics::par(mar = c(6, 5, 3, 1))
bp <- barplot(means[, "SensW"], col = c("grey60", "#db2777"), names.arg = rownames(means),
              las = 1, ylab = "Sensitivity at 95% specificity",
              main = "Detection under irregular visit spacing",
              ylim = c(0, max(means[, "SensW"] + sds[, "SensW"] / sqrt(length(acc)), na.rm = TRUE) * 1.25),
              cex.lab = 1.8, cex.axis = 1.6, cex.main = 1.7, cex.names = 1.6)
ses <- sds / sqrt(length(acc))          # Monte Carlo standard error of the mean
arrows(bp, means[, "SensW"] - ses[, "SensW"], bp, means[, "SensW"] + ses[, "SensW"],
       angle = 90, code = 3, length = 0.06, lwd = 1.6)   # error bars: +/- 1 MC SE
text(bp, means[, "SensW"] + ses[, "SensW"], sprintf("%.2f", means[, "SensW"]), pos = 3, cex = 1.6)
dev.off()
cat("Wrote reproduce_S4_sensitivity.png\n")
