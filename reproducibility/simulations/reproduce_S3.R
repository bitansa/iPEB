# Scenario 3 (longitudinal drift + serial correlation): value of the layer.
# Reproduces the S3 rows of the mechanism table and the sensitivity figure:
# PEB (intercept, i.i.d.), iPEB (intercept, i.i.d.), iPEB (slope + AR(1)/OU).
# Run:  source("reproduce_S3.R")

source("sim_helpers.R")

gen  <- list(n_ctrl = 600, n_case = 200, K = 8, dmax = 1080, dmin = 90,
             mu = 10, tau = 2, sd_noise = 1, trend_sd = 1.5, ar_phi = 0.6,
             hi_base_frac = 0.20, hi_base_shift = 8,
             pA = 0.60, delta1 = 6, delta2 = 6, leadA = 365, leadB = 365)
R    <- 200
SPEC <- 0.95

fit_ipeb <- function(df, tr, slope, innovation) {
  fit <- ipeb(df[df$IDvar %in% tr, ], markers = c("M1", "M2"),
              id = "IDvar", case = "D", time = "t", time_to_dx = "draw_to_dx_days",
              objective = "sensitivity", alpha = SPEC, slope = slope, innovation = innovation)
  ip <- predict(fit, df[!(df$IDvar %in% tr), ])
  eval_score(ip, df$D[!(df$IDvar %in% tr)], df$IDvar[!(df$IDvar %in% tr)],
             df$draw_to_dx_days[!(df$IDvar %in% tr)], fit$train_ctrl_scores, SPEC)
}

one_cohort <- function(r) {
  df  <- do.call(simulate_cohort, c(gen, seed = r))
  tr  <- split_ids(df, train_frac = 0.7, seed = 100 + r)
  te  <- !(df$IDvar %in% tr)
  ctr <- df$IDvar %in% tr & df$D == 0

  peb <- ipeb_innovations(df, "composite", reference_ids = tr, id = "IDvar", case = "D",
                          time = "t", slope = "off", innovation = "iid")[, 1]
  peb_row <- eval_score(peb[te], df$D[te], df$IDvar[te], df$draw_to_dx_days[te], peb[ctr], SPEC)

  rbind(`PEB (intercept, i.i.d.)`  = peb_row,
        `iPEB (intercept, i.i.d.)` = fit_ipeb(df, tr, "off", "iid"),
        `iPEB (slope + OU)`        = fit_ipeb(df, tr, "on",  "ar1"))
}

acc <- par_lapply(seq_len(R), function(r) tryCatch(one_cohort(r), error = function(e) NULL))
acc <- acc[!vapply(acc, is.null, logical(1))]
arr <- simplify2array(acc)
means <- apply(arr, c(1, 2), mean, na.rm = TRUE)
sds   <- apply(arr, c(1, 2), sd,   na.rm = TRUE)

cat(sprintf("\nScenario 3 (longitudinal trend): mean (sd) over %d cohorts\n", length(acc)))
print(matrix(sprintf("%.3f (%.3f)", means, sds), nrow(means), dimnames = dimnames(means)), quote = FALSE)
cat("\nPaper (Table 1): 0.880/0.529/0.548 ; 0.877/0.518/0.567 ; 0.974/0.846/0.632 (AUC / Sens / Lead)\n")

png("reproduce_S3_sensitivity.png", width = 2600, height = 2000, res = 300)
graphics::par(mar = c(9, 5, 2, 1))
bp <- barplot(means[, "SensW"], col = c("grey60", "#93c5fd", "#db2777"),
              names.arg = rownames(means), las = 2,
              ylab = "Sensitivity at 95% specificity", ylim = c(0, max(means[, "SensW"]) * 1.2))
text(bp, means[, "SensW"], sprintf("%.2f", means[, "SensW"]), pos = 3)
dev.off()
cat("Wrote reproduce_S3_sensitivity.png\n")
