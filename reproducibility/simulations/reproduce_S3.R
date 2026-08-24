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
cat("\nPaper (Table 1): 0.880/0.529/0.548 ; 0.877/0.518/0.563 ; 0.974/0.846/0.632 (AUC / Sens / Lead)\n")

# Persist the exact mean (sd) numbers behind the S3 rows of Table 1, mirroring the
# CSVs written for S4 and the factorial, so the table is fully reproducible.
s3_summary <- data.frame(
  method    = rownames(means),
  AUC_mean  = means[, "AUC"],   AUC_sd  = sds[, "AUC"],
  Sens_mean = means[, "SensW"], Sens_sd = sds[, "SensW"],
  Lead_mean = means[, "LT"],    Lead_sd = sds[, "LT"],
  row.names = NULL)
utils::write.csv(s3_summary, "reproduce_S3_means.csv", row.names = FALSE)
cat("Wrote reproduce_S3_means.csv\n")

png("reproduce_S3_sensitivity.png", width = 2600, height = 2100, res = 300)
graphics::par(mar = c(7, 5, 2, 1), cex.lab = 1.8, cex.axis = 1.6, cex.main = 1.7)
ses <- sds / sqrt(length(acc))          # Monte Carlo standard error of the mean
bp <- barplot(means[, "SensW"], col = c("grey60", "#93c5fd", "#db2777"),
              names.arg = rep("", nrow(means)),
              ylab = "Sensitivity at 95% specificity",
              ylim = c(0, max(means[, "SensW"] + ses[, "SensW"]) * 1.2))
arrows(bp, means[, "SensW"] - ses[, "SensW"], bp, means[, "SensW"] + ses[, "SensW"],
       angle = 90, code = 3, length = 0.06, lwd = 1.6)   # error bars: +/- 1 MC SE
text(bp, means[, "SensW"] + ses[, "SensW"], sprintf("%.2f", means[, "SensW"]), pos = 3, cex = 1.6)
# compact 3-line horizontal labels under each bar (upright, no overlap)
labs <- c("PEB\nintercept\ni.i.d.", "iPEB\nintercept\ni.i.d.", "iPEB\nslope\nOU")
text(x = bp, y = par("usr")[3] - 0.03 * diff(par("usr")[3:4]), labels = labs,
     xpd = TRUE, adj = c(0.5, 1), cex = 1.4)
dev.off()
cat("Wrote reproduce_S3_sensitivity.png\n")

# Healthy-trajectory figure (reproduces the S3 trend panel): example healthy
# subjects whose marker levels drift over time under the random slope + AR(1)/OU.
dfx   <- do.call(simulate_cohort, c(gen, seed = 1))
d2max <- max(dfx$draw_to_dx_days)
set.seed(1); ex <- sample(unique(dfx$IDvar[dfx$D == 0]), 8)
sub   <- dfx[dfx$IDvar %in% ex, ]
png("reproduce_S3_trajectories.png", width = 2600, height = 2200, res = 300)
graphics::par(mar = c(5, 5, 3, 1))
plot(NA, xlim = c(0, d2max / 365.25), ylim = range(sub$M1),
     xlab = "Time from first visit (years)", ylab = "Marker M1 (healthy subjects)",
     main = "Example healthy trajectories", cex.lab = 1.8, cex.axis = 1.6, cex.main = 1.7)
cols <- grDevices::hcl.colors(length(ex), "Dark 3")
for (k in seq_along(ex)) {
  s <- dfx[dfx$IDvar == ex[k], ]; o <- order(s$draw_to_dx_days, decreasing = TRUE)
  lines((d2max - s$draw_to_dx_days[o]) / 365.25, s$M1[o], type = "b", pch = 16,
        lwd = 1.8, col = cols[k])
}
dev.off()
cat("Wrote reproduce_S3_trajectories.png\n")
