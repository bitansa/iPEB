# Scenario 1 (marker heterogeneity): sensitivity / AUC. Reproduces the S1 row of
# the paper's mechanism-simulation table (PEB vs iPEB) and the mean-ROC figure.
# Run from this folder:  source("reproduce_S1.R")

source("sim_helpers.R")

gen  <- list(n_ctrl = 600, n_case = 200, K = 8, dmax = 1080, dmin = 90,
             hi_base_frac = 0.20, hi_base_shift = 8,
             pA = 0.60, delta1 = 6, delta2 = 6, leadA = 365, leadB = 365)
R      <- 200        # Monte-Carlo cohorts (reduce for a quick check)
SPEC   <- 0.95
SLOPE  <- "on"       # >= 5 visits here -> random slope; AR(1)/OU off (no serial correlation)
INNOV  <- "iid"

one_cohort <- function(r) {
  df  <- do.call(simulate_cohort, c(gen, seed = r))
  tr  <- split_ids(df, train_frac = 0.7, seed = 100 + r)
  te  <- !(df$IDvar %in% tr)
  ctr <- df$IDvar %in% tr & df$D == 0

  # PEB: the fixed logistic composite scored through the time-gap-aware layer.
  peb <- ipeb_innovations(df, markers = "composite", reference_ids = tr,
                          id = "IDvar", case = "D", time = "t",
                          slope = SLOPE, innovation = INNOV)[, 1]
  peb_row <- eval_score(peb[te], df$D[te], df$IDvar[te], df$draw_to_dx_days[te], peb[ctr], SPEC)

  # iPEB: learned combination of the M1/M2 innovations (sensitivity objective).
  fit <- ipeb(df[df$IDvar %in% tr, ], markers = c("M1", "M2"),
              id = "IDvar", case = "D", time = "t", time_to_dx = "draw_to_dx_days",
              objective = "sensitivity", alpha = SPEC, slope = SLOPE, innovation = INNOV)
  ip  <- predict(fit, df[te, ])
  ipeb_row <- eval_score(ip, df$D[te], df$IDvar[te], df$draw_to_dx_days[te],
                         fit$train_ctrl_scores, SPEC)

  list(tab = rbind(PEB = peb_row, iPEB = ipeb_row),
       roc_peb = roc_curve(peb[te], df$D[te], df$IDvar[te]),
       roc_ipeb = roc_curve(ip, df$D[te], df$IDvar[te]))
}

acc <- par_lapply(seq_len(R), function(r) tryCatch(one_cohort(r), error = function(e) NULL))
acc <- acc[!vapply(acc, is.null, logical(1))]

arr   <- simplify2array(lapply(acc, `[[`, "tab"))
means <- apply(arr, c(1, 2), mean, na.rm = TRUE)
sds   <- apply(arr, c(1, 2), sd,   na.rm = TRUE)

cat(sprintf("\nScenario 1 (marker heterogeneity): mean (sd) over %d cohorts\n", length(acc)))
print(matrix(sprintf("%.3f (%.3f)", means, sds), nrow(means), dimnames = dimnames(means)), quote = FALSE)
cat("\nPaper (Table 1): PEB 0.972 / 0.829 / 0.633 ; iPEB 0.984 / 0.923 / 0.619  (AUC / Sens / Lead)\n")

# Mean test ROC (PEB vs iPEB) with AUCs -- reproduces Web Figure 1.
g   <- seq(0, 1, 0.02)
mpb <- Reduce(`+`, lapply(acc, `[[`, "roc_peb"))  / length(acc)
mip <- Reduce(`+`, lapply(acc, `[[`, "roc_ipeb")) / length(acc)
png("reproduce_S1_ROC.png", width = 2400, height = 2000, res = 300)
plot(g, mip, type = "l", lwd = 3, col = "#db2777", xlab = "1 - specificity",
     ylab = "Sensitivity", main = "Simulation 1: mean ROC")
lines(g, mpb, lwd = 3, col = "#2563eb"); abline(0, 1, lty = 3, col = "grey60")
legend("bottomright", bty = "n", lwd = 3, col = c("#2563eb", "#db2777"),
       legend = sprintf("%-5s AUC %.3f", c("PEB", "iPEB"),
                        c(auc_trap(g, mpb), auc_trap(g, mip))))
dev.off()
cat("Wrote reproduce_S1_ROC.png\n")
