# Scenario 2 (differential timing): lead time. Reproduces the S2 row of the
# mechanism table (PEB vs iPEB) and the lead-time figure (PEB / PEB M1 / PEB M2 /
# iPEB, with each method's sensitivity annotated). Run:  source("reproduce_S2.R")

source("sim_helpers.R")

gen  <- list(n_ctrl = 600, n_case = 200, K = 10, dmax = 1080, dmin = 30, sd_noise = 1.5,
             hi_base_frac = 0.20, hi_base_shift = 8, both_markers = TRUE,
             delta1 = 6, delta2 = 6, leadA = 730, leadB = 90, dev_panel = TRUE)
R      <- 200
SPEC   <- 0.95
SLOPE  <- "on"
INNOV  <- "iid"

one_cohort <- function(r) {
  df  <- do.call(simulate_cohort, c(gen, seed = r))
  tr  <- split_ids(df, train_frac = 0.7, seed = 100 + r)
  te  <- !(df$IDvar %in% tr)
  ctr <- df$IDvar %in% tr & df$D == 0
  ev  <- function(s) eval_score(s[te], df$D[te], df$IDvar[te], df$draw_to_dx_days[te], s[ctr], SPEC)

  Rm  <- ipeb_innovations(df, c("M1", "M2"), reference_ids = tr, id = "IDvar",
                          case = "D", time = "t", slope = SLOPE, innovation = INNOV)
  peb <- ipeb_innovations(df, "composite", reference_ids = tr, id = "IDvar",
                          case = "D", time = "t", slope = SLOPE, innovation = INNOV)[, 1]

  fit <- ipeb(df[df$IDvar %in% tr, ], markers = c("M1", "M2"),
              id = "IDvar", case = "D", time = "t", time_to_dx = "draw_to_dx_days",
              objective = "leadtime", alpha = SPEC, slope = SLOPE, innovation = INNOV)
  ip  <- predict(fit, df[te, ])
  ipeb_row <- eval_score(ip, df$D[te], df$IDvar[te], df$draw_to_dx_days[te],
                         fit$train_ctrl_scores, SPEC)

  rbind(PEB = ev(peb), `PEB M1` = ev(Rm[, 1]), `PEB M2` = ev(Rm[, 2]), iPEB = ipeb_row)
}

acc <- par_lapply(seq_len(R), function(r) tryCatch(one_cohort(r), error = function(e) NULL))
acc <- acc[!vapply(acc, is.null, logical(1))]
arr <- simplify2array(acc)
means <- apply(arr, c(1, 2), mean, na.rm = TRUE)
sds   <- apply(arr, c(1, 2), sd,   na.rm = TRUE)

cat(sprintf("\nScenario 2 (differential timing): mean (sd) over %d cohorts\n", length(acc)))
print(matrix(sprintf("%.3f (%.3f)", means, sds), nrow(means), dimnames = dimnames(means)), quote = FALSE)
cat("\nPaper (Table 1): PEB 0.973 / 0.844 / 0.252 ; iPEB 0.995 / 0.998 / 1.194 (AUC / Sens / Lead)\n")

# Persist the exact mean (sd) numbers behind the S2 rows of Table 1, mirroring the
# CSVs written for S3, S4, and the factorial, so the table is fully reproducible.
s2_summary <- data.frame(
  method    = rownames(means),
  AUC_mean  = means[, "AUC"],   AUC_sd  = sds[, "AUC"],
  Sens_mean = means[, "SensW"], Sens_sd = sds[, "SensW"],
  Lead_mean = means[, "LT"],    Lead_sd = sds[, "LT"],
  row.names = NULL)
utils::write.csv(s2_summary, "reproduce_S2_means.csv", row.names = FALSE)
cat("Wrote reproduce_S2_means.csv\n")

# Lead-time figure (reproduces the S2 bar panel).
ord <- c("PEB", "PEB M1", "PEB M2", "iPEB")
png("reproduce_S2_leadtime.png", width = 2600, height = 2200, res = 300)
graphics::par(mar = c(8, 5, 2, 1), cex.lab = 1.8, cex.axis = 1.6, cex.main = 1.7)
bp <- barplot(means[ord, "LT"], col = c("#2563eb", "#60a5fa", "#a5b4fc", "#db2777"),
              names.arg = ord, las = 2, cex.names = 1.6, ylab = "Median lead time (yr)",
              ylim = c(0, max(means[, "LT"]) * 1.2))
text(bp, means[ord, "LT"], sprintf("%.2f", means[ord, "LT"]), pos = 3, cex = 1.6)
mtext(sprintf("Sens:  %s", paste(sprintf("%s %.2f", ord, means[ord, "SensW"]), collapse = "   ")),
      side = 1, line = 6.2, cex = 1.4)
dev.off()
cat("Wrote reproduce_S2_leadtime.png\n")

# Marker-trajectory figure (reproduces the S2 trajectory panel): mean case
# trajectories, showing M1 rising early and M2 rising late before diagnosis.
dfx <- do.call(simulate_cohort, c(gen, seed = 1))
cs  <- dfx[dfx$D == 1, ]
m1  <- tapply(cs$M1, cs$draw_to_dx_days, mean)
m2  <- tapply(cs$M2, cs$draw_to_dx_days, mean)
xx  <- as.numeric(names(m1)) / 365.25            # years before diagnosis
png("reproduce_S2_trajectories.png", width = 2600, height = 2200, res = 300)
graphics::par(mar = c(5, 5, 3, 1))
plot(xx, m1, type = "b", pch = 16, lwd = 2.5, col = "#2563eb", xlim = rev(range(xx)),
     ylim = range(m1, m2), xlab = "Years before diagnosis",
     ylab = "Mean marker level (cases)", main = "Mean case marker trajectories",
     cex.lab = 1.8, cex.axis = 1.6, cex.main = 1.7)
lines(xx, m2, type = "b", pch = 17, lwd = 2.5, col = "#db2777")
legend("topright", c("M1 (early-rising)", "M2 (late-rising)"), col = c("#2563eb", "#db2777"),
       pch = c(16, 17), lwd = 2.5, bty = "n", cex = 1.6)
dev.off()
cat("Wrote reproduce_S2_trajectories.png\n")
