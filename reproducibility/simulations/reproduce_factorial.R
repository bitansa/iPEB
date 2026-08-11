# Factorial study (48 cells): iPEB vs best-marginal PEB / equal-weight / PCA-1.
# Reproduces Web Tables 1-2 (per-cell means) and the summary figures. All scores
# go through library(iPEB): comparators via ipeb_innovations(), iPEB via ipeb().
#
# NOTE: this is the heavy simulation. With R_MC = 100 (the paper's setting) and
# p up to 10 it can take of order an hour even in parallel. For a quick sanity
# check set R_MC <- 20; use 100 for the reported numbers.
# Run:  source("reproduce_factorial.R")

source("sim_helpers.R")

## ---- general p-marker data generator (identical to the paper) --------------
simulate_p <- function(n_ctrl = 800, n_case = 200, p = 5, K = 6, dmax = 1080, dmin = 90,
                       tau = NULL, sigma = NULL, corr = "low", spar = 0.4, delta = 0.4,
                       W0_months = 12, trend_sd = 0, ar_phi = 0, misspec = "none",
                       irregular = TRUE, seed = 1) {
  set.seed(seed)
  if (is.null(tau))   tau   <- runif(p, 0.3, 0.8)
  if (is.null(sigma)) sigma <- runif(p, 0.5, 1.0)
  Rho <- switch(corr,
    low   = 0.2^abs(outer(1:p, 1:p, "-")),
    high  = 0.7^abs(outer(1:p, 1:p, "-")),
    block = { h <- ceiling(p / 2); m <- matrix(0.2, p, p)
              m[1:h, 1:h] <- 0.7; m[(h + 1):p, (h + 1):p] <- 0.7; diag(m) <- 1; m },
    diag(p))
  Sig <- diag(sigma) %*% Rho %*% diag(sigma) + diag(1e-6, p); Lc <- chol(Sig)
  nD  <- max(1, round(spar * p)); Dset <- sort(sample(p, nD))
  del <- numeric(p); del[Dset] <- runif(nD, 0.6 * delta, 1.4 * delta) * sigma[Dset]
  W0d <- W0_months * 30.4375
  reg_d2 <- round(seq(dmax, dmin, length.out = K))
  make_sched <- function() {
    if (!irregular) return(reg_d2)
    g <- sample(c(183, 365), K - 1, replace = TRUE)
    cum <- cumsum(c(0, g)); round(dmax - cum * (dmax - dmin) / cum[length(cum)])
  }
  resid <- function(d2v) {
    Z  <- if (misspec == "t5") matrix(rt(K * p, 5) / sqrt(5 / 3), K, p) else matrix(rnorm(K * p), K, p)
    gy <- c(0, abs(diff(d2v)) / 365.25)
    if (misspec == "ar1") for (j in 1:p) { x <- Z[, j]
      for (k in 2:K) { r <- 0.4^gy[k]; x[k] <- r * x[k - 1] + sqrt(1 - r^2) * Z[k, j] }; Z[, j] <- x }
    if (ar_phi > 0)       for (j in 1:p) { x <- Z[, j]
      for (k in 2:K) { r <- ar_phi^gy[k]; x[k] <- r * x[k - 1] + sqrt(1 - r^2) * Z[k, j] }; Z[, j] <- x }
    E <- Z %*% Lc
    if (misspec == "lognormal") { E <- exp(E); E <- sweep(E, 2, colMeans(E), "-")
      E <- sweep(E, 2, pmax(apply(E, 2, sd), 1e-8), "/") %*% diag(sigma) }
    E
  }
  mk <- function(id, D) {
    b <- rnorm(p, 0, tau); s <- rnorm(p, 0, trend_sd)
    d2 <- make_sched(); yr <- (max(d2) - d2) / 365.25
    rr <- pmax(0, 1 - d2 / W0d)
    Y <- sweep(resid(d2), 2, b, "+") + outer(yr, s)
    if (D == 1) Y <- Y + outer(rr, del)
    data.frame(IDvar = id, D = D, draw_to_dx_days = d2, t = -d2,
               stats::setNames(as.data.frame(Y), paste0("M", 1:p)))
  }
  do.call(rbind, c(lapply(1:n_ctrl,           function(i) mk(sprintf("c%05d", i), 0)),
                   lapply((1:n_case) + n_ctrl, function(i) mk(sprintf("k%05d", i), 1))))
}

## ---- all scores on p markers, through the package --------------------------
score_all_p <- function(df, tr, markers, spec, ar1, use_slope) {
  slope <- if (use_slope) "on" else "off"; innov <- if (ar1) "ar1" else "iid"
  p <- length(markers)
  R <- ipeb_innovations(df, markers, reference_ids = tr, id = "IDvar", case = "D",
                        time = "t", slope = slope, innovation = innov)
  ctr <- df$IDvar %in% tr & df$D == 0
  eqw <- as.numeric(R %*% rep(1, p)) / sqrt(p)
  Rc  <- R[ctr, , drop = FALSE]; Rc <- Rc[stats::complete.cases(Rc), , drop = FALSE]
  v   <- tryCatch(abs(prcomp(Rc, center = TRUE, scale. = FALSE)$rotation[, 1]),
                  error = function(e) rep(1, p))
  v   <- v / sqrt(sum(v^2)); pca1 <- as.numeric(R %*% v)
  train_df <- df[df$IDvar %in% tr, ]
  fit_s <- ipeb(train_df, markers, id = "IDvar", case = "D", time = "t", time_to_dx = "draw_to_dx_days",
                objective = "sensitivity", alpha = spec, slope = slope, innovation = innov)
  fit_l <- ipeb(train_df, markers, id = "IDvar", case = "D", time = "t", time_to_dx = "draw_to_dx_days",
                objective = "leadtime", alpha = spec, slope = slope, innovation = innov)
  list(marg = lapply(1:p, function(j) R[, j]), EqualW = eqw, PCA1 = pca1,
       iPEB_sens = predict(fit_s, df), iPEB_lead = predict(fit_l, df))
}

## ---- one replicate (75/25 split) -------------------------------------------
run_one_p <- function(df, markers, spec, ar1, use_slope, seed) {
  set.seed(seed); isC <- tapply(df$D, df$IDvar, max)
  cs <- names(isC)[isC == 1]; ct <- names(isC)[isC == 0]
  tr <- c(sample(cs, round(.75 * length(cs))), sample(ct, round(.75 * length(ct))))
  it <- !(df$IDvar %in% tr); ctr <- df$IDvar %in% tr & df$D == 0; itr <- df$IDvar %in% tr
  S  <- score_all_p(df, tr, markers, spec, ar1, use_slope)
  m  <- function(s, idx, which) unname(metrics_subj(s[idx], df$D[idx], df$IDvar[idx],
                                                    df$draw_to_dx_days[idx], s[ctr], spec)[which])
  jS <- which.max(vapply(S$marg, function(s) m(s, itr, "SensW"), numeric(1)))
  jL <- which.max(vapply(S$marg, function(s) m(s, itr, "LT"),    numeric(1)))
  list(
    sens = c(`Best-marginal PEB` = m(S$marg[[jS]], it, "SensW"), `Equal-weights` = m(S$EqualW, it, "SensW"),
             `PCA-1` = m(S$PCA1, it, "SensW"), iPEB = m(S$iPEB_sens, it, "SensW")),
    lead = c(`Best-marginal PEB` = m(S$marg[[jL]], it, "LT"), `Equal-weights` = m(S$EqualW, it, "LT"),
             `PCA-1` = m(S$PCA1, it, "LT"), iPEB = m(S$iPEB_lead, it, "LT")))
}

run_cell <- function(R, gen_args, markers, spec, ar1, use_slope) {
  acc <- par_lapply(1:R, function(r) {
    df <- do.call(simulate_p, c(gen_args, seed = r))
    tryCatch(run_one_p(df, markers, spec, ar1, use_slope, seed = 100 + r), error = function(e) NULL)
  })
  acc <- acc[vapply(acc, is.list, logical(1))]
  sensM <- sapply(acc, function(o) o$sens); leadM <- sapply(acc, function(o) o$lead)
  list(sensMean = rowMeans(sensM, na.rm = TRUE), leadMean = rowMeans(leadM, na.rm = TRUE))
}

## ---- factorial grid --------------------------------------------------------
R_MC <- 100     # paper setting (100 reps); use a small value (e.g. 20) for a quick check
SPEC <- 0.95
grid <- expand.grid(corr = c("low", "high", "block"), spar = c(0.4, 0.6),
                    misspec = c("none", "t5", "ar1", "lognormal"), stringsAsFactors = FALSE)
rows <- c("Best-marginal PEB", "Equal-weights", "PCA-1", "iPEB")

sens_cells <- list(); lead_cells <- list()
for (P in c(5, 10)) {
  MK <- paste0("M", 1:P)
  for (i in seq_len(nrow(grid))) {
    g  <- grid[i, ]
    ga <- list(n_ctrl = 800, n_case = 200, p = P, K = 6, spar = g$spar, delta = 0.4,
               corr = g$corr, misspec = g$misspec, W0_months = 12)
    res <- run_cell(R_MC, ga, MK, SPEC, ar1 = (g$misspec == "ar1"), use_slope = TRUE)
    tag <- sprintf("p=%d corr=%s spar=%.1f misspec=%s", P, g$corr, g$spar, g$misspec)
    cat(sprintf("[%s]  iPEB Sens %.3f  (best-marg %.3f, eq %.3f, pca %.3f)\n",
                tag, res$sensMean["iPEB"], res$sensMean["Best-marginal PEB"],
                res$sensMean["Equal-weights"], res$sensMean["PCA-1"]))
    sens_cells[[tag]] <- res$sensMean[rows]; lead_cells[[tag]] <- res$leadMean[rows]
  }
}

sens_tab <- do.call(rbind, sens_cells); lead_tab <- do.call(rbind, lead_cells)
write.csv(round(sens_tab, 4), "reproduce_factorial_SensW_by_cell.csv")
write.csv(round(lead_tab, 4), "reproduce_factorial_LT_by_cell.csv")

cat("\n=== mean Sens_W by method across all 48 cells ===\n")
print(round(colMeans(sens_tab), 3))
cat("Paper: iPEB 0.343 ; best-marginal 0.316 ; equal-weight 0.308 ; PCA-1 0.305\n")
cat(sprintf("iPEB beat: best-marginal in %d/48, equal-weight in %d/48, PCA-1 in %d/48 cells\n",
            sum(sens_tab[, "iPEB"] > sens_tab[, "Best-marginal PEB"]),
            sum(sens_tab[, "iPEB"] > sens_tab[, "Equal-weights"]),
            sum(sens_tab[, "iPEB"] > sens_tab[, "PCA-1"])))
cat("Paper: 39/48, 46/48, 47/48\n")

png("reproduce_factorial_summary.png", width = 2600, height = 2600, res = 300)
graphics::par(mar = c(16, 5, 1.5, 1), cex.lab = 1.8, cex.axis = 1.6, cex.main = 1.7)
m <- colMeans(sens_tab)
bp <- barplot(m, col = c("grey65", "#60a5fa", "#a78bfa", "#db2777"), las = 2, cex.names = 1.6,
              ylab = "Sensitivity at 95% specificity", ylim = c(0, max(m) * 1.25))
text(bp, m, sprintf("%.3f", m), pos = 3, cex = 1.6)
dev.off()
cat("Wrote reproduce_factorial_summary.png and the per-cell CSVs\n")
