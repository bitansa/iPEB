# Generates a synthetic six-marker longitudinal cohort that mimics the STRUCTURE
# of the PLCO lung analysis (centers, covariates, irregular visits, four
# informative markers plus two weaker ones). It is entirely simulated: the real
# PLCO data are controlled-access via NCI CDAS and are NOT distributed here, so
# the numbers this produces will NOT match the paper -- it exists only so the
# reproduce_lung.R workflow runs end-to-end on shareable data.
# Run:  source("make_synthetic_lung.R")   -> writes synthetic_lung.csv

set.seed(20260501)

markers <- c("CA125", "CEA", "CYFRA21", "PSFTPB", "OPN", "HE4")  # 4MP = first four
n_subj  <- 2000
p_case  <- 0.16
centers <- 1:10

make_subject <- function(idx) {
  is_case <- stats::rbinom(1, 1, p_case)
  center  <- sample(centers, 1)
  k       <- sample(2:6, 1)
  gaps    <- pmax(round(stats::rnorm(k - 1, 365, 45)), 120)
  last    <- if (is_case) sample(30:700, 1) else sample(250:1200, 1)
  ttd     <- rev(cumsum(c(last, gaps)))                 # days-to-dx, far -> near
  time    <- (max(ttd) - ttd) / 365.25
  age     <- round(stats::rnorm(1, 62, 5))
  batch   <- sample(1:4, 1)
  yr      <- sample(1993:2001, 1)

  b0  <- stats::rnorm(6, mean = c(2.5, 2.0, 1.5, 2.0, 3.0, 2.5), sd = 0.8)
  eps <- sapply(1:6, function(j) as.numeric(stats::arima.sim(list(ar = 0.4), k)) * 0.5)
  Y   <- sweep(eps, 2, b0, "+") + (age - 62) * 0.01     # mild covariate effect
  if (is_case) {
    ramp <- pmax(0, (730 - ttd) / 730)                  # rise over ~2 years pre-dx
    Y <- Y + outer(ramp, c(1.5, 1.8, 1.2, 1.6, 0.5, 0.6))  # 4MP strong, OPN/HE4 weak
  }
  out <- data.frame(id = idx, case = is_case, center = center, time = time,
                    time_to_dx = ttd, age = age, batch = batch, study_year = yr)
  out[markers] <- Y
  out
}

lung <- do.call(rbind, lapply(seq_len(n_subj), make_subject))
write.csv(lung, "synthetic_lung.csv", row.names = FALSE)
message(sprintf("Wrote synthetic_lung.csv: %d rows, %d subjects (%d cases).",
                nrow(lung), n_subj, sum(tapply(lung$case, lung$id, max))))
