# Generates the synthetic longitudinal example dataset shipped with the package.
# Run from the package root:  source("data-raw/make_example.R")
# It writes data/ipeb_example.rda. The data are entirely simulated; no real
# patient data are used (PLCO data are controlled-access via NCI CDAS).

set.seed(20260101)

n_ctrl <- 160L
n_case <- 80L
markers <- c("m1", "m2", "m3")

# One subject's longitudinal record.
make_subject <- function(idx, is_case) {
  k <- sample(4:8, 1)                                   # number of visits
  gaps <- round(stats::rnorm(k - 1, mean = 365, sd = 45))  # ~annual, irregular
  gaps <- pmax(gaps, 120)
  # days-to-diagnosis at each visit (decreasing toward the last visit)
  last_gap <- if (is_case) sample(30:200, 1) else sample(200:600, 1)
  ttd <- rev(cumsum(c(last_gap, gaps)))                 # visit 1 furthest from dx
  time <- (max(ttd) - ttd) / 365.25                     # years since first visit

  b0 <- stats::rnorm(3, mean = c(10, 9, 8), sd = c(2.0, 1.8, 1.6))  # random intercepts
  eps <- sapply(1:3, function(j) as.numeric(stats::arima.sim(list(ar = 0.5), n = k)) * 0.8)
  Y <- sweep(eps, 2, b0, "+")

  if (is_case) {
    # m1 rises from ~2 years out, m2 rises late (~6 months); m3 uninformative
    ramp1 <- pmax(0, (730 - ttd) / 730) * 5.0
    ramp2 <- pmax(0, (180 - ttd) / 180) * 6.0
    Y[, 1] <- Y[, 1] + ramp1
    Y[, 2] <- Y[, 2] + ramp2
  }
  data.frame(id = idx, case = as.integer(is_case), time = time,
             time_to_dx = ttd, m1 = Y[, 1], m2 = Y[, 2], m3 = Y[, 3])
}

recs <- vector("list", n_ctrl + n_case)
for (i in seq_len(n_ctrl)) recs[[i]] <- make_subject(i, FALSE)
for (i in seq_len(n_case)) recs[[n_ctrl + i]] <- make_subject(n_ctrl + i, TRUE)
ipeb_example <- do.call(rbind, recs)

# Deterministic train/test split (~70/30) by subject.
ids <- unique(ipeb_example$id)
train_ids <- sample(ids, round(0.70 * length(ids)))
ipeb_example$split <- ifelse(ipeb_example$id %in% train_ids, "train", "test")
rownames(ipeb_example) <- NULL

dir.create("data", showWarnings = FALSE)
save(ipeb_example, file = "data/ipeb_example.rda", compress = "xz")
message("Wrote data/ipeb_example.rda: ", nrow(ipeb_example), " rows, ",
        length(ids), " subjects.")
