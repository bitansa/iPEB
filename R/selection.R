# Objective-driven feature selection (internal). A candidate marker subset is
# scored by the objective's own value (= -loss) on a held-out validation split
# of the training data, using a scalar-combiner fit. Directions: "forward"
# grows from the empty set, "backward" prunes from the full set, "best" is an
# exhaustive search for small panels.
.select_markers <- function(dat, train_ids, marker_cols, objective, cov_cols,
                            alpha, window, use_slope, ar1, time_col, ttd_col,
                            val_frac, lam, direction, tol = 0, max_keep = Inf) {
  if (length(marker_cols) < 2L) return(marker_cols)
  if (is.null(lam)) lam <- .lambda_profile(objective)

  isc <- tapply(dat$D, dat$IDvar, max)
  trc <- intersect(train_ids, names(isc)[isc == 1])
  trk <- intersect(train_ids, names(isc)[isc == 0])
  f <- c(sample(trc, round((1 - val_frac) * length(trc))),
         sample(trk, round((1 - val_frac) * length(trk))))
  iv <- dat$IDvar %in% setdiff(train_ids, f)

  val_obj <- function(cols) {
    if (length(cols) == 0L) return(-Inf)
    s <- tryCatch(.fit_layer_and_combiner(dat, fit_ids = f, tune_ids = f,
            marker_cols = cols, variant = "scalar", objective = objective,
            cov_cols = cov_cols, alpha = alpha, window = window,
            use_slope = use_slope, ar1 = ar1, time_col = time_col,
            ttd_col = ttd_col, lam = lam)$score,
          error = function(e) NULL)
    if (is.null(s)) return(-Inf)
    -1 * .loss(s[iv], dat$D[iv], dat$IDvar[iv], dat[[ttd_col]][iv],
               s[iv & dat$D == 0], alpha, window, lam["lam1"], lam["lam2"], lam["lam3"])
  }

  forward <- function() {
    chosen <- character(0); pool <- marker_cols; best <- -Inf
    while (length(pool) > 0L) {
      sc <- vapply(pool, function(m) val_obj(c(chosen, m)), numeric(1))
      j <- which.max(sc)
      if (sc[j] > best + tol) {
        chosen <- c(chosen, pool[j]); pool <- pool[-j]; best <- sc[j]
      } else break
    }
    if (length(chosen) == 0L) marker_cols else chosen
  }
  backward <- function() {
    chosen <- marker_cols; best <- val_obj(chosen)
    while (length(chosen) > 1L) {
      sc <- vapply(seq_along(chosen), function(k) val_obj(chosen[-k]), numeric(1))
      j <- which.max(sc)
      if (length(chosen) > max_keep) { chosen <- chosen[-j]; best <- sc[j] }
      else if (is.infinite(max_keep) && sc[j] >= best - tol) { chosen <- chosen[-j]; best <- sc[j] }
      else break
    }
    chosen
  }
  bestsub <- function() {
    ns <- length(marker_cols)
    if (2^ns - 1 > 1023) { warning("best-subset too large (>10 markers); using forward selection"); return(forward()) }
    subs <- unlist(lapply(seq_len(ns), function(k) utils::combn(marker_cols, k, simplify = FALSE)),
                   recursive = FALSE)
    sc <- vapply(subs, val_obj, numeric(1))
    subs[[which.max(sc)]]
  }
  switch(direction, forward = forward(), backward = backward(), best = bestsub())
}
