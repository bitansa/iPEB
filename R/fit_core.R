# Internal engine: fit the covariate + innovation layer and the combiner for a
# single variant, returning both the in-sample scores and the parameters needed
# to score new data. `fit_ids` are the healthy controls used to estimate the
# layer; `tune_ids` are the subjects used to learn the combiner.
.fit_layer_and_combiner <- function(dat, fit_ids, tune_ids, marker_cols,
                                    variant, objective, cov_cols, alpha, window,
                                    use_slope, ar1, time_col, ttd_col, lam) {
  if (is.null(lam)) lam <- .lambda_profile(objective)

  cov_models <- .cov_fit(dat, fit_ids, marker_cols, cov_cols)
  d <- .cov_apply(dat, cov_models)

  uni_params <- NULL; mv_params <- NULL
  if (variant == "scalar") {
    uni_params <- lapply(marker_cols, function(col)
      .uni_fit(d, fit_ids, col, time_col = time_col, h_col = time_col,
               use_slope = use_slope, ar1 = ar1))
    names(uni_params) <- marker_cols
    X <- do.call(cbind, lapply(marker_cols, function(col)
      .uni_apply_marker(d, uni_params[[col]], col, time_col = time_col, h_col = time_col)))
  } else {
    Y <- as.matrix(d[, marker_cols, drop = FALSE])
    tt <- d[[time_col]]; ids <- d$IDvar
    fc <- ids %in% fit_ids & d$D == 0
    mv_params <- .mv_fit(Y[fc, , drop = FALSE], ids[fc], tt[fc], ar1 = ar1)
    X <- .mv_apply(d, mv_params, marker_cols, time_col = time_col)
  }

  it <- d$IDvar %in% tune_ids
  w <- .learn_combiner(X[it, , drop = FALSE], d$D[it], d$IDvar[it], d[[ttd_col]][it],
                       alpha, window, lam["lam1"], lam["lam2"], lam["lam3"],
                       whitened = (variant == "mv"))

  list(score = as.numeric(X %*% w), w = w, variant = variant,
       cov_models = cov_models, uni_params = uni_params, mv_params = mv_params,
       marker_cols = marker_cols)
}

# Score new data from stored layer + combiner parameters (used by predict()).
.score_newdata <- function(fit, newdata) {
  d <- .cov_apply(newdata, fit$cov_models)
  if (fit$variant == "scalar") {
    X <- do.call(cbind, lapply(fit$markers, function(col)
      .uni_apply_marker(d, fit$uni_params[[col]], col,
                        time_col = fit$time_col, h_col = fit$time_col)))
  } else {
    X <- .mv_apply(d, fit$mv_params, fit$markers, time_col = fit$time_col)
  }
  as.numeric(X %*% fit$weights)
}
