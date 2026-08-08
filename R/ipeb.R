# Map user column names onto the internal layout and validate inputs.
.standardize_data <- function(data, id, case, time, time_to_dx, markers, covariates) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  req <- c(id, case, time, time_to_dx, markers, covariates)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  if (length(markers) < 1L) stop("At least one marker column is required.", call. = FALSE)
  d <- as.data.frame(data)
  d$IDvar <- d[[id]]
  cs <- d[[case]]
  if (is.logical(cs)) cs <- as.integer(cs)
  if (!all(stats::na.omit(cs) %in% c(0, 1))) stop("`case` column must be 0/1 (or logical).", call. = FALSE)
  d$D <- as.integer(cs)
  d[[time]] <- as.numeric(d[[time]])
  d[[time_to_dx]] <- as.numeric(d[[time_to_dx]])
  for (m in markers) d[[m]] <- as.numeric(d[[m]])
  d
}

#' Fit an improved Parametric Empirical Bayes (iPEB) model
#'
#' Fits the iPEB model on longitudinal biomarker data: a time-gap-aware
#' standardization layer, optional covariate adjustment, objective-driven
#' multi-marker weighting, optional feature selection, and an automatically
#' chosen scalar or multivariate combiner. The returned object can score and
#' evaluate new subjects with \code{\link{predict.ipeb}} and
#' \code{\link{evaluate}}.
#'
#' @param data A long-format data frame with one row per subject-visit.
#' @param markers Character vector of biomarker column names.
#' @param id,case,time,time_to_dx Column names for the subject identifier, the
#'   case indicator (1 = case, 0 = control), the visit time, and each visit's
#'   time before diagnosis (in days; used to measure lead time).
#' @param covariates Character vector of covariate column names to adjust for
#'   (default none).
#' @param objective Clinical objective: \code{"sensitivity"} (default),
#'   \code{"leadtime"}, or \code{"combined"}. Each fixes an internal penalty
#'   profile rather than requiring manual tuning.
#' @param alpha Operating specificity at which the combiner is optimized
#'   (default \code{0.95}).
#' @param window Optional detection window in months restricting the objective
#'   to a fixed pre-diagnostic horizon. The default, \code{Inf}, uses the whole
#'   pre-diagnostic trajectory (no window).
#' @param slope Random-slope layer: \code{"auto"} (default; on when subjects
#'   have enough visits), \code{"on"}, or \code{"off"}.
#' @param innovation Residual model: \code{"auto"} (default) or \code{"ar1"}
#'   use gap-scaled AR(1)/OU innovations; \code{"iid"} uses independent
#'   innovations.
#' @param select Feature selection by objective-driven backward elimination:
#'   \code{"none"} (default) or \code{"backward"}.
#' @param n_markers Optional target panel size when \code{select} is not
#'   \code{"none"}.
#' @param validation_frac Fraction of training subjects held out as an internal
#'   validation slice to choose the combiner variant and drive selection
#'   (default \code{0.25}).
#'
#' @return An object of class \code{"ipeb"}: a list with the fitted weights,
#'   selected markers, chosen variant, layer parameters, and the training
#'   control scores used to calibrate operating thresholds.
#'
#' @examples
#' data(ipeb_example)
#' train <- subset(ipeb_example, split == "train")
#' \donttest{
#' fit <- ipeb(train, markers = c("m1", "m2", "m3"),
#'             id = "id", case = "case", time = "time",
#'             time_to_dx = "time_to_dx", objective = "sensitivity",
#'             innovation = "iid", slope = "off")
#' fit
#' }
#'
#' @seealso \code{\link{predict.ipeb}}, \code{\link{evaluate}}, \code{\link{ipeb_run}}
#' @export
ipeb <- function(data, markers, id = "id", case = "case", time = "time",
                 time_to_dx = "time_to_dx", covariates = character(0),
                 objective = c("sensitivity", "leadtime", "combined"),
                 alpha = 0.95, window = Inf,
                 slope = c("auto", "on", "off"),
                 innovation = c("auto", "ar1", "iid"),
                 select = c("none", "backward"),
                 n_markers = NULL, validation_frac = 0.25) {
  objective <- match.arg(objective)
  slope <- match.arg(slope)
  innovation <- match.arg(innovation)
  select <- match.arg(select)
  obj_internal <- if (objective == "sensitivity") "sens" else objective
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be strictly between 0 and 1.", call. = FALSE)

  dat <- .standardize_data(data, id, case, time, time_to_dx, markers, covariates)
  time_col <- time; ttd_col <- time_to_dx
  train_ids <- unique(dat$IDvar)

  visits <- as.numeric(table(dat$IDvar))
  use_slope <- switch(slope, on = TRUE, off = FALSE,
                      auto = stats::median(visits) >= 5)
  ar1 <- switch(innovation, ar1 = TRUE, iid = FALSE, auto = TRUE)
  lam <- .lambda_profile(obj_internal)
  max_keep <- if (is.null(n_markers)) Inf else n_markers

  marker_cols <- markers
  if (select != "none") {
    marker_cols <- .select_markers(dat, train_ids, marker_cols, obj_internal,
      covariates, alpha, window, use_slope, ar1, time_col, ttd_col,
      validation_frac, lam, direction = select, max_keep = max_keep)
  }

  # Choose the scalar/multivariate variant on an internal validation slice.
  isc <- tapply(dat$D, dat$IDvar, max)
  trc <- intersect(train_ids, names(isc)[isc == 1])
  trk <- intersect(train_ids, names(isc)[isc == 0])
  f <- c(sample(trc, round((1 - validation_frac) * length(trc))),
         sample(trk, round((1 - validation_frac) * length(trk))))
  iv <- dat$IDvar %in% setdiff(train_ids, f)
  variants <- if (length(marker_cols) >= 2L) c("scalar", "mv") else "scalar"
  vscore <- sapply(variants, function(v) {
    s <- .fit_layer_and_combiner(dat, fit_ids = f, tune_ids = f,
           marker_cols = marker_cols, variant = v, objective = obj_internal,
           cov_cols = covariates, alpha = alpha, window = window,
           use_slope = use_slope, ar1 = ar1, time_col = time_col,
           ttd_col = ttd_col, lam = lam)$score
    .metrics(s[iv], dat$D[iv], dat$IDvar[iv], dat[[ttd_col]][iv],
             s[iv & dat$D == 0], alpha, window)["SensW"]
  })
  best <- variants[which.max(vscore)]

  # Refit the winning variant on the full training set and store parameters.
  fit <- .fit_layer_and_combiner(dat, fit_ids = train_ids, tune_ids = train_ids,
    marker_cols = marker_cols, variant = best, objective = obj_internal,
    cov_cols = covariates, alpha = alpha, window = window,
    use_slope = use_slope, ar1 = ar1, time_col = time_col, ttd_col = ttd_col,
    lam = lam)

  structure(list(
    call = match.call(),
    objective = objective, objective_internal = obj_internal,
    alpha = alpha, window = window, lambda = lam,
    slope = use_slope, ar1 = ar1, select = select,
    variant = best, weights = stats::setNames(fit$w, marker_cols),
    markers = marker_cols, markers_input = markers, covariates = covariates,
    cov_models = fit$cov_models, uni_params = fit$uni_params,
    mv_params = fit$mv_params,
    id = id, case = case, time = time, time_to_dx = time_to_dx,
    time_col = time_col, ttd_col = ttd_col,
    train_ctrl_scores = fit$score[dat$D == 0],
    n_train_subjects = length(train_ids),
    n_train_cases = sum(isc == 1), n_train_controls = sum(isc == 0)
  ), class = "ipeb")
}
