#' Evaluate a fitted iPEB model
#'
#' Computes per-patient sensitivity and median lead time with per-visit
#' specificity at one or more operating points. Operating thresholds are
#' calibrated on the training control visits stored in the fitted object and
#' applied unchanged to \code{newdata}, so every reported number is
#' out-of-sample when \code{newdata} are held-out subjects.
#'
#' @param object A fitted \code{"ipeb"} object.
#' @param ... Passed to methods.
#'
#' @return A data frame with one row per requested specificity and columns
#'   \code{specificity} (target), \code{sensitivity}, \code{lead_time} (median,
#'   years), \code{realized_specificity} (achieved on \code{newdata}), and
#'   \code{auc} (subject-level, specificity-independent).
#'
#' @examples
#' data(ipeb_example)
#' train <- subset(ipeb_example, split == "train")
#' test  <- subset(ipeb_example, split == "test")
#' \donttest{
#' fit <- ipeb(train, markers = c("m1", "m2", "m3"),
#'             innovation = "iid", slope = "off")
#' evaluate(fit, test, specificities = c(0.90, 0.95, 0.99))
#' }
#'
#' @export
evaluate <- function(object, ...) UseMethod("evaluate")

#' @param newdata A data frame with the fitting schema, including the case
#'   indicator (needed to score sensitivity and lead time).
#' @param specificities Numeric vector of target specificities. Defaults to the
#'   specificity the model was optimized at.
#' @param window Optional evaluation detection window in months. If \code{NULL}
#'   (default), the window used at fitting is reused; supply a value (or
#'   \code{Inf} for the whole trajectory) to report at a different window than
#'   training.
#' @rdname evaluate
#' @method evaluate ipeb
#' @export
evaluate.ipeb <- function(object, newdata, specificities = NULL, window = NULL, ...) {
  if (is.null(specificities)) specificities <- object$alpha
  if (!is.numeric(specificities) || any(specificities <= 0 | specificities >= 1))
    stop("`specificities` must be strictly between 0 and 1.", call. = FALSE)
  win <- if (is.null(window)) object$window else window
  if (!is.numeric(win) || win <= 0)
    stop("`window` must be a positive number of months (or Inf).", call. = FALSE)
  dat <- .prep_newdata(object, newdata, need_labels = TRUE)
  score <- .score_newdata(object, dat)
  auc <- .auc_subject(score, dat$D, dat$IDvar)
  rows <- lapply(specificities, function(s) {
    mt <- .metrics(score, dat$D, dat$IDvar, dat[[object$ttd_col]],
                   object$train_ctrl_scores, s, win)
    data.frame(specificity = s, sensitivity = unname(mt["SensW"]),
               lead_time = unname(mt["LT"]),
               realized_specificity = unname(mt["Spec"]), auc = auc)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
