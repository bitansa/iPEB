#' Fit and evaluate iPEB in one call
#'
#' Convenience wrapper that fits \code{\link{ipeb}} on \code{train} and
#' evaluates the fitted model on \code{test}.
#'
#' @inheritParams ipeb
#' @param train Training data frame.
#' @param test Test data frame (same schema as \code{train}).
#' @param specificities Numeric vector of target specificities for evaluation.
#'   Defaults to \code{alpha}.
#' @param eval_window Optional evaluation detection window in months (passed to
#'   \code{\link{evaluate}}); \code{NULL} reuses the training window.
#'
#' @return A list with two elements: \code{fit} (the fitted \code{"ipeb"}
#'   object) and \code{evaluation} (the data frame returned by
#'   \code{\link{evaluate}} on \code{test}).
#'
#' @examples
#' data(ipeb_example)
#' train <- subset(ipeb_example, split == "train")
#' test  <- subset(ipeb_example, split == "test")
#' \donttest{
#' res <- ipeb_run(train, test, markers = c("m1", "m2", "m3"),
#'                 objective = "sensitivity", innovation = "iid", slope = "off",
#'                 specificities = c(0.90, 0.95, 0.99))
#' res$evaluation
#' }
#'
#' @seealso \code{\link{ipeb}}, \code{\link{evaluate}}
#' @param seed Optional single number passed to \code{\link{ipeb}} to fix its
#'   internal validation split, making the fit exactly reproducible.
#'
#' @export
ipeb_run <- function(train, test, markers, id = "id", case = "case",
                     time = "time", time_to_dx = "time_to_dx",
                     covariates = character(0),
                     objective = c("sensitivity", "leadtime", "combined"),
                     alpha = 0.95, window = Inf,
                     slope = c("auto", "on", "off"),
                     innovation = c("auto", "ar1", "iid"),
                     select = c("none", "backward"),
                     n_markers = NULL, validation_frac = 0.25,
                     specificities = NULL, eval_window = NULL, seed = NULL) {
  fit <- ipeb(train, markers = markers, id = id, case = case, time = time,
              time_to_dx = time_to_dx, covariates = covariates,
              objective = objective, alpha = alpha, window = window,
              slope = slope, innovation = innovation, select = select,
              n_markers = n_markers, validation_frac = validation_frac, seed = seed)
  ev <- evaluate(fit, test, specificities = specificities, window = eval_window)
  list(fit = fit, evaluation = ev)
}
