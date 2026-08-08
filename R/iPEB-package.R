#' iPEB: Improved Parametric Empirical Bayes for Longitudinal Biomarker Analysis
#'
#' The iPEB package fits the improved Parametric Empirical Bayes model, which
#' extends parametric empirical Bayes with a time-gap-aware standardization
#' layer, covariate adjustment, and objective-driven multi-marker weighting.
#'
#' The main entry point is \code{\link{ipeb}}, which fits the model on training
#' data. Fitted objects have \code{\link{predict.ipeb}} and
#' \code{\link{evaluate}} methods for scoring and evaluating new subjects, and
#' \code{\link{ipeb_run}} provides a one-call fit-and-evaluate wrapper.
#'
#' @section Data format:
#' Functions expect a long-format \code{data.frame} with one row per
#' (subject, visit): a subject identifier column, a case/control indicator
#' (1 = case, 0 = control), one column per biomarker, a visit-time column, and a
#' column giving each visit's time before diagnosis (used to score lead time).
#' Column names are supplied through the fitting arguments.
#'
#' @keywords internal
#' @importFrom nlme lme corCAR1 lmeControl fixef VarCorr getVarCov
#' @importFrom stats as.formula ave coef cov lm median optim predict quantile var
#' @importFrom graphics abline axis barplot legend lines mtext par text
#' @importFrom utils combn
"_PACKAGE"
