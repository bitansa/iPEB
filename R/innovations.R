#' Time-gap-aware standardized innovations
#'
#' Returns the innovation-standardized residuals from the iPEB layer: for each
#' marker, a random intercept (and optional slope) model with autocorrelated,
#' gap-scaled residuals is estimated on the healthy reference subjects, and every
#' visit is scored by sequential one-step prediction (history to innovation,
#' divided by its prediction standard deviation). The prediction variance grows
#' with the gap since the previous visit, so the innovations are comparable
#' across irregular schedules.
#'
#' This is the layer that underlies \code{\link{ipeb}}. It is exposed so that
#' history-adjusted baselines (for example, a fixed-panel PEB score, a
#' single-marker PEB score, or a principal component of the innovations) can be
#' built directly from the same layer.
#'
#' @param data A long-format data frame, one row per subject-visit.
#' @param markers Character vector of biomarker column names to standardize.
#' @param reference_ids Optional vector of subject ids whose healthy (control)
#'   visits define the reference for estimating the layer. Defaults to all
#'   subjects in \code{data}.
#' @param id,case,time Column names for the subject id, the case indicator
#'   (1 = case, 0 = control), and the visit time.
#' @param covariates Character vector of covariate columns to adjust for.
#' @param slope Random slope: \code{"off"} (default) or \code{"on"}.
#' @param innovation Residual model: \code{"iid"} (default) or \code{"ar1"} for
#'   gap-scaled AR(1)/OU residuals.
#'
#' @return A numeric matrix with one row per row of \code{data} and one column
#'   per marker (named by marker), containing the standardized innovations.
#'
#' @examples
#' data(ipeb_example)
#' train <- subset(ipeb_example, split == "train")
#' \donttest{
#' R <- ipeb_innovations(train, markers = c("m1", "m2", "m3"))
#' head(R)
#' }
#'
#' @seealso \code{\link{ipeb}}
#' @export
ipeb_innovations <- function(data, markers, reference_ids = NULL,
                             id = "id", case = "case", time = "time",
                             covariates = character(0),
                             slope = c("off", "on"),
                             innovation = c("iid", "ar1")) {
  slope <- match.arg(slope)
  innovation <- match.arg(innovation)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  req <- c(id, case, time, markers, covariates)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)

  d <- as.data.frame(data)
  d$IDvar <- d[[id]]
  d$D <- as.integer(d[[case]])
  d[[time]] <- as.numeric(d[[time]])
  for (m in markers) d[[m]] <- as.numeric(d[[m]])
  if (is.null(reference_ids)) reference_ids <- unique(d$IDvar)

  use_slope <- slope == "on"
  ar1 <- innovation == "ar1"

  cov_models <- .cov_fit(d, reference_ids, markers, covariates)
  d <- .cov_apply(d, cov_models)
  X <- do.call(cbind, lapply(markers, function(col) {
    p <- .uni_fit(d, reference_ids, col, time_col = time, h_col = time,
                  use_slope = use_slope, ar1 = ar1)
    .uni_apply_marker(d, p, col, time_col = time, h_col = time)
  }))
  colnames(X) <- markers
  X
}
