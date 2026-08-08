# Prepare new data for scoring/evaluation using a fitted model's schema.
.prep_newdata <- function(object, newdata, need_labels) {
  if (!is.data.frame(newdata)) stop("`newdata` must be a data.frame.", call. = FALSE)
  req <- c(object$id, object$time, object$time_to_dx, object$markers, object$covariates)
  if (need_labels) req <- c(req, object$case)
  miss <- setdiff(req, names(newdata))
  if (length(miss)) stop("Columns not found in `newdata`: ", paste(miss, collapse = ", "), call. = FALSE)
  d <- as.data.frame(newdata)
  d$IDvar <- d[[object$id]]
  d[[object$time_col]] <- as.numeric(d[[object$time_col]])
  d[[object$ttd_col]] <- as.numeric(d[[object$ttd_col]])
  for (m in object$markers) d[[m]] <- as.numeric(d[[m]])
  if (object$case %in% names(d)) d$D <- as.integer(d[[object$case]]) else d$D <- NA_integer_
  d
}

#' Score new subjects with a fitted iPEB model
#'
#' Computes the iPEB composite score for each visit in \code{newdata} using the
#' layer parameters and combiner weights stored in a fitted model. A subject is
#' scored from its own trajectory against the training healthy reference; higher
#' scores indicate stronger evidence of disease.
#'
#' @param object A fitted \code{"ipeb"} object from \code{\link{ipeb}}.
#' @param newdata A data frame with the same columns used at fitting (the
#'   selected markers, covariates, id, time, and time-to-diagnosis; the case
#'   column is optional here).
#' @param ... Unused.
#'
#' @return A numeric vector of iPEB scores, one per row of \code{newdata}.
#'
#' @examples
#' data(ipeb_example)
#' train <- subset(ipeb_example, split == "train")
#' test  <- subset(ipeb_example, split == "test")
#' \donttest{
#' fit <- ipeb(train, markers = c("m1", "m2", "m3"),
#'             innovation = "iid", slope = "off")
#' scores <- predict(fit, test)
#' }
#'
#' @seealso \code{\link{ipeb}}, \code{\link{evaluate}}
#' @method predict ipeb
#' @export
predict.ipeb <- function(object, newdata, ...) {
  dat <- .prep_newdata(object, newdata, need_labels = FALSE)
  .score_newdata(object, dat)
}
