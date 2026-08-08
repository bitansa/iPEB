#' @export
print.ipeb <- function(x, ...) {
  variant <- if (x$variant == "scalar") "scalar (iPEB-S)" else "multivariate (iPEB-M)"
  win <- if (is.infinite(x$window)) "whole trajectory" else paste(x$window, "months")
  layer <- paste0(if (x$slope) "random slope" else "intercept only", ", ",
                  if (x$ar1) "AR(1)/OU innovations" else "i.i.d. innovations")
  cat("iPEB model\n")
  cat("  Objective:          ", x$objective, "\n", sep = "")
  cat("  Combiner:           ", variant, "\n", sep = "")
  cat("  Operating spec (a): ", x$alpha, "\n", sep = "")
  cat("  Detection window:   ", win, "\n", sep = "")
  cat("  Layer:              ", layer, "\n", sep = "")
  if (length(x$covariates)) cat("  Covariates:         ", paste(x$covariates, collapse = ", "), "\n", sep = "")
  cat("  Markers (", length(x$markers), "): ", paste(x$markers, collapse = ", "), "\n", sep = "")
  w <- round(x$weights, 3)
  cat("  Weights:\n")
  for (m in names(w)) cat("    ", format(m, width = max(nchar(names(w)))), "  ", formatC(w[m], format = "f", digits = 3), "\n", sep = "")
  cat("  Training subjects:  ", x$n_train_subjects,
      " (", x$n_train_cases, " cases / ", x$n_train_controls, " controls)\n", sep = "")
  invisible(x)
}

#' @export
summary.ipeb <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Plot the fitted iPEB marker weights
#'
#' @param x A fitted \code{"ipeb"} object.
#' @param ... Passed to \code{\link[graphics]{barplot}}.
#' @return The fitted object, invisibly.
#' @method plot ipeb
#' @export
plot.ipeb <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(7, 4.5, 3, 1))
  graphics::barplot(x$weights, las = 2, ylab = "iPEB weight",
                    main = paste0("iPEB weights (", x$objective, " objective)"), ...)
  graphics::abline(h = 0, col = "grey60")
  invisible(x)
}
