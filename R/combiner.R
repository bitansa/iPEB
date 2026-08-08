# Learn the marker combiner (unit-norm weight vector) by the manuscript's two
# routes. Internal.
#
# ROUTE 1 (sensitivity objective, lam2 = lam3 = 0): the optimal combiner is the
#   closed form  w* proportional to Sigma_R^{-1} mu_Delta, used directly.
# ROUTE 2 (lead-time / combined objective): minimize the scalarized loss by
#   gradient-free search, started at the closed form and kept only if it lowers
#   the tuning loss. For the multivariate variant the innovations are already
#   whitened, so Sigma_R = I.
.learn_combiner <- function(X, D, ids, ttd, alpha, window,
                            lam1, lam2, lam3, whitened = FALSE) {
  ctrl <- D == 0
  Sig <- if (whitened) diag(ncol(X)) else stats::cov(X[ctrl, , drop = FALSE])
  mu <- colMeans(X[D == 1, , drop = FALSE], na.rm = TRUE) -
        colMeans(X[ctrl, , drop = FALSE], na.rm = TRUE)
  w0 <- tryCatch(as.numeric(solve(Sig + diag(1e-6, ncol(X)), mu)),
                 error = function(e) mu)
  w0 <- w0 / sqrt(sum(w0^2))
  if (lam2 == 0 && lam3 == 0) return(w0)                 # ROUTE 1: closed form
  obj <- function(w) {                                   # ROUTE 2: search
    wn <- w / sqrt(sum(w^2))
    .loss(as.numeric(X %*% wn), D, ids, ttd,
          as.numeric(X[ctrl, , drop = FALSE] %*% wn),
          alpha, window, lam1, lam2, lam3)
  }
  fit <- tryCatch(stats::optim(w0, obj, method = "Nelder-Mead",
                               control = list(maxit = 300)),
                  error = function(e) NULL)
  w <- if (!is.null(fit) && fit$value <= obj(w0)) fit$par else w0
  w / sqrt(sum(w^2))
}
