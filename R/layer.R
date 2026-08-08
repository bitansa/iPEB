# Time-gap-aware PEB layer: covariate adjustment and innovation standardization.
# All functions here are internal. They are split into a "fit" stage (estimate
# parameters from healthy training controls) and an "apply" stage (score any
# subject's trajectory), so that a fitted model can score genuinely new data.

# Small Moore-Penrose pseudo-inverse (avoids a hard dependency on MASS).
.ginv <- function(X) {
  s <- svd(X)
  d <- ifelse(s$d > 1e-10 * max(s$d), 1 / s$d, 0)
  s$v %*% (d * t(s$u))
}

# --- Covariate adjustment (manuscript eq. 2 term x_it' beta_j) ---------------
# Fit one linear model per marker on the healthy fit controls; store the model
# and the healthy grand mean so the same adjustment can be applied to new data.
.cov_fit <- function(dat, fit_ids, marker_cols, cov_cols) {
  cov_cols <- intersect(cov_cols, names(dat))
  if (length(cov_cols) == 0L) return(NULL)
  hc <- dat$IDvar %in% fit_ids & dat$D == 0
  models <- list()
  for (m in marker_cols) {
    fml <- stats::as.formula(paste(m, "~", paste(cov_cols, collapse = " + ")))
    fit <- tryCatch(stats::lm(fml, data = dat[hc, , drop = FALSE]),
                    error = function(e) NULL)
    if (is.null(fit)) { models[[m]] <- NULL; next }
    models[[m]] <- list(model = fit, gm = mean(dat[[m]][hc], na.rm = TRUE))
  }
  if (length(models) == 0L) NULL else models
}

# Apply stored covariate models to (new) data.
.cov_apply <- function(dat, cov_models) {
  if (is.null(cov_models)) return(dat)
  for (m in names(cov_models)) {
    cm <- cov_models[[m]]
    if (is.null(cm)) next
    pred <- tryCatch(as.numeric(stats::predict(cm$model, newdata = dat)),
                     error = function(e) rep(NA_real_, nrow(dat)))
    adj <- dat[[m]] - pred + cm$gm
    bad <- !is.finite(adj)
    adj[bad] <- dat[[m]][bad]
    dat[[m]] <- adj
  }
  dat
}

# --- Univariate innovation layer (one marker) --------------------------------
# Fit: random intercept [+ optional slope on the time index] with AR(1)/OU
# residuals, estimated from healthy fit controls. Returns the parameters needed
# to score any subject by sequential one-step prediction.
.uni_fit <- function(dat, fit_ids, col, time_col = "t", h_col = "t",
                     use_slope = FALSE, ar1 = TRUE) {
  d <- data.frame(ipy = dat[[col]], ipt = dat[[time_col]], iph = dat[[h_col]],
                  IDvar = dat$IDvar, D = dat$D)
  fc <- d[d$IDvar %in% fit_ids & d$D == 0, ]
  hbar <- mean(fc$iph, na.rm = TRUE)
  fc$iph <- fc$iph - hbar

  ranform <- if (use_slope) ~ 1 + iph | IDvar else ~ 1 | IDvar
  m <- NULL
  if (ar1) {
    m <- tryCatch(nlme::lme(ipy ~ 1, random = ranform,
                            correlation = nlme::corCAR1(form = ~ ipt | IDvar),
                            data = fc, control = nlme::lmeControl(opt = "optim")),
                  error = function(e) NULL)
  }
  if (is.null(m)) {
    m <- tryCatch(nlme::lme(ipy ~ 1, random = ranform, data = fc,
                            control = nlme::lmeControl(opt = "optim")),
                  error = function(e) NULL)
  }
  if (is.null(m) && use_slope) {
    m <- tryCatch(nlme::lme(ipy ~ 1, random = ~ 1 | IDvar, data = fc,
                            control = nlme::lmeControl(opt = "optim")),
                  error = function(e) NULL)
  }

  if (is.null(m)) {                         # moment fallback
    mu <- mean(fc$ipy, na.rm = TRUE); sigma2 <- stats::var(fc$ipy, na.rm = TRUE)
    D <- matrix(0, 1, 1); phi <- 0; slope <- FALSE
  } else {
    mu <- as.numeric(nlme::fixef(m))
    sigma2 <- m$sigma^2
    D <- tryCatch(as.matrix(nlme::getVarCov(m, type = "random.effects")),
                  error = function(e) matrix(as.numeric(nlme::VarCorr(m)[1, 1]), 1, 1))
    slope <- use_slope && nrow(D) == 2
    phi <- 0
    if (ar1 && !is.null(m$modelStruct$corStruct)) {
      phi <- as.numeric(stats::coef(m$modelStruct$corStruct, unconstrained = FALSE))[1]
    }
    if (!is.finite(phi)) phi <- 0
    if (any(!is.finite(D))) D <- matrix(sigma2 / 2, nrow(D), ncol(D))
  }
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- 1
  list(mu = mu, sigma2 = sigma2, D = D, phi = phi, slope = slope, hbar = hbar)
}

# Apply: sequential one-step innovation R_ijt = (e_j - conditional mean) /
# sqrt(prediction variance) for one subject's trajectory of a single marker.
.uni_score_subject <- function(y, tt, h, p) {
  o <- order(tt)
  y <- y[o]; tt <- tt[o]; hh <- h[o] - p$hbar
  k <- length(y)
  Z <- if (p$slope) cbind(1, hh) else cbind(rep(1, k))
  Vr <- Z %*% p$D %*% t(Z)
  Dst <- abs(outer(tt, tt, "-"))
  Va <- p$sigma2 * (p$phi^Dst); diag(Va) <- p$sigma2
  V <- Vr + Va
  e <- y - p$mu
  R <- rep(NA_real_, k)
  for (j in seq_len(k)) {
    if (j == 1L) { R[1] <- e[1] / sqrt(max(V[1, 1], 1e-8)); next }
    S11 <- V[1:(j - 1), 1:(j - 1), drop = FALSE] + diag(1e-8, j - 1)
    s1 <- V[j, 1:(j - 1)]
    Si <- tryCatch(solve(S11), error = function(z) .ginv(S11))
    cm <- as.numeric(s1 %*% Si %*% e[1:(j - 1)])
    cv <- V[j, j] - as.numeric(s1 %*% Si %*% s1)
    R[j] <- (e[j] - cm) / sqrt(max(cv, 1e-8))
  }
  out <- rep(NA_real_, k); out[o] <- R; out
}

# Score every subject/row of `dat` for one marker using fitted parameters `p`.
.uni_apply_marker <- function(dat, p, col, time_col = "t", h_col = "t") {
  R <- rep(NA_real_, nrow(dat))
  for (id in unique(dat$IDvar)) {
    ix <- which(dat$IDvar == id)
    R[ix] <- .uni_score_subject(dat[[col]][ix], dat[[time_col]][ix],
                                dat[[h_col]][ix], p)
  }
  R
}

# --- Multivariate whitened innovation layer ----------------------------------
# Joint AR(1)/OU model on the marker vector; one-step predictive covariance
# whitened by its inverse square root. Estimation (.mv_fit) and per-subject
# whitening (.mv_whiten_subject) are already separable.
.mv_fit <- function(Y, ids, tt, ar1 = TRUE) {
  p <- ncol(Y); mu <- colMeans(Y, na.rm = TRUE)
  sm <- t(sapply(split(seq_len(nrow(Y)), ids),
                 function(ix) colMeans(Y[ix, , drop = FALSE], na.rm = TRUE)))
  Rres <- apply(Y, 2, function(c) c - stats::ave(c, ids, FUN = function(z) mean(z, na.rm = TRUE)))
  vmk <- apply(Y, 2, stats::var, na.rm = TRUE); vmk[!is.finite(vmk) | vmk <= 0] <- 1
  Tb <- suppressWarnings(stats::cov(sm));   if (any(!is.finite(Tb))) Tb <- diag(vmk / 2, p)
  Sw <- suppressWarnings(stats::cov(Rres)); if (any(!is.finite(Sw))) Sw <- diag(vmk / 2, p)
  Tb <- (Tb + t(Tb)) / 2 + diag(1e-8, p); Sw <- (Sw + t(Sw)) / 2 + diag(1e-8, p)
  phi <- 0
  if (ar1) {
    ph <- sapply(seq_len(p), function(j) {
      m <- tryCatch(nlme::lme(y ~ 1, random = ~ 1 | ID,
                              correlation = nlme::corCAR1(form = ~ tt | ID),
                              data = data.frame(y = Y[, j], ID = ids, tt = tt),
                              control = nlme::lmeControl(opt = "optim")),
                    error = function(e) NULL)
      if (is.null(m) || is.null(m$modelStruct$corStruct)) NA_real_
      else as.numeric(stats::coef(m$modelStruct$corStruct, unconstrained = FALSE))[1]
    })
    phi <- suppressWarnings(stats::median(ph, na.rm = TRUE)); if (!is.finite(phi)) phi <- 0
  }
  list(mu = mu, Tb = Tb, Sw = Sw, phi = phi, p = p)
}

.mv_whiten_subject <- function(Yi, tti, par) {
  p <- par$p; o <- order(tti); Y <- Yi[o, , drop = FALSE]; tt <- tti[o]; k <- nrow(Y)
  Tb <- par$Tb; Sw <- par$Sw; phi <- par$phi; mu <- par$mu
  big <- matrix(0, k * p, k * p)
  for (s in 1:k) for (u in 1:k) {
    blk <- if (s == u) Tb + Sw else Tb + Sw * (phi^abs(tt[s] - tt[u]))
    big[((s - 1) * p + 1):(s * p), ((u - 1) * p + 1):(u * p)] <- blk
  }
  Rw <- matrix(NA_real_, k, p)
  for (j in 1:k) {
    ij <- ((j - 1) * p + 1):(j * p)
    if (j == 1L) { nu <- Y[1, ] - mu; S <- big[ij, ij, drop = FALSE]
    } else {
      ih <- 1:((j - 1) * p); Shh <- big[ih, ih, drop = FALSE]
      Shh <- (Shh + t(Shh)) / 2 + diag(1e-6 * mean(diag(Shh)) + 1e-8, length(ih))
      Sjh <- big[ij, ih, drop = FALSE]
      yh <- as.numeric(t(Y[1:(j - 1), , drop = FALSE])) - rep(mu, j - 1)
      Ay <- tryCatch(solve(Shh, yh), error = function(e) rep(0, length(yh)))
      As <- tryCatch(solve(Shh, t(Sjh)), error = function(e) matrix(0, length(ih), p))
      nu <- Y[j, ] - (mu + as.numeric(Sjh %*% Ay)); S <- big[ij, ij, drop = FALSE] - Sjh %*% As
    }
    if (any(!is.finite(S)) || any(!is.finite(nu))) { S <- Tb + Sw; nu <- Y[j, ] - mu }
    S <- (S + t(S)) / 2; eg <- eigen(S, symmetric = TRUE)
    val <- pmax(eg$values, 1e-6 * max(eg$values, na.rm = TRUE))
    Rw[j, ] <- as.numeric(eg$vectors %*% (t(eg$vectors) / sqrt(val)) %*% nu)
  }
  out <- matrix(NA_real_, k, p); out[o, ] <- Rw; out
}

.mv_apply <- function(dat, par, marker_cols, time_col = "t") {
  Y <- as.matrix(dat[, marker_cols, drop = FALSE]); tt <- dat[[time_col]]; ids <- dat$IDvar
  R <- matrix(NA_real_, nrow(dat), par$p)
  for (id in unique(ids)) {
    ix <- which(ids == id)
    R[ix, ] <- .mv_whiten_subject(Y[ix, , drop = FALSE], tt[ix], par)
  }
  R
}
