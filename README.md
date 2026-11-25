# iPEB

<!-- badges: start -->
[![R-CMD-check](https://github.com/bitansa/iPEB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bitansa/iPEB/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/iPEB)](https://CRAN.R-project.org/package=iPEB)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/iPEB)](https://CRAN.R-project.org/package=iPEB)
<!-- badges: end -->

**Improved Parametric Empirical Bayes for Longitudinal Biomarker Analysis**

iPEB extends parametric empirical Bayes (PEB) for longitudinal biomarker
screening with three ingredients:

1. a **time-gap-aware standardization layer** — a random intercept (and optional
   slope) with autocorrelated, gap-scaled residuals, so prediction uncertainty
   grows with the time between visits and per-visit specificity is preserved
   under irregular sampling;
2. **objective-driven multi-marker weighting** — weights are learned to optimize
   a user-selected clinical objective rather than to tune opaque penalties; and
3. optional **covariate adjustment** and **feature selection**, with an
   automatically chosen scalar or multivariate combiner.

## Try it in your browser

A hosted, browser-based version of the iPEB explorer is available at
**https://bitansa.github.io/iPEB/** — no R installation required. The app runs R
entirely inside your own browser (WebAssembly), so any CSV you load stays on
your machine and is never uploaded to a server.

## Installation

iPEB is on CRAN:

```r
install.packages("iPEB")
```

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("bitansa/iPEB")
```

## Quick start

```r
library(iPEB)

data(ipeb_example)
train <- subset(ipeb_example, split == "train")
test  <- subset(ipeb_example, split == "test")

# Fit for the sensitivity objective at 95% specificity
fit <- ipeb(train, markers = c("m1", "m2", "m3"),
            id = "id", case = "case", time = "time", time_to_dx = "time_to_dx",
            objective = "sensitivity", alpha = 0.95)
fit

# Score and evaluate held-out subjects
scores <- predict(fit, test)
evaluate(fit, test, specificities = c(0.90, 0.95, 0.99))

# One-call fit-and-evaluate
ipeb_run(train, test, markers = c("m1", "m2", "m3"),
         objective = "leadtime", specificities = c(0.90, 0.95, 0.99))
```

## Recovering conventional PEB

Conventional (marginal) PEB is a special case of iPEB: an intercept-only,
i.i.d. layer on a single marker (no gap-scaling, no multi-marker weighting).
Running it through the same package makes iPEB-versus-PEB comparisons exactly
apples-to-apples -- identical layer and thresholding, only the modeling options
differ. To score a plain PEB baseline for one marker:

```r
# Conventional PEB: one marker, i.i.d. innovations, intercept-only mean model.
peb <- ipeb(train, markers = c("m1"),
            id = "id", case = "case", time = "time", time_to_dx = "time_to_dx",
            objective = "sensitivity", alpha = 0.95,
            slope = "off", innovation = "iid")
evaluate(peb, test, specificities = c(0.90, 0.95, 0.99))
```

Turning on `innovation = "ar1"` and/or `slope = "on"` for the same data recovers
the time-gap-aware iPEB layer, so PEB and iPEB can be compared directly (this is
exactly how the baselines are computed in the reproducibility scripts). A fixed
multi-marker PEB composite (e.g. a published panel) can be scored the same way by
passing that composite as a single marker.

## Key options

| Argument | What it controls |
|----------|------------------|
| `objective` | `"sensitivity"`, `"leadtime"`, or `"combined"` (each fixes a penalty profile) |
| `alpha` | operating specificity for weight optimization |
| `window` | optional detection window in months; default `Inf` uses the whole trajectory |
| `slope` | random slope `"auto"`/`"on"`/`"off"` |
| `innovation` | `"auto"`/`"ar1"` (gap-scaled AR(1)/OU) or `"iid"` |
| `covariates` | columns to adjust for |
| `select`, `n_markers` | optional feature selection to a target panel size |

Sensitivity and lead time are scored per patient over the whole pre-diagnostic
trajectory; specificity is scored per visit and calibrated on the training
controls, then applied unchanged to new data.

## Shiny app

A point-and-click interface to the same engine ships with the package:

```r
iPEB::run_app()
```

It exposes the objective, operating specificity, slope and innovation switches,
covariates, optional feature selection, and evaluation specificities, and runs
`ipeb()`/`evaluate()` under the hood.

For users without R, the same app is hosted at
**https://bitansa.github.io/iPEB/**. It is built with
[shinylive](https://posit-dev.github.io/r-shinylive/), which compiles the app to
WebAssembly: R runs in the visitor's browser rather than on a server, so uploaded
data never leaves their machine. The first load downloads the R runtime and takes
roughly 30-60 seconds; subsequent visits are cached and fast.

## Data format

A long data frame with one row per subject-visit: a subject id, a case
indicator (`1` = case, `0` = control), a visit-time column, a column giving each
visit's days-to-diagnosis, and one column per biomarker. See `?ipeb_example`.

## Data availability

The lung analysis in the accompanying paper uses PLCO data, which are
controlled-access through the NCI Cancer Data Access System (CDAS) under its
standard data-use agreement and are not distributed with this package. The
bundled `ipeb_example` dataset is entirely synthetic.

## Citation

If you use iPEB, please cite the accompanying paper:

Sarkar, B., Kenney, A. M., Long, J. P., Fahrmann, J. F., Hanash, S., Do, K.-A., and
Irajizad, E. (2026). Multi-Objective Composite Longitudinal Biomarker Scores for Improved
Cancer Risk Assessment. arXiv preprint
[arXiv:2609.08225](https://arxiv.org/abs/2609.08225).

```bibtex
@article{sarkar2026ipeb,
  title   = {Multi-Objective Composite Longitudinal Biomarker Scores for Improved Cancer Risk Assessment},
  author  = {Sarkar, Bitan and Kenney, Ana Maria and Long, James P. and Fahrmann, Johannes F. and
             Hanash, Samir and Do, Kim-Anh and Irajizad, Ehsan},
  journal = {arXiv preprint arXiv:2609.08225},
  year    = {2026},
  doi     = {10.48550/arXiv.2609.08225}
}
```

## License

MIT © iPEB authors.
