# iPEB

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

## Installation

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

## Data format

A long data frame with one row per subject-visit: a subject id, a case
indicator (`1` = case, `0` = control), a visit-time column, a column giving each
visit's days-to-diagnosis, and one column per biomarker. See `?ipeb_example`.

## Data availability

The lung analysis in the accompanying paper uses PLCO data, which are
controlled-access through the NCI Cancer Data Access System (CDAS) under its
standard data-use agreement and are not distributed with this package. The
bundled `ipeb_example` dataset is entirely synthetic.

## License

MIT © iPEB authors. A manuscript describing the method is in preparation.
