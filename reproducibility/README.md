# Reproducibility

Scripts that reproduce the analyses in the iPEB paper using the released
package. They are **not** part of the CRAN package (excluded via
`.Rbuildignore`). Install the package first:

```r
remotes::install_github("bitansa/iPEB")
```

## Simulations

`simulations/reproduce_leadtime.R` is fully self-contained: it simulates a
differential-timing cohort and, using only `library(iPEB)`, shows that an
early-rising marker alone gives long lead but lower sensitivity, a late-rising
marker gives high sensitivity but short lead, and iPEB combines both — with the
lead-time objective extending lead further. Run:

```r
source("simulations/reproduce_leadtime.R")
```

## Lung cohort

The PLCO data are controlled-access through the NCI Cancer Data Access System
(CDAS) under its standard data-use agreement and are **not** included here, so
this folder ships a **synthetic** six-marker longitudinal cohort with the same
structure (centers, covariates, irregular visits, four informative markers plus
two weaker ones). Build it, then run the analysis:

```r
source("make_synthetic_lung.R")   # writes synthetic_lung.csv (entirely simulated)
source("reproduce_lung.R")        # four- / six- / selected-four x three objectives
```

`reproduce_lung.R` runs the whole workflow through `library(iPEB)` — the frozen
logistic panels via `ipeb_innovations()`, the iPEB arm via `ipeb()` with backward
selection — mirroring the structure of Table 2 and Web Tables 3-6. On the
synthetic data the numbers are illustrative only; to reproduce the paper's actual
numbers, point `lung` at your own approved CDAS extract with the same columns.
