# Reproducibility

Scripts that reproduce the analyses in the iPEB paper using the released
package. They are **not** part of the CRAN package (excluded via
`.Rbuildignore`). Install the package first:

```r
remotes::install_github("bitansa/iPEB")
```

## Simulations

Five self-contained scripts in `simulations/` reproduce the paper's simulation
results using only `library(iPEB)` (shared scaffolding in `sim_helpers.R`).
Each writes its summary numbers to a `*_means.csv` and its figure(s) to
300-dpi PNGs alongside the script:

- `reproduce_S1.R` — Scenario 1 (marker heterogeneity): sensitivity/AUC rows
  of Table 1 plus the mean-ROC figure (Web Figure 1).
- `reproduce_S2.R` — Scenario 2 (differential timing): an early-rising versus a
  late-rising marker; iPEB under the lead-time objective combines long lead
  with high sensitivity (Table 1 rows and the lead-time/trajectory figures).
- `reproduce_S3.R` — Scenario 3 (longitudinal drift + serial correlation):
  value of the time-gap-aware layer (intercept/i.i.d. versus slope + OU).
- `reproduce_S4.R` — Scenario 4 (irregular visit spacing): the S3 mechanism on
  an irregular 3-month-to-2-year visit grid; gap-ignorant PEB versus gap-aware
  iPEB (Web Figure 2).
- `reproduce_factorial.R` — the 48-cell factorial study (Web Tables 1-2 and
  the factorial summary figure).

```r
setwd("simulations")
source("reproduce_S1.R")   # likewise S2, S3, S4, factorial
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
