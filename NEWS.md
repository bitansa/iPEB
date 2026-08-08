# iPEB 0.1.0

* First release.
* `ipeb()` fits the improved Parametric Empirical Bayes model on training data:
  a time-gap-aware standardization layer (random intercept, optional random
  slope, optional AR(1)/OU residual autocorrelation), optional covariate
  adjustment, objective-driven weighting (sensitivity, lead-time, or combined
  objective), optional feature selection, and an automatically chosen scalar or
  multivariate combiner.
* `predict()` scores new subjects, and `evaluate()` reports per-patient
  sensitivity and lead time with per-visit specificity at user-chosen operating
  points.
* `ipeb_run()` provides a one-call fit-and-evaluate wrapper.
* `ipeb_innovations()` exposes the time-gap-aware standardization layer, so
  history-adjusted baselines can be built from the same layer as iPEB.
* `print()`, `summary()`, and `plot()` methods for fitted `ipeb` objects.
* Ships a small synthetic longitudinal example dataset, `ipeb_example`.
