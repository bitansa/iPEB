test_that("ipeb fits, predicts, and returns a well-formed object", {
  skip_on_cran()
  data(ipeb_example)
  tr <- subset(ipeb_example, split == "train")
  te <- subset(ipeb_example, split == "test")
  fit <- ipeb(tr, markers = c("m1", "m2", "m3"), objective = "sensitivity",
              innovation = "iid", slope = "off")
  expect_s3_class(fit, "ipeb")
  expect_length(fit$weights, 3)
  expect_true(fit$variant %in% c("scalar", "mv"))
  sc <- predict(fit, te)
  expect_length(sc, nrow(te))
  expect_true(all(is.finite(sc)))
})

test_that("evaluate returns one row per specificity with valid ranges", {
  skip_on_cran()
  data(ipeb_example)
  tr <- subset(ipeb_example, split == "train")
  te <- subset(ipeb_example, split == "test")
  fit <- ipeb(tr, markers = c("m1", "m2", "m3"), objective = "sensitivity",
              innovation = "iid", slope = "off")
  ev <- evaluate(fit, te, specificities = c(0.90, 0.95))
  expect_equal(nrow(ev), 2L)
  expect_true(all(c("specificity", "sensitivity", "lead_time",
                    "realized_specificity", "auc") %in% names(ev)))
  expect_true(all(ev$sensitivity >= 0 & ev$sensitivity <= 1))
})

test_that("ipeb_run wraps fit and evaluate", {
  skip_on_cran()
  data(ipeb_example)
  tr <- subset(ipeb_example, split == "train")
  te <- subset(ipeb_example, split == "test")
  res <- ipeb_run(tr, te, markers = c("m1", "m2", "m3"),
                  objective = "sensitivity", innovation = "iid", slope = "off",
                  specificities = 0.95)
  expect_named(res, c("fit", "evaluation"))
  expect_s3_class(res$fit, "ipeb")
  expect_equal(nrow(res$evaluation), 1L)
})

test_that("feature selection reduces the panel to the target size", {
  skip_on_cran()
  data(ipeb_example)
  tr <- subset(ipeb_example, split == "train")
  fit <- ipeb(tr, markers = c("m1", "m2", "m3"), objective = "sensitivity",
              innovation = "iid", slope = "off", select = "backward", n_markers = 2)
  expect_lte(length(fit$markers), 2L)
})
