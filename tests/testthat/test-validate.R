test_that("missing marker columns raise an informative error", {
  data(ipeb_example)
  expect_error(ipeb(ipeb_example, markers = "not_a_marker"), "not found")
})

test_that("alpha outside (0, 1) is rejected", {
  data(ipeb_example)
  expect_error(
    ipeb(ipeb_example, markers = "m1", alpha = 1.5),
    "between 0 and 1")
})

test_that("a non-0/1 case column is rejected", {
  data(ipeb_example)
  bad <- ipeb_example
  bad$case <- bad$case + 1L
  expect_error(ipeb(bad, markers = "m1"), "0/1")
})

test_that("evaluate rejects invalid specificities", {
  skip_on_cran()
  data(ipeb_example)
  tr <- subset(ipeb_example, split == "train")
  te <- subset(ipeb_example, split == "test")
  fit <- ipeb(tr, markers = c("m1", "m2"), objective = "sensitivity",
              innovation = "iid", slope = "off")
  expect_error(evaluate(fit, te, specificities = 1.2), "between 0 and 1")
})
