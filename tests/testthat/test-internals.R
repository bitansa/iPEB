test_that(".auc is 1 for perfect ordering, 0 for reversed, 0.5 at chance", {
  expect_equal(iPEB:::.auc(c(1, 2, 3, 4), c(0, 0, 1, 1)), 1)     # cases score highest
  expect_equal(iPEB:::.auc(c(1, 2, 3, 4), c(1, 1, 0, 0)), 0)     # cases score lowest
  expect_equal(iPEB:::.auc(c(1, 2, 3, 4), c(1, 0, 0, 1)), 0.5)   # interleaved -> chance
})

test_that(".metrics returns Sens/LT/Spec/Ref and honors the threshold", {
  score <- c(0, 0, 3, 3)
  D     <- c(0, 0, 1, 1)
  ids   <- c(1, 2, 3, 4)
  ttd   <- c(100, 100, 200, 300)
  ctrl  <- score[D == 0]
  mt <- iPEB:::.metrics(score, D, ids, ttd, ctrl, alpha = 0.95)
  expect_named(mt, c("SensW", "LT", "Spec", "Ref"))
  expect_equal(unname(mt["SensW"]), 1)      # both cases exceed the control threshold
  expect_true(mt["Spec"] >= 0 && mt["Spec"] <= 1)
})

test_that(".lambda_profile matches the manuscript profiles", {
  expect_equal(unname(iPEB:::.lambda_profile("sens")), c(50, 0, 0))
  expect_equal(unname(iPEB:::.lambda_profile("leadtime")), c(50, 0, 0.7))
  expect_equal(unname(iPEB:::.lambda_profile("combined")), c(50, 1, 0.5))
})
