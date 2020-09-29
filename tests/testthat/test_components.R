
test_that("calc.bounds for one stage returns fixed-sample test",{
  x <- list(n.stages=1, rE.seq=0.5, sig.level=0.05, power=0.8, power.efficacy=0.8, futility.type="none")
  xx <- calc.bounds(x, alpha.seq=0.05)
  expect_equal(xx$lower, -20)
  expect_equal(xx$upper, qnorm(0.95))
  expect_equal(xx$info, 24.7302289280791)
  expect_equal(xx$spending, 0.05)
 })

test_that("ztest.I is consistent with asbio::power.z.test",{
  i1 <- ztest.I(delta = 0.5, sig.level = 0.05, power=0.8)
  expect_equal(i1, 24.7302289280791)
  })

test_that("No data is needed if sig.level = power",{
  i2 <- ztest.I(delta = 0.5, sig.level = 0.1, power=0.1)
  expect_equal(i2, 0)
  })
