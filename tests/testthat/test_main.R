
test_that("Invalid rE-sequence inputs give an input error",{
  expect_error(gsDesignOC(n.stages = 2, rE.seq=c(3, 2, 1)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=c(3, 2)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=c(0.5, 1)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(2, 0.5)), "^Invalid input:")
})
test_that("Invalid rF-sequence inputs give an input error",{
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, rF.seq=c(-2, -1, 0)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, rF.seq=c(-2, -1)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, rF.seq=c(1, 0)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(3,2), rF.seq=c(-1, -2)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(3,2), rF.seq=NULL, futility.type="binding"), "^Invalid input:")
})
test_that("Integers and rounding in r-sequence do NOT give an input error",{
  res1 <- tryCatch(gsDesignOC(n.stages = 1, rE.seq = 1L), error=function(e)e)
  expect_true(!is(res1, "error") || !grepl(res1$message, "^Invalid input:"))

  res2 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = c(2L,1L)), error=function(e)e)
  expect_true(!is(res2, "error") || !grepl(res2$message, "^Invalid input:"))

  res3 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = c(2, 0.99999999999)), error=function(e)e)
  expect_true(!is(res3, "error") || !grepl(res3$message, "^Invalid input:"))

  res4 <- tryCatch(gsDesignOC(n.stages = 1, rE.seq = 1, rF.seq = 0L), error=function(e)e)
  expect_true(!is(res4, "error") || !grepl(res4$message, "^Invalid input:"))

  res5 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = 2, rF.seq = c(-1L,0L)), error=function(e)e)
  expect_true(!is(res5, "error") || !grepl(res5$message, "^Invalid input:"))

  res6 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = 2, rF.seq = c(-1, 0.000000001)), error=function(e)e)
  expect_true(!is(res6, "error") || !grepl(res6$message, "^Invalid input:"))
})
test_that("Invalid r_EN weights give an input error",{
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, r_EN=1, r_EN.w = c(1,2)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, r_EN=c(0, 1), r_EN.w = 1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, r_EN=c(0, 1), r_EN.w = c(-1, 1)), "^Invalid input:")
})
test_that("Invalid power inputs give an input error",{
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= 1.2), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, sig.level = 0.1, power= 0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= -0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= 0.8, power.efficacy = 0.9), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= 0.8, power.efficacy = -0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(2,1.5), power= 0.8, power.efficacy = c(0.8,0.5)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, sig.level=0.1, power.futility = 0.95), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, sig.level=0.1, power.futility = -0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(2,1.5), rF.seq=c(-1,-0.5), power= 0.8, power.futility = c(0.8,0.5)), "^Invalid input:")})

