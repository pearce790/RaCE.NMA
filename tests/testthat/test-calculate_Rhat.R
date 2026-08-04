test_that("ensure proper output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_equal(class(calculate_Rhat(res)),"gelman.diag")
  expect_equal(class(calculate_Rhat(res,names=c("1","2"))),"gelman.diag")
  expect_equal(class(calculate_Rhat(res,level=0.9)),"gelman.diag")
  expect_equal(class(calculate_Rhat(res,multivariate = T)),"gelman.diag")
})

test_that("output correctness", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  rhat <- calculate_Rhat(res)

  expected_psrf <- matrix(c(1.027629, 1.037689,
                            1.0954229, 1.0714156),
                          nrow = 2, dimnames = list(c("mu1","mu2"), c("Point est.","Upper C.I.")))
  expect_equal(unclass(rhat$psrf), expected_psrf, tolerance = 1e-6)
})
