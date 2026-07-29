test_that("ensure proper output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_equal(class(calculate_Rhat(res)),"gelman.diag")
  expect_equal(class(calculate_Rhat(res,names=c("1","2"))),"gelman.diag")
  expect_equal(class(calculate_Rhat(res,level=0.9)),"gelman.diag")
  expect_equal(class(calculate_Rhat(res,multivariate = T)),"gelman.diag")
})
