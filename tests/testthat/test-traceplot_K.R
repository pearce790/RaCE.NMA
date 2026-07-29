library(ggplot2)
test_that("ensure proper ggplot output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p <- traceplot_K(res)

  expect_true(is_ggplot(p))
  expect_equal(p$labels$colour,"Chain")
})
