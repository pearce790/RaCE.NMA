library(ggplot2)
test_that("ensure proper ggplot output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p <- forestplot_muhat(mcmc=res)
  q <- forestplot_muhat(data=toy_data)

  expect_true(is_ggplot(p))
  expect_true(is_ggplot(q))
})

test_that("error if mcmc supplied for data and vice versa",{
  expect_error(forestplot_muhat(mcmc=toy_data))
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_error(forestplot_muhat(data=res))
})


