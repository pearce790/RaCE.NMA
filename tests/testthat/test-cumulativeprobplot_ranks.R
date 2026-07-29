library(ggplot2)
test_that("ensure proper ggplot output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p1 <- cumulativeprobplot_ranks(mcmc=res)
  p2 <- cumulativeprobplot_ranks(mcmc=res,names=c("1","2"))

  data(toy_data)
  q1 <- cumulativeprobplot_ranks(data=toy_data)
  q2 <- cumulativeprobplot_ranks(data=toy_data,names=c("1","2","3","4"))

  expect_true(is_ggplot(p1))
  expect_true(is_ggplot(p2))
  expect_true(is_ggplot(q1))
  expect_true(is_ggplot(q2))
})

test_that("error if names argument of incorrect length", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_error(cumulativeprobplot_ranks(mcmc=res,names=c("1","2","3")))

  expect_error(cumulativeprobplot_ranks(data=toy_data,names=c("1","2","3")))
})

test_that("error if mcmc supplied for data argument and vice versa", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_error(cumulativeprobplot_ranks(data=res))
  expect_error(cumulativeprobplot_ranks(mcmc=toy_data))
})
