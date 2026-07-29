library(ggplot2)
test_that("ensure proper ggplot output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p <- traceplot_mu(res)
  q <- traceplot_mu(res,names = c("1","2"))

  expect_true(is_ggplot(p))
  expect_equal(p$labels$colour,"Chain")
  expect_equal(c(get_strip_labels(p)$facets$variable),c("Treatment 1","Treatment 2"))
  expect_equal(c(get_strip_labels(q)$facets$variable),c("1","2"))
})

test_that("error if names argument of incorrect length", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_error(traceplot_mu(res,names=c("1")))
  expect_error(traceplot_mu(res,names=c("1","2","3")))
})
