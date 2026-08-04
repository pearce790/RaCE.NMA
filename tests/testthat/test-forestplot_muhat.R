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

test_that("names argument is applied to the plotted treatment labels", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p_named <- forestplot_muhat(mcmc=res,names=c("Drug A","Drug B"))
  expect_true(is_ggplot(p_named))
  expect_setequal(levels(p_named$data$name), c("Drug A","Drug B"))
})

test_that("order_by_average=FALSE keeps treatments in original column order", {
  p_unordered <- forestplot_muhat(data=toy_data,order_by_average=FALSE)
  expect_true(is_ggplot(p_unordered))
  expect_equal(levels(p_unordered$data$name), as.character(1:ncol(toy_data)))
})

# 5 draws each for two treatments, chosen so the mean and the default 95% quantile
# interval can be hand-computed exactly:
#   A = 1:5      -> mean = 3,  95% CI = quantile(1:5, c(.025,.975)) = (1.1, 4.9)
#   B = 10*(1:5) -> mean = 30, 95% CI = (11, 49)  [B is just A scaled by 10]
test_that("forestplot_muhat computes the correct posterior mean and credible interval", {
  data_mat <- cbind(A = 1:5, B = 10 * (1:5))
  p <- forestplot_muhat(data = data_mat, names = c("A", "B"))
  pd <- p$data
  pd <- pd[match(c("A", "B"), as.character(pd$name)), ]

  expect_equal(pd$mean, c(3, 30))
  expect_equal(pd$lower_CI, c(1.1, 11), tolerance = 1e-8)
  expect_equal(pd$upper_CI, c(4.9, 49), tolerance = 1e-8)
})
