test_that("ensure proper output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1,2),s=c(.1,.1,.1),seed=1,chains=2,iter=100,verbose = F)
  expect_true(is.data.frame(calculate_SUCRA_MNBT(mcmc=res)))
  expect_true(is.data.frame(calculate_SUCRA_MNBT(mcmc=res,level=0.8)))
  expect_true(is.data.frame(calculate_SUCRA_MNBT(mcmc=res,names=c("1","2","3"))))
  expect_true(is.data.frame(calculate_SUCRA_MNBT(data=toy_data)))
  expect_true(is.data.frame(calculate_SUCRA_MNBT(data=toy_data,level=0.8)))
  expect_true(is.data.frame(calculate_SUCRA_MNBT(data=toy_data,names=c("1","2","3","4"))))
})

test_that("error if mcmc supplied for data and vice versa",{
  expect_error(calculate_SUCRA_MNBT(mcmc=toy_data))
  res <- mcmc_raceNMA(mu_hat=c(0,1,2),s=c(.1,.1,.1),seed=1,chains=2,iter=100,verbose = F)
  expect_error(calculate_SUCRA_MNBT(data=res))
})
