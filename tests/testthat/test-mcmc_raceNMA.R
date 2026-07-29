test_that("correct input warning messages", {
  data("toy_data")
  expect_warning(
    mcmc_raceNMA(posterior=toy_data, mu_hat=c(0,1),
                 seed=1,chains = 1, iter = 10, verbose = F)
  )
  expect_warning(
    mcmc_raceNMA(posterior=toy_data, cov=diag(c(.1,.1)),
                 seed=1,chains = 1, iter = 10, verbose = F)
  )
  expect_warning(
    mcmc_raceNMA(posterior=toy_data, s=c(.1,.1),
                 seed=1,chains = 1, iter = 10, verbose = F)
  )
  expect_warning(
    mcmc_raceNMA(mu_hat=c(0,1), cov=diag(c(.1,.1)),s=c(.1,.1),
                 seed=1,chains = 1, iter = 10, verbose = F)
  )
})

test_that("errors when incomplete/incongruous inputs given", {
  expect_error(
    mcmc_raceNMA(mu_hat=c(0,1),s=c(.1),
                 seed=1,chains = 1, iter = 10, verbose = F)
  )
  expect_error(
    mcmc_raceNMA(mu_hat=c(0,1,2),cov=diag(c(1,1)),
                 seed=1,chains = 1, iter = 10, verbose = F)
  )
})

test_that("correct output format", {
  expect_true({
    res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),
                        seed=1,chains = 1, iter = 10, verbose = F)
    is.data.frame(res)
  })
  expect_true({
    res <- mcmc_raceNMA(posterior=toy_data,
                        seed=1,chains = 1, iter = 10, verbose = F)
    is.data.frame(res)
  })
  expect_true({
    res <- mcmc_raceNMA(posterior=as.matrix(toy_data),
                        seed=1,chains = 1, iter = 10, verbose = F)
    is.data.frame(res)
  })
  expect_all_true({
    res <- mcmc_raceNMA(posterior=toy_data,seed=1,verbose = F,
                        chains=4, iter=3, nu_iter=4, burn=0.5)
    dim(res)==c(4*3*4/2,3*ncol(toy_data)+3)
  })
})

test_that("output correctness", {
  expect_equal({
    res <- mcmc_raceNMA(posterior=toy_data,seed=1,verbose = F,
                        chains=1,iter=50)
    res[1,4]
  },-0.9047409,tolerance = 1e-5)
})
