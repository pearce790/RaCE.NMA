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
  res <- mcmc_raceNMA(posterior=toy_data,seed=1,verbose = F,
                      chains=1,iter=50)

  expect_equal({res[1,4]},-0.9847569,tolerance = 1e-6)

  mu_means <- colMeans(res[,paste0("mu",1:4)])
  expected_means <- c(mu1=-0.999430254, mu2=0.007434866,
                      mu3=1.013860323, mu4=-0.128857193)
  expect_equal(mu_means, expected_means, tolerance = 1e-6)

  expect_equal(mean(res$K), 4, tolerance = 1e-6)

  last_row <- unname(unlist(res[nrow(res),paste0("mu",1:4)]))
  expect_equal(last_row, c(-1.01160135,-0.03766098,0.85879089,1.13997159),
               tolerance = 1e-6)
})

test_that("identical outputs across cores", {
  res_1 <- mcmc_raceNMA(posterior=toy_data,seed=1,verbose = F,
                      chains = 4, iter = 50,cores = 1)
  res_2 <- mcmc_raceNMA(posterior=toy_data,seed=1,verbose = F,
                        chains = 4, iter = 50,cores = 2)
  expect_true(identical(res_1,res_2))
})

test_that("non-identical outputs across cores when seed unspecified", {
  res_1_noseed <- mcmc_raceNMA(posterior=toy_data,verbose = F,
                        chains = 4, iter = 50,cores = 1)
  res_2_noseed <- mcmc_raceNMA(posterior=toy_data,verbose = F,
                        chains = 4, iter = 50,cores = 2)
  expect_false(identical(res_1_noseed,res_2_noseed))
})
