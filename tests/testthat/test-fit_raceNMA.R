test_that("output structure is correct for the independence path (mu_hat + s)", {
  res <- fit_raceNMA(mu_hat = c(0, 0, 1, 1), s = c(.1, .1, .1, .1),
                      mu0 = 0.5, sigma0 = 5, tau = 0.5, nu0 = NULL, iter = 20, nu_iter = 2)
  expect_type(res, "list")
  expect_named(res, c("mu", "nu", "g", "K"))
  expect_equal(dim(res$mu), c(20 * 2, 4))
  expect_equal(dim(res$g), c(20 * 2, 4))
  expect_length(res$K, 20 * 2)
  expect_true(all(res$K >= 1 & res$K <= 4))
})

test_that("output structure is correct for the correlation path (mu_hat + cov)", {
  cov_mat <- diag(rep(.01, 4))
  res <- fit_raceNMA(mu_hat = c(0, 0, 1, 1), cov = cov_mat,
                      mu0 = 0.5, sigma0 = 5, tau = 0.5, nu0 = NULL, iter = 20, nu_iter = 2)
  expect_equal(dim(res$mu), c(20 * 2, 4))
  expect_true(all(res$K >= 1 & res$K <= 4))
})

test_that("output structure is correct when posterior draws are supplied directly", {
  res <- fit_raceNMA(posterior = toy_data, mu0 = 0, sigma0 = 5, tau = 0.5,
                      nu0 = NULL, iter = 10, nu_iter = 2)
  expect_equal(dim(res$mu), c(10 * 2, ncol(toy_data)))
})

test_that("nu0 can be used to specify a custom starting partition", {
  res <- fit_raceNMA(mu_hat = c(0, 0.05, 2), s = c(.1, .1, .1),
                      mu0 = 0.5, sigma0 = 5, tau = 0.5, nu0 = c(1, 1, 2), iter = 5, nu_iter = 1)
  expect_equal(dim(res$mu), c(5, 3))
  expect_true(all(res$K >= 1 & res$K <= 3))
})

test_that("warns when posterior is supplied alongside mu_hat/s/cov, and ignores them", {
  expect_warning(
    fit_raceNMA(posterior = toy_data, mu_hat = c(0, 0, 0, 0),
                tau = 0.5, nu0 = NULL, iter = 3, nu_iter = 1)
  )
  expect_warning(
    fit_raceNMA(posterior = toy_data, s = rep(.1, 4),
                tau = 0.5, nu0 = NULL, iter = 3, nu_iter = 1)
  )
  expect_warning(
    fit_raceNMA(posterior = toy_data, cov = diag(rep(.1, 4)),
                tau = 0.5, nu0 = NULL, iter = 3, nu_iter = 1)
  )
})

test_that("errors on missing or incongruous inputs", {
  # no posterior and no mu_hat
  expect_error(
    fit_raceNMA(s = c(.1, .1), mu0 = 0, sigma0 = 1, tau = .1, nu0 = NULL, iter = 5, nu_iter = 1)
  )
  # mu_hat/s length mismatch
  expect_error(
    fit_raceNMA(mu_hat = c(0, 1, 2), s = c(.1, .1),
                mu0 = 0, sigma0 = 1, tau = .1, nu0 = NULL, iter = 5, nu_iter = 1)
  )
  # mu_hat/cov dimension mismatch
  expect_error(
    fit_raceNMA(mu_hat = c(0, 1), cov = diag(rep(1, 3)),
                mu0 = 0, sigma0 = 1, tau = .1, nu0 = NULL, iter = 5, nu_iter = 1)
  )
  # neither cov nor s supplied
  expect_error(
    fit_raceNMA(mu_hat = c(0, 1), mu0 = 0, sigma0 = 1, tau = .1, nu0 = NULL, iter = 5, nu_iter = 1)
  )
})

test_that("results are reproducible given the same seed", {
  set.seed(42)
  res1 <- fit_raceNMA(mu_hat = c(0, 1), s = c(.1, .1),
                       mu0 = 0.5, sigma0 = 1, tau = 0.5, nu0 = NULL, iter = 15, nu_iter = 1)
  set.seed(42)
  res2 <- fit_raceNMA(mu_hat = c(0, 1), s = c(.1, .1),
                       mu0 = 0.5, sigma0 = 1, tau = 0.5, nu0 = NULL, iter = 15, nu_iter = 1)
  expect_identical(res1, res2)
})
