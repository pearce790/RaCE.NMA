# These test the RJMCMC partition update used internally by fit_raceNMA.
# Instead of testing exact numeric draws, check invariants that must hold
# regardless, i.e., valid cluster labels, correct lengths, jump styles, etc.

test_that("sample_partition_correlation returns a valid partition structure", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  cov_mat <- diag(rep(0.05, J))

  set.seed(1)
  out <- sample_partition_correlation(mu_hat = mu_hat, J = J, nu = mu_hat, g = 1:J, K = J,
                                       mu0 = 1, sigma0 = 5, cov = cov_mat, tau = 0.3)
  expect_type(out, "list")
  expect_named(out, c("g", "nu", "K"))
  expect_length(out$g, J)
  expect_true(out$K >= 1 && out$K <= J)
  expect_length(out$nu, out$K)
  expect_setequal(unique(out$g), 1:out$K)
})

test_that("sample_partition_correlation can only merge (not split) from an all-singleton start", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  cov_mat <- diag(rep(0.05, J))

  set.seed(2)
  out <- sample_partition_correlation(mu_hat = mu_hat, J = J, nu = mu_hat, g = 1:J, K = J,
                                       mu0 = 1, sigma0 = 5, cov = cov_mat, tau = 0.3)
  expect_true(out$K <= J)
})

test_that("sample_partition_correlation can only split (not merge) from a single-cluster start", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  cov_mat <- diag(rep(0.05, J))

  set.seed(3)
  out <- sample_partition_correlation(mu_hat = mu_hat, J = J, nu = mean(mu_hat), g = rep(1, J), K = 1,
                                       mu0 = 1, sigma0 = 5, cov = cov_mat, tau = 0.3)
  expect_true(out$K %in% c(1, 2))
})

test_that("sample_partition_independence returns a valid partition structure", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  s <- rep(0.1, J)

  set.seed(1)
  out <- sample_partition_independence(mu_hat = mu_hat, J = J, nu = mu_hat, g = 1:J, K = J,
                                        mu0 = 1, sigma0 = 5, s = s, tau = 0.3)
  expect_type(out, "list")
  expect_named(out, c("g", "nu", "K"))
  expect_length(out$g, J)
  expect_true(out$K >= 1 && out$K <= J)
  expect_setequal(unique(out$g), 1:out$K)
})

test_that("sample_partition_independence can only merge from an all-singleton start", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  s <- rep(0.1, J)

  set.seed(2)
  out <- sample_partition_independence(mu_hat = mu_hat, J = J, nu = mu_hat, g = 1:J, K = J,
                                        mu0 = 1, sigma0 = 5, s = s, tau = 0.3)
  expect_true(out$K <= J)
})

test_that("sample_partition_independence can only split from a single-cluster start", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  s <- rep(0.1, J)

  set.seed(3)
  out <- sample_partition_independence(mu_hat = mu_hat, J = J, nu = mean(mu_hat), g = rep(1, J), K = 1,
                                        mu0 = 1, sigma0 = 5, s = s, tau = 0.3)
  expect_true(out$K %in% c(1, 2))
})

test_that("partition updates are reproducible given the same seed", {
  mu_hat <- c(0, 0.1, 2, 2.1)
  J <- 4
  s <- rep(0.1, J)

  set.seed(99)
  out1 <- sample_partition_independence(mu_hat = mu_hat, J = J, nu = mu_hat, g = 1:J, K = J,
                                         mu0 = 1, sigma0 = 5, s = s, tau = 0.3)
  set.seed(99)
  out2 <- sample_partition_independence(mu_hat = mu_hat, J = J, nu = mu_hat, g = 1:J, K = J,
                                         mu0 = 1, sigma0 = 5, s = s, tau = 0.3)
  expect_identical(out1, out2)
})
