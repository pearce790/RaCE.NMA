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
# The tests below use a small, fully deterministic set of draws so that SUCRA and
# MNBT can be hand-computed and checked exactly (rather than only checking types).
#
# 4 draws, 2 treatments. Treatment A is rank 1 in 3 of 4 draws, Treatment B
# is rank 1 in remaining 1 draw. Thus,
#   SUCRA_A = mean(cumulative rank-1 probability across ranks 1:(J-1)) = P(A ranked 1st) = 3/4
#   SUCRA_B = P(B ranked 1st) = 1/4
#   MNBT_A = median = 0, 25th pctile = 0, 75th pctile = 0.25
#   MNBT_B = median = 1, 25th pctile = 0.75, 75th pctile = 1
test_that("calculate_SUCRA_MNBT gives exactly correct values for a hand-computable case (data input)", {
  deterministic_data <- matrix(c(1, 2,
                                 1, 2,
                                 2, 1,
                                 1, 2), ncol = 2, byrow = TRUE)
  res <- calculate_SUCRA_MNBT(data = deterministic_data, names = c("A", "B"))
  res <- res[match(c("A", "B"), res$Treatment), ]

  expect_equal(res$SUCRA, c(0.75, 0.25))
  expect_equal(res[["MNBT (50% CI)"]], c("0 (0, 0.25)", "1 (0.75, 1)"))
})

test_that("calculate_SUCRA_MNBT gives identical, correct values for the equivalent mcmc-format input", {
  # Same draws as above, reshaped into the chain/iteration/K + mu1..muJ format
  # produced by mcmc_raceNMA, but constructed directly (no stochastic sampling)
  # so the expected numbers are exactly known.
  mock_mcmc <- data.frame(chain = factor(rep(1, 4)), iteration = 1:4, K = 2,
                          mu1 = c(1, 1, 2, 1), mu2 = c(2, 2, 1, 2))
  res <- calculate_SUCRA_MNBT(mcmc = mock_mcmc, names = c("A", "B"))
  res <- res[match(c("A", "B"), res$Treatment), ]

  expect_equal(res$SUCRA, c(0.75, 0.25))
  expect_equal(res[["MNBT (50% CI)"]], c("0 (0, 0.25)", "1 (0.75, 1)"))
})

test_that("calculate_SUCRA_MNBT handles tied values correctly (rank ties.method = 'min')", {
  # A single draw where two treatments are exactly tied for best: both should
  # receive rank 1 (ties.method = "min"), giving both an equal, perfect SUCRA.
  tied_data <- matrix(c(1, 1, 2), nrow = 1)
  res <- calculate_SUCRA_MNBT(data = tied_data, names = c("A", "B", "C"))
  res <- res[match(c("A", "B", "C"), res$Treatment), ]

  expect_equal(res$SUCRA, c(1, 1, 0))
})
