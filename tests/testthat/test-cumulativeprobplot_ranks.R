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

# Same hand-computable draws used in test-calculate_SUCRA_MNBT.R:
# A is ranked 1st in 3/4 draws, B in 1/4. Cumulative P(rank<=1), P(rank<=2):
#   A: (0.75, 1), B: (0.25, 1)
test_that("cumulativeprobplot_ranks computes the correct cumulative rank probabilities", {
  deterministic_data <- matrix(c(1, 2,
                                 1, 2,
                                 2, 1,
                                 1, 2), ncol = 2, byrow = TRUE)
  p <- cumulativeprobplot_ranks(data = deterministic_data, names = c("A", "B"))
  pd <- p$data

  a_vals <- pd$value[pd$Var2 == "A"][order(pd$Var1[pd$Var2 == "A"])]
  b_vals <- pd$value[pd$Var2 == "B"][order(pd$Var1[pd$Var2 == "B"])]

  expect_equal(a_vals, c(0.75, 1))
  expect_equal(b_vals, c(0.25, 1))
})
