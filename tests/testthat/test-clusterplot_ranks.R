library(ggplot2)
test_that("ensure proper ggplot output", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p <- clusterplot_ranks(mcmc=res)
  q <- clusterplot_ranks(data=toy_data)

  expect_true(is_ggplot(p))
  expect_true(is_ggplot(q))
})

test_that("error if mcmc supplied for data and vice versa",{
  expect_error(clusterplot_ranks(mcmc=toy_data))
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  expect_error(clusterplot_ranks(data=res))
})

test_that("label_ranks argument adds a text layer without erroring", {
  res <- mcmc_raceNMA(mu_hat=c(0,1),s=c(.1,.1),seed=1,chains=2,iter=50,verbose = F)
  p_no_labels <- clusterplot_ranks(mcmc=res)
  p_with_labels <- clusterplot_ranks(mcmc=res,label_ranks=1:2)
  q_with_labels <- clusterplot_ranks(data=toy_data,label_ranks=1:2)

  expect_true(is_ggplot(p_with_labels))
  expect_true(is_ggplot(q_with_labels))
  expect_equal(length(p_no_labels$layers), 1)
  expect_equal(length(p_with_labels$layers), 2)
})

# Same hand-computable draws used in test-calculate_SUCRA_MNBT.R:
# A is ranked 1st in 3/4 draws and 2nd in 1/4; B is the mirror image.
test_that("clusterplot_ranks computes the correct posterior rank probabilities", {
  deterministic_data <- matrix(c(1, 2,
                                 1, 2,
                                 2, 1,
                                 1, 2), ncol = 2, byrow = TRUE)
  p <- clusterplot_ranks(data = deterministic_data, names = c("A", "B"))
  pd <- p$data

  a_vals <- pd$value[pd$variable == "A"][order(pd$rank[pd$variable == "A"])]
  b_vals <- pd$value[pd$variable == "B"][order(pd$rank[pd$variable == "B"])]

  expect_equal(a_vals, c(0.75, 0.25))
  expect_equal(b_vals, c(0.25, 0.75))
})
