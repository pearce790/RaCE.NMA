test_that("toy_data has the documented structure", {
  data("toy_data")
  expect_true(is.data.frame(toy_data))
  expect_equal(ncol(toy_data), 4)     # 4 treatments
  expect_equal(nrow(toy_data), 10000) # 10k posterior draws
  expect_true(all(sapply(toy_data, is.numeric)))
  expect_false(anyNA(toy_data))
})

test_that("wang_posterior has the documented structure", {
  data("wang_posterior")
  expect_true(is.data.frame(wang_posterior))
  expect_equal(ncol(wang_posterior), 10)    # 10 baseline treatments
  expect_equal(nrow(wang_posterior), 60000) # 60k posterior draws
  expect_true(all(sapply(wang_posterior, is.numeric)))
  expect_false(anyNA(wang_posterior))
})
