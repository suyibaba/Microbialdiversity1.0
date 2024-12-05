test_that("Shannon diversity calculations are correct", {
  equal_abundances <- c(1, 1, 1, 1)
  expect_equal(calculate_shannon(equal_abundances), log(4), tolerance = 1e-4)
  
  unequal_abundances <- c(10, 20, 30, 40)
  expect_equal(calculate_shannon(unequal_abundances), 1.279, tolerance = 1e-3)
  
  expect_error(calculate_shannon(c(-1, 2, 3)))
})

test_that("Simpson diversity calculations are correct", {
  equal_abundances <- c(1, 1, 1, 1)
  expect_equal(calculate_simpson(equal_abundances), 0.75, tolerance = 1e-4)
  
  unequal_abundances <- c(10, 20, 30, 40)
  expect_equal(calculate_simpson(unequal_abundances), 0.70, tolerance = 1e-2)
  
  expect_error(calculate_simpson(c(-1, 2, 3)))
})

test_that("Inverse Simpson diversity calculations are correct", {
  equal_abundances <- c(1, 1, 1, 1)
  expect_equal(calculate_inverse_simpson(equal_abundances), 4.0, tolerance = 1e-3)
  
  unequal_abundances <- c(10, 20, 30, 40)
  expect_equal(calculate_inverse_simpson(unequal_abundances), 3.33, tolerance = 1e-2)
  
  expect_error(calculate_inverse_simpson(c(-1, 2, 3)))
})

test_that("Observed OTUs calculations are correct", {
  all_present <- c(1, 2, 3, 4)
  expect_equal(calculate_observed_otus(all_present), 4)
  
  some_zeros <- c(0, 2, 0, 4, 5)
  expect_equal(calculate_observed_otus(some_zeros), 3)
  
  expect_error(calculate_observed_otus(c(-1, 2, 3)))
})