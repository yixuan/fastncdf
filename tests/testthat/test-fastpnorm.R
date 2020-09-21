context("testing fastpnorm")

test_that("fastpnorm gives the correct result", {
  x <- seq(-6, 6, by = 1e-6)
  res <- fastpnorm(x)
  expect_equal(res, pnorm(x), tolerance = 1e-7)

  res <- fastpnorm(x, use_cubic = TRUE)
  expect_equal(res, pnorm(x), tolerance = 1e-7)
})

test_that("fastpnorm gives the correct result", {
  x <- seq(-6, 6, by = 1e-6)
  res <- rep(0., length(x))
  fastpnorm_preallocated(x, p = res)
  expect_equal(res, pnorm(x), tolerance = 1e-7)

  res <- rep(0., length(x))
  fastpnorm_preallocated(x, p = res, use_cubic = TRUE)
  expect_equal(res, pnorm(x), tolerance = 1e-7)
})
