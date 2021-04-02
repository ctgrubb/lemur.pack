test_that("lemur.control throws error with non-positive nchains", {
  expect_error(lemur.control(nchains = -1))
  expect_error(lemur.control(nchains = 0))
  expect_error(lemur.control(nchains = ""))
  expect_error(lemur.control(nchains = factor("")))
})

test_that("lemur.control throws error with non-positive nsteps", {
  expect_error(lemur.control(nsteps = -1))
  expect_error(lemur.control(nsteps = 0))
  expect_error(lemur.control(nsteps = ""))
  expect_error(lemur.control(nsteps = factor("")))
})

test_that("lemur.control throws error with non-positive burnin", {
  expect_error(lemur.control(burnin = -1))
  expect_error(lemur.control(burnin = 0))
  expect_error(lemur.control(burnin = ""))
  expect_error(lemur.control(burnin = factor("")))
})

test_that("lemur.control throws error with non-positive thin", {
  expect_error(lemur.control(thin = -1))
  expect_error(lemur.control(thin = 0))
  expect_error(lemur.control(thin = ""))
  expect_error(lemur.control(thin = factor("")))
})
