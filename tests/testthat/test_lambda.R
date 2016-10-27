library(grpSLOPE)

context("lambdaMCGroupSLOPE()")

n.obs <- 50
n.group <- 10
group <- c(10, 4, 7,10, 6, 9, 8, 2, 5, 1, 6, 2, 6, 1, 2, 1,
           2, 1, 5, 9, 5, 1, 3,10, 4, 5, 3, 7, 6, 9)
A <- as.matrix(read.table("./test_data/test_mat.txt"))
fdr <- 0.1

test_that("method 'BH' computes the sequence correctly", {
  # sol is obtained via the function create_lambda in package SLOPE
  sol <- c(2.575829,2.326348,2.170090,2.053749,1.959964,
           1.880794,1.811911,1.750686,1.695398,1.6448540)
  lambda.BH <- lambdaBH(fdr=fdr, n.group=n.group)
  expect_equal(lambda.BH, sol, tolerance=1e-6)
})

test_that("method 'gaussianMC' with n.MC=5 produces a non-increasing sequence of the correct length", {
  gaussian.MC <- lambdaMC(fdr=fdr, group=group, A=A,
                          method="gaussianMC", n.MC=5, MC.reps=50)

  # check length
  expect_equal(length(gaussian.MC), n.group)
  # check error messages and warnings
  expect_warning(lambdaMC(fdr=fdr, group=group, A=A, method="gaussianMC", 
                          n.MC = 12, MC.reps=50),
                 "n.MC is not allowed to exceed the number of groups.")
  expect_error(lambdaMC(fdr=fdr, group=group, A=A, method="gaussianMC", 
                        n.MC = 1.1, MC.reps=-50),
               "n.MC and MC.reps have to be positive integers.")
  expect_error(lambdaMC(fdr=fdr, group=group, A=A, method="gaussianMC", 
                        n.MC = 5, MC.reps=5.05),
               "n.MC and MC.reps have to be positive integers.")
  # check for NA and NaN
  expect_true(!anyNA(gaussian.MC))
  # check non-increasing
  expect_true(all(diff(gaussian.MC) <= 0))
})

test_that("method 'gaussianMC' with n.MC=10 produces a non-increasing sequence of the correct length", {
  gaussian.MC <- lambdaMC(fdr=fdr, group=group, A=A,
                          method="gaussianMC", n.MC=10, MC.reps=50)

  # check length
  expect_equal(length(gaussian.MC), n.group)
  # check error messages and warnings
  expect_warning(lambdaMC(fdr=fdr, group=group, A=A, method="gaussianMC", 
                          n.MC=12, MC.reps=50),
                 "n.MC is not allowed to exceed the number of groups.")
  # check for NA and NaN
  expect_true(!anyNA(gaussian.MC))
  # check non-increasing
  expect_true(all(diff(gaussian.MC) <= 0))
})

test_that("method 'chiMC' with n.MC=5 produces a non-increasing sequence of the correct length", {
  y  <- A %*% c(rep(3, 10), rep(0, 20))
  wt <- sqrt(sapply(grpSLOPE::getGroupID(group), length))
  chi.MC <- lambdaMC(fdr=fdr, group=group, A=A, y=y, wt=wt,
                     method="chiMC", n.MC=5, MC.reps=50)

  # check length
  expect_equal(length(chi.MC), n.group)
  # check error messages and warnings
  expect_warning(lambdaMC(fdr=fdr, group=group, A=A, y=y, wt=wt, 
                          method="chiMC", n.MC = 112, MC.reps=50),
                 "n.MC is not allowed to exceed the number of groups.")
  expect_error(lambdaMC(fdr=fdr, group=group, A=A, wt=wt, 
                          method="chiMC", n.MC=5, MC.reps=50),
               "y and wt need to be passed as arguments when method is 'chiMC'.")
  expect_error(lambdaMC(fdr=fdr, group=group, A=A, y = y, 
                          method="chiMC", n.MC=5, MC.reps=50),
               "y and wt need to be passed as arguments when method is 'chiMC'.")
  # check for NA and NaN
  expect_true(!anyNA(chi.MC))
  # check non-increasing
  expect_true(all(diff(chi.MC) <= 0))
})

test_that("method 'chiMC' with n.MC=10 produces a non-increasing sequence of the correct length", {
  y  <- A %*% c(rep(3, 10), rep(0, 20))
  wt <- sqrt(sapply(grpSLOPE::getGroupID(group), length))
  chi.MC <- lambdaMC(fdr=fdr, group=group, A=A, y=y, wt=wt,
                     method="chiMC", n.MC=10, MC.reps=50)

  # check length
  expect_equal(length(chi.MC), n.group)
  # check error messages and warnings
  expect_warning(lambdaMC(fdr=fdr, group=group, A=A, y=y, wt=wt, 
                          method="chiMC", n.MC=11, MC.reps=50),
                 "n.MC is not allowed to exceed the number of groups.")
  expect_error(lambdaMC(fdr=fdr, group=group, A=A, wt=wt, 
                          method="chiMC", n.MC=10, MC.reps=50),
               "y and wt need to be passed as arguments when method is 'chiMC'.")
  expect_error(lambdaMC(fdr=fdr, group=group, A=A, y = y, 
                          method="chiMC", n.MC=10, MC.reps=50),
               "y and wt need to be passed as arguments when method is 'chiMC'.")
  # check for NA and NaN
  expect_true(!anyNA(chi.MC))
  # check non-increasing
  expect_true(all(diff(chi.MC) <= 0))
})
