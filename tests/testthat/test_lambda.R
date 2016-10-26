library(grpSLOPE)

context("lambdaMCGroupSLOPE()")

n.obs <- 50
n.group <- 10
group <- c(10, 4, 7,10, 6, 9, 9, 2, 5, 1, 6, 2, 6, 1, 2, 1,
           2, 1, 5, 9, 5, 1, 3,10, 4, 5, 3, 7, 6, 9)
wt <- sapply(getGroupID(group), length)
A <- as.matrix(read.table("./test_data/gaussianMC_test_mat.txt"))
fdr <- 0.1

test_that("method 'gaussianMC' produces a non-increasing sequence", {
  lambda.MC <- lambdaMCGroupSLOPE(fdr=fdr, group=group, A=A, n.obs=n.obs,
                                  method="gaussianMC", MC.reps=50)
  # check for NA and NaN
  expect_true(!anyNA(lambda.MC))
  # check non-increasing
  expect_true(all(diff(lambda.MC) <= 0))
})
