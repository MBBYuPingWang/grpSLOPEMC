###############################################################################
#
#    grpSLOPEMC
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#' @useDynLib grpSLOPEMC
#' @importFrom Rcpp sourceCpp
#' @importFrom stats qnorm pchisq qchisq uniroot
NULL
#> NULL

#' Regularizing sequence for Group SLOPE
#'
#' Generate the regularizing sequence \code{lambda} for the Group SLOPE
#' problem according to one of multiple methods.
#'
#' Multiple methods are available to generate the regularizing sequence \code{lambda}:
#' \itemize{
#'   \item "BH" -- method of Theorem 1.1 in Bogdan et. al. (2015)
#'   \item "gaussian" -- method of Section 3.2.2 in Bogdan et. al. (2015)
#'   \item "gaussianMC" -- method introduced in Gossmann et. al. (2015)
#'   \item "chiOrthoMax" -- lambdas as in Theorem 2.5 in Brzyski et. al. (2015)
#'   \item "chiOrthoMean" -- lambdas of equation (2.14) in Brzyski et. al. (2015)
#'   \item "chiEqual" -- Procedure 1 in Brzyski et. al. (2015)
#'   \item "chiMean" -- Procedure 2 in Brzyski et. al. (2015)
#'   \item "chiMC" -- (experimental) A Monte Carlo lambda selection method based on equation (2.25)
#'            in Brzyski et. al. (2015). Requires that rank(\code{A}) is greater than
#'            the sum of the number of elements in any \code{n.MC} groups. 
#' }
#' When \code{method} is "gaussianMC" or "chiMC", the corrections of the entries of lambda will be 
#' computed up to the index given by \code{n.MC} only. \code{n.MC} should be
#' less than or equal to \code{n.group}. Since lambda sequences obtained via MC tend to
#' flatten out quickly, it is reasonable to choose \code{n.MC} to be much smaller than the
#' number of groups.
#'
#' @param fdr Target false discovery rate
#' @param n.group Number of groups
#' @param group A vector describing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#' @param A The model matrix
#' @param y The response variable
#' @param wt A named vector of weights, one weight per group of predictors (named according to names as in vector \code{group})
#' @param n.obs Number of observations (i.e., number of rows in \code{A})
#' @param method Possible values are "BH", "gaussian", "gaussianMC",
#'    "chiOrthoMax", "chiOrthoMean",  "chiEqual", "chiMean", "chiMC". See under Details.
#' @param n.MC When \code{method} is "gaussianMC" or "chiMC", the corrections of the entries of lambda will be 
#'    computed up to the index given by \code{n.MC} only. See details.
#' @param MC.reps The number of repetitions of the Monte Carlo procedure
#'
#' @examples
#' fdr     <- 0.1
#' n.obs   <- 700
#' n.group <- 90
#' # generate 90 groups of sizes 5, 10, and 20
#' group   <- vector()
#' for (i in 1:30) {
#'   tmp <- rep((i-1)*3+c(1,2,3), c(5,10,20))
#'   group <- c(group, tmp)
#' }
#' # set the weight for each group to the square root of the group's size
#' wt <- rep(c(sqrt(5), sqrt(10), sqrt(20)), 30)
#' names(wt) <- names(getGroupID(group))
#' # compute different lambda sequences
#' lambda.BH <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, method="BH")
#' lambda.G <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, n.obs=n.obs, method="gaussian")
#' lambda.max <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, method="chiOrthoMax") 
#' lambda.mean <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, method="chiOrthoMean") 
#' lambda.chi <- lambdaGroupSLOPE(fdr=fdr, n.obs=n.obs, group=group, wt=wt, method="chiMean")
#'
#' @references D. Brzyski, W. Su, M. Bogdan (2015), \emph{Group SLOPE -- adaptive selection of groups of predictors}, \url{http://arxiv.org/abs/1511.09078}
#' @references A. Gossmann, S. Cao, Y.-P. Wang (2015), \emph{Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE}, \url{http://dx.doi.org/10.1145/2808719.2808743}
#' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE -- Adaptive variable selection via convex optimization}, Annals of Applied Statistics
#'
#' @export
lambdaMC <- function(fdr=0.1, group, A, y=NULL, wt=NULL, method,
                     n.MC, MC.reps=5000)
{
  # Prepare grouping information
  group.id <- grpSLOPE::getGroupID(group)
  n.group  <- length(group.id)

  if (n.MC > n.group) {
    warning("n.MC is not allowed to exceed the number of groups.")
    n.MC <- n.group
  }

  if (method=="gaussianMC") {
    return(lambdaGaussianMC(fdr=fdr, n.group=n.group, group.id=group.id,
                            A=A, n.MC=n.MC, MC.reps=MC.reps))
    
  } else if (method == "chiMC") {
    if (is.null(y) || is.null(wt)) {
      stop("y and wt need to be passed as arguments when method is 'chiMC'.")
    }

    return(lambdaChiMC(fdr=fdr, X=A, y=y, group.id=group.id, wt=wt,
                       n.MC=n.MC, MC.reps=MC.reps))

  } else {
    stop(paste(method, "is not a valid method."))
  }
}
