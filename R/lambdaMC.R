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

#' Monte Carlo based regularizing sequences for Group SLOPE
#'
#' Generate the regularizing sequence \code{lambda} for the Group SLOPE
#' problem according to one of multiple Monte Carlo based methods.
#'
#' Two Monte Carlo based methods are currently available to generate the regularizing sequence \code{lambda}:
#' \itemize{
#'   \item "gaussianMC" -- method introduced in Gossmann et. al. (2015)
#'   \item "chiMC" -- A Monte Carlo lambda selection method based on equation (G.10)
#'            in Brzyski et. al. (2016). Requires that rank(\code{A}) is greater than
#'            the sum of the number of elements in any \code{n.MC} groups. 
#' }
#' The corrections of the entries of lambda will be computed up to the index given 
#' by \code{n.MC} only (the remaining tail of the returned sequence will be flat). 
#' \code{n.MC} should be less than or equal to \code{n.group}. 
#' Since lambda sequences obtained via MC tend to flatten out quickly, 
#' it is reasonable to choose \code{n.MC} to be much smaller than the
#' number of groups.
#'
#' @param fdr Target gFDR (group false discovery rate)
#' @param group A vector describing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#' @param A The model matrix
#' @param y The response variable
#' @param wt A named vector of weights, one weight per group of predictors (named according to names as in vector \code{group})
#' @param method Possible values are "gaussianMC" and "chiMC". See details.
#' @param n.MC The corrections of the entries of lambda will be 
#'    computed up to the index given by \code{n.MC} only. See details.
#' @param MC.reps The number of repetitions of the Monte Carlo procedure
#'
#' @examples
#' # generate some data
#' A <- matrix(rnorm(20000), 200, 100)
#' b <- c(rep(1, 10), rep(0, 90))
#' y <- A %*% b
#' # randomly divide 100 predictors into 30 groups
#' group <- sample(1:30, 100, repl = TRUE)
#' # set the weight for each group to the square root of the group's size
#' wt <- c(sqrt(table(group)))
#' # compute different lambda sequences
#' gaussianMC <- lambdaMC(fdr=0.1, group=group, A=A, method="gaussianMC",
#'                        n.MC=10, MC.reps=100)
#' chiMC <- lambdaMC(fdr=0.2, group=group, A=A, method="chiMC", 
#'                   y=y, wt=wt, n.MC=10, MC.reps=100)
#'
#' @references D. Brzyski, A. Gossmann, W. Su, M. Bogdan (2016), \emph{Group SLOPE - adaptive selection of groups of predictors}, \url{https://arxiv.org/abs/1610.04960}
#' @references A. Gossmann, S. Cao, Y.-P. Wang (2015), \emph{Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE}, \url{http://dx.doi.org/10.1145/2808719.2808743}
#'
#' @export
lambdaMC <- function(fdr=0.1, group, A, y=NULL, wt=NULL, method,
                     n.MC, MC.reps=5000)
{
  # Prepare grouping information
  group.id <- grpSLOPE::getGroupID(group)
  n.group  <- length(group.id)

  is.pos.int <- function(x) { abs(x) - round(x) < .Machine$double.eps^0.5 }
  if (!(is.pos.int(n.MC) && is.pos.int(MC.reps))) {
    stop("n.MC and MC.reps have to be positive integers.")
  }
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
