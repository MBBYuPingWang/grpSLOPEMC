##############################################################################
#
#    grpSLOPEMC
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and#or modify
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
#    along with this program.  If not, see <http:##www.gnu.org#licenses#>.
#
##############################################################################

# method of Theorem 1.1 in Bogdan et. al. (2015)
lambdaBH <- function(fdr, n.group) {
  lambda.BH <- rep(NA,n.group)
  for (i in 1:n.group){
    lambda.BH[i] <- qnorm(1-(i*fdr)/(2*n.group))
  }

  return(lambda.BH)
}

# method introduced in Gossmann et. al. (2015)
lambdaGaussianMC <- function(fdr, n.group, group.id, A, n.MC, MC.reps) {
  lambda.BH <- lambdaBH(fdr=fdr, n.group=n.group)

  mA <- matrix(NA, nrow(A), n.group)
  for (i in 1:n.group) {
    mA[ , i] <- apply(as.matrix(A[ , group.id[[i]] ]), 1, mean)
    mA[ , i] <- mA[ , i] / sqrt(sum(mA[ , i]^2))
  }

  # Monte Carlo corrections for lambda.BH
  lambda.MC <- lambdaGaussianMCRcpp(lambda.BH, mA, n.MC, MC.reps)
  lambda.MC <- c(lambda.MC, rep(lambda.MC[n.MC], n.group-n.MC))

  return(lambda.MC)
}

# approximates the variance of (G.10) in Brzyski et. al. (2016) via Monte Carlo
lambdaChiMC <- function(fdr, X, y, group.id, wt, n.MC, MC.reps) {
  n.group     <- length(group.id)
  group.sizes <- sapply(group.id, FUN=length)
  lambda.MC   <- vector()

  # make sure weights are in the same order as group.id
  wt <- wt[names(group.id)]

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  # get upper and lower bounds for lambda.MC[1]
  qchi.seq <- rep(NA, n.group)
  for (j in 1:n.group) {
    qchi.seq[j] <- sqrt(qchisq(1 - fdr/n.group, df=group.sizes[j])) / wt[j]
  }
  upperchi <- max(qchi.seq)
  lowerchi <- min(qchi.seq)

  if (upperchi == lowerchi) {
    lambda.MC[1] <- upperchi
  } else {
    lambda.MC[1] <- uniroot(function(y) (cdfMean(y) - (1-fdr/n.group)),
                            lower = lowerchi, upper = upperchi, extendInt="yes")$root
  }

  # get lambda.MC[2:n.MC]
  for (i in 2:n.MC) {
    s <- lambdaChiMCAdjustmentRcpp(y=y, X=X, group_id=group.id, lambda=lambda.MC,
                                   w=wt, number_of_drawings=MC.reps)

    cdfMean <- function(x) {
      pchi.seq <- rep(NA, n.group)
      for (j in 1:n.group) {
        pchi.seq[j] <- pchisq((wt[j]/s * x)^2, df=group.sizes[j])
      }
      return(mean(pchi.seq))
    }

    cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*i/n.group)),
                           lower = 0, upper = upperchi, extendInt="upX")$root

    if (cdfMean.inv <= lambda.MC[i-1]) {
      lambda.MC[i] <- cdfMean.inv 
    } else {
      lambda.MC[i] <- lambda.MC[i-1] 
    }
  }

  # get lambda.MC[n.MC:n.group]
  lambda.MC <- c(lambda.MC, rep(lambda.MC[n.MC], n.group-n.MC))
  
  return(lambda.MC)
}
