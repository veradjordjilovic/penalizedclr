#' Data adaptive candidate vector of penalty factors for L1 penalty in conditional logistic regression
#'
#' Performs cross validation to determine reasonable default values for L1 penalty
#' in a conditional logistic regression
#'
#'
#' @inheritParams penalized.clr
#' @param X A matrix of covariates, with the number of rows equaling the number
#' of observations.
#' @param Y A binary response variable.
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @param type.step1 Should the tentative model be fit on all covariates jointly or to each block separately.
#' @return A numeric vector of length 1 to 3 (depending on the problem)
#' giving L1 penalties
#'
#'




default.pf <- function(X, Y, stratum, nfolds = 10, alpha = 1, standardize = TRUE,
                           type.step1, p){
  if (alpha == 1 & missing(type.step1)) {
    warning("type.step1 is set to sep")
    type.step1 <- "sep"
  }
  if (alpha == 0 & missing(type.step1)) {
    warning("type.step1 is set to comb")
    type.step1 <- "comb"
  }
  if (alpha > 0 & alpha < 1 & missing(type.step1)) {
    stop("type.step1 should be sep or comb")
  }

  if (standardize) X <- apply(X, 2, function(x)
    x/sqrt((length(x)-1)/length(x)*var(x)))

  blocks <- rep(1:length(p),  p)

  if (type.step1 == "comb"){
    fit_t <- clogitL1::clogitL1(X, Y, stratum, alpha = alpha)
    cvfit_t <- clogitL1::cv.clogitL1(fit_t, numFolds = nfolds)
    ind <- which(cvfit_t$lambda == cvfit_t$minCV1se_lambda)
    beta_est_abs <- abs(fit_t$beta[ind, ])
    means <- as.numeric(tapply(beta_est_abs, blocks, mean))

  }

  if (type.step1 == "sep") {
    means <- rep(0, length(p))
    for (i in 1:length(p)) {
      fit_t <- clogitL1::clogitL1(X[, blocks == i], Y, stratum, alpha = alpha)
      cvfit_t <- clogitL1::cv.clogitL1(fit_t, numFolds = nfolds)
      ind <- which(cvfit_t$lambda == cvfit_t$minCV1se_lambda)
      beta_est_abs <- abs(fit_t$beta[ind, ])
      means[i] <- mean(beta_est_abs)
      }
  }
  exc <- NULL
  if (all(means!=0)) {pf <- 1/means} else{
    exc <- which(means==0)
    pf <- 1/means[means!=0]
  }
  if (is.null(exc)){return(list("The data adaptive vector of penalty factors is ", pf = pf))}else{
    return(list("The data adaptive vector of penalty factors is ", pf = pf, "and the following modality(ies) have been excluded ", exc = exc))
  }

}
