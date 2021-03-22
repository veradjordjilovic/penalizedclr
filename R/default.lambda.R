#' Default values for L1 penalty in conditional logistic regression
#'
#' Performs cross validation to determine reasonable default values for L1 penalty
#' in a conditional logistic regression
#'
#'
#' @inheritParams penalized.clr
#' @param X A matrix of covariates, with the number of rows equalling the number
#' of observations.
#' @param Y A binary response variable.
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @return A numeric vector of length 1 to 3 (depending on the problem)
#' giving L1 penalties
#'
#'




default.lambda <- function(X, Y, stratum, nfolds = 10, alpha = 1){
  fit <- clogitL1::clogitL1(X, Y, stratum, alpha = alpha)
  cvfit <- clogitL1::cv.clogitL1(fit, numFolds = nfolds)
  d_lambda <- unique(c(
                exp((3*cvfit$minCV_lambda - cvfit$minCV1se_lambda)/2),
                exp(cvfit$minCV_lambda),
                exp((cvfit$minCV_lambda + cvfit$minCV1se_lambda)/2),
                exp(cvfit$minCV1se_lambda)))
  return(d_lambda)
}
