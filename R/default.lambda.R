#' Default values for L1 penalty in conditional logistic regression
#'
#' Internal function that performs cross validation to determine reasonable default values for L1 penalty
#' in a conditional logistic regression
#'
#'
#' @inheritParams penalized.clr
#' @param X A matrix of covariates, with the number of rows equaling the number
#' of observations.
#' @param Y A binary response variable.
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @return A numeric value for `lambda` minimizing cross validated deviance.
#'
#'




default.lambda <- function(X, Y, stratum, nfolds = 10, alpha = 1){
  fit <- clogitL1::clogitL1(X, Y, stratum, alpha = alpha)
  cvfit <- clogitL1::cv.clogitL1(fit, numFolds = nfolds)
  d_lambda <- exp(cvfit$minCV_lambda)

  return(d_lambda)
}
