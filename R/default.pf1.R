#' Data adaptive candidate vector of penalty factors for L1 penalty in conditional logistic regression with covariates divided in blocks
#'
#' Computes a data adaptive vector of penalty factors for  blocks of covariates by fitting
#' a tentative penalized conditional logistic regression model. The penalty for the `i`th block is obtained
#' as the inverse of the arithmetic mean of coefficient estimates for its covariates.
#'
#'
#' @inheritParams penalized.clr
#' @param X A matrix of covariates, with the number of rows equaling the number
#' of observations.
#' @param Y A binary response variable.
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @param type.step1 Should the tentative model be fit on all covariates jointly or to each block separately.
#' @return The function returns a list containing the  vector of penalty factors correspondng to different blocks.
#'
#' @references  Schulze G. (2017) Clinical Outcome Prediction based on Multi-Omics Data: Extension of IPF-LASSO. Master Thesis.
#'
#' @export
#' @details Blocks that contain covariates with large estimated coefficients will be penalized less in the second step.
#' If all estimated coefficients pertaining to a block are zero, the functions returns a message.
#' A tentative conditional logistic regression model is fit either to each covariates block separately (`type.step1 = "sep"`) or  jointly to all blocks (`type.step1 = "comb"`).
#' @seealso [find.default.lambda]





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

 # if (standardize) X <- apply(X, 2, function(x)
 #   x/sqrt((length(x)-1)/length(x)*var(x)))

  if (standardize) X <- scale(X)

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
  if (is.null(exc)){print("The data adaptive vector of penalty factors is ")
    return(list(pf = pf))}else{
   print("The data adaptive vector of penalty factors and the  modality(ies) that have been excluded are ")
    return(list(pf = pf, exc = exc))
  }

}
