#' Data adaptive candidate vector of penalty factors for L1 penalty in conditional logistic regression with covariates divided in blocks
#'
#' Computes a data adaptive vector of penalty factors for  blocks of covariates by fitting
#' a tentative penalized conditional logistic regression model. The penalty for the `i`th block is obtained
#' as the inverse of the arithmetic mean of coefficient estimates for its covariates.
#'
#'
#' @inheritParams penalized.clr
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @param type.step1 Should the tentative model be fit on all covariates jointly (`comb`) or to each block separately (`sep`).
#' @param verbose Logical. Should the message about the obtained penalty factors be printed?
#' @return The function returns a list containing the  vector of penalty factors correspondng to different blocks.
#'
#' @references  Schulze G. (2017) Clinical Outcome Prediction based on Multi-Omics Data: Extension of IPF-LASSO. Master Thesis.
#'
#' @export
#' @details Blocks that contain covariates with large estimated coefficients will obtain a smaller penalty.
#' If all estimated coefficients pertaining to a block are zero, the function returns a message.
#' A tentative conditional logistic regression model is fit either to each covariates block separately (`type.step1 = "sep"`) or  jointly to all blocks (`type.step1 = "comb"`).
#' Note that `unpenalized = NULL` is the only implemented option  in this function as of now.
#' @seealso [find.default.lambda]





default.pf <- function(response,
                       stratum,
                       penalized,
                       unpenalized = NULL,
                       alpha = 1,
                       p = NULL,
                       standardize = TRUE,
                       event,
                       nfolds = 10,
                       type.step1,
                       verbose = FALSE){
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


  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (standardize) {X <- scale(penalized)}else {X <- penalized}

  blocks <- rep(1:length(p),  p)

  if (type.step1 == "comb"){
    fit_t <- clogitL1::clogitL1(X, response, stratum, alpha = alpha)
    cvfit_t <- clogitL1::cv.clogitL1(fit_t, numFolds = nfolds)
    ind <- which(cvfit_t$lambda == cvfit_t$minCV1se_lambda)
    beta_est_abs <- abs(fit_t$beta[ind, ])
    means <- as.numeric(tapply(beta_est_abs, blocks, mean))

  }

  if (type.step1 == "sep") {
    means <- rep(0, length(p))
    for (i in 1:length(p)) {
      fit_t <- clogitL1::clogitL1(X[, blocks == i], response, stratum, alpha = alpha)
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

  if (is.null(exc)){if (verbose)
    {message(paste0("The data adaptive vector of penalty factors is ", pf, collapse = ""))}
    return(list(pf = pf))}
  else{if(verbose){
   message(paste0("The data adaptive vector of penalty factors is ", pf," and the modality(ies) that have been excluded are ", exc, collapse =""))}
    return(list(pf = pf, exc = exc))
  }

}
