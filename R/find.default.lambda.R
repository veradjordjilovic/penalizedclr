#' Default values for L1 penalty in conditional logistic regression
#'
#' Performs cross validation to determine reasonable values for L1 penalty
#' in a conditional logistic regression.
#'
#'
#' @inheritParams penalized.clr
#'
#'
#' @param pf.list List of vectors of penalty factors.
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @return A single numeric value if \code{p} and \code{pf.list} are missing, or a list of numeric values
#' with L1 penalties for each vector of penalty factors supplied.
#'
#' @export
#'
#' @details
#'
#' The function is based on cross-validation implemented in the `clogitL1` package and returns
#' the value of `lambda` that minimizes cross validated deviance.
#' In the presence of blocks of covariates, a user specifies a list of
#' candidate vectors of penalty factors. For each candidate vector of penalty factors a
#' single `lambda` value is obtained.  Note that
#' cross-validation includes random data splitting, meaning
#' that obtained values can vary significantly between different runs.
#'
#'@seealso [default.pf]
#'@examples
#' set.seed(123)
#' # simulate covariates (pure noise in two blocks of 20 and 80 variables)
#' X <- cbind(matrix(rnorm(4000, 0, 1), ncol = 20), matrix(rnorm(16000, 2, 0.6), ncol = 80))
#' p <- c(20,80)
#' pf.list <- list(c(0.5, 1), c(2, 0.9))
#' # stratum membership
#' stratum <- sort(rep(1:100, 2))
#'
#' # the response
#' Y <- rep(c(1, 0), 100)
#'
#' # obtain a list with vectors of penalty factors
#' \donttest{
#' lambda.list <- find.default.lambda(response = Y,
#'                                    penalized = X, stratum = stratum, p = p, pf.list = pf.list)}
#'
#' # when `p` and `pf.list` are not provided all covariates are treated as a single block
#' \donttest{
#' lambda <- find.default.lambda(response = Y,
#'                                    penalized = X, stratum = stratum)}




find.default.lambda <- function(response, stratum, penalized,
                                unpenalized = NULL,
                                alpha = 1,
                                p = NULL,
                                standardize = TRUE,
                                event,
                                pf.list = NULL,
                                nfolds = 10){

  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized)) {X <- cbind(unpenalized, penalized)
  nc <- ncol(unpenalized)} else{
    X <- penalized
    nc <- 0
  }

  #if(standardize == TRUE) X <- apply(X, 2, function(x)
  #  x/sqrt((length(x)-1)/length(x)*var(x)))

  if (standardize == TRUE) X <- scale(X)

  if(missing(p) | is.null(p)){
    lambda.seq <- default.lambda(X, response, stratum, alpha, nfolds = nfolds)
  } else{
    if (missing(pf.list)) stop("Argument pf.list is missing with no default.")
    lambda.seq <- list()
    for(i in 1:length(pf.list)){
      w <- as.numeric(pf.list[[i]])
      weights <- c(rep(1, nc),
                   rep(1/w,  p))
      X_w <- t(apply(X, 1, function(x)  x* weights))
      lambda.seq[[i]] <- default.lambda(X = X_w,
                                   response,
                                   stratum,
                                   alpha,
                                   nfolds = nfolds)
    }
  }

  return("lambda.seq" = lambda.seq)
}
