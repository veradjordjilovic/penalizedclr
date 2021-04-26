#' Default values for L1 penalty in conditional logistic regression
#'
#' Performs cross validation to determine reasonable default values for L1 penalty
#' in a conditional logistic regression
#'
#'
#' @inheritParams penalized.clr
#'
#' @param nfolds The number of folds used in cross-validation. Default is 10.
#' @return A numeric vector
#' giving a  sequence of L1 penalties if \code{p} is missing, or a list
#' of per-block penalty sequences otherwise.
#'
#' @export
#'
#' @details The function `find.default.lambda` is used to obtain a sequence of reasonable values
#' of the parameter `lambda`. In the presence of blocks of covariates, a separate
#' sequence is obtained for each block and the function returns a list. The function is
#' based on cross-validation implemented in the `clogitL1` package. For each block,
#' a separate conditional logistic model is fit (together with unpenalized covariates if provided)
#' and  two `lambda` values a) `minCV_lambda` minimizing cross-validated deviance  and b) `minCV1se_lambda`
#' satisfying 1-SE rule, are obtained. Two additional values of `lambda` are obtained
#' as a midpoint between the two `(minCV_lambda + minCV1se_lambda)/2` and symmetrically `(3*minCV_lambda - minCV1se_lambda)/2)`
#' (on a log-scale). If `minCV_lambda = minCV1se_lambda`, a shorter sequnce is returned. Note that
#' cross-validation includes random data splitting, meaning
#' that obtained values can vary significantly between different runs.
#'
#'
#'@examples
#' set.seed(123)
#' # simulate covariates (pure noise in two blocks of 20 and 80 variables)
#' X <- cbind(matrix(rnorm(4000, 0, 1), ncol = 20), matrix(rnorm(16000, 2, 0.6), ncol = 80))
#' p <- c(20,80)
#' # stratum membership
#' stratum <- sort(rep(1:100, 2))
#'
#' # the response
#' Y <- rep(c(1, 0), 100)
#'
#' # obtain a list with a separate sequence for each block
#' \donttest{
#' lambda.list <- find.default.lambda(response = Y,
#'                                    penalized = X, stratum = stratum, p = p)}
#'
#' # obtain a single sequence
#' \donttest{
#' lambda.seq <- find.default.lambda(response = Y,
#'                                    penalized = X, stratum = stratum)}




find.default.lambda <- function(response, stratum, penalized,
                                unpenalized = NULL,
                                alpha = 1,
                                p = NULL,
                                standardize = FALSE,
                                event,
                                nfolds = 10){

  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized)) {X <- cbind(unpenalized, penalized)
  nc <- ncol(unpenalized)} else{
    X <- penalized
    nc <- 0
  }

  if(standardize == T) X <- apply(X, 2, function(x)
    x/sqrt((length(x)-1)/length(x)*var(x)))

  if(missing(p) | is.null(p)){
    lambda.seq <- default.lambda(X, response, stratum, alpha, nfolds = nfolds)
  } else{
    g <- length(p)
    le <- c(1, p[-g]+1) + nc
    ue <- cumsum(p) + nc
    lambda.seq <- list()
    for(i in 1:g){
      lambda.seq[[i]] <- default.lambda(X = X[, c((nc != 0)*(1:nc),le[i]:ue[i])],
                                   response,
                                   stratum,
                                   alpha,
                                   nfolds = nfolds)
    }
  }

  return("lambda.seq" = lambda.seq)
}
