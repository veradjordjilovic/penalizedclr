#' Stability selection based on penalized conditional logistic regression
#'
#' Performs stability selection for conditional logistic regression models with
#' L1 and L2 penalty.
#'
#'
#' @inheritParams penalized.clr
#' @inheritParams subsample.clr
#' @param lambda.seq a sequence of non-negative values to be used as tuning
#'    parameters for L1
#'
#' @return A list with a  numeric vector \code{Pilambda}
#'         giving selection probabilities for each penalized covariate, and
#'         a sequence \code{lambda.seq} used.
#'
#'
#' @export
#' @examples
#'  set.seed(123)
#'
#' # simulate covariates (pure noise in two blocks of 20 and 80 variables)
#' X <- cbind(matrix(rnorm(4000, 0, 1), ncol = 20), matrix(rnorm(16000, 2, 0.6), ncol = 80))
#'
#' # stratum membership
#' stratum <- sort(rep(1:100, 2))
#'
#' # the response
#' Y <- rep(c(1, 0), 100)
#'
#' # sequence of L1 penalties
#' lambda.seq <- find.default.lambda(response = Y,
#'                                    penalized = X,
#'                                    stratum = stratum)
#'
#' # perform stability selection
#' \donttest{
#' stable1 <- stable.clr(response = Y, penalized = X, stratum = stratum,
#'                          lambda.seq = lambda.seq)}
#'
#' # when lambda.seq is not provided,
#' # it is computed within the function (slightly different results might occur due to the
#' # randomness inherent to cross-validation)
#'
#'\donttest{
#' stable2 <- stable.clr.g(response = Y, penalized = X, stratum = stratum)}
#'
#'
#' @seealso  \code{\link{stable.clr.g}} for stability selection
#' in penalized conditional logistic regression with multiple penalties for block structured covariates.


stable.clr <- function(response,
                       stratum,
                       penalized,
                       unpenalized = NULL,
                       lambda.seq = NULL,
                       alpha = 1,
                       B = 100,
                       parallel = TRUE,
                       standardize = FALSE,
                       event) {
  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized) && !is.numeric(dim(unpenalized))) {
    unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
  }

  if (is.null(lambda.seq)){
    lambda.seq <- find.default.lambda(response,
                                      stratum,
                                      penalized,
                                      unpenalized,
                                      alpha,
                                      p = NULL,
                                      standardize,
                                      event)
  }

  fit <- subsample.clr(
    response = response,
    stratum = stratum,
    penalized = penalized,
    unpenalized = unpenalized,
    lambda = lambda.seq[1],
    alpha = alpha,
    B = B,
    matB = NULL,
    return.matB = TRUE,
    parallel = parallel,
    standardize = standardize
  )
  matB <- fit$matB
  if (parallel) {
    cl <- parallel::makeCluster(getOption("cl.cores", 2), setup_timeout = 0.5)
    parallel::clusterExport(cl, varlist = c("penalized.clr"))

    P <- expand.grid("B" = 1:nrow(matB), "lambda" = lambda.seq[-1])

    res <- t(parallel::parApply(cl,
      P,
      1,
      FUN = function(x,
                     response,
                     stratum,
                     penalized,
                     unpenalized,
                     matB,
                     alpha,
                     standardize) {
        #require("penalized")
        ind <- stratum %in% matB[x[1], ]

        (penalized.clr(
          response = response[ind],
          stratum = stratum[ind],
          penalized = penalized[ind, ],
          unpenalized = unpenalized[ind, ],
          lambda = x[2],
          alpha = alpha,
          standardize = standardize
        )$penalized != 0) * 1
      },
      response,
      stratum,
      penalized,
      unpenalized,
      matB,
      alpha,
      standardize
    ))
    res1 <- as.data.frame(cbind(P, res))
    res2 <- stats::aggregate(
      res,
      list(lambda = res1$lambda),
      mean
    )
    res <- t(rbind(fit$Pilambda, res2[, -c(1)]))
    parallel::stopCluster(cl)
  } else {
    res <- subsample.clr.v(
      response = response,
      stratum = stratum,
      penalized = penalized,
      unpenalized = unpenalized,
      lambda = lambda.seq[-1],
      alpha = alpha,
      B = B,
      matB = fit$matB,
      parallel = FALSE,
      standardize = standardize
    )
    res <- cbind(fit$Pilambda, res)
  }

  Pistab <- apply(res, 1, max)
  names(Pistab) <- names(fit$Pilambda)
  return(list(Pistab = Pistab, lambda.seq = lambda.seq))
}
