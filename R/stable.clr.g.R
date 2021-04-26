#' Stability selection based on penalized conditional logistic regression
#'
#' Performs stability selection for conditional logistic regression models with L1 and L2 penalty
#' allowing for different penalties for different blocks
#' (groups) of covariates (different data sources).
#'
#'
#'
#'
#' @inheritParams penalized.clr
#' @inheritParams subsample.clr
#' @inheritParams stable.clr
#' @param lambda.list A list with per-block sequences of L1 penalties.
#'                    If NULL, `find.default.lambda` function is called.
#'
#' @return A list containing a numeric vector \code{Pilambda},
#'         giving selection probabilities for all penalized covariates and
#'         `lambda.list`.
#'
#'@seealso
#'    \code{\link{find.default.lambda}} for obtaining default sequences of L1 penalties.
#'
#'
#' @references  1. Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473.
#'
#' 2. Shah, R. D., & Samworth, R. J. (2013). Variable selection with error control:
#' another look at stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(1), 55-80.
#' @export
#' @details This function implements stability selection (Meinshausen and Bühlmann, 2010) in
#' a conditional logistic regression. The implementation is based on the modification of Shah and
#' Samworth (2013) featuring complementary subsamples. Note that this means that the number
#' of subsamples will be `2B` instead of `B`. Subsampling procedure is repeated
#' `2B` times for each combination of per-block penalties resulting each time in a vector of
#' selection frequences (frequency of non-zero coefficient estimate of each covariate).
#' The final selection probability `Pilambda` is obtained by taking the maximum over
#' all considered values of penalties.
#'
#'@examples
#' set.seed(123)
#'
#' # simulate covariates (pure noise in two blocks of 20 and 80 variables)
#' X <- cbind(matrix(rnorm(4000, 0, 1), ncol = 20), matrix(rnorm(16000, 2, 0.6), ncol = 80))
#' p <- c(20,80)
#'
#' # stratum membership
#' stratum <- sort(rep(1:100, 2))
#'
#' # the response
#' Y <- rep(c(1, 0), 100)
#'
#' # list of L1 penalties
#' lambda.list <- find.default.lambda(response = Y,
#'                                    penalized = X, stratum = stratum, p = p)
#'
#' # perform stability selection
#' \donttest{
#' stable.g1 <- stable.clr.g(response = Y, penalized = X, stratum = stratum,
#'                          p = p, lambda.list = lambda.list)}
#'
#' # when lambda.list is not provided,
#' # it is computed within the function (slightly different results might occur due to the
#' # randomness inherent to cross-validation)
#'\donttest{
#' stable.g2 <- stable.clr.g(response = Y, penalized = X, stratum = stratum,
#'                           p = p)}
#'
#' # if p is not provided, all covariates are penalized equally
#' \donttest{
#' stable.g3 <- stable.clr.g(response = Y, penalized = X, stratum = stratum)}
#'


stable.clr.g <- function(response,
                         stratum,
                         penalized,
                         unpenalized = NULL,
                         p = NULL,
                         lambda.list = NULL,
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

  if (missing(p)| is.null(p) | length(p) == 1){
    warning("valid p is not provided:
            all covariates are penalized equally.")
    if (is.null(lambda.list)) {

    temp <- stable.clr(response, stratum, penalized,
                       unpenalized, lambda.seq = NULL,
                       alpha, B, parallel, standardize, event)}else{
    temp <- stable.clr(response, stratum, penalized,
                                            unpenalized, lambda.seq = unlist(lambda.list),
                                            alpha, B, parallel, standardize, event)
                       }
    Pistab <- temp$Pistab
    lambda.list <- temp$lambda.seq
    }else{

  g <- length(p) # the number of groups

  ind.pair <- unique(stratum)
  b <- length(ind.pair)
  subsample.size <- ceiling(b / 2)

  matB <- matrix(0, nrow = 2 * B, ncol = subsample.size)
  for (i in 1:B) {
    matB[i, ] <- sample(ind.pair,
      size = subsample.size,
      replace = FALSE
    )
    matB[B + i, 1:(b - subsample.size)] <- setdiff(ind.pair, matB[i, ])
  }

  if (is.null(lambda.list) | missing(lambda.list)) {
    lambda.list <- find.default.lambda(response,
                                       stratum,
                                       penalized,
                                       unpenalized,
                                       alpha,
                                       p,
                                       standardize)
  }
  if (parallel) {
    cl <- parallel::makeCluster(getOption("cl.cores", 2), setup_timeout = 0.5)
    parallel::clusterExport(cl, varlist = c("penalized.clr"))

    subslist <- c(list(B = 1:(2 * B)), lambda.list)
    P <- expand.grid(subslist)

    res <- t(parallel::parApply(cl,
      P,
      1,
      FUN = function(x,
                     response, stratum, penalized, unpenalized,
                     p, matB, alpha, standardize) {
            ind <- stratum %in% matB[x[1], ]
            (penalized.clr(
              response = response[ind],
              stratum = stratum[ind],
              penalized = penalized[ind, ],
              unpenalized = unpenalized[ind, ],
              p = p,
              lambda = as.numeric(x[-1]),
              alpha = alpha,
              standardize = standardize)$penalized != 0) * 1
      },
      response,
      stratum,
      penalized,
      unpenalized,
      p,
      matB,
      alpha,
      standardize))

    res1 <- as.data.frame(cbind(P, res))
    res2 <- stats::aggregate(res, by = as.list(res1[, 2:(g + 1)]), mean)
    res <- t(res2[, -c(1:g)])
    parallel::stopCluster(cl)
  } else {
    for (i in 1:nrow(P)) {
      ind <- stratum %in% matB[P[i, 1], ]
      fit <- penalized.clr(response[ind],
        stratum[ind],
        penalized = penalized[ind, ],
        unpenalized = unpenalized[ind, ],
        lambda = as.numeric(P[i, -1]),
        alpha = alpha,
        standardize = standardize )
      selB[i, ] <- (fit$penalized != 0) * 1
      colnames(selB) <- names(fit$penalized)
    }
  }

  Pistab <- apply(res, 1, max)
  names(Pistab) <- colnames(penalized)}

  return(list(Pistab = Pistab, lambda.list = lambda.list))
}
