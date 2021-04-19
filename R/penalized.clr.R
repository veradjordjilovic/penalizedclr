#' Penalized conditional logistic regression
#'
#' Fits conditional logistic regression models with L1 and L2 penalty allowing
#' for different penalties for different blocks of covariates.
#'
#' @param response The response variable, either a 0/1 vector or a factor with two levels.
#' @param stratum A numeric vector with stratum membership of each observation.
#' @param penalized A matrix of penalized covariates.
#' @param unpenalized A matrix of additional unpenalized covariates.
#' @param lambda The tuning parameter for L1. Either a single non-negative number,
#'  or a numeric vector of the length equal to the number of blocks. If NULL, function `find.default.lambda` is called.  See p below.
#' @param alpha The elastic net mixing parameter, a number between 0 and 1.
#'              alpha=0 would give pure ridge; alpha=1 gives lasso. Pure ridge penalty is never obtained in this implementation since alpha must be positive.
#' @param p The sizes of blocks of covariates,
#'         a numerical vector of the length equal to the number of blocks,
#'         and with the sum equal to the number of penalized covariates.
#'         If missing, all covariates are treated the same and a single penalty is applied.
#' @param standardize Should the covariates be standardized, a logical value.
#' @param event  If response is a factor, the level that
#'              should be considered a success in the logistic regression.
#' @return A list with the following elements:
#'              \itemize{
#'                   \item \code{penalized} - Regression coefficients for the penalized covariates.
#'                   \item \code{unpenalized} - Regression coefficients for the unpenalized covariates.
#'                   \item \code{converged} - Whether the fitting process was judged to be converged.
#'                   \item \code{lambda} - The tuning parameter for L1 used.
#'                   \item \code{alpha} - The elastic net mixing parameter used.
#'                       }
#'
#' @seealso \code{\link{stable.clr}} and \code{\link{stable.clr.g}} for variable selection through stability selection
#'    in penalized conditional logistic regression with a single penalty factor and multiple penalty factors, respectively.
#'
#' @examples
#' set.seed(123)
#' # simulate covariates (pure noise in two blocks of 20 and 80 variables)
#' X <- cbind(matrix(rnorm(4000, 0, 1), ncol = 20), matrix(rnorm(16000, 2, 0.6), ncol = 80))
#'
#' # stratum membership
#' stratum <- sort(rep(1:100, 2))
#'
#' # the response
#' Y <- rep(c(1, 0), 100)
#'
#' fit <- penalized.clr( response = Y, stratum = stratum,
#'   penalized = X, lambda = c(1, 0.3),
#'   p = c(20, 80), standardize = TRUE)
#' fit$penalized
#' fit$converged
#' fit$lambda
#' @export
#' @details The `penalized.clr` function fits a conditional logistic regression
#'  model for a given combination of L1 (`lambda`) and L2 penalties. L2 penalty is
#'  obtained from `lambda` and `alpha` as `lambda*(1-alpha)/(2*alpha)`.
#'  Note that `lambda` is a single number if all covariates are to be penalized
#'  equally, and a vector of penatlies, if predictors are divided in blocks (of sizes provided in
#'  `p`) that are to be penalized differently.  If `lambda` is not provided by the user,
#'  a default value is computed by the `find.default.lambda`
#'  function (which slows down the computation). The `penalized.clr` function
#'  is based on the Cox model routine available in the
#'  `penalized` package.
#' @importFrom survival strata
#' @importFrom stats var
#'
#'


penalized.clr <- function(response,
                          stratum,
                          penalized,
                          unpenalized = NULL,
                          lambda = NULL,
                          alpha = 1,
                          p = NULL,
                          standardize = FALSE,
                          event) {

  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized) && !is.numeric(dim(unpenalized))) {
    unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
  }
  if (length(p) > 0 && sum(p) != ncol(penalized)) stop("elements of p must sum to the number of penalized covariates.")



  if(is.null(lambda)){
    lambda <- find.default.lambda(response,
                                  stratum,
                                  penalized,
                                  unpenalized,
                                  alpha,
                                  p,
                                  standardize,
                                  event)
    if (is.numeric(lambda)) {lambda <- lambda[1]} else{
      lambda <- unlist(lapply(lambda, function(x) x[1]))
    }
    }else{
    if (length(lambda) > 1 && (missing(p) | is.null(p))) stop("multiple penalties are supplied in lambda, but p is missing.")
    if (length(p) != length(lambda) && !missing(p) && !is.null(p)) stop("lambda and p are not of the same length.")
    if (sum(lambda < 0) > 0) stop("lambda must be non-negative.")
  }


  if (missing(p) | is.null(p)) {
    lambda1 <- lambda
  } else {
    lambda1 <- rep(lambda, p)
  }

  if (alpha <= 0 | alpha > 1) stop("alpha must be between zero and one.")
  lambda2 <- lambda1 * (1 - alpha) / (2 * alpha)


  Y <- survival::Surv(rep(1, length(response)),
    event = (response == 1)
  )



  if (is.null(unpenalized)) {
    fit <- penalized::penalized(Y ~ strata(stratum),
      penalized = penalized,
      lambda1 = lambda1,
      lambda2 = lambda2,
      standardize = standardize
    )
  }
  else {
    fit <- penalized::penalized(Y ~ strata(stratum) + unpenalized,
      penalized = penalized,
      lambda1 = lambda1,
      lambda2 = lambda2,
      standardize = standardize
    )
  }
  return(list(
    penalized = fit@penalized,
    unpenalized = fit@unpenalized,
    converged = fit@converged,
    lambda = unique(lambda1),
    alpha = alpha,
     Y=Y))
}
