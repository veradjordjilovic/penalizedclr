#' Stability selection for penalized conditional logistic regression
#'
#' Internal function used by `stable.clr` and `stable.clr.g`.
#'
#' @inheritParams penalized.clr
#' @param B A single positive number for the number of subsamples.
#' @param matB  A 2B x ceiling(unique(stratum)/2) matrix with index set of selected strata in each of 2B subsamples
#' @param return.matB Logical. Should the matrix matB be returned?
#' @param parallel Logical. Should the computation be parallelized?
#'
#' @return If `return.matB` is TRUE, a list with two elements, a numeric vector `Pistab`,
#'         giving selection probabilities for each covariate and a matrix `matB`;
#'         otheriwise only `Pistab`.
#'
#'



subsample.clr <- function(response,
                          stratum,
                          penalized,
                          unpenalized = NULL,
                          lambda,
                          alpha = 1,
                          B = 100,
                          matB = NULL,
                          return.matB = FALSE,
                          parallel = TRUE,
                          standardize = TRUE) {
  ind.pair <- unique(stratum)
  b <- length(ind.pair)
  subsample.size <- ceiling(b / 2)

  if (is.null(matB)) {
    matB <- matrix(0, nrow = 2 * B, ncol = subsample.size)
    for (i in 1:B) {
      matB[i, ] <- sample(ind.pair,
        size = subsample.size,
        replace = FALSE
      )
      matB[B + i, 1:(b - subsample.size)] <- setdiff(ind.pair, matB[i, ])
    }
  }

  if (is.null(unpenalized)) {
    selB <- matrix(0, ncol = ncol(penalized), nrow = 2 * B)
  } else {
    if (!is.numeric(dim(unpenalized))) {
      unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
    }
    selB <- matrix(0, ncol = ncol(penalized), nrow = 2 * B)
  }

  if (parallel) {
    cl <- parallel::makeCluster(getOption("cl.cores", 2), setup_timeout = 0.5)
    parallel::clusterExport(cl, varlist = c("penalized.clr"))
    selB <- t(parallel::parSapply(cl,
      X = 1:(2 * B),
      FUN = function(x,
                     response,
                     stratum,
                     penalized,
                     unpenalized,
                     lambda, alpha, matB,
                     standardize) {
        #require(penalized)
        ind <- stratum %in% matB[x, ]
        (penalized.clr(response[ind],
          stratum[ind],
          penalized = penalized[ind, ],
          unpenalized = unpenalized[ind, ],
          lambda = lambda,
          alpha = alpha,
          standardize = standardize
        )$penalized != 0) * 1
      },
      response, stratum, penalized,
      unpenalized, lambda, alpha, matB,
      standardize
    ))
    parallel::stopCluster(cl)
  } else {
    for (i in 1:(2 * B)) {
      ind <- stratum %in% matB[i, ]
      fit <- penalized.clr(response[ind],
        stratum[ind],
        penalized = penalized[ind, ],
        unpenalized = unpenalized[ind, ],
        lambda = lambda,
        alpha = alpha,
        standardize = standardize
      )
      selB[i, ] <- (fit$penalized != 0) * 1
      colnames(selB) <- names(fit$penalized)
    }
  }
  Pistab <- colMeans(selB)
  if (return.matB) {
    return(list("Pistab" = Pistab, "matB" = matB))
  } else {
    return(list("Pistab" = Pistab))
  }
}
subsample.clr.v <- Vectorize(subsample.clr, vectorize.args = "lambda")
