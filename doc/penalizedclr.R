## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,eval =TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("penalizedclr")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  install_github("veradjordjilovic/penalizedclr")

## ----setup--------------------------------------------------------------------
library(penalizedclr)

## ----message=FALSE, warning=FALSE---------------------------------------------
set.seed(123)
require(tidyverse)

## -----------------------------------------------------------------------------
# two groups of predictors
p <- c(80, 20)

# percentage and number of non-null variables
p_nz <- c(0.2, 0.8)
m_nz <- round(p*p_nz, 0)

# number of different strata (case-control pairs)
K <- 125

# number of cases and controls in each stratum (not necessarily 1:1 matching, other designs are also allowed)
n_cases <- 1
n_ctrl <- 1


# generating covariates
X = cbind(matrix(rnorm(p[1] * K * (n_cases + n_ctrl), 0, 1), ncol = p[1]),
          matrix(rnorm(p[2] * K * (n_cases + n_ctrl), 0, 2), ncol = p[2]))

# coefficients
beta <- as.matrix(c(rnorm(m_nz[1], 0, 0.8),
                    rep(0, p[1] - m_nz[1]),
                    rnorm(m_nz[2], 0.1, 0.4),
                    rep(0, p[2] - m_nz[2])), ncol = 1)



beta_stratum <- rep(rnorm(K, 0, 2), each= n_cases+n_ctrl)

# stratum membership
stratum <- rep(1:K, each= n_cases+n_ctrl)

# linear predictor
lin_pred <- beta_stratum + X %*% beta

prob_case <- exp(lin_pred) / (1 + exp(lin_pred))


# generate the response

Y <- rep(0, length(stratum))

data_sim <- as_tibble(data.frame(stratum = stratum,
                                 probability = prob_case,
                                 obs_id = 1 : length(stratum)))
data_sim_cases <- data_sim %>%
  group_by(stratum)%>%
  sample_n(n_cases, weight = probability)

Y[data_sim_cases$obs_id] <- 1

## ----results = 'hide'---------------------------------------------------------
fit1 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      lambda = c(6,7), p = p, standardize = TRUE)


## -----------------------------------------------------------------------------
nonzero_index <- (beta != 0) * 1 #index of nonzero coefficients
table(fitted = (fit1$penalized != 0) * 1, nonzero_index)

## ----results = 'hide'---------------------------------------------------------
fit2 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      p = p,
                      standardize = TRUE)


## -----------------------------------------------------------------------------
fit2$lambda

## ----results = 'hide'---------------------------------------------------------
fit3 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      standardize = TRUE)

## -----------------------------------------------------------------------------
fit3$lambda

## ----results = 'hide'---------------------------------------------------------
X1 <- X[, 1:p[1]]
X2 <- X[, (p[1]+1):(p[1]+p[2])]
fit4 <- penalized.clr(response = Y, penalized = X1, unpenalized = X2, 
                      stratum = stratum, p = p[1], 
                      standardize = TRUE)

## ----results = 'hide'---------------------------------------------------------
fit2 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      p = p,
                      standardize = TRUE, alpha = 0.6)

## -----------------------------------------------------------------------------
stable1 <- stable.clr(response = Y,
                      penalized = X, stratum = stratum,
                      lambda.seq = c(10,20))

## -----------------------------------------------------------------------------
which(stable1$P>0.6)


## -----------------------------------------------------------------------------
lambda.seq <- find.default.lambda(response = Y, 
                                  stratum = stratum, 
                                  penalized = X, 
                                  alpha=1,
                                  nfolds = 10)

lambda.seq

## -----------------------------------------------------------------------------
stable2 <- stable.clr(response = Y,
                      penalized = X, stratum = stratum,
                      lambda.seq = lambda.seq)

## -----------------------------------------------------------------------------
which(stable2$P>0.6)


## -----------------------------------------------------------------------------
stable3 <- stable.clr(response = Y,
                      penalized = X, 
                      stratum = stratum)

## -----------------------------------------------------------------------------
stable3$lambda.seq

## -----------------------------------------------------------------------------
which(stable3$P>0.6)


## -----------------------------------------------------------------------------

lambda.list <- find.default.lambda(response = Y,
                    penalized = X, stratum = stratum,
                    p = p)

lambda.list

## -----------------------------------------------------------------------------
stable.g1 <- stable.clr.g(response = Y,
                          penalized = X,
                          stratum = stratum,
                          p = p,
                          lambda.list = lambda.list)

## -----------------------------------------------------------------------------

which(stable.g1$P>0.6)

