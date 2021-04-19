# penalizedclr



The R package <tt>penalizedclr</tt> provides an implementation of the penalized logistic regression model that can be used in the analysis of matched case-control studies. The implementation allows to apply different penalties to different blocks of covariates, and  is therefore particularly useful in the presence of multi-source heterogenous data, such as multiple layers of omics measurements. Both L1 and L2 penalties are implemented. 
Additionally, the package implements stability selection to allow for  variable selection in the considered regression model. 

## Installation

You can install the released version of <tt>penalizedclr</tt> from [CRAN](https://CRAN.R-project.org) with:


```{r eval = FALSE}
install.packages("penalizedclr")
```


And the development version from [GitHub](https://github.com/) with:

```{r eval = FALSE}
library(devtools)
install_github("veradjordjilovic/penalizedclr")
```

Load the package with: 

```{r setup}
library(penalizedclr)
```



## Examples
Some examples of how to fit a penalized conditional regression model with source-specific penalty parameters and how to perform variable selection with <tt>penalizedclr</tt>. 

Initial settings and libraries to be loaded:

```{r message=FALSE, warning=FALSE}
set.seed(123)
require(tidyverse)
```

### Data generation

We generate a simple multi-source data set, with two groups of covariates containing 80 and 20 variables, respectively. As specified above, each case is matched to a control, and the probability of being a case in each stratum (case-control pair) is obtained from the linear predictor. An intercept is generated independently for each stratum.

```{r}
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
```

### Penalized conditional logistic regression: model fit

The function <tt>penalized.clr</tt>  fits a penalized conditional logistic regression model with different penalties for different blocks of covariates. The L1 penalty parameter `lambda` (a numeric vector of the length equal to the number of blocks) can be specified by the user or computed internally. In the latter case, cross-validation is performed internally for each data layer separately. It is also possible to apply elastic net penalty through parameter `alpha`. 


#### Case 1: Penalty parameters provided by the user

This code illustrates how to fit <tt>penalized.clr</tt> with penalty parameters specified by the user. Here <tt>Y</tt> is the response vector, <tt>X</tt> is the multi-source matrix of covariates, <tt>stratum</tt> is the vector of ids of the matching pairs, and <tt>p</tt> is the vector of block dimensions. It has  the same dimension as the vector of penalty parameters <tt>lambda</tt>. Penalty parameters are thus specified for each data source, i.e., a block of covariates. It is possible to standardize the variables by setting <tt>standardize = TRUE</tt>  (<tt>FALSE</tt> by default).

```{r results = 'hide'}
fit1 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      lambda = c(6,7), p = p, standardize = TRUE)

```
The function returns a list with elements:

1. <tt>penalized</tt> - Regression coefficients for the penalized covariates.

2.  <tt>unpenalized</tt> - Regression coefficients for the unpenalized covariates.

3. <tt>converged</tt> - Whether the fitting process was judged to be converged.

4. <tt>lambda</tt> - The tuning parameter for L1 used.

5. <tt>alpha</tt> - The elastic net mixing parameter used.


For instance, `fit1$penalized` is a numeric vector containing the regression coefficients for the penalized covariates.  We can compare estimated coefficients with the true coefficients (the ones that generated our data), by constructing a `2 x 2` contingency table cross-tabulating true non-zero and estimated non-zero coefficients:

```{r}
nonzero_index <- (beta != 0) * 1 #index of nonzero coefficients
table(fitted = (fit1$penalized != 0) * 1, nonzero_index)
```



#### Case 2: Penalty parameters computed internally

When the user does not specify the penalty parameters, they are computed internally for each data source separately based on cross-validation. 

```{r  results = 'hide'}
fit2 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      p = p,
                      standardize = TRUE)

```

The selected penalty coefficients are: 

```{r}
fit2$lambda
```
The package uses fast cross-validation implemented in the R package <tt>clogitL1</tt> (see the function <tt>find.default.lambda</tt> described below).
We recommend to inspect manually the obtained `lambda` parameters. It is also recommended to investigate different ratios of data-source specific  penalties (see, for instance, Boulesteix et al. 2017 and examples of penalty grids therein). Note that cross validation depends on random data splits, and different runs will return different values of `lambda`. 





#### Case 3: Blocks not provided

This code shows what happens when no information is given about the blocks of predictors (`p`). In this case, all covariates are considered to form a single block and are penalized equally. 

```{r  results = 'hide'}
fit3 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      standardize = TRUE)
```

The selected penalty coefficient is: 

```{r}
fit3$lambda
```






#### Case 4: Unpenalized covariates

Somethimes,  a subset of covariates should be excluded from penalization. This can be achieved by specifying the `unpenalized` argument. In what follows, we penalize  the first block of covariates, and leave the remaining block unpenalized. 

```{r  results = 'hide'}
X1 <- X[, 1:p[1]]
X2 <- X[, (p[1]+1):(p[1]+p[2])]
fit4 <- penalized.clr(response = Y, penalized = X1, unpenalized = X2, 
                      stratum = stratum, p = p[1], 
                      standardize = TRUE)
```

This can be particularly useful when performing variable selection on omics variables (penalized) while adjusting for clinical covariates that should not be penalized. 

#### Case 5: Elastic net penalty

The package <tt>penalizedclr</tt>  is not limited to L1, i.e. lasso, penalty. While L1 penalty is more suited for variable selection, in the presence of highly correlated covariates, it can be useful to add small amount of  L2 penalty.  The two can be  combined in an elastic net framework by specifying  the mixing parameter `alpha` that can assume values between 0 and 1. Default `alpha=1` gives lasso. 

```{r  results = 'hide'}
fit2 <- penalized.clr(response = Y, penalized = X, stratum = stratum,
                      p = p,
                      standardize = TRUE, alpha = 0.6)
```


### Penalized conditional logistic regression: stability selection

The function <tt>stable.clr</tt> performs stability selection (Meinshausen and Bühlmann 2010 and  Shah and  Samworth 2013) for variable selection in penalized conditional logistic regression. For details on stability selection, we refer to the original publications. Briefly, a set of L1 penalties is considered, and for each considered value  of the penalty, $2B$ random subsamples of  $\lfloor K/2 \rfloor$ matched pairs are taken and a penalized model is estimated. For each covariate and penalty value, selection probability is estimated as the proportion of estimated models in which the associated coefficient estimate is different from zero. Finally, the overall selection probability of a variable is obtained by taking the maximum selection probability over all considered penalty values. 

Parameter $B$ is set to 100 by default, but can be changed by the user. Note that this choice will have an impact on the computation time, and higher values of $B$ will lead to a slower computation.  

The function returns a list containing the selection probabilities for each covariate, i.e. the proportion of estimated models in which the associated coefficient estimate is different from zero. 
The user can then set a threshold  for  selection probability (values ranging from 0.55 to 0.9 are recommended) and obtain the set of selected covariates. 

#### Case 1: Penalty parameters specified by the user

The following code performs stability selection when all covariates are considered as a single block (a single data source) and a sequence of L1 penalties contains only 2 values.

```{r}
stable1 <- stable.clr(response = Y,
                      penalized = X, stratum = stratum,
                      lambda.seq = c(10,20))
```
The function returns a list with two elements:

1. <tt>Pilambda</tt> a numeric vector giving estimates of selection probabilities for each penalized covariate.

2. <tt>lambda.seq</tt> a sequence of L1 penalties used.  

To inspect the results, we can, for instance, print the covariates with selection probability higher than 0.6:

```{r}
which(stable1$P>0.6)

```



#### Case 2: Penalty parameters computed externally

It is possible to obtain the sequence of L1 penalty parameters via the function <tt>find.default.lambda</tt>. It relies on the `cv.clogitL1` function of the `clogitL1` package to perform cross-validation to determine a suitable `lambda` sequence. Note that it considers each data source separately and  returns a sequence of four values per source (the optimal penalty along with three additional values, see package documentation for details). The number of folds is set to 10 by default but can also be specified by the user.

```{r}
lambda.seq <- find.default.lambda(response = Y, 
                                  stratum = stratum, 
                                  penalized = X, 
                                  alpha=1,
                                  nfolds = 10)

lambda.seq
```

```{r}
stable2 <- stable.clr(response = Y,
                      penalized = X, stratum = stratum,
                      lambda.seq = lambda.seq)
```

Covariates with selection probability higher than 0.6:

```{r}
which(stable2$P>0.6)

```

#### Case 3: Penalty parameters computed internally

If we do not specify  the penalty parameters as in case 2, <tt>stable.clr</tt> will compute them by default by using <tt>find.default.lambda</tt> with default options. We only need to run: 

```{r}
stable3 <- stable.clr(response = Y,
                      penalized = X, 
                      stratum = stratum)
```

The selected sequence of `lambda` is: 

```{r}
stable3$lambda.seq
```


Covariates with selection probability higher than 0.6:

```{r}
which(stable3$P>0.6)

```



#### Case 4: Covariates divided into blocks

This code implements stability selection while taking into account the block structure of covariates ($p$). To achieve this, we use the function <tt>stable.clr.g</tt> which extends <tt>stable.clr</tt> to allow for blocks (groups) of covariates. In this case, the function takes as an argument `lambda.list`,  a list of  sequences of L1 penalties, one for each block.  Note that the computation time depends on the length of the block-specific sequences: for each combination of penalties, $2B$ subsamples are taken and used to estimate a penalized conditional logistic regression model. 

```{r}

lambda.list <- find.default.lambda(response = Y,
                    penalized = X, stratum = stratum,
                    p = p)

lambda.list
```


```{r}
stable.g1 <- stable.clr.g(response = Y,
                          penalized = X,
                          stratum = stratum,
                          p = p,
                          lambda.list = lambda.list)
```


```{r}

which(stable.g1$P>0.6)
```
Note that if <tt>p</tt> is not specified, the package will automatically run <tt>stable.clr</tt> and apply equal penalties for all covariates. 

# References


1. Boulesteix, A. L., De Bin, R., Jiang, X., & Fuchs, M. (2017). IPF-LASSO: Integrative-penalized regression with penalty factors for prediction based on multi-omics data. Computational and mathematical methods in medicine, 2017.

2. Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473.

3.  Shah, R. D., & Samworth, R. J. (2013). Variable selection with error control:
another look at stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(1), 55-80.




