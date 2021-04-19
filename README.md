# penalizedclr

## Introduction

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


