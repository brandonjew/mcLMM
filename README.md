# mcLMM

An ultra-fast linear mixed model (LMM) framework for association studies across **m**ultiple **c**ontexts.

### What is a multiple-context association study?

This term refers to a study design where multiple measurements are taken from the same individuals. For example, measuring BMI at several time points or measuring gene expression in different tissues from the same individuals. If you wish to model these response variables (e.g. BMI, gene expression) in an LMM framework with features that remain constant throughout each measurement (such as genetics) and random effects that account for within-individual variation, mcLMM can efficiently optimize this model.

### How fast is ultra-fast?

The mcLMM optimization algorithm is specifically designed for LMMs with structured design matrices. We are able to estimate the optimal MLE and REMLE variance components in linear time with respect to the number of individuals. A naive LMM optimization algorithm is typically cubic in runtime. This means that we can process even the largest datasets in seconds rather than hours (or days sometimes). In the general case (possible missing response measurements), the algorithm is iterative, with each iteration cubic in time complexity with respect to the number of contexts. When there is no missing data, the algorithm is exact and linear with respect to the number of samples and contexts.

## Installation

mcLMM is an R package and can be installed as follows:

```r
devtools::install_github("brandonjew/mcLMM")
```

## Getting Started

```r
library(mcLMM)
sim.data <- mcLMM::simulate_data(ni=100, context.weights=c(1, 0.9, 0.6),
                                 var.e=0.2, var.g=0.4)
```

mcLMM can provide the optimal MLE or REMLE variance components for a model as follows:

```r
mle.mdl <- mcLMM::mc_mle(Y=sim.data$Y.mc, X=sim.data$X.mc)
```

```r
remle.mdl <- mcLMM::mc_remle(Y=sim.data$Y.mc, X=sim.data$X.mc)
```

where Y is the response variable encoded as a matrix with individuals as rows and contexts as columns, X is the covariate matrix with individuals as rows and covariates as columns. These functions return a list of the LMM parameters (variance components and estimated covariate coefficients and coefficient correlations)

mcLMM also implements the Meta-Tissue meta-analytic framework. It can be called as follows

```r
meta.mdl <- mcLMM:meta_tissue(expr=sim.data$Y.mc, geno=sim.data$X.mc[,"Genotype"])
```
