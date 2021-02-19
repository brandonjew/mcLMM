#' file MASS/R/mvrnorm.R
#' copyright (C) 1994-2015 W. N. Venables and B. D. Ripley
#'
#'  This program is free software; you can redistribute it and/or modify
#'  it under the terms of the GNU General Public License as published by
#'  the Free Software Foundation; either version 2 or 3 of the License
#'  (at your option).
#'
#'  This program is distributed in the hope that it will be useful,
#'  but WITHOUT ANY WARRANTY; without even the implied warranty of
#'  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'  GNU General Public License for more details.
#'
#'  A copy of the GNU General Public License is available at
#'  http://www.r-project.org/Licenses/
#'
#'  Slightly edited code to remove some input checks.
#'  
#'  @param n Number of samples
#'  @param mu Mean vector
#'  @param Sigma Covariance matrix
#'  
#'  @return MVN sample in \code{n} by \code{length(mu)} matrix
mcmvrnorm <-
  function(n, mu, Sigma)
  {
    tol <- 1e-6
    p <- length(mu)
    if(!all(dim(Sigma) == c(p,p))) stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if(!all(ev >= -tol*abs(ev[1L]))) stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    if(n == 1) drop(X) else t(X)
  }

#' Faster mvrnorm for Sigma=1*vg + I*ve
#'
#' Avoid eigen() since eigenvalues and eigenvectors are known
#' However, eigenvectors still need to be orthonormalized.
#' I'm assuming gramSchmidt() is faster than eigen()... hopefully.
#'
#' @return Transpose of what the normal mvrnorm returns.
fast_mcmvrnorm <- 
  function(n, t, ve, vg)
  {
    eigen.vectors <- diag(1,t)
    eigen.vectors[1,] <- -1
    eigen.vectors[,1] <- 1
    eigen.vectors <- pracma::gramSchmidt(eigen.vectors)$Q
    eigen.values <- c(t*vg + ve, rep(ve, t-1))
    X <- matrix(rnorm(t * n), n)
    X <- eigen.vectors %*% diag(sqrt(eigen.values)) %*% t(X)
    if(n == 1) drop(X) else X
  }

#' Simulate test data for mcLMM
#' 
#' Simulates response Y under the mcLMM model with the given parameters.
#' Univariate LMM with 012 genotype value as covariate (sampled from 
#' binomial with 0.5 MAF)
#' 
#' @param ni Number of individuals
#' @param context.weights Sampling rates for each context. Will append
#'                        a 1 to these rates to ensure each sample has
#'                        at least 1 context. Total number of contexts
#'                        will be \code{length(context.weights) + 1}.
#' @param var.e Residual variance
#' @param var.g Within-individual residual variance + covariance
#' @param random.seed Seed for reproducibility
#' 
#' @return A list with slots \code{X}, \code{K}, and \code{Y} serving as
#'         input for EMMA R package (design matrix, kinship matrix, and 
#'         response vector, respectively) and slots \code{X.mc} and \code{Y.mc}
#'         as input for the mcLMM functions. See \code{mc_mle} or 
#'         \code{mc_remle} documentation for more details on mcLMM input.
#' @export
simulate_data <- function(ni, context.weights=NULL, var.e, var.g, 
                          random.seed=NULL){
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
  n.individuals <- ni
  n.measurements <- length(context.weights) + 1
  context.weights <- c(1, context.weights)
  measurement.identities <- do.call(cbind, lapply(context.weights, function(p){
    np <- round(ni*p)
    x <- rep(FALSE, ni)
    x.indices <- sample(1:ni, np)
    x[x.indices] <- TRUE
    return(x)
    #rbinom(ni, 1, p) == 1
  }))
  measurement.identities <- lapply(1:n.individuals, function(i){
    which(measurement.identities[i,,drop=T])
  })
  measurement.sizes <- sapply(measurement.identities, function(x){
    length(x)
  })
  measurement.ordering <- order(measurement.sizes)
  measurement.identities <- measurement.identities[measurement.ordering]
  measurement.sizes <- measurement.sizes[measurement.ordering]
  group.sizes <- table(measurement.sizes)
  group.sizes <- group.sizes[as.character(sort(unique(measurement.sizes)))]
  n.samples <- sum(measurement.sizes)
  covs <- rbinom(n.individuals, 2, 0.5)
  x <- cbind(rep(1,n.individuals), covs)
  colnames(x) <- c("Intercept", "Genotype")
  K = Matrix::bdiag(lapply(measurement.sizes, function(ti) {
    matrix(1,nrow=ti,ncol=ti)
  }))
  X <- matrix(0,nrow=n.samples,ncol=n.measurements)
  cur.row <- 1
  cur.col <- 1
  for (i in 1:n.individuals){
    for (j in measurement.identities[[i]]){
      X[cur.row,j] <- 1
      cur.row <- cur.row + 1
    }
  }
  X <- cbind(X, (X*do.call(c, lapply(1:n.individuals, function(i) {
    rep((covs[i]), measurement.sizes[i])
  }))))
  beta <- rnorm(ncol(X))
  Y <- X%*%beta
  residuals <- do.call(c, lapply(names(group.sizes), function(group.size){
    ti <- as.numeric(group.size)
    sub.sigma <- (matrix(1,nrow=ti,ncol=ti) * var.g) + diag(var.e,ti)
    as.vector(t(mcmvrnorm(n=group.sizes[group.size], mu=rep(0,ti),
                        Sigma=sub.sigma)))
  }))
  Y <- Y + residuals
  tis <- unlist(measurement.identities)
  ind <- unlist(lapply(1:n.individuals, function(x){
    rep(x,length(measurement.identities[[x]]))
  })) 
  Y.mc <- matrix(NA, nrow=n.individuals, ncol=n.measurements)
  for (i in 1:n.samples){
    Y.mc[ind[i],tis[i]] <- Y[i]
  }
  colnames(Y.mc) <- paste0("Tissue_",1:n.measurements)
  rownames(Y.mc) <- paste0("Sample_",1:n.individuals)
  return(list(Y=Y, Y.mc=Y.mc, X=X, X.mc=x, K=K,
              ve=var.e, vg=var.g, beta=beta))
}



#' Simulate null test data for mcLMM
#' 
#' Simulates response Y under the mcLMM model with the given parameters.
#' Samples from MVN with Sigma = K*v_g + I*v_e
#' 
#' @param ni Number of individuals
#' @param nc Number of contexts/tissues
#' @param var.e Residual variance
#' @param var.g Within-individual residual variance + covariance
#' @param random.seed Seed for reproducibility
#' 
#' @return A list with slots \code{K}, and \code{Y} serving as
#'         input for EMMA R package (kinship matrix and 
#'         response vector, respectively) and slot \code{Y.mc}
#'         as input for the mcLMM functions. See \code{mc_mle} or 
#'         \code{mc_remle} documentation for more details on mcLMM input.
#' @export
simulate_null_data <- function(ni, nc, var.e, var.g, 
                               random.seed=NULL){
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
  n.individuals <- ni
  n.measurements <- nc
  sub.sigma <- matrix(var.g,nrow=n.measurements,ncol=n.measurements) + diag(var.e,n.measurements)
  #Y <- as.vector(mcmvrnorm(n=n.individuals, mu=rep(0,n.measurements), Sigma=sub.sigma))
  Y <- as.vector(fast_mcmvrnorm(n.individuals, n.measurements, var.e, var.g))
  tis <- unlist(lapply(1:n.individuals, function(x){
    1:n.measurements
  }))
  ind <- unlist(lapply(1:n.individuals, function(x){
    rep(x,n.measurements)
  })) 
  Y.mc <- matrix(NA, nrow=n.individuals, ncol=n.measurements)
  n.samples <- n.individuals*n.measurements
  for (i in 1:n.samples){
    Y.mc[ind[i],tis[i]] <- Y[i]
  }
  colnames(Y.mc) <- paste0("Tissue_",1:n.measurements)
  rownames(Y.mc) <- paste0("Sample_",1:n.individuals)
  K = Matrix::bdiag(lapply(rep(n.measurements, n.individuals), function(ti) {
    matrix(1,nrow=ti,ncol=ti)
  }))
  return(list(Y=Y, Y.mc=Y.mc, K=K,
              ve=var.e, vg=var.g))
}
