#' Calculate MLE negative log likelihood given detla
#' 
#' @param d Value of delta to evaluate
#' @param const A named list of constants used to calculate likelihood as
#'              described in documentation of \code{get_constants}. These
#'              variables are loaded into the environment with \code{list2env}.
#' @return Negative log likelihood of MLE model given delta (float)
mle_nll <- function(d, const){
  list2env(const, env=environment())
  d.vec <- g.sizes + d
  B[upper.tri(B, diag=TRUE)] <- E - colSums(D/d.vec)
  # We assume chol() only uses upper.tri (things will break if it doesn't)
  tryCatch(expr = {B <- chol2inv(chol(B))},
           error = function(err.msg){
             stop("Caught error with inverse of t(X)%*%solve(H)%*%X. ",
                  "Likely due to linearly dependent columns in X.\n",
                  "Original error message: ", err.msg)
           })
  YHX <- sxy - colSums(gxsy/d.vec)
  YHY <- sy2 - sum(syi2g/d.vec)
  YHXBXHY <- sum((YHX^2)*diag(B)) + 2*(sum(sapply(1:((ntc)-1), function(i){
    sum(sapply((i+1):(ntc), function(j){
      YHX[i]*YHX[j]*B[i,j]
    }))
  }))) 
  R <- YHY - YHXBXHY
  return(ns*log(R) - ni*log(d) + sum(g.ind.sizes*log(d.vec)))
}

#' Calculate REMLE negative log likelihood given detla
#' 
#' @param d Value of delta to evaluate
#' @param const A named list of constants used to calculate likelihood as
#'              described in documentation of \code{get_constants}. These
#'              variables are loaded into the environment with \code{list2env}.
#' @return Negative log likelihood of REMLE model given delta (float)
remle_nll <- function(d, const){
  list2env(const, env=environment())
  d.vec <- g.sizes + d
  B[upper.tri(B, diag=TRUE)] <- E - colSums(D/d.vec)
  # We assume chol() only uses upper.tri (things will break if it doesn't)
  tryCatch(expr = {B <- chol(B)},
           error = function(err.msg){
             stop("Caught error with inverse of t(X)%*%solve(H)%*%X. ",
                  "Likely due to linearly dependent columns in X.\n",
                  "Original error message: ", err.msg)
           })
  S <- prod(diag(B)^2) # Determinant of t(X)%*%solve(H)%*%X
  B <- chol2inv(B)
  YHX <- sxy - colSums(gxsy/d.vec)
  YHY <- sy2 - sum(syi2g/d.vec)
  YHXBXHY <- sum((YHX^2)*diag(B)) + 2*(sum(sapply(1:((ntc)-1), function(i){
    sum(sapply((i+1):(ntc), function(j){
      YHX[i]*YHX[j]*B[i,j]
    }))
  }))) 
  R = YHY - YHXBXHY
  return((ns-ntc)*log(R) - ni*log(d) + sum(g.ind.sizes*log(d.vec)) + log(S))
}

#' Calculate optimal variance component MLE estimates given delta
#' 
#' @param d Value of delta to evaluate
#' @param const A named list of constants used to calculate parameters as
#'              described in documentation of \code{get_constants}. These
#'              variables are loaded into the environment with \code{list2env}.
#' @return A list with ve in slot \code{ve}, vg in slot \code{vg}, 
#'         estimated effect sizes and uncorrelated standard errors
#'         in slot \code{coef}.
mle_get_params <- function(d, const){
  list2env(const, env=environment())
  d.vec <- g.sizes + d
  B[upper.tri(B, diag=TRUE)] <- E - colSums(D/d.vec)
  # We assume chol() only uses upper.tri (things will break if it doesn't)
  B <- chol2inv(chol(B))
  YHX <- sxy - colSums(gxsy/d.vec)
  YHY <- sy2 - sum(syi2g/d.vec)
  YHXBXHY <- sum((YHX^2)*diag(B)) + 2*(sum(sapply(1:((ntc)-1), function(i){
    sum(sapply((i+1):(ntc), function(j){
      YHX[i]*YHX[j]*B[i,j]
    }))
  }))) 
  R <- YHY - YHXBXHY
  vg <- R/ns/d
  ve <- d * vg
  b.stderr <- sqrt(diag(B)*ve)
  b <- as.vector(B %*% YHX)
  beta <- cbind(b, b.stderr)
  rownames(B) <- b.names
  colnames(B) <- b.names
  rownames(beta) <- b.names
  colnames(beta) <- c("Estimate", "StandardError")
  return(list(ve=ve, vg=vg, coef=beta, b.cov=B))
}

#' Calculate optimal variance component REMLE estimates given delta
#' 
#' @param d Value of delta to evaluate
#' @param const A named list of constants used to calculate parameters as
#'              described in documentation of \code{get_constants}. These
#'              variables are loaded into the environment with \code{list2env}.
#' @return A list with ve in slot \code{ve}, vg in slot \code{vg}, 
#'         estimated effect sizes and uncorrelated standard errors
#'         in slot \code{coef}.
remle_get_params <- function(d, const){
  list2env(const, env=environment())
  d.vec <- g.sizes + d
  B[upper.tri(B, diag=TRUE)] <- E - colSums(D/d.vec)
  # We assume chol() only uses upper.tri (things will break if it doesn't)
  B <- chol2inv(chol(B))
  YHX <- sxy - colSums(gxsy/d.vec)
  YHY <- sy2 - sum(syi2g/d.vec)
  YHXBXHY <- sum((YHX^2)*diag(B)) + 2*(sum(sapply(1:((ntc)-1), function(i){
    sum(sapply((i+1):(ntc), function(j){
      YHX[i]*YHX[j]*B[i,j]
    }))
  }))) 
  R <- YHY - YHXBXHY
  vg <- R/(ns-ntc)/d
  ve <- d * vg
  b.stderr <- sqrt(diag(B)*ve)
  b <- as.vector(B %*% YHX)
  beta <- cbind(b, b.stderr)
  rownames(B) <- b.names
  colnames(B) <- b.names
  rownames(beta) <- b.names
  colnames(beta) <- c("Estimate", "StandardError")
  return(list(ve=ve, vg=vg, coef=beta, b.cor=cov2cor(B)))
}

#' Calculate optimal variance component REMLE estimates given delta using Lin
#' and Sullivan approach or standard for Meta-Tissue
#' 
#' @param d Value of delta to evaluate
#' @param const A named list of constants used to calculate parameters as
#'              described in documentation of \code{get_constants}. These
#'              variables are loaded into the environment with \code{list2env}.
#' @return A list with ve in slot \code{ve}, vg in slot \code{vg}, 
#'         estimated effect sizes and uncorrelated standard errors
#'         in slot \code{coef}.
remle_get_params_meta <- function(d, const, newRE=TRUE){
  list2env(const, env=environment())
  d.vec <- g.sizes + d
  B[upper.tri(B, diag=TRUE)] <- E - colSums(D/d.vec)
  # We assume chol() only uses upper.tri (things will break if it doesn't)
  B <- chol2inv(chol(B))
  YHX <- sxy - colSums(gxsy/d.vec)
  YHY <- sy2 - sum(syi2g/d.vec)
  YHXBXHY <- sum((YHX^2)*diag(B)) + 2*(sum(sapply(1:((ntc)-1), function(i){
    sum(sapply((i+1):(ntc), function(j){
      YHX[i]*YHX[j]*B[i,j]
    }))
  }))) 
  R <- YHY - YHXBXHY
  vg <- R/(ns-ntc)/d
  ve <- d * vg
  if (newRE){
    b.stderr <- sqrt(diag(B)*ve)[(ntc-nt+1):ntc]
  }
  else {
    C <- chol2inv(chol((B*ve)[(ntc-nt+1):ntc,(ntc-nt+1):ntc]))
    b.stderr <- 1/sqrt(rowSums(C))
  }
  b <- as.vector(B %*% YHX)
  beta <- cbind(b[(ntc-nt+1):ntc], b.stderr)
  rownames(B) <- b.names
  colnames(B) <- b.names
  rownames(beta) <- b.names[(ntc-nt+1):ntc]
  colnames(beta) <- c("Estimate", "StandardError")
  return(list(ve=ve, vg=vg, coef=beta, 
              b.cor=cov2cor(B)[(ntc-nt+1):ntc,(ntc-nt+1):ntc]))
}

#' Calculate constants independent of delta in MLE and REMLE likelihoods
#' 
#' @param Y Matrix with individuals as rows and measurements as columns
#'          Missing measurements must be NA.
#' @param X Matrix with \code{n} rows and \code{c} columns for \code{n}
#'          individuals and \code{c} covariates. 
#' @return A list of named constants. These constants are difficult to
#'         describe briefly. See Jew et al, 2020 for a more in-depth 
#'         discussion of these constants. \code{ns} is the total number of 
#'         measured responses. \code{nt} is the number of contexts.
#'         \code{ni} is the number of individuals. \code{nc} is the number
#'         of covariates. \code{ntc} is the number of measurement contexts 
#'         multiplied by the number of covariates. \code{g.sizes} is a vector
#'         of the unique number of measurements across all individuals.
#'         \code{g.ind.sizes} is a vector of the number of individuals
#'         that have a certain number of measurements. \code{B} is a dummy
#'         matrix that will store the inverse of XHX (beta hat covariance).
#'         \code{D} and \code{E} are vector representations of the upper 
#'         triangle of symmetric matrices that are used to calculate \code{B}.
#'         \code{sxy} is a vector of the sum of responses multiplied by each
#'         covariate within each context. \code{gxsy} is a matrix where each
#'         row represents a group and column represents the sum of the 
#'         covariate in each context multiplied by the sum of responses for
#'         each individual. \code{sy2} is the sum of squared responses.
#'         \code{syi2g} provides the sum of squared responses within each
#'         group of measurement numbers. 
get_constants <- function(Y,X){
  nt <- ncol(Y)
  ni <- nrow(Y)
  nc <- ncol(X)
  if (ni != nrow(X)){
    stop("X and Y should have the same number of rows (1 per individual)")
  }
  ntc <- nt*nc
  covariate.names <- colnames(X)
  context.names <- colnames(Y)
  if (is.null(context.names)){
    context.names <- paste0("T", 1:nt)
  }
  if (is.null(covariate.names)){
    covariate.names <- paste0("C", 1:nc)
  }
  covariate.names <- paste0('[',covariate.names,']')
  context.names <- paste0('[',context.names,']')
  b.names <- unlist(lapply(covariate.names, function(covariate.name){
    paste(covariate.name, context.names, sep=',')
  }))
  # given an index 1:nt, provides indices of individuals with this measurement
  t.indices <- lapply(1:nt, function(ti){
    ti.indices <- which(!is.na(Y[,ti]))
    if (length(ti.indices) == 0){
      stop("Detected measurement with no values ",
           sprintf("Column %i is all NA and should be removed", ti))
    }
    return(ti.indices)
  })
  ns <- sum(sapply(1:nt, function(ti){length(t.indices[[ti]])}))
  # sum of all responses squared
  sy2 <- 0
  # given index of individual provides sum of their responses
  syi <- rep(0,ni)
  # given index of group, provides sum of syi^2 for each individual in group
  syi2g <- rep(0,nt)
  # given an index 1:ng (defined after), provides indices of individuals with
  # this number of measurements
  g.indices <- lapply(1:nt, function(ti){c()})
  for (ind in 1:ni){
    g <- sum(!is.na(Y[ind,]))
    if (g == 0){
      stop("Detected individual with no measurements ",
           sprintf("(Row %i is all NA and should be removed).", ind))
    }
    sy2 <- sy2 + sum(Y[ind,]^2, na.rm=TRUE)
    syi[ind] <- sum(Y[ind,], na.rm=TRUE)
    syi2g[g] <- syi2g[g] + syi[ind]^2
    g.indices[[g]] <- c(g.indices[[g]], ind)
  }
  # vector of unique number of measurements in experiment
  g.sizes <- which(sapply(1:nt, function(g){length(g.indices[[g]]) != 0}))
  # vector of number of individuals in each group
  g.ind.sizes <- sapply(g.sizes, function(g){length(g.indices[[g]])})
  syi2g <- syi2g[g.sizes]
  ng <- length(g.sizes)
  # given row g in 1:ng and column i in 1:(nt*nc) , provides sum of covi*syi for 
  # each individual in group g
  gxsy <- do.call(rbind, lapply(g.sizes, function(g){
    sapply(1:(ntc), function(i){
      ti <- ((i-1) %%  nt) + 1
      ci <- ((i-1) %/% nt) + 1
      sub.ind.indices <- intersect(g.indices[[g]], t.indices[[ti]])
      sum(X[sub.ind.indices,ci]*syi[sub.ind.indices])
    })
  }))
  sxy <- sapply(1:(ntc), function(i){
    ti <- ((i-1) %%  nt) + 1
    ci <- ((i-1) %/% nt) + 1
    sum(X[t.indices[[ti]], ci]*Y[t.indices[[ti]],ti])
  })
  # we are storing upper.tri(M,diag=T) for E and D matrices
  # these are the matrices that form t(X)%*%solve(H)%*%X
  E <- unlist(lapply(1:(ntc), function(j){
    tj <- ((j-1) %%  nt) + 1
    cj <- ((j-1) %/% nt) + 1
    lapply(1:j, function(i){
      ti <- ((i-1) %% nt) + 1
      if (ti != tj){
        return(0)
      }
      ci <- ((i-1) %/% nt) + 1
      return(sum(sapply(t.indices[[ti]], function(ind){
        X[ind,ci]*X[ind,cj]
      })))
    })
  }))
  D <- do.call(rbind, lapply(g.sizes, function(g){
    unlist(lapply(1:(ntc), function(j){
      tj <- ((j-1) %%  nt) + 1
      cj <- ((j-1) %/% nt) + 1
      ind.indices <- intersect(t.indices[[tj]], g.indices[[g]])
      lapply(1:j, function(i){
        ti <- ((i-1) %% nt) + 1
        ci <- ((i-1) %/% nt) + 1
        sub.ind.indices <- intersect(ind.indices, t.indices[[ti]])
        if (length(sub.ind.indices) == 0){
          return(0)
        }
        return(sum(sapply(sub.ind.indices, function(ind){
          X[ind,ci]*X[ind,cj]
        })))
      })
    }))
  }))
  B <- matrix(nrow=ntc, ncol=ntc)
  return(list(ns=ns, ni=ni, nt=nt, ntc=ntc, nc=nc, g.sizes=g.sizes, 
              g.ind.sizes=g.ind.sizes, B=B, D=D, E=E, sxy=sxy, 
              gxsy=gxsy, sy2=sy2, syi2g=syi2g,
              b.names=b.names))
}

#' Ultra-fast MLE for multi-context LMM
#'
#' @param Y Matrix with individuals as rows and measurements as columns
#'          Missing measurements must be NA
#' @param X Matrix with \code{n} rows and \code{c} columns for \code{n}
#'          individuals and \code{c} covariates. 
#' @return List of MLE parameters in slots \code{delta}, \code{ve},
#'         and \code{vg}.
#' @export
mc_mle <- function(Y, X){
  const <- get_constants(Y,X)
  d <- stats::optimize(f=mle_nll, interval=c(exp(-10), exp(10)), const)$minimum
  est.params <- mle_get_params(d, const)
  return(est.params)
}

#' Ultra-fast REMLE for multi-context LMM
#'
#' @param Y Matrix with individuals as rows and measurements as columns
#'          Missing measurements must be NA
#' @param X Matrix with \code{n} rows and \code{c} columns for \code{n}
#'          individuals and \code{c} covariates. 
#' @return List of REMLE parameters in slots \code{delta}, \code{ve},
#'         and \code{vg}.
#' @export
mc_remle <- function(Y, X){
  const <- get_constants(Y,X)
  d <- stats::optimize(f=remle_nll, interval=c(exp(-10), exp(10)), const)$minimum
  est.params <- remle_get_params(d, const)
  return(est.params)
}

#' Ultra-fast meta-tissue algorithm
#' 
#' Provides identical results as meta-tissue in linear time with respect to
#' the number of samples rather than cubic. Slight differences may be 
#' due to different optimization algorithms of the same likelihood function.
#' 
#' @param expr Matrix with individuals as rows and tissues as columns. Missing
#'          gene expression must be NA
#' @param geno Vector of genotypes for each individual
#' @param covs Matrix with individuals as rows and covariates as columns.
#' @param heuristic Boolean. Uses heuristic for scaling standard errors. 
#'                  increases sensitivity at the cost of higher FPR.
#' @param newRE Boolean. Use new random effects model to perform meta
#'              analyses as discussed in Sul et al.
#' @return List of estimated coefficients \code{beta}, coefficient correlation
#'         \code{corr}, and \code{sigma_g}. 
#' @export
meta_tissue <- function(expr, geno, covs=NULL, heuristic=FALSE, newRE=TRUE){
  Y <- scale(expr)
  X <- rep(1,nrow(Y))
  if (!is.null(covs)){
    X <- cbind(X, covs)
  }
  X <- cbind(X, geno)
  const <- get_constants(Y,X)
  ntc <- const$ntc
  nt <- const$nt
  nc <- const$nc
  d <- stats::optimize(f=remle_nll, interval=c(exp(-10), exp(10)), const)$minimum
  est.params <- remle_get_params_meta(d, const, newRE)
  beta <- est.params$coef
  corr <- est.params$b.cor
  if (heuristic){
    # Run separate OLS
    separate.std <- sapply(1:ncol(Y), function(i){
      separate.Y <- Y[,i]
      indices <- !is.na(separate.Y)
      separate.Y <- separate.Y[indices]
      separate.X <- X[indices,,drop=F]
      sqrt(diag(vcov(lm(separate.Y~separate.X +0))))[nc]
    })
    # Run combined OLS (This can probably be sped up with some same tricks used for EMMA)
    combined.Y <- as.vector(t(Y))
    combined.Y <- combined.Y[!is.na(combined.Y)]
    combined.base.X <- matrix(0,nrow=length(combined.Y),ncol=ncol(Y))
    cur.row <- 1
    cur.col <- 1
    for (i in 1:nrow(X)){
      tissues <- which(!is.na(Y[i,]))
      for (j in tissues){
        combined.base.X[cur.row,j] <- 1
        cur.row <- cur.row + 1
      }
    }
    combined.X <- combined.base.X
    for (i in 2:ncol(X)){
      combined.X <- cbind(combined.X, combined.base.X*do.call(c,lapply(1:nrow(X), function(j) rep((X[j,i]),sum(!is.na(Y[j,]))))))
    }
    combined.std <- sqrt(diag(vcov(lm(combined.Y~combined.X +0))))[(ntc-nt+1):ntc]
    beta[,"StandardError"] <- beta[,"StandardError"]*separate.std/combined.std
  }
  zscores <- beta[,"Estimate"]/beta[,"StandardError"]
  metazscores <- sapply(1:ncol(Y), function(i){
    zscore <- zscores[i]
    pval <- pt(abs(zscore),
               df=(sum(!is.na(Y[,i])) - ncol(X)),
               lower.tail = FALSE)
    metazscore <- (zscore/abs(zscore))*abs(qnorm(pval))
    return(metazscore)
  })
  beta[,"StandardError"] <- beta[,"Estimate"]/metazscores
  return(list(beta=beta,corr=corr, sigmag=est.params$vg/(est.params$vg+est.params$ve)))
}