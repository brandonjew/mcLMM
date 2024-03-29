% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcLMM.R
\name{mc_remle}
\alias{mc_remle}
\title{Efficient REMLE for multi-context LMM}
\usage{
mc_remle(Y, X, force.iter = FALSE, verbose = TRUE)
}
\arguments{
\item{Y}{Matrix with individuals as rows and measurements as columns
Missing measurements must be NA}

\item{X}{Matrix with \code{n} rows and \code{c} columns for \code{n}
individuals and \code{c} covariates.}

\item{force.iter}{Boolean. If TRUE, force iterative method even when there
is no missing data. This is included for testing purposes
only. The optimal non-iterative method for no missing data
is exact and way faster.}

\item{verbose}{Boolean. Output logging info.}
}
\value{
List of REMLE parameters in slots \code{coef}, \code{ve},
        and \code{vg}.
}
\description{
Efficient REMLE for multi-context LMM
}
