% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss_funs.R
\name{get_res}
\alias{get_res}
\title{Fit a linear model on a subset of covariates and return the residuals.}
\usage{
get_res(X_tr, Y_tr, X_te, Y_te, S, intercept = TRUE)
}
\arguments{
\item{X_tr}{A matrix of dimension \verb{n_tr x d}.}

\item{Y_tr}{A vector of length \code{n_tr}.}

\item{X_te}{A matrix of dimension \verb{n_te x d}.}

\item{Y_te}{A vector of length \code{n_te}.}

\item{S}{A vector of either indices or names corresponding to the columns of
\code{X_tr} and \code{X_te}.}

\item{intercept}{Whether to include an intercept in the linear regressions.}
}
\value{
A list of two vectors of length \code{n_tr} and \code{n_te} respectively.
}
\description{
Fit a linear model on a subset of covariates and return the residuals.
}
\details{
Given two sets of data, the 1st set is used as the training set
and the 2nd is used as the test set. Both the training residual and the
test residuals are returned.
}
