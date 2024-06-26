% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedy_algs.R
\name{prune_candi}
\alias{prune_candi}
\title{Prune candidates and remove the ones that are not causal change points.}
\usage{
prune_candi(X, Y, candi, PS, alpha, test_name, ...)
}
\value{
A vector of estimated causal change points
}
\description{
Prune candidates and remove the ones that are not causal change points.
}
