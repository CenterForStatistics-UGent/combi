% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqM.R
\name{seqM}
\alias{seqM}
\title{A small auxiliary function for the indices of the lagrange multipliers}
\usage{
seqM(y, normal = TRUE, nLambda1s = 1)
}
\arguments{
\item{y}{an integer, the current dimension}

\item{normal}{a logical, is there a normalization restriction?}

\item{nLambda1s}{the number of centering restrictions}
}
\value{
a vector containing the ranks of the current lagrangian multipliers
}
\description{
A small auxiliary function for the indices of the lagrange multipliers
}
