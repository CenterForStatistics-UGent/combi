% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estMeanVarTrend.R
\name{estMeanVarTrend}
\alias{estMeanVarTrend}
\title{Estimate a column-wise mean-variance trend}
\usage{
estMeanVarTrend(
  data,
  meanMat,
  baseAbundances,
  libSizes,
  plot = FALSE,
  meanVarFit,
  degree = 2L,
  constraint = "none",
  ...
)
}
\arguments{
\item{data}{the data matrix with n rows}

\item{meanMat}{the estimated mean matrix}

\item{baseAbundances}{The baseline abundances}

\item{libSizes}{Library sizes}

\item{plot}{A boolean, should the trend be plotted?}

\item{meanVarFit}{A character string describing the type of trend to be
fitted: either "spline" or "cubic"}

\item{degree}{The degree of the spline}

\item{constraint}{Constraint to the spline}

\item{...}{additional arguments passed on to the plot() function}
}
\value{
A list with components
\item{meanVarTrend}{An smoothed trend function, that can map a mean on a variance}
\item{meanVarTrendDeriv}{A derivative function of this}
}
\description{
Estimate a column-wise mean-variance trend
}
