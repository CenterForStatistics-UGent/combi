% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildConfMat.R
\name{buildConfMat}
\alias{buildConfMat}
\title{Build confounder design matrices with and without intercepts}
\usage{
buildConfMat(confounders)
}
\arguments{
\item{confounders}{A dataframe of confounding variables
#' For the preliminary trimming, we do not include an intercept,
but we do include all the levels of the factors using contrasts=FALSE:
 we want to do the trimming in every subgroup, so no hidden reference levels
  For the filtering we just use a model with an intercept and
   treatment coding, here the interest is only in adjusting the offset}
}
\value{
a list with components
\item{confModelMatTrim}{A confounder matrix without intercept, with all
 levels of factors present. This will be used to trim out taxa that have
  zero abundances in any subgroup defined by confounders}
\item{confModelMat}{A confounder matrix with intercept,
and with reference levels for factors absent.
This will be used to fit the model to modify the independence model,
and may include continuous variables}
}
\description{
Build confounder design matrices with and without intercepts
}
