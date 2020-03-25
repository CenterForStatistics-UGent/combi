#' Microbiomes of mice that underwent Pulsed Antibiotic Treatment
#'  (PAT) and controls
#'
#' Microbiome of mice that underwent Pulsed Antibiotic Treatment
#'  (PAT) and controls. The data were extracted from the source
#'  \url{https://www.ibdmdb.org/}, and then only the samples matching between
#'  microbiome and metabolome were retained.
#'
#' @format A phyloseq object containing microbiome data
#' \describe{
#'   \item{zhangMicrobio}{The microbiome dataset pruned for matches with the metabolome object}
#'  }
#' @usage data(Zhang)
#' @source \url{https://www.ibdmdb.org/}
"zhangMicrobio"
#' Metabolomes of mice that underwent Pulsed Antibiotic Treatment
#'  (PAT) and controls
#'
#' Metabolome of mice that underwent Pulsed Antibiotic Treatment
#'  (PAT) and controls
#'
#' @format SummarizedExperiment with metabolome data
#' \describe{
#'   \item{zhangMetabo}{The metabolome data as a SummarizedExperiment object}
#' }
#' @usage data(Zhang)
#' @source \url{https://www.ibdmdb.org/}
"zhangMetabo"
#' Baseline sample variables of PAT and control mice
#'
#' Baseline covariates of PAT mice and healthy controls
#'
#' @format A dataframe with baseline sample variables
#' \describe{
#'   \item{zhangMetavars}{The metadata on the mice}
#' }
#' @usage data(Zhang)
#' @source \url{https://www.ibdmdb.org/}
"zhangMetavars"
