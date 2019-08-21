#' Microbiomes of IBD patients and healthy controls
#'
#' Microbiome sequencing data of IBD patients and healthy controls,
#' together with baseline covariates
#'
#' @format A phyloseq object containing microbiome data on IBD patients and healthy controls
#' \describe{
#'   \item{microPruneVir}{The microbiome dataset pruned for matches with the virome object}
#'  }
#' @source \url{https://www.ibdmdb.org/}
"microPruneVir"
#' Viromes of IBD patients and healthy controls
#'
#' Virome sequencing data of IBD patients and healthy controls
#'
#' @format A phyloseq object containing virome data on IBD patients and healthy controls
#' \describe{
#'   \item{virPrune}{The virome dataset pruned for matches with the microbiome object}
#' }
#' @source \url{https://www.ibdmdb.org/}
"virPrune"
#' Baseline sample variables of IBD patients
#'
#' Baseline covariates of IBD patients and healthy controls,
#'
#' @format A dataframe with sample variables
#' \describe{
#'   \item{hmp2samVar}{Sample variables matching the microbiome and virome datasets}
#' }
#' @source \url{https://www.ibdmdb.org/}
"hmp2samVar"