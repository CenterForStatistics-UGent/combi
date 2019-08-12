#' Microbiomes and viromes of IBD patients and healthy controls
#'
#' Microbiome and virome sequencing data of IBD patients and healthy controls,
#' together with baseline covariates
#'
#' @format Two phyloseq objects containing microbiome and virome data on IBD patients and healthy controls
#' \describe{
#'   \item{microPruneVir}{The microbiome dataset pruned for matches with the virome object}
#'   \item{virPrune}{The virome dataset pruned for matches with the microbiome object}
#' }
#' @source \url{https://www.ibdmdb.org/}
"microPruneVir"
"virPrune"