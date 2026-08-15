#' ancemb: Bayesian Ancestral Reconstruction of Protein Embeddings
#'
#' Provides functions to reconstruct ancestral protein embeddings from
#' extant (tip) embeddings and a phylogenetic tree. Models are fit with
#' pre-compiled Stan models using Brownian Motion or Ornstein-Uhlenbeck
#' evolutionary processes.
#'
#' @docType package
#' @name ancemb-package
#' @aliases ancemb
#' @keywords internal
"_PACKAGE"

#' @importFrom methods is
#' @importFrom stats rnorm sd median quantile
#' @importFrom ape vcv node.depth.edgelength dist.nodes reorder.phylo
NULL
