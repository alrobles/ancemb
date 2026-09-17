#' Assert that a phylo tree has usable branch lengths
#'
#' Brownian Motion and Ornstein-Uhlenbeck models are defined in terms of
#' branch lengths; several `ape`/`phytools` calls fail later with less
#' actionable messages when `edge.length` is missing or malformed.
#'
#' @param tree A `phylo` object.
#' @return `tree`, invisibly.
#' @keywords internal
.assert_branch_lengths <- function(tree) {
  if (is.null(tree$edge.length)) {
    stop("'tree' has no branch lengths; supply a tree with 'edge.length' ",
         "(e.g. ape::compute.brlen()) before phylogenetic reconstruction.",
         call. = FALSE)
  }
  if (length(tree$edge.length) != nrow(tree$edge)) {
    stop("'tree$edge.length' does not match the number of edges.",
         call. = FALSE)
  }
  if (any(!is.finite(tree$edge.length))) {
    stop("'tree' has non-finite branch lengths.", call. = FALSE)
  }
  if (any(tree$edge.length < 0)) {
    stop("'tree' has negative branch lengths.", call. = FALSE)
  }
  invisible(tree)
}
