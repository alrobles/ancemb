#' Prepare embedding data for Stan
#'
#' Converts observed tip embeddings and a phylogenetic tree into the named
#' list required by the Brownian Motion and Ornstein-Uhlenbeck Stan models.
#' The tree is stored as an attribute on the returned list so that
#' downstream functions (e.g. `extract_ancestral()`) can reconstruct internal
#' node states without needing the tree to be passed again.
#'
#' @param tree A `phylo` object. Tip labels must be present.
#' @param z_obs Numeric matrix (`S x D`) of observed tip embeddings. Rownames
#'   must match `tree$tip.label` exactly.
#' @param obs_sigma Non-negative observation noise sd (default `0.1`).
#' @param rate_scale Half-normal scale for BM rates / OU diffusion variance
#'   (default `1`).
#' @param alpha_scale Half-normal scale for OU selection strength (default `1`).
#'
#' @return A named list suitable for `fit_embedding_bm()` and
#'   `fit_embedding_ou()`. The `tree` argument is attached as a `tree`
#'   attribute.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   set.seed(1)
#'   tree <- ape::rtree(8)
#'   tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
#'   z_obs <- matrix(stats::rnorm(8 * 4), nrow = 8)
#'   rownames(z_obs) <- tree$tip.label
#'   stan_data <- prepare_embedding_data(tree, z_obs)
#'   str(stan_data)
#' }
#' }
prepare_embedding_data <- function(tree, z_obs,
                                   obs_sigma = 0.1,
                                   rate_scale = 1,
                                   alpha_scale = 1) {
  checkmate::assert_class(tree, "phylo")
  checkmate::assert_matrix(z_obs, mode = "numeric", any.missing = FALSE,
                           min.rows = 2L, min.cols = 1L)
  checkmate::assert_numeric(obs_sigma, lower = 0, finite = TRUE, len = 1L)
  checkmate::assert_numeric(rate_scale, lower = 0, finite = TRUE, len = 1L)
  checkmate::assert_numeric(alpha_scale, lower = 0, finite = TRUE, len = 1L)

  if (is.null(tree$tip.label)) {
    stop("'tree' must have tip labels.")
  }
  if (is.null(rownames(z_obs))) {
    stop("'z_obs' must have rownames matching 'tree$tip.label'.")
  }
  if (!setequal(rownames(z_obs), tree$tip.label)) {
    stop("Rownames of 'z_obs' do not match the tip labels of 'tree'.")
  }

  S <- nrow(z_obs)
  D <- ncol(z_obs)

  z_obs <- z_obs[tree$tip.label, , drop = FALSE]

  C <- ape::vcv(tree)
  C <- C[tree$tip.label, tree$tip.label, drop = FALSE]

  tip_heights <- ape::node.depth.edgelength(tree)[seq_len(S)]

  z_anc_prior_mean <- colMeans(z_obs)
  z_anc_prior_sd   <- apply(z_obs, 2L, stats::sd)
  z_anc_prior_sd   <- pmax(z_anc_prior_sd, 0.1)

  stan_data <- list(
    S                = S,
    D                = D,
    z_obs            = z_obs,
    C                = C,
    tip_heights      = as.numeric(tip_heights),
    obs_sigma        = as.numeric(obs_sigma),
    z_anc_prior_mean = as.numeric(z_anc_prior_mean),
    z_anc_prior_sd   = as.numeric(z_anc_prior_sd),
    rate_scale       = as.numeric(rate_scale),
    alpha_scale      = as.numeric(alpha_scale)
  )

  attr(stan_data, "tree") <- tree
  stan_data
}
