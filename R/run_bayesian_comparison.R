#' Compare Brownian Motion and Ornstein-Uhlenbeck embedding fits
#'
#' Fits both the BM and OU Stan models to the same prepared data and returns
#' the fits, ancestral summaries, and a simple comparison of the inferred root
#' embeddings.
#'
#' @param data Stan data list produced by `prepare_embedding_data()`.
#' @param tree A `phylo` object. Required only if `data` does not carry a
#'   `tree` attribute.
#' @param chains Number of MCMC chains per model (default `4`).
#' @param iter_warmup Number of warmup iterations (default `1000`).
#' @param iter_sampling Number of post-warmup iterations (default `1000`).
#' @param ... Additional arguments passed to `fit_embedding_bm()` and
#'   `fit_embedding_ou()`.
#'
#' @return A list with components `bm_fit`, `ou_fit`, `bm_summary`,
#'   `ou_summary`, and `root_comparison`.
#' @export
#'
#' @examples
#' \donttest{
#' if (instantiate::stan_cmdstan_exists()) {
#'   set.seed(1)
#'   tree <- ape::rtree(8)
#'   tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
#'   z_obs <- matrix(stats::rnorm(8 * 4), nrow = 8)
#'   rownames(z_obs) <- tree$tip.label
#'   data <- prepare_embedding_data(tree, z_obs)
#'   cmp <- run_bayesian_comparison(data, chains = 1, iter_warmup = 200,
#'                                  iter_sampling = 200)
#'   str(cmp$root_comparison)
#' }
#' }
run_bayesian_comparison <- function(data, tree = NULL,
                                    chains = 4, iter_warmup = 1000,
                                    iter_sampling = 1000, ...) {
  checkmate::assert_list(data, names = "named")
  if (is.null(tree)) {
    tree <- attr(data, "tree")
  }
  if (is.null(tree)) {
    stop("'tree' is required (or set as a 'tree' attribute on data).")
  }
  .assert_branch_lengths(tree)
  checkmate::assert_int(chains, lower = 1L)
  checkmate::assert_int(iter_warmup, lower = 1L)
  checkmate::assert_int(iter_sampling, lower = 1L)

  fit_bm <- fit_embedding_bm(
    data = data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )
  fit_ou <- fit_embedding_ou(
    data = data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )

  bm_summary <- extract_ancestral(fit_bm, data = data, tree = tree)
  ou_summary <- extract_ancestral(fit_ou, data = data, tree = tree)

  root_idx <- which(bm_summary$node == data$S + 1L)
  if (length(root_idx) == 0L) {
    root_idx <- 1L
  }

  bm_root <- bm_summary$mean[root_idx, , drop = TRUE]
  ou_root <- ou_summary$mean[root_idx, , drop = TRUE]

  rmse  <- sqrt(mean((bm_root - ou_root)^2))
  cosine <- cosine_sim_r(bm_root, ou_root)

  list(
    bm_fit       = fit_bm,
    ou_fit       = fit_ou,
    bm_summary   = bm_summary,
    ou_summary   = ou_summary,
    root_comparison = list(
      bm_root  = bm_root,
      ou_root  = ou_root,
      rmse     = rmse,
      cosine   = cosine
    )
  )
}
