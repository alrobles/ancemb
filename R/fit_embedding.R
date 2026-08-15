#' Fit a Brownian Motion embedding model
#'
#' Loads the `embedding_phylo_bm` Stan model (compiling it on first use if
#' necessary) and runs MCMC sampling. The returned `CmdStanMCMC` object
#' carries `ancemb_data` and `ancemb_family` attributes for use with
#' `extract_ancestral()`.
#'
#' @param data A named list produced by `prepare_embedding_data()`.
#' @param chains Number of MCMC chains (default `4`).
#' @param iter_warmup Number of warmup iterations per chain (default `1000`).
#' @param iter_sampling Number of post-warmup iterations per chain (default `1000`).
#' @param ... Additional arguments passed to `CmdStanModel$sample()`.
#'
#' @return A `CmdStanMCMC` fit object.
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
#'   fit <- fit_embedding_bm(data, chains = 1, iter_warmup = 200,
#'                           iter_sampling = 200)
#' }
#' }
fit_embedding_bm <- function(data, chains = 4, iter_warmup = 1000,
                           iter_sampling = 1000, ...) {
  .fit_embedding_model(
    data = data,
    family = "bm",
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )
}

#' Fit an Ornstein-Uhlenbeck embedding model
#'
#' Loads the `embedding_phylo_ou` Stan model (compiling it on first use if
#' necessary) and runs MCMC sampling.
#'
#' @inheritParams fit_embedding_bm
#' @return A `CmdStanMCMC` fit object.
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
#'   fit <- fit_embedding_ou(data, chains = 1, iter_warmup = 200,
#'                           iter_sampling = 200)
#' }
#' }
fit_embedding_ou <- function(data, chains = 4, iter_warmup = 1000,
                             iter_sampling = 1000, ...) {
  .fit_embedding_model(
    data = data,
    family = "ou",
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )
}

#' Internal helper that compiles and samples a Stan model
#'
#' @param data Stan data list.
#' @param family Either `"bm"` or `"ou"`.
#' @param chains,iter_warmup,iter_sampling Passed to `$sample()`.
#' @param ... Additional arguments passed to `$sample()`.
#'
#' @keywords internal
.fit_embedding_model <- function(data, family,
                                   chains = 4, iter_warmup = 1000,
                                   iter_sampling = 1000, ...) {
  checkmate::assert_list(data, names = "named")
  checkmate::assert_choice(family, c("bm", "ou"))
  checkmate::assert_int(chains, lower = 1L)
  checkmate::assert_int(iter_warmup, lower = 1L)
  checkmate::assert_int(iter_sampling, lower = 1L)

  required <- c("S", "D", "z_obs", "C", "tip_heights",
                "obs_sigma", "z_anc_prior_mean",
                "z_anc_prior_sd", "rate_scale")
  if (family == "ou") {
    required <- c(required, "alpha_scale")
  }
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing required Stan data fields: ", paste(missing, collapse = ", "))
  }

  if (!instantiate::stan_cmdstan_exists()) {
    stop(
      "CmdStan is not available. Install it with:\n",
      "  cmdstanr::install_cmdstan()\n",
      "and reload 'ancemb' to compile the Stan models."
    )
  }

  model_name <- paste0("embedding_phylo_", family)
  model <- instantiate::stan_package_model(
    name = model_name,
    package = "ancemb",
    compile = TRUE
  )

  fit <- model$sample(
    data = data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )

  attr(fit, "ancemb_data") <- data
  attr(fit, "ancemb_family") <- family
  fit
}
