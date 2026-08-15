#' Summarise posterior ancestral embeddings
#'
#' Takes a `CmdStanMCMC` fit object produced by `fit_embedding_bm()` or
#' `fit_embedding_ou()` and computes posterior summaries for the internal
#' (ancestral) nodes. Internal node estimates are obtained by applying the
#' conditional expectation formula to each posterior draw of the root state
#' and evolutionary parameters.
#'
#' @param fit A `CmdStanMCMC` fit object.
#' @param nodes Optional integer vector of ape node numbers to summarise.
#'   If `NULL`, all internal nodes are returned.
#' @param data Optional Stan data list. If omitted, `fit` must carry an
#'   `ancemb_data` attribute (set automatically by `fit_embedding_*()`).
#' @param tree Optional `phylo` object. If omitted, `fit` must carry an
#'   `ancemb_tree` attribute, or `data` must carry a `tree` attribute.
#' @param probs Quantile probabilities for the uncertainty interval
#'   (default `c(0.05, 0.95)`).
#' @param ... Currently unused.
#'
#' @return A list with `node` (ape node numbers), `mean`, `median`, `sd`,
#'   `q_lower`, `q_upper` (each an `n_nodes x D` matrix), and `draws`, a
#'   `posterior::draws_matrix` object with one column per node-dimension
#'   combination.
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
#'   anc <- extract_ancestral(fit, data = data, tree = tree)
#'   str(anc)
#' }
#' }
extract_ancestral <- function(fit, nodes = NULL, data = NULL, tree = NULL,
                            probs = c(0.05, 0.95), ...) {
  checkmate::assert_function(fit$draws, args = c("inc_warmup", "format"))
  checkmate::assert_numeric(nodes, min.len = 0L, null.ok = TRUE)
  checkmate::assert_numeric(probs, lower = 0, upper = 1, len = 2L,
                            any.missing = FALSE)

  if (is.null(data)) {
    data <- attr(fit, "ancemb_data")
  }
  if (is.null(data)) {
    stop("'data' is required (or set as an 'ancemb_data' attribute on fit).")
  }
  if (is.null(tree)) {
    tree <- attr(fit, "ancemb_tree")
    if (is.null(tree)) {
      tree <- attr(data, "tree")
    }
  }
  if (is.null(tree)) {
    stop("'tree' is required (or set as a 'tree' attribute on data/fit).")
  }

  family <- attr(fit, "ancemb_family")
  if (is.null(family)) {
    family <- .detect_family(fit)
  }
  checkmate::assert_choice(family, c("bm", "ou"))

  draws <- posterior::as_draws_matrix(fit$draws())
  n_draws <- nrow(draws)
  d <- data$D
  n_nodes <- tree$Nnode
  n_tips  <- data$S
  n_total <- n_tips + n_nodes

  if (is.null(nodes)) {
    node_idx <- seq_len(n_nodes)
    node_numbers <- (n_tips + 1L):n_total
  } else {
    node_numbers <- as.integer(nodes)
    if (any(node_numbers <= n_tips | node_numbers > n_total)) {
      stop("All 'nodes' must be internal node numbers between ",
           n_tips + 1L, " and ", n_total, ".")
    }
    node_idx <- node_numbers - n_tips
  }

  # Preallocate draws array: n_draws x n_nodes x d
  node_draws <- array(NA_real_, dim = c(n_draws, length(node_idx), d))
  colnames <- names(draws)

  if (family == "bm") {
    for (i in seq_len(n_draws)) {
      z_anc <- as.numeric(draws[i, paste0("z_anc[", seq_len(d), "]"), drop = FALSE])
      rate  <- as.numeric(draws[i, paste0("rate[",  seq_len(d), "]"), drop = FALSE])
      ref <- reconstruct_ancestral_bm_r(
        tree = tree,
        z_obs = data$z_obs,
        rate = rate,
        z_anc = z_anc,
        obs_sigma = data$obs_sigma
      )
      node_draws[i, , ] <- ref$node_estimates[node_idx, , drop = FALSE]
    }
  } else {
    for (i in seq_len(n_draws)) {
      z_anc  <- as.numeric(draws[i, paste0("z_anc[",  seq_len(d), "]"), drop = FALSE])
      sigma2 <- as.numeric(draws[i, paste0("sigma2[", seq_len(d), "]"), drop = FALSE])
      alpha  <- as.numeric(draws[i, paste0("alpha[",  seq_len(d), "]"), drop = FALSE])
      ref <- reconstruct_ancestral_ou_r(
        tree = tree,
        z_obs = data$z_obs,
        sigma2 = sigma2,
        alpha = alpha,
        z_anc = z_anc,
        obs_sigma = data$obs_sigma
      )
      node_draws[i, , ] <- ref$node_estimates[node_idx, , drop = FALSE]
    }
  }

  mean_mat  <- apply(node_draws, c(2L, 3L), mean)
  median_mat <- apply(node_draws, c(2L, 3L), stats::median)
  sd_mat    <- apply(node_draws, c(2L, 3L), stats::sd)
  q_lower   <- apply(node_draws, c(2L, 3L), stats::quantile, probs = min(probs))
  q_upper   <- apply(node_draws, c(2L, 3L), stats::quantile, probs = max(probs))

  rownames(mean_mat) <- rownames(median_mat) <- node_numbers
  rownames(sd_mat) <- rownames(q_lower) <- rownames(q_upper) <- node_numbers
  colnames(mean_mat) <- colnames(median_mat) <- colnames(sd_mat) <-
    colnames(q_lower) <- colnames(q_upper) <- paste0("dim", seq_len(d))

  # Flatten draws for posterior::draws_matrix
  flat <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = length(node_idx) * d
  )
  var_names <- character(length(node_idx) * d)
  for (j in seq_along(node_idx)) {
    for (k in seq_len(d)) {
      col <- (j - 1L) * d + k
      flat[, col] <- node_draws[, j, k]
      var_names[col] <- paste0("node", node_numbers[j], "[", k, "]")
    }
  }
  colnames(flat) <- var_names
  draws_out <- posterior::as_draws_matrix(flat)

  list(
    node    = node_numbers,
    mean    = mean_mat,
    median  = median_mat,
    sd      = sd_mat,
    q_lower = q_lower,
    q_upper = q_upper,
    draws   = draws_out
  )
}

#' Detect model family from a Stan fit
#'
#' Uses the presence of `rate[1]` (BM) or `sigma2[1]` (OU) in the posterior
#' draws.
#'
#' @param fit A `CmdStanMCMC` object.
#' @return Character `"bm"` or `"ou"`.
#' @keywords internal
.detect_family <- function(fit) {
  draws <- posterior::as_draws_matrix(fit$draws())
  nms <- colnames(draws)
  if (any(grepl("^rate\\[", nms))) {
    "bm"
  } else if (any(grepl("^sigma2\\[", nms))) {
    "ou"
  } else {
    stop("Cannot determine model family from fit draws.")
  }
}
