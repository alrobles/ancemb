#' Reference implementations for ancestral embedding reconstruction
#'
#' Internal pure-R functions that reproduce the Brownian Motion (BM) and
#' Ornstein-Uhlenbeck (OU) phylogenetic embedding math line-by-line. They
#' are not exported and have no input validation -- callers are responsible
#' for providing correctly formatted objects.
#'
#' @keywords internal
#' @name reference
NULL

#' Simulate BM evolution of an embedding along a phylogeny
#'
#' @param tree A `phylo` object.
#' @param z_root Numeric vector, root state (length `D`).
#' @param rates Numeric vector, BM rate per dimension (length `D`).
#'
#' @return List with `tip_values` (`n_tips` x `D`), `node_values`
#'   (`n_nodes` x `D`), `all_values` (`n_total` x `D`) and `root_value`.
#' @keywords internal
simulate_bm_embedding_r <- function(tree, z_root, rates) {
  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  d       <- length(z_root)

  stopifnot(length(rates) == d)

  values <- matrix(NA_real_, nrow = n_total, ncol = d)
  root_idx <- n_tips + 1L
  values[root_idx, ] <- z_root

  tree_reordered <- ape::reorder.phylo(tree, "cladewise")
  for (i in seq_len(nrow(tree_reordered$edge))) {
    parent <- tree_reordered$edge[i, 1]
    child  <- tree_reordered$edge[i, 2]
    bl     <- tree_reordered$edge.length[i]

    noise <- vapply(seq_len(d), function(k) {
      stats::rnorm(1, mean = 0, sd = sqrt(rates[k] * bl))
    }, numeric(1))

    values[child, ] <- values[parent, ] + noise
  }

  list(
    tip_values  = values[1:n_tips, , drop = FALSE],
    node_values = values[(n_tips + 1):n_total, , drop = FALSE],
    all_values  = values,
    root_value  = z_root
  )
}

#' Simulate OU evolution of an embedding along a phylogeny
#'
#' Uses Euler-Maruyama on each branch with step size `dt`.
#'
#' @param tree A `phylo` object.
#' @param z_root Numeric vector, root state.
#' @param theta Numeric vector, OU optimum. Defaults to `z_root`.
#' @param alpha Numeric vector, selection strength per dimension.
#' @param sigma2 Numeric vector, diffusion variance per dimension.
#' @param dt Euler-Maruyama step size as a fraction of branch length.
#'
#' @return Same structure as `simulate_bm_embedding_r`.
#' @keywords internal
simulate_ou_embedding_r <- function(tree, z_root, theta = z_root,
                                    alpha, sigma2, dt = 0.01) {
  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  d       <- length(z_root)

  stopifnot(
    length(alpha) == d,
    length(sigma2) == d,
    length(theta) == d
  )

  values <- matrix(NA_real_, nrow = n_total, ncol = d)
  root_idx <- n_tips + 1L
  values[root_idx, ] <- z_root

  tree_reordered <- ape::reorder.phylo(tree, "cladewise")
  for (i in seq_len(nrow(tree_reordered$edge))) {
    parent <- tree_reordered$edge[i, 1]
    child  <- tree_reordered$edge[i, 2]
    bl     <- tree_reordered$edge.length[i]

    n_steps <- max(1L, as.integer(ceiling(bl / dt)))
    actual_dt <- bl / n_steps

    state <- values[parent, ]
    for (step in seq_len(n_steps)) {
      drift <- alpha * (theta - state) * actual_dt
      noise <- vapply(seq_len(d), function(k) {
        stats::rnorm(1, mean = 0, sd = sqrt(sigma2[k] * actual_dt))
      }, numeric(1))
      state <- state + drift + noise
    }
    values[child, ] <- state
  }

  list(
    tip_values  = values[1:n_tips, , drop = FALSE],
    node_values = values[(n_tips + 1):n_total, , drop = FALSE],
    all_values  = values,
    root_value  = z_root
  )
}

#' Cosine similarity between two vectors
#'
#' @param a,b Numeric vectors.
#' @keywords internal
cosine_sim_r <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
}

#' Build the phylogenetic shared-height matrix for all nodes
#'
#' Returns the matrix `M` where `M[i,j]` is the height of the MRCA of nodes
#' `i` and `j`. This is the generalisation of `ape::vcv` to all nodes.
#'
#' @param tree A `phylo` object.
#' @return A symmetric matrix with `n_tips + n_nodes` rows/columns.
#' @keywords internal
build_all_node_covariance_r <- function(tree) {
  heights <- ape::node.depth.edgelength(tree)
  dist_mat <- ape::dist.nodes(tree)
  M <- (outer(heights, heights, "+") - dist_mat) / 2
  M[M < 0] <- 0
  M
}

#' Brownian Motion covariance for all node pairs
#'
#' @param tree A `phylo` object.
#' @param rate Numeric vector of length 1 or `D`.
#' @return Matrix or 3D array. If `rate` is length 1, a matrix is returned;
#'   otherwise an `n_total x n_total x D` array.
#' @keywords internal
bm_covariance_r <- function(tree, rate) {
  M <- build_all_node_covariance_r(tree)
  d <- length(rate)
  if (d == 1L) {
    return(rate * M)
  }
  n_total <- nrow(M)
  V <- array(NA_real_, dim = c(n_total, n_total, d))
  for (k in seq_len(d)) {
    V[, , k] <- rate[k] * M
  }
  V
}

#' Ornstein-Uhlenbeck covariance for all node pairs
#'
#' Uses the fixed-root OU covariance:
#'   V[i,j,k] = (sigma2[k] / (2*alpha[k])) *
#'              (exp(-alpha[k] * (t_i + t_j - 2*t_mrca(i,j)))
#'               - exp(-alpha[k] * (t_i + t_j)))
#'
#' @param tree A `phylo` object.
#' @param alpha Numeric vector of length 1 or `D`.
#' @param sigma2 Numeric vector of length 1 or `D`.
#' @return Matrix or 3D array as in `bm_covariance_r`.
#' @keywords internal
ou_covariance_r <- function(tree, alpha, sigma2) {
  heights <- ape::node.depth.edgelength(tree)
  dist_mat <- ape::dist.nodes(tree)
  M <- (outer(heights, heights, "+") - dist_mat) / 2
  M[M < 0] <- 0
  total <- outer(heights, heights, "+")
  d <- max(length(alpha), length(sigma2))
  n_total <- length(heights)
  if (d == 1L) {
    scale <- sigma2 / (2 * alpha)
    return(scale * (exp(-alpha * (total - 2 * M)) - exp(-alpha * total)))
  }
  V <- array(NA_real_, dim = c(n_total, n_total, d))
  for (k in seq_len(d)) {
    a <- alpha[min(k, length(alpha))]
    s <- sigma2[min(k, length(sigma2))]
    scale <- s / (2 * a)
    V[, , k] <- scale * (exp(-a * (total - 2 * M)) - exp(-a * total))
  }
  V
}

#' Reconstruct ancestral embeddings under Brownian Motion (pure-R reference)
#'
#' Estimates the root state with the GLS estimator and then computes the
#' conditional expectation of every internal node. With `obs_sigma = 0`
#' the resulting node estimates match `phytools::fastAnc`.
#'
#' @param tree A `phylo` object.
#' @param z_obs Numeric matrix (`S x D`) of observed tip embeddings.
#' @param rate Optional per-dimension BM rate. If `NULL`, estimated via ML
#'   under `obs_sigma = 0`.
#' @param z_anc Optional per-dimension root state. If `NULL`, estimated by GLS.
#' @param obs_sigma Observation noise sd added to the diagonal of the tip
#'   covariance matrix.
#'
#' @return List with `z_anc` (length `D`), `rate` (length `D`), and
#'   `node_estimates` (`n_nodes x D`). Row `u` of `node_estimates` corresponds
#'   to ape node number `n_tips + u`.
#' @keywords internal
reconstruct_ancestral_bm_r <- function(tree, z_obs,
                                       rate = NULL, z_anc = NULL,
                                       obs_sigma = 0) {
  if (!is.null(rownames(z_obs))) {
    z_obs <- z_obs[tree$tip.label, , drop = FALSE]
  }

  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  d       <- ncol(z_obs)

  M <- build_all_node_covariance_r(tree)
  C <- M[seq_len(n_tips), seq_len(n_tips), drop = FALSE]

  # ML / GLS estimates from the noise-free tip covariance C
  if (is.null(rate) || is.null(z_anc)) {
    invC <- solve(C)
    one  <- rep(1, n_tips)
    den  <- as.numeric(crossprod(one, invC %*% one))
  }

  if (is.null(z_anc)) {
    z_anc_est <- as.numeric(crossprod(one, invC %*% z_obs)) / den
  } else {
    z_anc_est <- rep_len(z_anc, d)
  }

  if (is.null(rate)) {
    rate_est <- vapply(seq_len(d), function(k) {
      r <- z_obs[, k] - z_anc_est[k]
      as.numeric(crossprod(r, invC %*% r)) / n_tips
    }, numeric(1))
  } else {
    rate_est <- rep_len(rate, d)
  }

  node_est <- matrix(NA_real_, nrow = n_nodes, ncol = d)
  colnames(node_est) <- colnames(z_obs)
  rownames(node_est) <- (n_tips + 1L):n_total

  for (k in seq_len(d)) {
    V <- rate_est[k] * C
    if (obs_sigma > 0) {
      diag(V) <- diag(V) + obs_sigma^2
    }
    V_inv <- solve(V)
    for (u in seq_len(n_nodes)) {
      node_idx <- n_tips + u
      c_u <- rate_est[k] * M[node_idx, seq_len(n_tips)]
      resid <- z_obs[, k] - z_anc_est[k]
      node_est[u, k] <- z_anc_est[k] +
        as.numeric(c_u %*% V_inv %*% resid)
    }
  }

  list(
    z_anc         = z_anc_est,
    rate          = rate_est,
    node_estimates = node_est
  )
}

#' Reconstruct ancestral embeddings under Ornstein-Uhlenbeck (pure-R reference)
#'
#' Computes the fixed-root OU covariance among all nodes and then evaluates
#' the conditional expectation of the internal nodes given the tip
#' observations. If `z_anc` is `NULL` it is estimated by GLS using the
#' provided (or default) OU parameters.
#'
#' @param tree A `phylo` object.
#' @param z_obs Numeric matrix (`S x D`) of observed tip embeddings.
#' @param sigma2 Diffusion variance per dimension. Defaults to 1.
#' @param alpha Selection strength per dimension. Defaults to 1.
#' @param z_anc Optional root state per dimension.
#' @param obs_sigma Observation noise sd.
#'
#' @return List with `z_anc`, `sigma2`, `alpha`, and `node_estimates`.
#' @keywords internal
reconstruct_ancestral_ou_r <- function(tree, z_obs,
                                     sigma2 = NULL, alpha = NULL,
                                     z_anc = NULL, obs_sigma = 0) {
  if (!is.null(rownames(z_obs))) {
    z_obs <- z_obs[tree$tip.label, , drop = FALSE]
  }

  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  d       <- ncol(z_obs)

  if (is.null(sigma2)) sigma2 <- rep(1, d)
  if (is.null(alpha))  alpha  <- rep(1, d)
  sigma2 <- rep_len(sigma2, d)
  alpha  <- rep_len(alpha, d)

  V_all <- ou_covariance_r(tree, alpha = alpha, sigma2 = sigma2)

  V_tt <- V_all[seq_len(n_tips), seq_len(n_tips), , drop = FALSE]
  V_it <- V_all[(n_tips + 1):n_total, seq_len(n_tips), , drop = FALSE]

  if (is.null(z_anc)) {
    z_anc_est <- vapply(seq_len(d), function(k) {
      V <- V_tt[, , k]
      if (obs_sigma > 0) diag(V) <- diag(V) + obs_sigma^2
      V_inv <- solve(V)
      one <- rep(1, n_tips)
      as.numeric(crossprod(one, V_inv %*% z_obs[, k])) /
        as.numeric(crossprod(one, V_inv %*% one))
    }, numeric(1))
  } else {
    z_anc_est <- rep_len(z_anc, d)
  }

  node_est <- matrix(NA_real_, nrow = n_nodes, ncol = d)
  colnames(node_est) <- colnames(z_obs)
  rownames(node_est) <- (n_tips + 1L):n_total

  for (k in seq_len(d)) {
    V <- V_tt[, , k]
    if (obs_sigma > 0) diag(V) <- diag(V) + obs_sigma^2
    V_inv <- solve(V)
    for (u in seq_len(n_nodes)) {
      c_u <- V_it[u, , k]
      resid <- z_obs[, k] - z_anc_est[k]
      node_est[u, k] <- z_anc_est[k] +
        as.numeric(c_u %*% V_inv %*% resid)
    }
  }

  list(
    z_anc         = z_anc_est,
    sigma2        = sigma2,
    alpha         = alpha,
    node_estimates = node_est
  )
}

#' Reconstruct ancestral embeddings with `phytools::fastAnc`
#'
#' Convenience wrapper around `phytools::fastAnc` that returns the full
#' matrix of ancestral node estimates (one column per embedding dimension).
#'
#' @param tree A `phylo` object.
#' @param tip_embeddings Numeric matrix (`n_tips x D`).
#' @return Numeric matrix (`n_nodes x D`). Row 1 is the root estimate.
#' @keywords internal
reconstruct_ancestral_ml_r <- function(tree, tip_embeddings) {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is required for 'reconstruct_ancestral_ml_r()' ",
         "but is not installed.", call. = FALSE)
  }
  d <- ncol(tip_embeddings)
  n_nodes <- tree$Nnode

  anc_matrix <- matrix(NA_real_, nrow = n_nodes, ncol = d)
  for (k in seq_len(d)) {
    tip_values <- tip_embeddings[, k]
    names(tip_values) <- tree$tip.label
    anc_matrix[, k] <- phytools::fastAnc(tree, tip_values)
  }
  rownames(anc_matrix) <- (length(tree$tip.label) + 1L):(length(tree$tip.label) + n_nodes)
  anc_matrix
}
