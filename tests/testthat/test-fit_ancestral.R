#' Test suite for ancestral embedding reconstruction

set.seed(42)

make_test_tree <- function(n_tips = 8) {
  tree <- ape::rtree(n_tips)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  tree
}

make_bm_embeddings <- function(tree, d = 3, rate = 0.5, root = NULL) {
  if (is.null(root)) root <- stats::rnorm(d)
  sim <- simulate_bm_embedding_r(tree, z_root = root, rates = rep(rate, d))
  Z <- sim$tip_values
  rownames(Z) <- tree$tip.label
  Z
}

test_that("pure-R BM reconstruction matches phytools::fastAnc", {
  tree <- make_test_tree(8)
  Z <- make_bm_embeddings(tree, d = 5, rate = 0.3)

  ref <- reconstruct_ancestral_bm_r(tree, Z)
  ml  <- reconstruct_ancestral_ml_r(tree, Z)

  expect_equal(ref$z_anc, ml[1L, ], tolerance = 1e-10)
  expect_equal(ref$node_estimates, ml, tolerance = 1e-10)
})

test_that("pure-R OU reconstruction recovers supplied root", {
  tree <- make_test_tree(8)
  d <- 4
  z_root <- stats::rnorm(d)
  alpha  <- runif(d, 0.5, 1.5)
  sigma2 <- runif(d, 0.2, 0.8)

  sim <- simulate_ou_embedding_r(tree, z_root, alpha = alpha, sigma2 = sigma2)
  Z <- sim$tip_values
  rownames(Z) <- tree$tip.label

  ref <- reconstruct_ancestral_ou_r(tree, Z, alpha = alpha, sigma2 = sigma2,
                                    z_anc = z_root)
  expect_equal(ref$node_estimates[1L, ], z_root, tolerance = 1e-10)
})

test_that("OU with small alpha approximates BM", {
  tree <- make_test_tree(8)
  d <- 3
  rate <- 0.4
  alpha <- rep(1e-4, d)
  sigma2 <- rate * 2 * alpha
  z_root <- stats::rnorm(d)

  sim <- simulate_ou_embedding_r(tree, z_root, alpha = alpha, sigma2 = sigma2)
  Z <- sim$tip_values
  rownames(Z) <- tree$tip.label

  ref_ou <- reconstruct_ancestral_ou_r(tree, Z, alpha = alpha, sigma2 = sigma2,
                                       z_anc = z_root)
  ref_bm <- reconstruct_ancestral_bm_r(tree, Z, rate = rate, z_anc = z_root)

  expect_equal(ref_ou$node_estimates, ref_bm$node_estimates, tolerance = 0.05)
})

test_that("prepare_embedding_data builds a valid Stan data list", {
  tree <- make_test_tree(8)
  Z <- make_bm_embeddings(tree, d = 4)

  data <- prepare_embedding_data(tree, Z, obs_sigma = 0.05)

  expect_type(data, "list")
  expect_equal(data$S, 8)
  expect_equal(data$D, 4)
  expect_equal(rownames(data$z_obs), tree$tip.label)
  expect_equal(dim(data$C), c(8, 8))
  expect_equal(length(data$tip_heights), 8)
  expect_type(data$tip_heights, "double")
  expect_true(!is.null(attr(data, "tree")))
})

test_that("prepare_embedding_data validates inputs", {
  tree <- make_test_tree(8)
  Z <- make_bm_embeddings(tree, d = 4)

  rownames(Z)[1] <- "wrong"
  expect_error(prepare_embedding_data(tree, Z))

  expect_error(prepare_embedding_data(tree, Z, obs_sigma = -1))
})

test_that("Stan BM fit and ancestral extraction agree with fastAnc", {
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  tree <- make_test_tree(8)
  Z <- make_bm_embeddings(tree, d = 3, rate = 0.3)

  data <- prepare_embedding_data(tree, Z, obs_sigma = 0.05)
  fit <- fit_embedding_bm(data, chains = 1, iter_warmup = 100,
                          iter_sampling = 100, refresh = 0)

  expect_s3_class(fit, "CmdStanMCMC")

  anc <- extract_ancestral(fit, data = data, tree = tree)
  root_node <- which(anc$node == data$S + 1L)

  ml_root <- unname(phytools::fastAnc(tree, Z[, 1])[1])

  expect_equal(unname(anc$mean[root_node, "dim1"]), ml_root, tolerance = 0.15)

  ref <- reconstruct_ancestral_bm_r(tree, Z)
  expect_equal(unname(anc$mean[root_node, ]), ref$z_anc, tolerance = 0.15)
})

test_that("Stan OU fit runs and extract_ancestral returns summaries", {
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  tree <- make_test_tree(8)
  d <- 2
  z_root <- stats::rnorm(d)
  alpha  <- runif(d, 0.5, 1.0)
  sigma2 <- runif(d, 0.2, 0.5)

  sim <- simulate_ou_embedding_r(tree, z_root, alpha = alpha, sigma2 = sigma2)
  Z <- sim$tip_values
  rownames(Z) <- tree$tip.label

  data <- prepare_embedding_data(tree, Z, obs_sigma = 0.05)
  fit <- fit_embedding_ou(data, chains = 1, iter_warmup = 100,
                          iter_sampling = 100, refresh = 0)

  expect_s3_class(fit, "CmdStanMCMC")

  anc <- extract_ancestral(fit, data = data, tree = tree)
  expect_equal(rownames(anc$mean), as.character((data$S + 1L):(data$S + tree$Nnode)))
  expect_true(all(anc$sd > 0))
})

test_that("run_bayesian_comparison returns a structured comparison", {
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  tree <- make_test_tree(8)
  Z <- make_bm_embeddings(tree, d = 2, rate = 0.4)
  data <- prepare_embedding_data(tree, Z, obs_sigma = 0.05)

  cmp <- run_bayesian_comparison(data, tree, chains = 1, iter_warmup = 100,
                                 iter_sampling = 100, refresh = 0)

  expect_named(cmp, c("bm_fit", "ou_fit", "bm_summary", "ou_summary",
                      "root_comparison"))
  expect_named(cmp$root_comparison, c("bm_root", "ou_root", "rmse", "cosine"))
  expect_length(cmp$root_comparison$rmse, 1L)
  expect_length(cmp$root_comparison$cosine, 1L)
})
