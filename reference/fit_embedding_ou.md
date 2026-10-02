# Fit an Ornstein-Uhlenbeck embedding model

Loads the \`embedding_phylo_ou\` Stan model (compiling it on first use
if necessary) and runs MCMC sampling.

## Usage

``` r
fit_embedding_ou(
  data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  ...
)
```

## Arguments

- data:

  A named list produced by \`prepare_embedding_data()\`.

- chains:

  Number of MCMC chains (default \`4\`).

- iter_warmup:

  Number of warmup iterations per chain (default \`1000\`).

- iter_sampling:

  Number of post-warmup iterations per chain (default \`1000\`).

- ...:

  Additional arguments passed to \`CmdStanModel\$sample()\`.

## Value

A \`CmdStanMCMC\` fit object.

## Examples

``` r
# \donttest{
if (instantiate::stan_cmdstan_exists()) {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  z_obs <- matrix(stats::rnorm(8 * 4), nrow = 8)
  rownames(z_obs) <- tree$tip.label
  data <- prepare_embedding_data(tree, z_obs)
  fit <- fit_embedding_ou(data, chains = 1, iter_warmup = 200,
                          iter_sampling = 200)
}
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 400 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: Matrix m is not positive definite (in '/tmp/Rtmp1wWU0u/model-1dc01388f553.stan', line 62, column 6 to column 58)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 400 [ 25%]  (Warmup) 
#> Chain 1 Iteration: 200 / 400 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 201 / 400 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 300 / 400 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 400 / 400 [100%]  (Sampling) 
#> Chain 1 finished in 1.5 seconds.
# }
```
