# Compare Brownian Motion and Ornstein-Uhlenbeck embedding fits

Fits both the BM and OU Stan models to the same prepared data and
returns the fits, ancestral summaries, and a simple comparison of the
inferred root embeddings.

## Usage

``` r
run_bayesian_comparison(
  data,
  tree = NULL,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  ...
)
```

## Arguments

- data:

  Stan data list produced by \`prepare_embedding_data()\`.

- tree:

  A \`phylo\` object. Required only if \`data\` does not carry a
  \`tree\` attribute.

- chains:

  Number of MCMC chains per model (default \`4\`).

- iter_warmup:

  Number of warmup iterations (default \`1000\`).

- iter_sampling:

  Number of post-warmup iterations (default \`1000\`).

- ...:

  Additional arguments passed to \`fit_embedding_bm()\` and
  \`fit_embedding_ou()\`.

## Value

A list with components \`bm_fit\`, \`ou_fit\`, \`bm_summary\`,
\`ou_summary\`, and \`root_comparison\`.

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
  cmp <- run_bayesian_comparison(data, chains = 1, iter_warmup = 200,
                                 iter_sampling = 200)
  str(cmp$root_comparison)
}
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 400 [  0%]  (Warmup) 
#> Chain 1 Iteration: 100 / 400 [ 25%]  (Warmup) 
#> Chain 1 Iteration: 200 / 400 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 201 / 400 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 300 / 400 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 400 / 400 [100%]  (Sampling) 
#> Chain 1 finished in 0.2 seconds.
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 400 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/Rtmp49Uj6K/model-1e511c7e7ff4.stan', line 62, column 6 to column 58)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 400 [ 25%]  (Warmup) 
#> Chain 1 Iteration: 200 / 400 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 201 / 400 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 300 / 400 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 400 / 400 [100%]  (Sampling) 
#> Chain 1 finished in 1.9 seconds.
#> List of 4
#>  $ bm_root: Named num [1:4] 0.888 -0.148 0.412 0.275
#>   ..- attr(*, "names")= chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ ou_root: Named num [1:4] 0.822 -0.198 0.355 0.294
#>   ..- attr(*, "names")= chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ rmse   : num 0.051
#>  $ cosine : num 0.997
# }
```
