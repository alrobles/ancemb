# Fit a Brownian Motion embedding model

Loads the \`embedding_phylo_bm\` Stan model (compiling it on first use
if necessary) and runs MCMC sampling. The returned \`CmdStanMCMC\`
object carries \`ancemb_data\` and \`ancemb_family\` attributes for use
with \`extract_ancestral()\`.

## Usage

``` r
fit_embedding_bm(
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
  fit <- fit_embedding_bm(data, chains = 1, iter_warmup = 200,
                          iter_sampling = 200)
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
# }
```
