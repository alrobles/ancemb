# Summarise posterior ancestral embeddings

Takes a \`CmdStanMCMC\` fit object produced by \`fit_embedding_bm()\` or
\`fit_embedding_ou()\` and computes posterior summaries for the internal
(ancestral) nodes. Internal node estimates are obtained by applying the
conditional expectation formula to each posterior draw of the root state
and evolutionary parameters.

## Usage

``` r
extract_ancestral(
  fit,
  nodes = NULL,
  data = NULL,
  tree = NULL,
  probs = c(0.05, 0.95),
  ...
)
```

## Arguments

- fit:

  A \`CmdStanMCMC\` fit object.

- nodes:

  Optional integer vector of ape node numbers to summarise. If \`NULL\`,
  all internal nodes are returned.

- data:

  Optional Stan data list. If omitted, \`fit\` must carry an
  \`ancemb_data\` attribute (set automatically by
  \`fit_embedding\_\*()\`).

- tree:

  Optional \`phylo\` object. If omitted, \`fit\` must carry an
  \`ancemb_tree\` attribute, or \`data\` must carry a \`tree\`
  attribute.

- probs:

  Quantile probabilities for the uncertainty interval (default \`c(0.05,
  0.95)\`).

- ...:

  Currently unused.

## Value

A list with \`node\` (ape node numbers), \`mean\`, \`median\`, \`sd\`,
\`q_lower\`, \`q_upper\` (each an \`n_nodes x D\` matrix), and
\`draws\`, a \`posterior::draws_matrix\` object with one column per
node-dimension combination.

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
  anc <- extract_ancestral(fit, data = data, tree = tree)
  str(anc)
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
#> List of 7
#>  $ node   : int [1:7] 9 10 11 12 13 14 15
#>  $ mean   : num [1:7, 1:4] 0.888 0.821 0.686 0.308 0.651 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ median : num [1:7, 1:4] 0.908 0.828 0.689 0.308 0.653 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ sd     : num [1:7, 1:4] 0.2274 0.0744 0.0288 0.0096 0.0118 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ q_lower: num [1:7, 1:4] 0.476 0.688 0.635 0.293 0.629 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ q_upper: num [1:7, 1:4] 1.244 0.938 0.732 0.324 0.666 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ draws  : 'draws_matrix' num [1:200, 1:28] 0.914 0.969 1.43 0.388 0.77 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ draw    : chr [1:200] "1" "2" "3" "4" ...
#>   .. ..$ variable: chr [1:28] "node9[1]" "node9[2]" "node9[3]" "node9[4]" ...
#>   ..- attr(*, "nchains")= int 1
# }
```
