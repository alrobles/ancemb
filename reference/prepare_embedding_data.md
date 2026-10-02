# Prepare embedding data for Stan

Converts observed tip embeddings and a phylogenetic tree into the named
list required by the Brownian Motion and Ornstein-Uhlenbeck Stan models.
The tree is stored as an attribute on the returned list so that
downstream functions (e.g. \`extract_ancestral()\`) can reconstruct
internal node states without needing the tree to be passed again.

## Usage

``` r
prepare_embedding_data(
  tree,
  z_obs,
  obs_sigma = 0.1,
  rate_scale = 1,
  alpha_scale = 1
)
```

## Arguments

- tree:

  A \`phylo\` object. Tip labels and branch lengths must be present.

- z_obs:

  Numeric matrix (\`S x D\`) of observed tip embeddings. Rownames must
  match \`tree\$tip.label\` exactly.

- obs_sigma:

  Non-negative observation noise sd (default \`0.1\`).

- rate_scale:

  Half-normal scale for BM rates / OU diffusion variance (default
  \`1\`).

- alpha_scale:

  Half-normal scale for OU selection strength (default \`1\`).

## Value

A named list suitable for \`fit_embedding_bm()\` and
\`fit_embedding_ou()\`. The \`tree\` argument is attached as a \`tree\`
attribute.

## Examples

``` r
# \donttest{
if (requireNamespace("ape", quietly = TRUE)) {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  z_obs <- matrix(stats::rnorm(8 * 4), nrow = 8)
  rownames(z_obs) <- tree$tip.label
  stan_data <- prepare_embedding_data(tree, z_obs)
  str(stan_data)
}
#> List of 10
#>  $ S               : int 8
#>  $ D               : int 4
#>  $ z_obs           : num [1:8, 1:4] 1.1249 -0.0449 -0.0162 0.9438 0.8212 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:8] "t5" "t7" "t3" "t2" ...
#>   .. ..$ : NULL
#>  $ C               : num [1:8, 1:8] 0.196 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:8] "t5" "t7" "t3" "t2" ...
#>   .. ..$ : chr [1:8] "t5" "t7" "t3" "t2" ...
#>  $ tip_heights     : num [1:8] 0.196 0.66 1 0.816 0.476 ...
#>  $ obs_sigma       : num 0.1
#>  $ z_anc_prior_mean: num [1:4] 0.64 -0.38 -0.082 0.163
#>  $ z_anc_prior_sd  : num [1:4] 0.441 0.909 0.773 0.701
#>  $ rate_scale      : num 1
#>  $ alpha_scale     : num 1
#>  - attr(*, "tree")=List of 4
#>   ..$ edge       : int [1:14, 1:2] 9 9 10 11 12 12 13 13 11 10 ...
#>   ..$ tip.label  : chr [1:8] "t5" "t7" "t3" "t2" ...
#>   ..$ Nnode      : int 7
#>   ..$ edge.length: num [1:14] 0.196 0.127 0.183 0.253 0.097 ...
#>   ..- attr(*, "class")= chr "phylo"
#>   ..- attr(*, "order")= chr "cladewise"
# }
```
