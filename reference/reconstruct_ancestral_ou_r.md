# Reconstruct ancestral embeddings under Ornstein-Uhlenbeck (pure-R reference)

Computes the fixed-root OU covariance among all nodes and then evaluates
the conditional expectation of the internal nodes given the tip
observations. If \`z_anc\` is \`NULL\` it is estimated by GLS using the
provided (or default) OU parameters.

## Usage

``` r
reconstruct_ancestral_ou_r(
  tree,
  z_obs,
  sigma2 = NULL,
  alpha = NULL,
  z_anc = NULL,
  obs_sigma = 0
)
```

## Arguments

- tree:

  A \`phylo\` object.

- z_obs:

  Numeric matrix (\`S x D\`) of observed tip embeddings.

- sigma2:

  Diffusion variance per dimension. Defaults to 1.

- alpha:

  Selection strength per dimension. Defaults to 1.

- z_anc:

  Optional root state per dimension.

- obs_sigma:

  Observation noise sd.

## Value

List with \`z_anc\`, \`sigma2\`, \`alpha\`, and \`node_estimates\`.
