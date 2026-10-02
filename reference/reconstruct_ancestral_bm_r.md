# Reconstruct ancestral embeddings under Brownian Motion (pure-R reference)

Estimates the root state with the GLS estimator and then computes the
conditional expectation of every internal node. With \`obs_sigma = 0\`
the resulting node estimates match \`phytools::fastAnc\`.

## Usage

``` r
reconstruct_ancestral_bm_r(
  tree,
  z_obs,
  rate = NULL,
  z_anc = NULL,
  obs_sigma = 0
)
```

## Arguments

- tree:

  A \`phylo\` object.

- z_obs:

  Numeric matrix (\`S x D\`) of observed tip embeddings.

- rate:

  Optional per-dimension BM rate. If \`NULL\`, estimated via ML under
  \`obs_sigma = 0\`.

- z_anc:

  Optional per-dimension root state. If \`NULL\`, estimated by GLS.

- obs_sigma:

  Observation noise sd added to the diagonal of the tip covariance
  matrix.

## Value

List with \`z_anc\` (length \`D\`), \`rate\` (length \`D\`), and
\`node_estimates\` (\`n_nodes x D\`). Row \`u\` of \`node_estimates\`
corresponds to ape node number \`n_tips + u\`.
