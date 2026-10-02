# Simulate OU evolution of an embedding along a phylogeny

Uses Euler-Maruyama on each branch with step size \`dt\`.

## Usage

``` r
simulate_ou_embedding_r(tree, z_root, theta = z_root, alpha, sigma2, dt = 0.01)
```

## Arguments

- tree:

  A \`phylo\` object.

- z_root:

  Numeric vector, root state.

- theta:

  Numeric vector, OU optimum. Defaults to \`z_root\`.

- alpha:

  Numeric vector, selection strength per dimension.

- sigma2:

  Numeric vector, diffusion variance per dimension.

- dt:

  Euler-Maruyama step size as a fraction of branch length.

## Value

Same structure as \`simulate_bm_embedding_r\`.
