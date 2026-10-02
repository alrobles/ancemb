# Brownian Motion covariance for all node pairs

Brownian Motion covariance for all node pairs

## Usage

``` r
bm_covariance_r(tree, rate)
```

## Arguments

- tree:

  A \`phylo\` object.

- rate:

  Numeric vector of length 1 or \`D\`.

## Value

Matrix or 3D array. If \`rate\` is length 1, a matrix is returned;
otherwise an \`n_total x n_total x D\` array.
