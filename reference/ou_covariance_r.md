# Ornstein-Uhlenbeck covariance for all node pairs

Uses the fixed-root OU covariance: V\[i,j,k\] = (sigma2\[k\] /
(2\*alpha\[k\])) \* (exp(-alpha\[k\] \* (t_i + t_j - 2\*t_mrca(i,j))) -
exp(-alpha\[k\] \* (t_i + t_j)))

## Usage

``` r
ou_covariance_r(tree, alpha, sigma2)
```

## Arguments

- tree:

  A \`phylo\` object.

- alpha:

  Numeric vector of length 1 or \`D\`.

- sigma2:

  Numeric vector of length 1 or \`D\`.

## Value

Matrix or 3D array as in \`bm_covariance_r\`.
