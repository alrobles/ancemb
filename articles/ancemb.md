# ancemb: Bayesian Ancestral Protein Embeddings

`ancemb` provides a small set of functions to reconstruct ancestral
protein embeddings from extant (tip) embeddings and a phylogenetic tree.
The package uses pre-compiled Stan models for Brownian Motion (BM) and
Ornstein-Uhlenbeck (OU) evolution of continuous embedding dimensions.

## Setup

``` r

library(ancemb)
set.seed(42)
```

## Simulated data

We create a small tree and a synthetic set of tip embeddings:

``` r

tree <- ape::rtree(8)
tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))

z_obs <- matrix(stats::rnorm(8 * 4, sd = 0.5), nrow = 8)
rownames(z_obs) <- tree$tip.label
```

## Fitting the BM model

``` r

stan_data <- prepare_embedding_data(tree, z_obs, obs_sigma = 0.1)
fit <- fit_embedding_bm(stan_data, chains = 1, iter_warmup = 100,
                        iter_sampling = 100, refresh = 0)
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10 
#> Chain 1 finished in 0.1 seconds.
```

## Extracting ancestral states

[`extract_ancestral()`](https://alrobles.github.io/ancemb/reference/extract_ancestral.md)
computes posterior summaries for all internal nodes:

``` r

anc <- extract_ancestral(fit, data = stan_data, tree = tree)
str(anc)
#> List of 7
#>  $ node   : int [1:7] 9 10 11 12 13 14 15
#>  $ mean   : num [1:7, 1:4] 0.432 0.213 0.305 0.296 0.106 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ median : num [1:7, 1:4] 0.435 0.219 0.305 0.296 0.107 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ sd     : num [1:7, 1:4] 0.17104 0.06296 0.01587 0.01196 0.00736 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ q_lower: num [1:7, 1:4] 0.1049 0.0933 0.2773 0.2762 0.0933 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ q_upper: num [1:7, 1:4] 0.673 0.306 0.328 0.313 0.117 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "9" "10" "11" "12" ...
#>   .. ..$ : chr [1:4] "dim1" "dim2" "dim3" "dim4"
#>  $ draws  : 'draws_matrix' num [1:100, 1:28] -0.0164 0.6737 0.6598 0.6612 0.6677 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ draw    : chr [1:100] "1" "2" "3" "4" ...
#>   .. ..$ variable: chr [1:28] "node9[1]" "node9[2]" "node9[3]" "node9[4]" ...
#>   ..- attr(*, "nchains")= int 1
```

## Comparing model families

``` r

cmp <- run_bayesian_comparison(stan_data, tree, chains = 1, iter_warmup = 100,
                               iter_sampling = 100, refresh = 0)
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: normal_lpdf: Location parameter[1] is inf, but must be finite! (in '/tmp/Rtmp1wWU0u/model-1dc07a94686b.stan', line 56, column 4 to column 39)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1
#> Chain 1 finished in 0.1 seconds.
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: Matrix m is not positive definite (in '/tmp/Rtmp1wWU0u/model-1dc01388f553.stan', line 62, column 6 to column 58)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1
#> Chain 1 finished in 0.7 seconds.
cmp$root_comparison
#> $bm_root
#>        dim1        dim2        dim3        dim4 
#> 0.415469476 0.626433717 0.001769575 0.007975022 
#> 
#> $ou_root
#>        dim1        dim2        dim3        dim4 
#>  0.35081586  0.53592974 -0.03307235 -0.01349669 
#> 
#> $rmse
#> [1] 0.05925802
#> 
#> $cosine
#> [1] 0.9980271
```

## Pure-R reference

Every exported function is paired with a pure-R reference implementation
(`*_r`) that follows the math line-by-line. For example,
[`reconstruct_ancestral_bm_r()`](https://alrobles.github.io/ancemb/reference/reconstruct_ancestral_bm_r.md)
matches
[`phytools::fastAnc()`](https://rdrr.io/pkg/phytools/man/fastAnc.html)
for noise-free BM data, which is used in the package test suite as the
cross-language QA benchmark.
