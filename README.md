
# ancemb

[![R-CMD-check](https://github.com/alrobles/ancemb-devel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alrobles/ancemb-devel/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/alrobles/ancemb-devel/actions/workflows/pkgdown.yaml/badge.svg)](https://alrobles.github.io/ancemb/)

`ancemb` implements Bayesian phylogenetic reconstruction of ancestral
protein embeddings. It treats each embedding dimension as a continuous
trait evolving along a phylogeny under Brownian Motion (BM) or an
Ornstein-Uhlenbeck (OU) process, and fits pre-compiled Stan models via
the [`instantiate`](https://wlandau.github.io/instantiate/) package.

## Installation

`ancemb` requires a working
[CmdStan](https://mc-stan.org/users/interfaces/cmdstan) installation.
Install CmdStan from R with:

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

Then install `ancemb` from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("alrobles/ancemb-devel")
```

## Quick start

``` r
library(ancemb)
set.seed(42)
```

Prepare a small tree and some observed tip embeddings:

``` r
tree <- ape::rtree(8)
tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))

z_obs <- matrix(stats::rnorm(8 * 4, sd = 0.5), nrow = 8)
rownames(z_obs) <- tree$tip.label
```

Build the Stan data list and fit the Brownian Motion model:

``` r
stan_data <- prepare_embedding_data(tree, z_obs, obs_sigma = 0.1)
fit <- fit_embedding_bm(stan_data, chains = 2, iter_warmup = 300,
                        iter_sampling = 300, refresh = 0)
#> Running MCMC with 2 sequential chains...
#> 
#> Chain 1 finished in 0.2 seconds.
#> Chain 2 finished in 0.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 0.5 seconds.
```

Extract ancestral states and compare the inferred root with
`phytools::fastAnc`:

``` r
anc <- extract_ancestral(fit, data = stan_data, tree = tree)
root_node <- which(anc$node == stan_data$S + 1)

ml_root <- phytools::fastAnc(tree, z_obs[, 1])[1]
cat("ML root (dim 1):", ml_root, "\n")
#> ML root (dim 1): 0.5041557
cat("Bayes root mean (dim 1):", anc$mean[root_node, "dim1"], "\n")
#> Bayes root mean (dim 1): 0.4410044
```

## Reference implementation

The package follows the SSDLC software-factory pattern: every exported
Stan wrapper is paired with a pure-R reference function (`*_r`) that
reproduces the same math line-by-line. For example,
`reconstruct_ancestral_bm_r()` gives the same conditional expectations
as `phytools::fastAnc()` for noise-free data.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/code_of_conduct.html).
By participating in this project you agree to abide by its terms.
