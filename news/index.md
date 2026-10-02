# Changelog

## ancemb 0.1.0

- First public release under `alrobles/ancemb`.
- Optional reticulate backend adapter for NumPy/NPZ embedding arrays
  ([`configure_python_backend()`](https://alrobles.github.io/ancemb/reference/configure_python_backend.md),
  [`load_embedding_array()`](https://alrobles.github.io/ancemb/reference/load_embedding_array.md)).
- Trees without branch lengths now fail fast with an explicit error.
- Exported pure-R simulation and reconstruction references
  (`simulate_*_embedding_r()`, `reconstruct_ancestral_*_r()`).

## ancemb 0.0.1

- Initial CRAN-ready package skeleton.
- Exported
  [`prepare_embedding_data()`](https://alrobles.github.io/ancemb/reference/prepare_embedding_data.md),
  [`fit_embedding_bm()`](https://alrobles.github.io/ancemb/reference/fit_embedding_bm.md),
  [`fit_embedding_ou()`](https://alrobles.github.io/ancemb/reference/fit_embedding_ou.md),
  [`extract_ancestral()`](https://alrobles.github.io/ancemb/reference/extract_ancestral.md)
  and
  [`run_bayesian_comparison()`](https://alrobles.github.io/ancemb/reference/run_bayesian_comparison.md).
- Added pure-R reference implementations (`*_r`) for BM/OU simulation
  and ancestral reconstruction.
- Added cross-method `testthat` tests comparing Stan posterior means to
  [`phytools::fastAnc`](https://rdrr.io/pkg/phytools/man/fastAnc.html)
  and the pure-R reference.
- Set up `pkgdown` documentation and GitHub Actions for `R-CMD-check`
  and `pkgdown` deployment.
