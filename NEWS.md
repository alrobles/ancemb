# ancemb 0.0.1

* Initial CRAN-ready package skeleton.
* Exported `prepare_embedding_data()`, `fit_embedding_bm()`,
  `fit_embedding_ou()`, `extract_ancestral()` and `run_bayesian_comparison()`.
* Added pure-R reference implementations (`*_r`) for BM/OU simulation and
  ancestral reconstruction.
* Added cross-method `testthat` tests comparing Stan posterior means to
  `phytools::fastAnc` and the pure-R reference.
* Set up `pkgdown` documentation and GitHub Actions for `R-CMD-check` and
  `pkgdown` deployment.
