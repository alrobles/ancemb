## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... [2s] NOTE
  New submission.
  
  Suggests or Enhances not in mainstream repositories: cmdstanr
  Availability using Additional_repositories specification: cmdstanr yes https://mc-stan.org/r-packages/

  The package requires CmdStan and uses the `cmdstanr` interface.
  `cmdstanr` is available via the Stan r-universe listed in
  `Additional_repositories`.

  Several URLs reported 404 because the public pkgdown site and GitHub
  Pages are not yet live. These will resolve once the repository is made
  public and `pkgdown` deploys to `https://alrobles.github.io/ancemb/`.

* checking for future file timestamps ... NOTE
  unable to verify current time

  This is an environment-specific NOTE and not a package issue.

## Test environments

- Local: Ubuntu 24.04, R 4.1.2, CmdStan 2.39.0
- GitHub Actions: ubuntu-latest (R release) via `r-lib/actions`
