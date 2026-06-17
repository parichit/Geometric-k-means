## Submission

This is a new submission of **geokmeans** (version 0.1.0). It provides fast C++
implementations of several k-means clustering algorithms (including the
Geometric-k-means method of Sharma et al., 2026,
<doi:10.1007/s10994-025-06891-1>) behind a uniform R interface, via 'Rcpp' and
'RcppEigen'.

## Resubmission

This is a resubmission addressing the points raised in the previous review:

* Removed the duplicated `Description:` label at the start of the `Description`
  field in DESCRIPTION.
* Wrapped the software name `C++` in single quotes in the `Description` field
  (`'Rcpp'` and `'RcppEigen'` were already quoted).
* No longer modify the global environment (`.GlobalEnv`): `R/kmeans-algos.R` no
  longer reads, assigns, or removes `.Random.seed` in `globalenv()`. The `seed`
  argument now defaults to `NULL` and only calls `set.seed()` when the user
  supplies a seed, leaving the RNG untouched otherwise.

## Test environments

* local macOS, R release
* win-builder, R-devel and R-release (`devtools::check_win_devel()`,
  `devtools::check_win_release()`)
* R CMD check run with `--as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission, so the "New submission" note is expected.

## Notes for the reviewer

* The package contains compiled C++. Cluster initialisation uses R's own random
  number generator (`unif_rand()`); the package does not call the C library
  `rand()`/`srand()`. Reproducibility is controlled from R via the `seed`
  argument and `set.seed()`.
* All console output from compiled code is routed through R's streams
  (`Rcpp::Rcout` / `Rcpp::Rcerr`); nothing is written to stdout/stderr directly.
* If a NOTE about installed size appears on some platforms, it comes from
  linking the header-only Eigen template library through 'RcppEigen'; the source
  package itself is small.
* Ball k-means is adapted from <https://github.com/syxiaa/ball-k-means>, credited
  in the source header; all code is redistributed under GPL-3.

## Downstream dependencies

There are currently no downstream dependencies (new submission).
