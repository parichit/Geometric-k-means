# geokmeans 0.1.0

* Initial CRAN release.

* Provides seven fast *k*-means clustering algorithms behind a single, uniform
  interface, wrapping high-performance C++ implementations via 'Rcpp' and
  'RcppEigen':
  * `geo_kmeans()` — Geometric-*k*-means, the bound-free method of Sharma et al.
    (2026) <doi:10.1007/s10994-025-06891-1>.
  * `ball_kmeans()` — Ball *k*-means++.
  * `lloyd_kmeans()`, `elkan_kmeans()`, `hamerly_kmeans()`, `annulus_kmeans()`,
    and `exponion_kmeans()`.
  * `kmeans_dc()` — dispatcher to select any of the above by name.

* `centers` accepts either the number of clusters or a matrix of initial
  centroids (mirroring `stats::kmeans()`); initialisation can be `"random"` or
  `"sequential"`.

* Returns a `geokmeans` object with the final centroids, per-point cluster
  assignments, iteration count, and number of distance computations, along with
  a `print()` method.

* Random initialisation uses R's random number generator and is reproducible via
  `set.seed()` or the `seed` argument, without disturbing the caller's RNG
  stream.

* Safeguards for degenerate input: an informative error when more clusters are
  requested than there are distinct observations, and optional removal of empty
  clusters via `drop_empty`.

* Ships two example datasets in `inst/extdata` (`Breastcancer.csv` and
  `CreditRisk.csv`) and a "Getting started with geokmeans" vignette.
