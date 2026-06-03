# R Package for Geometric-K-means (geo-k-means)

Fast, eco-friendly k-means clustering for R. The package wraps high-performance
C++ implementations of seven k-means variants behind a single, uniform R
interface:

* **Geometric-k-means** — the bound-free method of Sharma et al. (2025)
* Lloyd's algorithm
* Elkan
* Hamerly
* Annulus
* Exponion
* Ball k-means++

All distance-heavy work runs in C++ (via `Rcpp` / `RcppEigen`), so you keep
native speed while calling everything from R.

## Installation

```r
# from a local source tarball
install.packages("geokmeans_0.1.0.tar.gz", repos = NULL, type = "source")

# or from GitHub once published
# remotes::install_github("parichit/Geometric-k-means", subdir = "...")
```

A C++17 compiler is required (R's standard toolchains qualify).

## Usage

```r
library(geokmeans)

set.seed(1)
X <- rbind(matrix(rnorm(200, 0), ncol = 2),
           matrix(rnorm(200, 6), ncol = 2))

fit <- geo_kmeans(X, centers = 2)
fit                      # print method: iterations, distances, centroids
fit$cluster              # per-point cluster ids
fit$centroids            # final centres

# pick the algorithm by name
kmeans_dc(X, centers = 2, method = "elkan")

# supply your own starting centroids
geo_kmeans(X, centers = X[c(1, 201), ])
```

Every function shares the same signature:

```r
geo_kmeans(data, centers, iter_max = 100, threshold = 1e-3,
           init = c("random", "sequential"), seed = 0,
           with_labels = TRUE, verbose = FALSE)
```

`centers` is either the number of clusters `k` or a matrix of initial
centroids (mirroring `stats::kmeans()`).

## Example data

Two small datasets ship with the package:

```r
bc <- as.matrix(read.csv(
  system.file("extdata", "Breastcancer.csv", package = "geokmeans"),
  header = FALSE))
geo_kmeans(bc, centers = 2)
```

## Reference

Sharma, P., Stanislaw, M., Kurban, H., Kulekci, O., & Dalkilic, M. (2025).
*Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means.*
<https://doi.org/10.1007/s10994-025-06891-1>
