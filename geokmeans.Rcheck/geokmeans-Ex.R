pkgname <- "geokmeans"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('geokmeans')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("kmeans_algorithms")
### * kmeans_algorithms

flush(stderr()); flush(stdout())

### Name: kmeans_algorithms
### Title: k-Means clustering algorithms
### Aliases: kmeans_algorithms geo_kmeans lloyd_kmeans elkan_kmeans
###   hamerly_kmeans annulus_kmeans exponion_kmeans ball_kmeans

### ** Examples

set.seed(1)
X <- rbind(matrix(rnorm(100, 0), ncol = 2),
           matrix(rnorm(100, 5), ncol = 2))
fit <- geo_kmeans(X, centers = 2)
fit$centroids
table(fit$cluster)

# Supplying explicit starting centroids:
geo_kmeans(X, centers = X[c(1, 51), ])




cleanEx()
nameEx("kmeans_dc")
### * kmeans_dc

flush(stderr()); flush(stdout())

### Name: kmeans_dc
### Title: Run a k-means variant by name
### Aliases: kmeans_dc

### ** Examples

set.seed(1)
X <- rbind(matrix(rnorm(100, 0), ncol = 2),
           matrix(rnorm(100, 5), ncol = 2))
kmeans_dc(X, centers = 2, method = "elkan")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
