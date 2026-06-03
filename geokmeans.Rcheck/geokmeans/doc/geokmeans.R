## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)
set.seed(1)

## ----library------------------------------------------------------------------
library(geokmeans)

## ----first--------------------------------------------------------------------
X <- rbind(
  matrix(rnorm(200, mean = 0), ncol = 2),
  matrix(rnorm(200, mean = 6), ncol = 2)
)

fit <- geo_kmeans(X, centers = 2)
fit

## ----structure----------------------------------------------------------------
str(fit)

## ----plot, fig.alt = "Two clusters coloured by assignment with centroids marked"----
plot(X, col = fit$cluster, pch = 19, cex = 0.6,
     xlab = "x1", ylab = "x2", main = "geo_kmeans result")
points(fit$centroids, pch = 8, cex = 2, lwd = 2)

## ----dispatch-----------------------------------------------------------------
kmeans_dc(X, centers = 2, method = "elkan")$centroids

## ----compare------------------------------------------------------------------
set.seed(42)
Y <- do.call(rbind, lapply(seq_len(6), function(i)
  matrix(rnorm(500, mean = 4 * i), ncol = 2)))

methods <- c("lloyd", "hamerly", "annulus", "exponion",
             "ball", "geokmeans")

comparison <- data.frame(
  method = methods,
  distance_calcs = vapply(methods, function(m) {
    kmeans_dc(Y, centers = 6, method = m, seed = 1)$distance_calculations
  }, numeric(1)),
  row.names = NULL
)
comparison[order(comparison$distance_calcs), ]

## ----seed---------------------------------------------------------------------
a <- geo_kmeans(X, centers = 2, seed = 7)
b <- geo_kmeans(X, centers = 2, seed = 7)
identical(a$centroids, b$centroids)

## ----custom-init--------------------------------------------------------------
init <- X[c(1, 101), ]
geo_kmeans(X, centers = init)$centroids

## ----realdata-----------------------------------------------------------------
path <- system.file("extdata", "Breastcancer.csv", package = "geokmeans")
bc <- as.matrix(read.csv(path, header = FALSE))
dim(bc)

bc_fit <- geo_kmeans(bc, centers = 2, seed = 1)
table(bc_fit$cluster)

## ----toomany, error = TRUE----------------------------------------------------
try({
D <- rbind(matrix(0.1, 50, 2), matrix(9, 50, 2))  # only 2 distinct rows
geo_kmeans(D, centers = 3)
})

## ----cite, eval = FALSE-------------------------------------------------------
# citation("geokmeans")

