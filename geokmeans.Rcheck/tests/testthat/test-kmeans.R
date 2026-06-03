test_that("all algorithms recover two well-separated clusters", {
  set.seed(42)
  X <- rbind(matrix(rnorm(200, mean = 0), ncol = 2),
             matrix(rnorm(200, mean = 8), ncol = 2))
  algos <- list(geo_kmeans, lloyd_kmeans, elkan_kmeans, hamerly_kmeans,
                annulus_kmeans, exponion_kmeans, ball_kmeans)

  for (fn in algos) {
    fit <- fn(X, centers = 2, seed = 1)
    expect_s3_class(fit, "geokmeans")
    expect_equal(dim(fit$centroids), c(2L, 2L))
    expect_length(fit$cluster, nrow(X))
    expect_true(all(fit$cluster %in% 1:2))
    expect_gte(fit$iterations, 1L)
    # the two recovered centres should be far apart
    expect_gt(stats::dist(fit$centroids)[1], 5)
  }
})

test_that("matrix 'centers' are used as initial centroids", {
  set.seed(1)
  X <- rbind(matrix(rnorm(100, 0), ncol = 2),
             matrix(rnorm(100, 5), ncol = 2))
  init <- X[c(1, 51), ]
  fit <- geo_kmeans(X, centers = init)
  expect_equal(fit$k, 2L)
  expect_equal(dim(fit$centroids), c(2L, 2L))
})

test_that("dispatcher matches the direct call", {
  set.seed(7)
  X <- matrix(rnorm(60), ncol = 3)
  a <- kmeans_dc(X, centers = 2, method = "elkan", seed = 3)
  b <- elkan_kmeans(X, centers = 2, seed = 3)
  expect_equal(a$centroids, b$centroids)
})

test_that("input validation works", {
  X <- matrix(rnorm(20), ncol = 2)
  expect_error(geo_kmeans(X, centers = 0))
  expect_error(geo_kmeans(X, centers = 1000))
  Xna <- X; Xna[1, 1] <- NA
  expect_error(geo_kmeans(Xna, centers = 2))
})

test_that("with_labels = FALSE omits cluster assignments", {
  X <- matrix(rnorm(40), ncol = 2)
  fit <- geo_kmeans(X, centers = 2, with_labels = FALSE)
  expect_null(fit$cluster)
})

test_that("requesting more clusters than distinct rows errors clearly", {
  X <- rbind(matrix(0.1, 50, 2), matrix(9, 50, 2))  # only 2 distinct rows
  expect_error(geo_kmeans(X, centers = 3), "distinct")
  expect_error(lloyd_kmeans(X, centers = 5), "distinct")
})

test_that("degenerate duplicate-heavy data does not crash and stays consistent", {
  X <- rbind(matrix(0.1, 50, 2), matrix(9, 50, 2))  # 2 distinct rows, k = 2 allowed
  for (m in c("geokmeans", "lloyd", "elkan", "hamerly", "annulus",
              "exponion", "ball")) {
    fit <- suppressMessages(kmeans_dc(X, centers = 2, method = m, seed = 1))
    expect_s3_class(fit, "geokmeans")
    # returned object is self-consistent: labels match the number of centroids
    expect_equal(nrow(fit$centroids), fit$k)
    expect_setequal(sort(unique(fit$cluster)), seq_len(fit$k))
    expect_true(fit$k >= 1L && fit$k <= 2L)
  }
})

test_that("drop_empty = FALSE keeps the requested number of centroids", {
  X <- rbind(matrix(0.1, 50, 2), matrix(9, 50, 2))
  fit <- geo_kmeans(X, centers = 2, seed = 1, drop_empty = FALSE)
  expect_equal(nrow(fit$centroids), 2L)
})
