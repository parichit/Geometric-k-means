# Internal helpers -----------------------------------------------------------

# Validate and coerce the input data to a numeric double matrix.
.prep_data <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    stop("'data' must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(data))
    stop("'data' must not contain missing values.", call. = FALSE)
  storage.mode(data) <- "double"
  data
}

# Resolve the 'centers' argument (mirrors stats::kmeans):
#   * a single number             -> that many clusters, initialised via 'init'
#   * a matrix / data frame       -> used as the initial centroids
# Returns a list(k, init_type, cleanup) where 'init_type' is either a keyword
# ("random"/"sequential") or a path to a temporary CSV holding the centroids.
.resolve_centers <- function(centers, data, init) {
  if (is.matrix(centers) || is.data.frame(centers)) {
    centers <- as.matrix(centers)
    storage.mode(centers) <- "double"
    if (ncol(centers) != ncol(data))
      stop("'centers' matrix must have the same number of columns as 'data'.",
           call. = FALSE)
    f <- tempfile(fileext = ".csv")
    utils::write.table(centers, f, sep = ",", row.names = FALSE,
                       col.names = FALSE)
    return(list(k = nrow(centers), init_type = f, cleanup = f))
  }
  if (length(centers) != 1L || centers < 1)
    stop("'centers' must be a single positive integer or a matrix of centroids.",
         call. = FALSE)
  k <- as.integer(centers)
  if (k > nrow(data))
    stop("'centers' cannot exceed the number of observations.", call. = FALSE)
  list(k = k, init_type = match.arg(init, c("random", "sequential")),
       cleanup = NULL)
}

# Common driver shared by every algorithm wrapper.
.run_kmeans <- function(cpp_fun, method, data, centers, iter_max, threshold,
                        init, seed, with_labels, verbose, drop_empty) {
  data <- .prep_data(data)
  res  <- .resolve_centers(centers, data, init)
  if (!is.null(res$cleanup)) on.exit(unlink(res$cleanup), add = TRUE)

  # Guard: you cannot form more non-empty clusters than there are distinct
  # observations. This is the genuinely ill-posed case (e.g. duplicate-heavy
  # data with k larger than the number of unique rows).
  n_distinct <- nrow(unique(data))
  if (res$k > n_distinct) {
    stop(sprintf(paste0(
      "Requested %d clusters, but 'data' has only %d distinct row(s).\n",
      "k-means cannot form more non-empty clusters than there are distinct ",
      "observations.\n",
      "Set 'centers' to at most %d, or provide data with at least %d distinct rows."),
      res$k, n_distinct, n_distinct, res$k), call. = FALSE)
  }

  # Random initialisation draws from R's RNG. A non-NULL 'seed' makes results
  # reproducible while leaving the caller's RNG stream untouched on return;
  # 'seed = NULL' uses the ambient RNG (i.e. honours a prior set.seed()).
  if (!is.null(seed)) {
    genv <- globalenv()
    if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = genv), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = genv)), add = TRUE)
    }
    set.seed(as.integer(seed))
  }

  seed_arg <- if (is.null(seed)) 0L else as.integer(seed)
  call_cpp <- function()
    cpp_fun(data, res$k, as.numeric(threshold), as.integer(iter_max),
            res$init_type, seed_arg, isTRUE(with_labels))

  # The C++ routines emit a short convergence message; capture it unless asked.
  out <- if (isTRUE(verbose)) call_cpp()
         else { invisible(utils::capture.output(out <- call_cpp())); out }

  out$method <- method
  out$k      <- res$k

  # Empty clusters: a centroid that captured no observations in the final
  # partition. They are not errors, but the returned object is more useful
  # without them. With drop_empty = TRUE we keep only non-empty clusters and
  # renumber the labels; the caller is told via a message.
  if (isTRUE(drop_empty) && !is.null(out$cluster)) {
    sizes <- tabulate(out$cluster, nbins = res$k)
    empty <- which(sizes == 0L)
    if (length(empty) > 0L) {
      keep  <- setdiff(seq_len(res$k), empty)
      remap <- integer(res$k)
      remap[keep]   <- seq_along(keep)
      out$cluster   <- remap[out$cluster]
      out$centroids <- out$centroids[keep, , drop = FALSE]
      out$k         <- length(keep)
      message(sprintf(
        paste0("%d of %d requested clusters were empty and have been dropped; ",
               "returning %d non-empty cluster(s).\n",
               "This usually means the data has limited distinct structure for ",
               "the requested 'centers'; try a smaller 'centers' or a different 'seed'."),
        length(empty), res$k, length(keep)))
    }
  }

  class(out) <- "geokmeans"
  out
}

# Algorithm wrappers ----------------------------------------------------------

#' k-Means clustering algorithms
#'
#' Run one of the bundled k-means variants on a numeric data matrix. All
#' functions share the same interface and return value; they differ only in the
#' acceleration strategy used internally. [geo_kmeans()] runs the bound-free
#' Geometric-k-means method.
#'
#' @param data A numeric matrix or data frame with observations in rows and
#'   features in columns. Missing values are not allowed.
#' @param centers Either a single positive integer giving the number of clusters
#'   `k`, or a numeric matrix of initial cluster centres (one centroid per row,
#'   with `ncol(centers) == ncol(data)`).
#' @param iter_max Maximum number of iterations.
#' @param threshold Convergence threshold on centroid movement.
#' @param init Initialisation strategy when `centers` is a number: `"random"`
#'   (random observations) or `"sequential"` (the first `k` observations).
#'   Ignored when `centers` is a matrix.
#' @param seed Integer seed for the random initialisation, or `NULL`.
#'   Initialisation uses R's random number generator: a non-`NULL` `seed` gives
#'   reproducible results (and the caller's RNG stream is restored afterwards),
#'   while `NULL` uses the ambient RNG so a preceding [set.seed()] applies.
#' @param with_labels Logical; if `TRUE` (default) the result includes a
#'   per-observation cluster assignment computed from the final centroids.
#' @param verbose Logical; if `TRUE`, print the algorithm's convergence message.
#' @param drop_empty Logical; if `TRUE` (default), clusters that end up with no
#'   assigned observations are removed from the result and the remaining cluster
#'   labels are renumbered, with a message. Requesting more clusters than the
#'   number of distinct rows in `data` is always an error.
#'
#' @return An object of class `"geokmeans"`: a list with components
#'   \describe{
#'     \item{centroids}{A `k x ncol(data)` matrix of final cluster centres.}
#'     \item{cluster}{Integer vector of cluster ids (1-based), if
#'       `with_labels = TRUE`.}
#'     \item{iterations}{Number of iterations performed.}
#'     \item{distance_calculations}{Total number of point-to-centroid distance
#'       computations.}
#'     \item{method}{The algorithm used.}
#'     \item{k}{The number of clusters.}
#'   }
#'
#' @references
#' Sharma, P., Stanislaw, M., Kurban, H., Kulekci, O., and Dalkilic, M. (2026).
#' Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means.
#' \doi{10.1007/s10994-025-06891-1}
#'
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(100, 0), ncol = 2),
#'            matrix(rnorm(100, 5), ncol = 2))
#' fit <- geo_kmeans(X, centers = 2)
#' fit$centroids
#' table(fit$cluster)
#'
#' # Supplying explicit starting centroids:
#' geo_kmeans(X, centers = X[c(1, 51), ])
#'
#' @name kmeans_algorithms
NULL

#' @rdname kmeans_algorithms
#' @export
geo_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                       init = c("random", "sequential"), seed = 0L,
                       with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_geo_kmeans, "geokmeans", data, centers, iter_max, threshold,
              init, seed, with_labels, verbose, drop_empty)
}

#' @rdname kmeans_algorithms
#' @export
lloyd_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                         init = c("random", "sequential"), seed = 0L,
                         with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_lloyd_kmeans, "lloyd", data, centers, iter_max, threshold,
              init, seed, with_labels, verbose, drop_empty)
}

#' @rdname kmeans_algorithms
#' @export
elkan_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                         init = c("random", "sequential"), seed = 0L,
                         with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_elkan_kmeans, "elkan", data, centers, iter_max, threshold,
              init, seed, with_labels, verbose, drop_empty)
}

#' @rdname kmeans_algorithms
#' @export
hamerly_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                           init = c("random", "sequential"), seed = 0L,
                           with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_hamerly_kmeans, "hamerly", data, centers, iter_max, threshold,
              init, seed, with_labels, verbose, drop_empty)
}

#' @rdname kmeans_algorithms
#' @export
annulus_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                           init = c("random", "sequential"), seed = 0L,
                           with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_annulus_kmeans, "annulus", data, centers, iter_max, threshold,
              init, seed, with_labels, verbose, drop_empty)
}

#' @rdname kmeans_algorithms
#' @export
exponion_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                            init = c("random", "sequential"), seed = 0L,
                            with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_exponion_kmeans, "exponion", data, centers, iter_max,
              threshold, init, seed, with_labels, verbose, drop_empty)
}

#' @rdname kmeans_algorithms
#' @export
ball_kmeans <- function(data, centers, iter_max = 100L, threshold = 1e-3,
                        init = c("random", "sequential"), seed = 0L,
                        with_labels = TRUE, verbose = FALSE, drop_empty = TRUE) {
  .run_kmeans(cpp_ball_kmeans, "ball", data, centers, iter_max, threshold,
              init, seed, with_labels, verbose, drop_empty)
}

#' Run a k-means variant by name
#'
#' A thin dispatcher over the individual algorithm functions.
#'
#' @inheritParams kmeans_algorithms
#' @param method The algorithm to use. One of `"geokmeans"`, `"lloyd"`,
#'   `"elkan"`, `"hamerly"`, `"annulus"`, `"exponion"`, `"ball"`.
#' @param ... Further arguments passed to the chosen algorithm.
#'
#' @return An object of class `"geokmeans"`; see [geo_kmeans()].
#'
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(100, 0), ncol = 2),
#'            matrix(rnorm(100, 5), ncol = 2))
#' kmeans_dc(X, centers = 2, method = "elkan")
#'
#' @export
kmeans_dc <- function(data, centers,
                      method = c("geokmeans", "lloyd", "elkan", "hamerly",
                                 "annulus", "exponion", "ball"), ...) {
  method <- match.arg(method)
  fun <- switch(method,
                geokmeans = geo_kmeans,
                lloyd     = lloyd_kmeans,
                elkan     = elkan_kmeans,
                hamerly   = hamerly_kmeans,
                annulus   = annulus_kmeans,
                exponion  = exponion_kmeans,
                ball      = ball_kmeans)
  fun(data, centers, ...)
}

#' @export
print.geokmeans <- function(x, ...) {
  cat(sprintf("k-means clustering (%s) with %d clusters\n", x$method, x$k))
  cat(sprintf("Iterations: %d | Distance computations: %.0f\n",
              x$iterations, x$distance_calculations))
  if (!is.null(x$cluster))
    cat("Cluster sizes: ",
        paste(tabulate(x$cluster, nbins = x$k), collapse = ", "), "\n", sep = "")
  cat("Centroids:\n")
  print(round(x$centroids, 4))
  invisible(x)
}
