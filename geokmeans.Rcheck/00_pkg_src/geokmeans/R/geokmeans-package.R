#' geokmeans: Fast and Eco-Friendly k-Means Clustering Algorithms
#'
#' Fast C++ implementations of several k-means clustering algorithms exposed to
#' R through a uniform interface: Lloyd's algorithm, Elkan, Hamerly, Annulus,
#' Exponion, Ball k-means, and the bound-free Geometric-k-means method.
#'
#' The main entry points are [geo_kmeans()], [lloyd_kmeans()], [elkan_kmeans()],
#' [hamerly_kmeans()], [annulus_kmeans()], [exponion_kmeans()], [ball_kmeans()],
#' and the dispatcher [kmeans_dc()].
#'
#' @references
#' Sharma, P., Stanislaw, M., Kurban, H., Kulekci, O., and Dalkilic, M. (2026).
#' Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means.
#' \doi{10.1007/s10994-025-06891-1}
#'
#' @keywords internal
#' @useDynLib geokmeans, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils write.table capture.output
"_PACKAGE"
