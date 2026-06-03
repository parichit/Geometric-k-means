// Rcpp wrappers for the Geometric-k-means C++ library.
//
// One [[Rcpp::export]] function per algorithm. Each takes an R numeric matrix
// (rows = observations, columns = features) plus parameters, calls the
// corresponding C++ routine, and returns an R list with the final centroids,
// the number of iterations, and the number of distance computations.
//
// CRAN note: the bundled headers print progress via the unqualified stream
// names `cout`/`cerr`. We pre-include every standard header the library uses
// (so <iostream>'s own `std::cout` declaration is processed first and guarded),
// then redirect the unqualified names to R's own output streams. This keeps all
// compiled output going through R instead of to stdout/stderr directly.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// Pull in the standard headers the library relies on *before* the redirect,
// so the macros below only affect the library's bare `cout`/`cerr` usages.
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <locale>
#include <map>
#include <numeric>
#include <algorithm>
#include <random>
#include <utility>

// Redirect library console output to R's streams (CRAN requirement).
#define cout Rcpp::Rcout
#define cerr Rcpp::Rcerr

#include "data_holder.h"
#include "IOutils.h"
#include "algo_utils.h"
#include "misc_utils.h"
#include "lloyd_kmeans.h"
#include "geokmeans.h"
#include "elkan.h"
#include "hamerly.h"
#include "annulus.h"
#include "exponion.h"
#include "ball_kmeans_xf.h"

using namespace Rcpp;

// ------------------------------------------------------------------ helpers

// R numeric matrix -> vector<vector<float> >
static std::vector<std::vector<float> > matrix_to_vv(const NumericMatrix& m) {
  const int n = m.nrow(), d = m.ncol();
  std::vector<std::vector<float> > out(n, std::vector<float>(d));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < d; ++j)
      out[i][j] = static_cast<float>(m(i, j));
  return out;
}

// vector<vector<float> > -> R numeric matrix
static NumericMatrix vv_to_matrix(const std::vector<std::vector<float> >& v) {
  const int k = static_cast<int>(v.size());
  const int d = k > 0 ? static_cast<int>(v[0].size()) : 0;
  NumericMatrix m(k, d);
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < d; ++j)
      m(i, j) = v[i][j];
  return m;
}

// Eigen::MatrixXf -> R numeric matrix
static NumericMatrix eigen_to_matrix(const Eigen::MatrixXf& E) {
  NumericMatrix m(E.rows(), E.cols());
  for (int i = 0; i < E.rows(); ++i)
    for (int j = 0; j < E.cols(); ++j)
      m(i, j) = E(i, j);
  return m;
}

// nearest-centroid assignment; returns 1-based cluster ids for R
static IntegerVector assign_clusters(const NumericMatrix& data,
                                     const NumericMatrix& centers) {
  const int n = data.nrow(), d = data.ncol(), k = centers.nrow();
  IntegerVector cl(n);
  for (int i = 0; i < n; ++i) {
    double best = R_PosInf;
    int bestc = 0;
    for (int c = 0; c < k; ++c) {
      double s = 0.0;
      for (int j = 0; j < d; ++j) {
        const double t = data(i, j) - centers(c, j);
        s += t * t;
      }
      if (s < best) { best = s; bestc = c; }
    }
    cl[i] = bestc + 1;
  }
  return cl;
}

static List make_result(const NumericMatrix& data, const NumericMatrix& centroids,
                        int iterations, double n_dist, bool with_labels) {
  List res = List::create(
    _["centroids"]             = centroids,
    _["iterations"]            = iterations,
    _["distance_calculations"] = n_dist);
  if (with_labels) res["cluster"] = assign_clusters(data, centroids);
  return res;
}

// Brackets a call so the library's R-based RNG (unif_rand) is valid; restores
// R's RNG state on exit. R's seed is set on the R side (see set.seed handling).
struct RNGGuard {
  RNGGuard()  { GetRNGstate(); }
  ~RNGGuard() { PutRNGstate(); }
};

// --------------------------------------------------- exported entry points

// [[Rcpp::export]]
List cpp_lloyd_kmeans(NumericMatrix data, int num_clusters, double threshold,
                      int num_iterations, std::string init_type, int seed,
                      bool with_labels) {
  RNGGuard rng_;  std::vector<std::vector<float> > X = matrix_to_vv(data);
  output_data r = lloyd_kmeans(X, num_clusters, static_cast<float>(threshold),
                               num_iterations, data.ncol(), init_type, seed);
  return make_result(data, vv_to_matrix(r.centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}

// [[Rcpp::export]]
List cpp_elkan_kmeans(NumericMatrix data, int num_clusters, double threshold,
                      int num_iterations, std::string init_type, int seed,
                      bool with_labels) {
  RNGGuard rng_;  std::vector<std::vector<float> > X = matrix_to_vv(data);
  output_data r = elkan_kmeans(X, num_clusters, static_cast<float>(threshold),
                               num_iterations, data.ncol(), init_type, seed);
  return make_result(data, vv_to_matrix(r.centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}

// [[Rcpp::export]]
List cpp_hamerly_kmeans(NumericMatrix data, int num_clusters, double threshold,
                        int num_iterations, std::string init_type, int seed,
                        bool with_labels) {
  RNGGuard rng_;  std::vector<std::vector<float> > X = matrix_to_vv(data);
  output_data r = hamerly_kmeans(X, num_clusters, static_cast<float>(threshold),
                                 num_iterations, data.ncol(), init_type, seed);
  return make_result(data, vv_to_matrix(r.centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}

// [[Rcpp::export]]
List cpp_annulus_kmeans(NumericMatrix data, int num_clusters, double threshold,
                        int num_iterations, std::string init_type, int seed,
                        bool with_labels) {
  RNGGuard rng_;  std::vector<std::vector<float> > X = matrix_to_vv(data);
  output_data r = annulus(X, num_clusters, static_cast<float>(threshold),
                          num_iterations, data.ncol(), init_type, seed);
  return make_result(data, vv_to_matrix(r.centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}

// [[Rcpp::export]]
List cpp_exponion_kmeans(NumericMatrix data, int num_clusters, double threshold,
                         int num_iterations, std::string init_type, int seed,
                         bool with_labels) {
  RNGGuard rng_;  std::vector<std::vector<float> > X = matrix_to_vv(data);
  output_data r = exponion(X, num_clusters, static_cast<float>(threshold),
                           num_iterations, data.ncol(), init_type, seed);
  return make_result(data, vv_to_matrix(r.centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}

// [[Rcpp::export]]
List cpp_geo_kmeans(NumericMatrix data, int num_clusters, double threshold,
                    int num_iterations, std::string init_type, int seed,
                    bool with_labels) {
  RNGGuard rng_;  std::vector<std::vector<float> > X = matrix_to_vv(data);
  output_data r = geokmeans(X, num_clusters, static_cast<float>(threshold),
                            num_iterations, data.ncol(), init_type, seed);
  return make_result(data, vv_to_matrix(r.centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}

// [[Rcpp::export]]
List cpp_ball_kmeans(NumericMatrix data, int num_clusters, double threshold,
                     int num_iterations, std::string init_type, int seed,
                     bool with_labels) {
  RNGGuard rng_;  const int n = data.nrow(), d = data.ncol();
  Eigen::MatrixXf D(n, d);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < d; ++j)
      D(i, j) = static_cast<float>(data(i, j));
  // 'detail = false' suppresses the routine's own progress printing.
  output_data r = ball_k_means_Ring(D, false, num_clusters, threshold,
                                    num_iterations, init_type, seed);
  return make_result(data, eigen_to_matrix(r.ballkm_centroids), r.loop_counter,
                     static_cast<double>(r.num_dists), with_labels);
}
