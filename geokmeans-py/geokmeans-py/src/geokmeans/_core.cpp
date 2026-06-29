// pybind11 bindings for the Geometric-k-means C++ library.
//
// This is the Python analogue of the package's Rcpp wrapper. A single
// dispatch function `_run` selects one of the seven algorithms, converts the
// incoming NumPy array to the format each routine expects, runs the compute
// with the GIL released, and returns centroids / iterations / distance count /
// labels back to Python.
//
// Two design points worth noting:
//
//  * Output handling. The C++ library prints progress via the unqualified
//    stream names `cout`/`cerr`. With `using namespace std;` below, those
//    resolve to std::cout/std::cerr, whose stream buffers we redirect at the
//    C++ level: to a null sink by default (quiet), or to a captured buffer
//    when verbose=true (flushed to Python's stdout afterwards). This is pure
//    C++ and therefore safe to do while the GIL is released.
//
//  * The GIL. Each compute call runs inside `py::gil_scoped_release`, so a
//    long fit on a large dataset does not freeze the interpreter: other Python
//    threads keep running, the Jupyter kernel stays responsive, and heavy jobs
//    can be dispatched to a background process/thread without blocking.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

// Standard headers the library relies on, included before the algorithm
// headers (mirrors the Rcpp wrapper's ordering).
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
#include <streambuf>

using namespace std;  // so the library's bare `cout`/`cerr` resolve to std::

#include "data_holder.hpp"
#include "IOutils.hpp"
#include "algo_utils.hpp"
#include "misc_utils.hpp"
#include "lloyd_kmeans.hpp"
#include "geokmeans.hpp"
#include "elkan.hpp"
#include "hamerly.hpp"
#include "annulus.hpp"
#include "exponion.hpp"
#include "ball_kmeans_xf.hpp"

namespace py = pybind11;

// --------------------------------------------------------------- output sink

// A streambuf that swallows everything written to it (quiet mode).
struct NullBuffer : std::streambuf {
  int overflow(int c) override { return c; }
};
static NullBuffer g_null_buffer;

// RAII redirect of std::cout/std::cerr. target == nullptr -> discard.
struct StreamRedirect {
  std::streambuf *old_out, *old_err;
  explicit StreamRedirect(std::streambuf *target) {
    old_out = std::cout.rdbuf(target);
    old_err = std::cerr.rdbuf(target);
  }
  ~StreamRedirect() {
    std::cout.rdbuf(old_out);
    std::cerr.rdbuf(old_err);
  }
};

// --------------------------------------------------------------- conversions

// NumPy (n x d, row-major double) -> vector<vector<float> >
static std::vector<std::vector<float> >
numpy_to_vv(const py::array_t<double, py::array::c_style | py::array::forcecast> &arr) {
  auto buf = arr.request();
  const ssize_t n = buf.shape[0], d = buf.shape[1];
  const double *p = static_cast<double *>(buf.ptr);
  std::vector<std::vector<float> > out(n, std::vector<float>(d));
  for (ssize_t i = 0; i < n; ++i)
    for (ssize_t j = 0; j < d; ++j)
      out[i][j] = static_cast<float>(p[i * d + j]);
  return out;
}

// NumPy -> Eigen::MatrixXf (Ball k-means takes Eigen).
static Eigen::MatrixXf
numpy_to_eigen(const py::array_t<double, py::array::c_style | py::array::forcecast> &arr) {
  auto buf = arr.request();
  const ssize_t n = buf.shape[0], d = buf.shape[1];
  const double *p = static_cast<double *>(buf.ptr);
  Eigen::MatrixXf E(n, d);
  for (ssize_t i = 0; i < n; ++i)
    for (ssize_t j = 0; j < d; ++j)
      E(i, j) = static_cast<float>(p[i * d + j]);
  return E;
}

// vector<vector<float> > -> NumPy (k x d, double).
static py::array_t<double> vv_to_numpy(const std::vector<std::vector<float> > &v) {
  const ssize_t k = v.size();
  const ssize_t d = k > 0 ? static_cast<ssize_t>(v[0].size()) : 0;
  py::array_t<double> out({k, d});
  auto buf = out.request();
  double *p = static_cast<double *>(buf.ptr);
  for (ssize_t i = 0; i < k; ++i)
    for (ssize_t j = 0; j < d; ++j)
      p[i * d + j] = v[i][j];
  return out;
}

// Eigen::MatrixXf -> NumPy (k x d, double).
static py::array_t<double> eigen_to_numpy(const Eigen::MatrixXf &E) {
  py::array_t<double> out({static_cast<ssize_t>(E.rows()),
                           static_cast<ssize_t>(E.cols())});
  auto buf = out.request();
  double *p = static_cast<double *>(buf.ptr);
  for (ssize_t i = 0; i < E.rows(); ++i)
    for (ssize_t j = 0; j < E.cols(); ++j)
      p[i * E.cols() + j] = E(i, j);
  return out;
}

// Nearest-centroid assignment. The algorithms return centroids but not point
// labels, so (as in the Rcpp wrapper) we recompute them. 0-based for Python.
static py::array_t<int>
assign_labels(const py::array_t<double, py::array::c_style | py::array::forcecast> &data,
              const std::vector<std::vector<float> > &centers) {
  auto buf = data.request();
  const ssize_t n = buf.shape[0], d = buf.shape[1];
  const ssize_t k = centers.size();
  const double *X = static_cast<double *>(buf.ptr);

  py::array_t<int> labels(n);
  int *out = static_cast<int *>(labels.request().ptr);

  for (ssize_t i = 0; i < n; ++i) {
    double best = std::numeric_limits<double>::infinity();
    int bestc = 0;
    for (ssize_t c = 0; c < k; ++c) {
      double s = 0.0;
      for (ssize_t j = 0; j < d; ++j) {
        const double t = X[i * d + j] - static_cast<double>(centers[c][j]);
        s += t * t;
      }
      if (s < best) { best = s; bestc = static_cast<int>(c); }
    }
    out[i] = bestc;
  }
  return labels;
}

// ------------------------------------------------------------------ dispatch

static py::dict run_algorithm(const std::string &algo,
                              py::array_t<double, py::array::c_style | py::array::forcecast> data,
                              int num_clusters, double threshold,
                              int num_iterations, int seed, bool verbose) {
  auto buf = data.request();
  if (buf.ndim != 2)
    throw std::invalid_argument("data must be a 2-D array (n_samples x n_features)");
  const int n_cols = static_cast<int>(buf.shape[1]);
  const std::string init_type = "random";

  output_data r;
  std::ostringstream captured;
  std::streambuf *sink = verbose ? captured.rdbuf()
                                 : static_cast<std::streambuf *>(&g_null_buffer);

  // Ball k-means needs Eigen; the other six take vector<vector<float> >.
  const bool is_ball = (algo == "ball");
  std::vector<std::vector<float> > X;
  Eigen::MatrixXf D;
  if (is_ball) D = numpy_to_eigen(data);
  else         X = numpy_to_vv(data);

  {
    StreamRedirect redirect(sink);
    py::gil_scoped_release release;  // keep Python responsive during compute
    const float thr = static_cast<float>(threshold);
    if      (algo == "lloyd")    r = lloyd_kmeans(X, num_clusters, thr, num_iterations, n_cols, init_type, seed);
    else if (algo == "elkan")    r = elkan_kmeans(X, num_clusters, thr, num_iterations, n_cols, init_type, seed);
    else if (algo == "hamerly")  r = hamerly_kmeans(X, num_clusters, thr, num_iterations, n_cols, init_type, seed);
    else if (algo == "annulus")  r = annulus(X, num_clusters, thr, num_iterations, n_cols, init_type, seed);
    else if (algo == "exponion") r = exponion(X, num_clusters, thr, num_iterations, n_cols, init_type, seed);
    else if (algo == "geo")      r = geokmeans(X, num_clusters, thr, num_iterations, n_cols, init_type, seed);
    else if (algo == "ball")     r = ball_k_means_Ring(D, false, num_clusters, threshold, num_iterations, init_type, seed);
    else throw std::invalid_argument("unknown algorithm: " + algo);
  }

  if (verbose) {
    const std::string s = captured.str();
    if (!s.empty()) py::print(s, py::arg("end") = "");
  }

  // Centroids: Ball returns them in a separate Eigen field.
  std::vector<std::vector<float> > centroid_vv;
  py::array_t<double> centroids;
  if (is_ball) {
    centroids = eigen_to_numpy(r.ballkm_centroids);
    const Eigen::MatrixXf &E = r.ballkm_centroids;
    centroid_vv.assign(E.rows(), std::vector<float>(E.cols()));
    for (int i = 0; i < E.rows(); ++i)
      for (int j = 0; j < E.cols(); ++j)
        centroid_vv[i][j] = E(i, j);
  } else {
    centroids = vv_to_numpy(r.centroids);
    centroid_vv = r.centroids;
  }

  py::array_t<int> labels = assign_labels(data, centroid_vv);

  py::dict out;
  out["centroids"]             = centroids;
  out["labels"]                = labels;
  out["n_iter"]                = r.loop_counter;
  out["distance_calculations"] = static_cast<unsigned long long>(r.num_dists);
  return out;
}

PYBIND11_MODULE(_core, m) {
  m.doc() = "Low-level C++ bindings for geokmeans (internal).";
  m.def("run", &run_algorithm,
        py::arg("algorithm"), py::arg("data"), py::arg("num_clusters"),
        py::arg("threshold"), py::arg("num_iterations"), py::arg("seed"),
        py::arg("verbose") = false,
        "Run one k-means variant and return a dict with centroids, labels, "
        "n_iter and distance_calculations.");
  m.attr("ALGORITHMS") = py::make_tuple(
      "geo", "lloyd", "elkan", "hamerly", "annulus", "exponion", "ball");
}
