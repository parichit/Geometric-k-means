"""BLAS-vectorised Geometric-k-means (the `geo` algorithm only).

Why this exists
---------------
The C++ core computes distances in scalar loops. Measured on sensor, that costs
~146 ns per distance for `geo` against ~2.1 ns for scikit-learn's chunked BLAS
GEMM -- a ~60x per-distance handicap that swamps geo's 24-34x reduction in the
*number* of distances. This module keeps geo's geometric pruning exactly as the
C++ implements it, but expresses every bulk numeric step as a matrix product so
it lands in BLAS.

Three kernels dominate and all three become GEMMs here:
  * initial assignment            -- full n x k distance matrix
  * inter-centroid distances      -- k x k
  * per-cluster candidate work    -- the affine/hyperplane test and the
                                     candidate distances, blocked per cluster

Semantics mirror `src/geokm_utils.hpp` exactly, including the details that
matter for equivalence:
  * distances are true Euclidean (not squared), as `calc_euclidean` returns
  * a cluster's radius is the running max distance *seen* this pass, which is
    an upper bound and is what the next pass's neighbour radius uses
  * a cluster that loses all its points gets the zero vector as its centroid,
    because `update_centroids` skips the division when the count is zero
  * convergence is `sqrt(S/k) <= tol` over the summed squared centroid shift

Distance accounting reports two numbers, because they differ here and the
difference is the whole point:
  * ``n_distance_calculations_``  -- distances the *algorithm* logically needs,
    i.e. what the C++ counter would report. Comparable across implementations.
  * ``n_distances_materialised_`` -- distances this implementation actually
    computes. Blocking by cluster means whole (points x neighbours) tiles are
    evaluated to keep the work in GEMM, so this is larger. The bet is that
    cheap-per-distance beats few-distances.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

__all__ = ["GeoKMeansBLAS"]


def _pairwise_distances(A: np.ndarray, B: np.ndarray,
                        A_sq: np.ndarray | None = None,
                        B_sq: np.ndarray | None = None) -> np.ndarray:
    """Euclidean distances between rows of A and rows of B, via one GEMM.

    Uses ||a-b||^2 = ||a||^2 - 2 a.b + ||b||^2 so the O(nmd) term is a single
    matrix product. The expansion can go slightly negative through rounding, so
    the result is clipped before the square root.
    """
    if A_sq is None:
        A_sq = np.einsum("ij,ij->i", A, A)
    if B_sq is None:
        B_sq = np.einsum("ij,ij->i", B, B)
    G = A @ B.T
    G *= -2.0
    G += A_sq[:, None]
    G += B_sq[None, :]
    np.maximum(G, 0.0, out=G)
    return np.sqrt(G, out=G)


def _segment_sum(X: np.ndarray, order: np.ndarray, starts: np.ndarray,
                 k: int) -> np.ndarray:
    """Per-label column sums of X, as one sparse-dense matrix product.

    `order` and `starts` (the label-sorted point indices and the k+1 group
    boundaries) are exactly a CSR indices/indptr pair for the k x n cluster
    indicator matrix, so the sums are `Indicator @ X` with no gather of X and
    no Python loop over features. `np.add.at` and a per-feature `np.bincount`
    were both profiled and are far slower.
    """
    n = X.shape[0]
    ind = sparse.csr_matrix(
        (np.ones(n, dtype=X.dtype), order, starts), shape=(k, n))
    return ind @ X


def _segment_max_sorted(values_sorted: np.ndarray, starts: np.ndarray,
                        k: int, dtype) -> np.ndarray:
    """Per-label maximum from label-sorted values, via `maximum.reduceat`.

    `starts` holds the k+1 group boundaries. Empty groups get 0. Avoids
    `np.maximum.at`, for the same reason as `_segment_sum`.
    """
    out = np.zeros(k, dtype=dtype)
    nonempty = np.flatnonzero(starts[1:] > starts[:-1])
    if nonempty.size:
        red = np.maximum.reduceat(values_sorted, starts[nonempty])
        out[nonempty] = red
    return out


class GeoKMeansBLAS:
    """Geometric-k-means with BLAS-backed inner kernels.

    Parameters
    ----------
    n_clusters : int
    init : ndarray of shape (n_clusters, n_features)
        Explicit initial centroids. Required -- this class exists for
        controlled comparisons, so there is no random initialisation path.
    tol : float, default 1e-4
        Convergence threshold on ``sqrt(S/k)``, matching the C++ core.
    max_iter : int, default 300
    dtype : np.dtype, default np.float64
        Working precision. The C++ core uses float32; float64 is the default
        here so timings are comparable with scikit-learn, which also uses
        float64.

    Attributes
    ----------
    cluster_centers_, labels_, n_iter_
    n_distance_calculations_ : int
        Logical distance count, comparable with the C++ counter.
    n_distances_materialised_ : int
        Distances actually evaluated, including those computed as part of a
        GEMM tile but discarded.
    """

    def __init__(self, n_clusters, *, init, tol=1e-4, max_iter=300,
                 dtype=np.float64):
        self.n_clusters = int(n_clusters)
        self.init = init
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.dtype = dtype

    def fit(self, X) -> "GeoKMeansBLAS":
        X = np.ascontiguousarray(X, dtype=self.dtype)
        n, d = X.shape
        k = self.n_clusters

        C = np.ascontiguousarray(self.init, dtype=self.dtype).copy()
        if C.shape != (k, d):
            raise ValueError(f"init must have shape ({k}, {d}), got {C.shape}")

        X_sq = np.einsum("ij,ij->i", X, X)

        n_logical = 0
        n_material = 0

        # --- initial assignment: full n x k distance matrix -----------------
        D = _pairwise_distances(X, C, X_sq)
        labels = D.argmin(axis=1).astype(np.int64)
        mind = D[np.arange(n), labels]
        n_logical += n * k
        n_material += n * k
        del D

        counts = np.bincount(labels, minlength=k)
        order = np.argsort(labels, kind="stable")
        starts = np.append(np.searchsorted(labels[order], np.arange(k)), n)
        radius = _segment_max_sorted(mind[order], starts, k, X.dtype)

        arange_k = np.arange(k)
        converged = False
        it = 0

        for it in range(1, self.max_iter + 1):
            # --- update centroids ------------------------------------------
            Cnew = _segment_sum(X, order, starts, k).astype(self.dtype, copy=False)
            nz = counts > 0
            Cnew[nz] /= counts[nz, None]
            Cnew[~nz] = 0.0
            # Clusters with no members keep the zero vector, mirroring
            # update_centroids, which skips the division when the count is 0.

            # --- convergence: sqrt(S/k) <= tol -----------------------------
            shift = Cnew - C
            S = float(np.einsum("ij,ij->", shift, shift))
            if np.sqrt(S / k) <= self.tol:
                C = Cnew
                converged = True
                break

            # --- find_neighbors --------------------------------------------
            Dc = _pairwise_distances(Cnew, Cnew)
            Dc *= 0.5
            np.fill_diagonal(Dc, np.inf)
            n_logical += k * (k - 1) // 2
            n_material += k * k
            closest = Dc.min(axis=1)

            # radius from the previous pass drives the neighbour search
            nbr_mask = Dc <= (radius + closest)[:, None]

            # --- determine_data_expression ---------------------------------
            diff = X - Cnew[labels]
            my_dist = np.sqrt(np.einsum("ij,ij->i", diff, diff))
            del diff
            n_logical += n
            n_material += n

            # order/starts were computed when `labels` last changed and serve
            # the centroid sum, the radius reduction and the grouping below.
            lab_ordered = labels[order]

            # radii are recomputed from scratch each pass, from the distance
            # of every point to its currently assigned centre
            radius = _segment_max_sorted(my_dist[order], starts, k, X.dtype)

            active = my_dist > closest[labels]
            active_sorted = order[active[order]]
            lab_sorted = lab_ordered[active[order]]
            bounds = np.append(np.searchsorted(lab_sorted, arange_k),
                               lab_sorted.size)

            new_labels = labels.copy()
            for c in range(k):
                idx = active_sorted[bounds[c]:bounds[c + 1]]
                if idx.size == 0:
                    continue
                N = np.flatnonzero(nbr_mask[c])
                if N.size == 0:
                    continue

                Xc = X[idx]
                CN = Cnew[N]
                # affine vector a = (c_j - c_i)/2, midpoint m = (c_i + c_j)/2
                A = (CN - Cnew[c]) * 0.5
                M = (CN + Cnew[c]) * 0.5

                # hyperplane test dot(x - m, a) > 0, as one GEMM
                side = Xc @ A.T
                side -= np.einsum("ij,ij->i", M, A)[None, :]

                cand = side > 0
                cand &= my_dist[idx][:, None] > Dc[c, N][None, :]

                rows = cand.any(axis=1)
                if not rows.any():
                    continue
                sub = np.flatnonzero(rows)
                Xs = Xc[sub]
                cand_s = cand[sub]

                Dblk = _pairwise_distances(Xs, CN, X_sq[idx[sub]])
                n_logical += int(cand_s.sum())
                n_material += Dblk.size

                # every examined neighbour updates its radius with the
                # distance seen, whether or not the point ends up moving.
                # Column-wise max over the candidate mask, in one pass.
                seen = np.where(cand_s, Dblk, -np.inf).max(axis=0)
                np.maximum(radius[N], seen, out=seen)
                radius[N] = seen

                Dm = np.where(cand_s, Dblk, np.inf)
                best = Dm.argmin(axis=1)
                bestd = Dm[np.arange(sub.size), best]
                moved = bestd < my_dist[idx[sub]]
                if moved.any():
                    new_labels[idx[sub][moved]] = N[best[moved]]

            labels = new_labels
            counts = np.bincount(labels, minlength=k)
            order = np.argsort(labels, kind="stable")
            starts = np.append(np.searchsorted(labels[order], arange_k), n)
            C = Cnew

        self.cluster_centers_ = np.asarray(C, dtype=np.float64)
        self.labels_ = labels.astype(np.int64)
        self.n_iter_ = it
        self.converged_ = converged
        self.n_distance_calculations_ = int(n_logical)
        self.n_distances_materialised_ = int(n_material)
        return self
